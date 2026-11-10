"""Check the physics of the interactive page web/crystallite-size.html against the Python package.

The page contains a block of pure JavaScript functions between the markers
"// ---- physics: begin" and "// ---- physics: end". This script runs that block in Node.js for
random parameter sets and compares every number with xrd_analysis:

- peak positions with Bragg's law and the JCPDS table of the package,
- Scherrer constants with get_k_for_hkl,
- the TCH width and eta, and its inverse, with xrd_analysis.fitting.pseudo_voigt,
- the instrument width, corrected sample width, size and validity flag with ScherrerCalculator.

Usage: python scripts/verify_web_page.py [--cases 2000]
Requires Node.js on the PATH.
"""

from __future__ import annotations

import argparse
import json
import math
import random
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from xrd_analysis.core.copper_crystal import CU_JCPDS_EXTENDED, get_k_for_hkl  # noqa: E402
from xrd_analysis.fitting.pseudo_voigt import (  # noqa: E402
    tch_components_from_eta,
    tch_eta_from_components,
)
from xrd_analysis.methods.scherrer import ScherrerCalculator  # noqa: E402

PAGE = ROOT / "web" / "crystallite-size.html"
PEAKS = [(1, 1, 1), (2, 0, 0), (2, 2, 0), (3, 1, 1), (2, 2, 2)]
METHODS = ["quadratic", "voigt", "none"]
TOLERANCE = 1e-9

NODE_HARNESS = r"""
const fs = require("fs"), vm = require("vm");
const page = fs.readFileSync(process.argv[1], "utf8");  // with -e the first argument is argv[1]
const after = page.split("// ---- physics: begin")[1], block = after.slice(after.indexOf("\n")).split("// ---- physics: end")[0];
const ctx = {}; vm.createContext(ctx);
vm.runInContext(block + "\n;this.api = { kFor, twoTheta, caglioti, tchEta, tchComponents, forward, scherrer, CU_KA1, CU_KA2, A_CU };", ctx);
const api = ctx.api, cases = JSON.parse(fs.readFileSync(0, "utf8")), out = [];
for (const c of cases) {
  const kTrue = api.kFor(c.hkl, c.shape === "cube");
  const fw = api.forward(c.hkl, c.D, c.eps, kTrue, c.U, c.V, c.W);
  const an = api.scherrer(fw.tt, fw.fObs, fw.fInst, c.hkl, c.method, fw.eta, c.kmode === "cube");
  out.push({ tt: fw.tt, fInst: fw.fInst, fObs: fw.fObs, eta: fw.eta, fL: fw.fL, fGeps: fw.fGeps,
             K: an.K, size: an.size, fSample: an.fSample, flag: an.flag,
             inv: api.tchComponents(c.f, c.etaIn) });
}
process.stdout.write(JSON.stringify(out));
"""


def random_cases(n, seed=0):
    rng = random.Random(seed)
    cases = []
    for _ in range(n):
        cases.append({
            "hkl": list(rng.choice(PEAKS)),
            "D": 10 ** rng.uniform(math.log10(2), math.log10(300)),
            "eps": rng.choice([0.0, rng.uniform(0, 0.005)]),
            "U": rng.choice([0.0, rng.uniform(0, 0.05)]),
            "V": rng.uniform(-0.03, 0.01),
            "W": rng.uniform(0, 0.02),
            "shape": rng.choice(["cube", "sphere"]),
            "kmode": rng.choice(["cube", "sphere"]),
            "method": rng.choice(METHODS),
            "f": rng.uniform(0.01, 2.0),
            "etaIn": rng.uniform(0, 1),
        })
    return cases


def reference(c):
    """The same quantities computed with xrd_analysis."""
    h, k, l = c["hkl"]
    lam = 1.540562
    d = 3.6150 / math.sqrt(h * h + k * k + l * l)
    tt = 2 * math.degrees(math.asin(lam / (2 * d)))
    theta = math.radians(tt / 2)
    k_true = get_k_for_hkl(h, k, l, use_cubic_habit=c["shape"] == "cube")
    f_l = math.degrees(k_true * lam / (c["D"] * 10 * math.cos(theta)))
    f_g_eps = math.degrees(4 * c["eps"] * math.tan(theta))
    calc = ScherrerCalculator(use_cubic_habit=c["kmode"] == "cube", caglioti_params=(c["U"], c["V"], c["W"]))
    f_inst = float(calc._calculate_caglioti(tt))
    f_obs, eta = tch_eta_from_components(math.hypot(f_g_eps, f_inst), f_l)
    if c["method"] == "none":
        res = calc.calculate(tt, f_obs, fwhm_instrumental=0.0, hkl=(h, k, l), correction_method="quadratic")
    else:
        res = calc.calculate(tt, f_obs, fwhm_instrumental=f_inst, hkl=(h, k, l), correction_method=c["method"],
                             eta_observed=eta, eta_instrumental=0.0)
    return {"tt": tt, "fInst": f_inst, "fObs": f_obs, "eta": eta, "fL": f_l, "fGeps": f_g_eps,
            "K": res.k_factor, "size": res.size_nm, "fSample": res.fwhm_sample, "flag": res.validity_flag.value.lower(),
            "inv": list(tch_components_from_eta(c["f"], c["etaIn"]))}


def rel(a, b):
    """Relative difference; NaN on both sides (JSON null from the page) counts as equal."""
    a_nan = a is None or math.isnan(a)
    b_nan = b is None or math.isnan(b)
    if a_nan or b_nan:
        return 0.0 if a_nan and b_nan else math.inf
    return abs(a - b) / max(abs(b), 1e-300) if b != 0 else abs(a)


def run(n_cases=2000, page=PAGE):
    node = shutil.which("node")
    if node is None:
        raise RuntimeError("Node.js is needed to run the page's JavaScript")
    cases = random_cases(n_cases)
    out = subprocess.run([node, "-e", NODE_HARNESS, str(page)], input=json.dumps(cases), capture_output=True,
                         text=True, check=False)
    if out.returncode != 0:
        raise RuntimeError("Node.js failed:\n" + out.stderr)
    got = json.loads(out.stdout)
    worst, flags_differ = {}, 0
    for c, js in zip(cases, got):
        ref = reference(c)
        for key in ("tt", "fInst", "fObs", "eta", "fL", "fGeps", "K", "size", "fSample"):
            worst[key] = max(worst.get(key, 0.0), rel(js[key], ref[key]))
        worst["inverse TCH"] = max(worst.get("inverse TCH", 0.0), rel(js["inv"][0], ref["inv"][0]),
                                   rel(js["inv"][1], ref["inv"][1]))
        flags_differ += js["flag"] != ref["flag"]
    # peak positions of the page against the JCPDS table of the package (rounded to 0.001 degree there)
    table = max(abs(reference({"hkl": list(p), "D": 30, "eps": 0, "U": 0, "V": 0, "W": 0, "shape": "cube",
                                "kmode": "cube", "method": "none", "f": 1, "etaIn": 0.5})["tt"]
                    - CU_JCPDS_EXTENDED[p]["two_theta"]) for p in PEAKS)
    return worst, flags_differ, table


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cases", type=int, default=2000)
    args = parser.parse_args()
    worst, flags_differ, table = run(args.cases)
    print(f"{args.cases} random parameter sets, largest relative difference page vs xrd_analysis:")
    for key, value in worst.items():
        print(f"  {key:12s} {value:.1e}")
    print(f"validity flags that differ: {flags_differ}")
    print(f"largest difference of the Bragg 2-theta from the package's JCPDS table: {table:.4f} deg")
    ok = max(worst.values()) <= TOLERANCE and flags_differ == 0 and table < 0.002
    print("OK" if ok else "MISMATCH")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
