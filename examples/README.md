# Example patterns

Simulated patterns for the Quick Start, written by `python scripts/make_example_data.py`. They are not measurements, so the true answer is known.

| File | Content | Range |
|---|---|---|
| `lab6_synthetic.txt` | LaB6 standard (a = 4.156826 Å, NIST SRM 660c), instrument broadening only | 20–80° 2θ |
| `cu_synthetic.txt` | Cu, 40 nm cube-shaped crystallites, no microstrain | 40–100° 2θ |

Both are measured on the same simulated instrument: a Gaussian width from the Caglioti equation with U = 0.01, V = −0.005, W = 0.004 deg² (0.06° at 43°, 0.10° at 95°). Every reflection is a Cu Kα1 + Kα2 pair (Bragg positions for each wavelength, intensity ratio 0.5) with an exact Voigt profile from `scipy.special.voigt_profile`, which is independent of the fitting code. The Lorentzian part is the Scherrer size width K λ / (D cos θ), with K for cube-shaped grains (Langford & Wilson 1978). A smooth background and Poisson counting noise (fixed seed) are added; the step is 0.02°.

Files are Bruker-style TXT: comment lines, then a `[Data]` section with `Angle,PSD` columns.

What the Quick Start returns for these files, and why, is described in the main [README](../README.md#run-the-example); `tests/test_examples.py` runs the same two commands.
