"""XRD Data Loader XRD 數據加載器
================================

Load and parse XRD data files from various instrument formats.
從各種儀器格式加載和解析 XRD 數據文件。
"""

import re
from pathlib import Path
from typing import Any, Dict, Tuple

import numpy as np


def load_bruker_txt(filepath: str) -> Tuple[np.ndarray, np.ndarray]:
    """Load Bruker TXT format XRD data.

    Returns:
        Tuple of (two_theta, intensity) arrays

    """
    two_theta = []
    intensity = []
    in_data_section = False

    with open(filepath, encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()

            if "[Data]" in line:
                in_data_section = True
                continue

            if in_data_section and line:
                # Skip header row
                if "Angle" in line or "PSD" in line:
                    continue

                # Parse data
                parts = line.replace(",", " ").split()
                if len(parts) >= 2:
                    try:
                        theta = float(parts[0])
                        counts = float(parts[1])
                        two_theta.append(theta)
                        intensity.append(counts)
                    except ValueError:
                        continue

    return np.array(two_theta), np.array(intensity)


def parse_filename(filepath: str) -> Dict[str, Any]:
    """Parse sample info from filename.

    Format: YYYYMMDD_Xml_Xh.txt or YYYYMMDD_Xml_Xh_Xmin.txt
    """
    name = Path(filepath).stem

    result: Dict[str, Any] = {
        "name": name,
        "concentration_ml": None,
        "time_hours": None,
    }

    # Extract concentration (e.g., "0ml", "4.5ml", "9ml", "18ml")
    conc_match = re.search(r"(\d+\.?\d*)ml", name)
    if conc_match:
        result["concentration_ml"] = float(conc_match.group(1))

    # Extract time (e.g., "2h", "0h_15min", "24h")
    time_hours = 0
    hour_match = re.search(r"(\d+)h", name)
    if hour_match:
        time_hours = int(hour_match.group(1))

    min_match = re.search(r"(\d+)min", name)
    if min_match:
        time_hours += int(min_match.group(1)) / 60

    result["time_hours"] = time_hours

    return result
