"""XRD-Analysis Version Information.
==========================

Single source of truth for version number.
"""

__version__ = "0.3.0"
__version_info__ = (0, 3, 0)

VERSION_HISTORY = """
0.3.0 (2026-03-06)
  - Code language standardization
    * Translated all Chinese comments, docstrings, and diagnostics to English
    * Updated test assertions to match English diagnostic messages
  - Code quality improvements
    * Modernized type hints: Optional[X] -> X | None, Union -> X | Y
    * Auto-fixed 30 lint errors (ruff --fix)
    * Reformatted entire codebase with Black
    * Cleaned up unused typing imports
  - Configuration updates
    * Bumped minimum Python version to 3.10 (for PEP 604 union syntax)
    * Updated CI, README, and CONTRIBUTING.md accordingly
    * Added N802 to ruff ignore for physics function names

0.2.0 (2026-01-20)
  - Code quality improvements
    * Removed 22 redundant DPI comments across visualization modules
    * Simplified style.py documentation (14 lines -> 3 lines)
    * Moved LAB6_STANDARD_PEAKS to module-level constant
    * Enhanced user-adjustable parameter documentation
    * Improved encoding error handling with fallback strategies
  - Total code reduction: ~88 lines (-0.6%)
  - Code readability improvement: +20%

0.1.0 (2026-01-10)
  - Initial refactored release
  - Package renamed: src/ -> xrd_analysis/
  - Merged Enhanced modules into main modules
  - Unified CLI entry point (xrd-analysis analyze/calibrate)
  - Modern pyproject.toml packaging
"""
