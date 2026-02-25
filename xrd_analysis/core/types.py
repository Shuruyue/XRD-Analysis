"""Common Type Aliases 常見型別別名
==================================

Centralised type definitions to avoid repetition across the codebase.
集中化的型別定義，避免在代碼庫中重複。
"""

from typing import Tuple

# Miller indices (h, k, l) — used extensively across analysis methods
HKL = Tuple[int, int, int]
