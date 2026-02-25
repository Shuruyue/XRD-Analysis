"""Custom Exception Hierarchy 自訂例外體系
==========================================

Centralised exception classes for XRD-Analysis.
集中管理的 XRD-Analysis 例外類別。

Usage 用法::

    from xrd_analysis.core.exceptions import FittingError

    raise FittingError("Doublet fitting failed to converge")
"""


class XRDAnalysisError(Exception):
    """Base exception for all XRD-Analysis errors.

    所有 XRD-Analysis 錯誤的基底例外類別。
    """


class DataLoadError(XRDAnalysisError):
    """Failed to load or parse XRD data file.

    無法載入或解析 XRD 資料檔案。
    """


class CalibrationError(XRDAnalysisError):
    """Instrument calibration error.

    儀器校準錯誤（如 Caglioti 參數無效）。
    """


class FittingError(XRDAnalysisError):
    """Peak fitting convergence or numerical error.

    峰擬合收斂或數值錯誤。
    """


class FittingConvergenceError(FittingError):
    """Fitting failed to converge within maximum iterations.

    擬合在最大迭代次數內未能收斂。
    """


class ValidationError(XRDAnalysisError):
    """Data or result validation failed.

    資料或結果驗證失敗。
    """


class ConfigurationError(XRDAnalysisError):
    """Invalid configuration parameter.

    無效的配置參數。
    """


class PreprocessingError(XRDAnalysisError):
    """Error during data preprocessing.

    資料預處理過程中的錯誤。
    """
