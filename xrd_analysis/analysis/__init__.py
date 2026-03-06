"""Analysis Pipeline Module.
=========================

High-level analysis orchestration for XRD data.
"""

from xrd_analysis.analysis.pipeline import (
    AnalysisConfig,
    PipelineResult,
    XRDAnalysisPipeline,
    batch_analyze,
    run_full_analysis,
)

# Alias for backward compatibility
AnalysisPipeline = XRDAnalysisPipeline

__all__ = [
    "XRDAnalysisPipeline",
    "AnalysisPipeline",
    "AnalysisConfig",
    "PipelineResult",
    "run_full_analysis",
    "batch_analyze",
]
