"""
Comprehensive Report Generator 綜合報告產生器
============================================

Final report aggregating all xrd_analysis analysis results (Phase 04-07).
彙總所有 xrd_analysis 分析結果 (Phase 04-07) 的最終報告。
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Tuple, List, Dict, Any
from datetime import datetime


# =============================================================================
# Comprehensive Analysis Result
# =============================================================================

@dataclass
class ComprehensiveResult:
    """
    Complete xrd_analysis analysis result aggregating all phases.
    
    Contains results from:
    - Phase 04: Scherrer crystallite size
    - Phase 05: Williamson-Hall strain analysis  
    - Phase 06: Texture analysis
    - Phase 07: Defect and stress analysis
    """
    # Sample metadata
    sample_name: str = "Unknown"
    sample_age_hours: Optional[float] = None
    analysis_date: str = field(default_factory=lambda: datetime.now().strftime("%Y-%m-%d %H:%M"))
    
    # Phase 04: Scherrer
    scherrer_size_nm: Optional[float] = None
    scherrer_validity: str = ""
    scherrer_k_used: Optional[float] = None
    
    # Phase 05: Williamson-Hall
    wh_size_nm: Optional[float] = None
    wh_strain: Optional[float] = None
    wh_r_squared: Optional[float] = None
    wh_quality: str = ""
    
    # Phase 06: Texture (DATA ONLY)
    dominant_orientation: Optional[Tuple[int, int, int]] = None
    dominant_tc: Optional[float] = None
    is_random_texture: bool = True
    
    # Phase 07: Defects
    peak_separation_deg: Optional[float] = None
    stacking_fault_alpha: Optional[float] = None
    stacking_fault_severity: str = ""
    lattice_constant: Optional[float] = None
    lattice_status: str = ""
    annealing_state: str = "unknown"
    
    # Overall
    recommendations: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)


# =============================================================================
# Report Generator
# =============================================================================

def generate_comprehensive_report(result: ComprehensiveResult) -> str:
    """
    Generate comprehensive xrd_analysis analysis report.
    產生 xrd_analysis 綜合分析報告。
    
    Args:
        result: ComprehensiveResult with all analysis data
        
    Returns:
        Formatted report string
    """
    lines = [
        "=" * 60,
        "xrd_analysis 電鍍銅微結構綜合分析報告",
        "=" * 60,
        "",
        f"樣品名稱：{result.sample_name}",
        f"分析日期：{result.analysis_date}",
    ]
    
    if result.sample_age_hours is not None:
        lines.append(f"樣品存放：{result.sample_age_hours:.1f} 小時")
    else:
        lines.append("樣品存放：未知")
    
    # Phase 04: Scherrer
    lines.extend([
        "",
        "-" * 40,
        "【Phase 04: Scherrer 晶粒尺寸】",
        "-" * 40,
    ])
    
    if result.scherrer_size_nm is not None:
        lines.append(f"  晶粒尺寸 D = {result.scherrer_size_nm:.1f} nm")
        if result.scherrer_k_used:
            lines.append(f"  K 值 = {result.scherrer_k_used:.3f}")
        lines.append(f"  狀態：{result.scherrer_validity or 'N/A'}")
    else:
        lines.append("  未執行 Scherrer 分析")
    
    # Phase 05: Williamson-Hall
    lines.extend([
        "",
        "-" * 40,
        "【Phase 05: Williamson-Hall 應變分析】",
        "-" * 40,
    ])
    
    if result.wh_size_nm is not None:
        lines.append(f"  晶粒尺寸 D = {result.wh_size_nm:.1f} nm")
        if result.wh_strain is not None:
            lines.append(f"  微觀應變 ε = {result.wh_strain:.2e}")
        if result.wh_r_squared is not None:
            lines.append(f"  R² = {result.wh_r_squared:.3f}")
        lines.append(f"  品質：{result.wh_quality or 'N/A'}")
    else:
        lines.append("  未執行 W-H 分析")
    
    # Phase 06: Texture
    lines.extend([
        "",
        "-" * 40,
        "【Phase 06: 織構分析】",
        "-" * 40,
    ])
    
    if result.dominant_orientation:
        hkl = result.dominant_orientation
        hkl_str = f"({hkl[0]}{hkl[1]}{hkl[2]})"
        lines.append(f"  主要取向：{hkl_str}")
        if result.dominant_tc:
            lines.append(f"  TC = {result.dominant_tc:.2f}")
        lines.append(f"  隨機性：{'是' if result.is_random_texture else '否'}")
    else:
        lines.append("  未執行織構分析")
    
    # Phase 07: Defects
    lines.extend([
        "",
        "-" * 40,
        "【Phase 07: 缺陷與應力分析】",
        "-" * 40,
    ])
    
    # Stacking faults
    if result.peak_separation_deg is not None:
        lines.append(f"  峰間距：{result.peak_separation_deg:.3f}°")
        if result.stacking_fault_alpha is not None:
            lines.append(f"  層錯機率 α ≈ {result.stacking_fault_alpha:.2f}%")
        lines.append(f"  層錯狀態：{result.stacking_fault_severity or 'N/A'}")
    
    # Lattice constant
    if result.lattice_constant is not None:
        lines.append(f"  晶格常數 a = {result.lattice_constant:.4f} Å")
        lines.append(f"  晶格狀態：{result.lattice_status or 'N/A'}")
    
    # Self-annealing
    lines.append(f"  自退火狀態：{result.annealing_state}")
    
    # Warnings
    if result.warnings:
        lines.extend([
            "",
            "-" * 40,
            "⚠️ 警告",
            "-" * 40,
        ])
        for w in result.warnings:
            lines.append(f"  • {w}")
    
    # Recommendations
    if result.recommendations:
        lines.extend([
            "",
            "-" * 40,
            "📋 建議",
            "-" * 40,
        ])
        for r in result.recommendations:
            lines.append(f"  • {r}")
    
    lines.extend([
        "",
        "=" * 60,
        "xrd_analysis - Advanced XRD Crystallite Size Analysis System",
        "=" * 60,
    ])
    
    return "\n".join(lines)


def generate_csv_summary(result: ComprehensiveResult) -> str:
    """
    Generate CSV-format summary for data export.
    
    Returns:
        CSV string with header and values
    """
    headers = [
        "Sample",
        "SampleAge_h",
        "Scherrer_nm",
        "WH_Size_nm",
        "WH_Strain",
        "WH_R2",
        "Dominant_hkl",
        "TC_dominant",
        "Peak_Sep_deg",
        "SF_alpha_pct",
        "Lattice_A",
        "Anneal_State",
    ]
    
    dom_hkl = ""
    if result.dominant_orientation:
        h, k, l = result.dominant_orientation
        dom_hkl = f"({h}{k}{l})"
    
    values = [
        result.sample_name,
        f"{result.sample_age_hours:.1f}" if result.sample_age_hours else "",
        f"{result.scherrer_size_nm:.1f}" if result.scherrer_size_nm else "",
        f"{result.wh_size_nm:.1f}" if result.wh_size_nm else "",
        f"{result.wh_strain:.2e}" if result.wh_strain else "",
        f"{result.wh_r_squared:.3f}" if result.wh_r_squared else "",
        dom_hkl,
        f"{result.dominant_tc:.2f}" if result.dominant_tc else "",
        f"{result.peak_separation_deg:.3f}" if result.peak_separation_deg else "",
        f"{result.stacking_fault_alpha:.2f}" if result.stacking_fault_alpha else "",
        f"{result.lattice_constant:.4f}" if result.lattice_constant else "",
        result.annealing_state,
    ]
    
    return ",".join(headers) + "\n" + ",".join(values)


def generate_process_recommendations(result: ComprehensiveResult) -> List[str]:
    """
    Generate process optimization recommendations based on analysis.
    根據分析結果產生製程優化建議。
    """
    recommendations = []
    
    # Size recommendations
    if result.scherrer_size_nm and result.scherrer_size_nm < 30:
        recommendations.append(
            "晶粒尺寸 < 30 nm：可能導致高電阻率，建議降低電流密度或增加溫度"
        )
    
    # Texture recommendations
    if result.dominant_orientation == (2, 2, 0) and result.dominant_tc and result.dominant_tc > 1.5:
        recommendations.append(
            "強(220)織構：可能表示高應力沉積，建議降低電流密度"
        )
    
    if result.dominant_orientation == (1, 1, 1) and result.dominant_tc and result.dominant_tc > 1.2:
        recommendations.append(
            "(111)擇優取向：有利於電遷移抗性，製程狀態良好"
        )
    
    # Stacking fault recommendations
    if result.peak_separation_deg and result.peak_separation_deg < 7.0:
        recommendations.append(
            "峰間距 < 7.0°：SPS 濃度可能過高，建議降低加速劑濃度"
        )
    
    # Lattice constant recommendations
    if result.lattice_constant and result.lattice_constant > 3.618:
        recommendations.append(
            "晶格常數 > 3.618 Å：嚴重雜質固溶，需檢查添加劑純度"
        )
    
    # Self-annealing recommendations
    if result.annealing_state == "as-deposited":
        recommendations.append(
            "鍍態樣品：建議 7 天後重測以獲得穩定結構數據"
        )
    
    return recommendations
