
import pandas as pd
from pathlib import Path

def generate_markdown_report():
    # Setup paths
    project_root = Path(__file__).resolve().parent.parent
    data_dir = project_root / "outputs/tables/fwhm_data"
    output_path = project_root / "docs/engineering_specs/06_Appendix_FWHM_Data.md"
    
    # Define files to process
    peaks = ["111", "200", "220"]
    approx_pos = {"111": "43.3", "200": "50.4", "220": "74.1"}
    
    md_content = [
        "# 附錄：FWHM 完整數據表 (Appendix: FWHM Complete Data Tables)\n",
        "**文件編號**: DOC-ENG-06",
        "**數據來源**: 原始 XRD 圖譜 (True Voigt Doublet Fitting)",
        "**描述**: 本文件收錄了所有 52 個樣品的擬合結果，包含 $K\\alpha_1$/$K\\alpha_2$ 峰位分離與誤差分析。\n",
        "---\n"
    ]
    
    for peak in peaks:
        csv_path = data_dir / f"fwhm_{peak}.csv"
        if not csv_path.exists():
            print(f"Warning: {csv_path} not found.")
            continue
            
        df = pd.read_csv(csv_path)
        
        df['Leveler'] = pd.to_numeric(df['Leveler (mL)'], errors='coerce')
        df['Time'] = pd.to_numeric(df['Time (h)'], errors='coerce')
        
        # DEBUG: Check types
        print(f"DEBUG: Peak {peak} dtypes:\n{df.dtypes}")
        
        df = df.sort_values(by=['Leveler', 'Time'], ascending=[True, True])
        
        # DEBUG: Print first 10 rows of sorted frame
        print(f"DEBUG: Peak {peak} sorted head:\n{df[['Sample', 'Leveler', 'Time']].head(10)}")
        
        # Drop helper sort columns if not needed for display, 
        # but requested format includes them.
        
        # Format columns
        # Filename, Ka1, Ka2, FWHM, Error
        
        md_content.append(f"## Table 6.{peaks.index(peak)+1}: ({peak}) Peak Data ($2\\theta \\approx {approx_pos[peak]}^\\circ$)\n")
        
        # Create display string
        # Sample Name | Ka1 | Ka2 | FWHM (Error)
        
        table_header = "| Sample Name | Kα1 (deg) | Kα2 (deg) | FWHM (deg) | Error (±) |\n| :--- | :--- | :--- | :--- | :--- |"
        md_content.append(table_header)
        
        for _, row in df.iterrows():
            # Format values
            ka1 = f"{row['Ka1 (deg)']:.3f}" if pd.notnull(row['Ka1 (deg)']) else "N/A"
            ka2 = f"{row['Ka2 (deg)']:.3f}" if pd.notnull(row['Ka2 (deg)']) else "N/A"
            fwhm = f"{row['FWHM (deg)']:.4f}" if pd.notnull(row['FWHM (deg)']) else "N/A"
            err = f"{row['FWHM Error']:.5f}" if pd.notnull(row['FWHM Error']) else "N/A"
            
            line = f"| {row['Sample']} | {ka1} | {ka2} | {fwhm} | {err} |"
            md_content.append(line)
            
        md_content.append("\n")
        md_content.append("### Raw Data for Copy-Paste (CSV)")
        md_content.append("```csv")
        md_content.append("Sample Name,Ka1 (deg),Ka2 (deg),FWHM (deg),Error (+/-)")
        
        for _, row in df.iterrows():
            ka1 = f"{row['Ka1 (deg)']:.3f}" if pd.notnull(row['Ka1 (deg)']) else ""
            ka2 = f"{row['Ka2 (deg)']:.3f}" if pd.notnull(row['Ka2 (deg)']) else ""
            fwhm = f"{row['FWHM (deg)']:.4f}" if pd.notnull(row['FWHM (deg)']) else ""
            err = f"{row['FWHM Error']:.5f}" if pd.notnull(row['FWHM Error']) else ""
            md_content.append(f"{row['Sample']},{ka1},{ka2},{fwhm},{err}")
            
        md_content.append("```")
        
        md_content.append("\n")
        md_content.append("### Simplified FWHM Data (For plotting)")
        md_content.append("```csv")
        md_content.append("Sample Name,FWHM (deg)")
        
        for _, row in df.iterrows():
            fwhm = f"{row['FWHM (deg)']:.4f}" if pd.notnull(row['FWHM (deg)']) else ""
            md_content.append(f"{row['Sample']},{fwhm}")
            
        md_content.append("```")
        md_content.append("\n---\n")
        
    # Write to file
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write("\n".join(md_content))
    
    print(f"Successfully generated markdown at: {output_path}")

if __name__ == "__main__":
    generate_markdown_report()
