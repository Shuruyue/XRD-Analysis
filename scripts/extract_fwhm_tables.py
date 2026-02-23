
import pandas as pd
from pathlib import Path
import numpy as np
from xrd_analysis.analysis.pipeline import load_bruker_txt, parse_filename
from xrd_analysis.fitting.peak_fitter import fit_peak_with_diagnosis
from xrd_analysis.core.copper_crystal import get_standard_peaks

def extract_fwhm_tables():
    # Setup paths
    project_root = Path(__file__).resolve().parent.parent
    data_dir = project_root / "data/raw/202511"
    output_dir = project_root / "outputs/tables/fwhm_data"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Get all sample files
    data_files = sorted(data_dir.glob("*.txt"))
    print(f"Found {len(data_files)} samples.")
    
    # Define peaks of interest
    peaks_map = {
        "111": 43.316,
        "200": 50.448,
        "220": 74.124
    }
    
    # Initialize storage
    results_111 = []
    results_200 = []
    results_220 = []
    
    print("Processing samples and fitting peaks...")
    
    for i, filepath in enumerate(data_files, 1):
        filename = filepath.stem
        try:
            # 1. Load Data
            two_theta, intensity = load_bruker_txt(str(filepath))
            
            # 2. Parse filename info (Leveler, Time)
            info = parse_filename(filename)
            leveler = info.get('concentration_ml', 0.0)
            time = info.get('time_hours', 0.0)
            
            # 3. Fit each peak
            # (111)
            fit111 = fit_peak_with_diagnosis(two_theta, intensity, peaks_map["111"], use_doublet=True)
            
            # (200)
            fit200 = fit_peak_with_diagnosis(two_theta, intensity, peaks_map["200"], use_doublet=True)
            
            # (220)
            fit220 = fit_peak_with_diagnosis(two_theta, intensity, peaks_map["220"], use_doublet=True)
            
            # Helper to extract data safely
            def get_data(fit_res):
                if fit_res.get('success', False):
                    return {
                        "Ka1 (deg)": fit_res.get('center', np.nan),
                        "Ka2 (deg)": fit_res.get('center_ka2', np.nan),
                        "FWHM (deg)": fit_res.get('fwhm', np.nan),
                        "FWHM Error": fit_res.get('fwhm_err', np.nan)
                    }
                else:
                    return {
                        "Ka1 (deg)": np.nan, "Ka2 (deg)": np.nan,
                        "FWHM (deg)": np.nan, "FWHM Error": np.nan
                    }

            # Store row data
            row_base = {
                "Sample": filename, 
                "Leveler (mL)": leveler, 
                "Time (h)": time
            }
            
            # Append to lists
            row111 = row_base.copy()
            row111.update(get_data(fit111))
            results_111.append(row111)
            
            row200 = row_base.copy()
            row200.update(get_data(fit200))
            results_200.append(row200)
            
            row220 = row_base.copy()
            row220.update(get_data(fit220))
            results_220.append(row220)
            
            if i % 10 == 0:
                print(f"  Processed {i}/{len(data_files)}...")
                
        except Exception as e:
            print(f"Error processing {filename}: {e}")
            
    # Convert to DataFrames
    df_111 = pd.DataFrame(results_111)
    df_200 = pd.DataFrame(results_200)
    df_220 = pd.DataFrame(results_220)
    
    # Save to CSV
    path_111 = output_dir / "fwhm_111.csv"
    path_200 = output_dir / "fwhm_200.csv"
    path_220 = output_dir / "fwhm_220.csv"
    
    df_111.to_csv(path_111, index=False)
    df_200.to_csv(path_200, index=False)
    df_220.to_csv(path_220, index=False)
    
    print("\n" + "="*50)
    print(" extraction Complete")
    print("="*50)
    print(f"Saved tables to: {output_dir}")
    print("\n--- (111) Table Preview ---")
    print(df_111.head().to_markdown(index=False))
    print("\n--- (200) Table Preview ---")
    print(df_200.head().to_markdown(index=False))
    print("\n--- (220) Table Preview ---")
    print(df_220.head().to_markdown(index=False))

if __name__ == "__main__":
    extract_fwhm_tables()
