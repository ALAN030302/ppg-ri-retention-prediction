#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
mzML peak extraction and compound matching program
:
1. Extract peak data from mzML files
2. Load compound m/z information from Excel tables
3. Match peak data by compound m/z
4. Save matching results to a new Excel table
Author: Qianlei Yao
Version: 6.0
Date: 2025
"""
import os
import sys
import threading
import traceback
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Any
import warnings

# Suppress warnings
warnings.filterwarnings('ignore')

try:
    import pandas as pd
    import numpy as np
except ImportError:
    print("Error: please install pandas and numpy first")
    print("Installation command: pip install pandas numpy")
    sys.exit(1)

# Try to import GUI libraries
try:
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox, scrolledtext
    from tkinter.font import Font

    GUI_AVAILABLE = True
except ImportError:
    GUI_AVAILABLE = False
    print("Warning: tkinter is not installed, so the GUI is unavailable")


class MzMLProcessor:
    """mzML file processor - core processing logic"""

    def __init__(self, progress_callback=None, log_callback=None):
        """Initialize the processor"""
        self.results_summary = []
        self.temp_files = []
        self.progress_callback = progress_callback
        self.log_callback = log_callback
        self.cancel_requested = False

    def log(self, message: str, level: str = "INFO"):
        """Record log messages"""
        timestamp = datetime.now().strftime("%H:%M:%S")
        log_message = f"[{timestamp}] {message}"

        if self.log_callback:
            self.log_callback(log_message, level)
        else:
            print(log_message)

    def update_progress(self, value: int, message: str = ""):
        """Update progress"""
        if self.progress_callback:
            self.progress_callback(value, message)

    def check_mzml_file(self, file_path: Path) -> Tuple[bool, str]:
        """Check whether the mzML file is valid"""
        try:
            if not file_path.exists():
                return False, f"file does not exist: {file_path}"

            file_size = file_path.stat().st_size
            if file_size == 0:
                return False, "file is empty"

            # Check file format
            with open(file_path, 'rb') as f:
                header = f.read(500)
                header_str = header.decode('utf-8', errors='ignore')

                if '<mzML' in header_str or '<indexedmzML' in header_str:
                    return True, f"file is valid ({file_size / 1024 / 1024:.2f} MB)"
                else:
                    return False, "not a valid mzML file"

        except Exception as e:
            return False, f"error while checking the file: {str(e)}"

    def extract_peaks_with_pyopenms(self, mzml_path: str, intensity_threshold: float = 1000) -> Tuple[
        Optional[List[Dict]], str]:
        """Use pyOpenMS for peak extraction"""
        try:
            from pyopenms import MSExperiment, MzMLFile

            self.log(f"Load file with pyOpenMS: {Path(mzml_path).name}")
            exp = MSExperiment()
            MzMLFile().load(mzml_path, exp)

            if exp.size() == 0:
                return None, "loaded 0 spectra"

            peaks = []
            total_spectra = exp.size()

            for i, spectrum in enumerate(exp):
                if self.cancel_requested:
                    return None, "user cancelled the operation"

                if spectrum.size() > 0:
                    mz_array, intensity_array = spectrum.get_peaks()
                    rt = spectrum.getRT()

                    for j in range(len(mz_array)):
                        if intensity_array[j] > intensity_threshold:
                            peaks.append({
                                'spectrum_index': i,
                                'retention_time(RT)': rt,
                                '(m/z)': float(mz_array[j]),
                                '': float(intensity_array[j]),
                                'MS': spectrum.getMSLevel(),
                                'ID': f"F{i:05d}_{j:04d}"
                            })

                # Update progress
                if i % 10 == 0 or i == total_spectra - 1:
                    progress = int((i + 1) / total_spectra * 50) # 50%load
                    self.update_progress(progress, f" {i + 1}/{total_spectra}")

            return peaks, f"pyOpenMS: {len(peaks)}peak, {exp.size()}"

        except ImportError:
            return None, "pyOpenMS"
        except Exception as e:
            return None, f"pyOpenMSfailed: {str(e)}"

    def extract_peaks_with_pymzml(self, mzml_path: str, intensity_threshold: float = 1000) -> Tuple[
        Optional[List[Dict]], str]:
        """UsepymzMLpeak"""
        try:
            import pymzml

            self.log(f"UsepymzMLloadfile: {Path(mzml_path).name}")
            peaks = []

            run = pymzml.run.Reader(mzml_path)
            spectrum_count = 0
            total_spectra = None

            try:
                # translated note
                total_spectra = len(run)
            except:
                pass

            for spectrum in run:
                if self.cancel_requested:
                    return None, "user cancelled the operation"

                if spectrum.ms_level == 1:
                    mz_array = spectrum.mz
                    intensity_array = spectrum.i

                    if len(mz_array) > 0:
                        for idx, (mz, intensity) in enumerate(zip(mz_array, intensity_array)):
                            if intensity > intensity_threshold:
                                peaks.append({
                                    'spectrum_index': spectrum_count,
                                    'retention_time(RT)': spectrum.scan_time[0] if hasattr(spectrum,
                                                                                     'scan_time') and spectrum.scan_time else spectrum_count * 0.5,
                                    '(m/z)': float(mz),
                                    '': float(intensity),
                                    'MS': 1,
                                    'ID': f"F{spectrum_count:05d}_{idx:04d}"
                                })

                        spectrum_count += 1

                # Update progress
                if spectrum_count % 10 == 0:
                    progress_msg = f" {spectrum_count}"
                    if total_spectra:
                        progress = int(spectrum_count / total_spectra * 50)
                        progress_msg += f"/{total_spectra}"
                    else:
                        progress = min(50, spectrum_count // 2)

                    self.update_progress(progress, progress_msg)

            if peaks:
                return peaks, f"pymzML: {len(peaks)}peak"
            else:
                return None, "pymzMLconditionpeak"

        except ImportError:
            return None, "pymzML"
        except Exception as e:
            return None, f"pymzMLfailed: {str(e)}"

    def save_results(self, peaks: List[Dict], output_path: Path, base_name: str, method_used: str) -> Dict[str, str]:
        """saveresultsCSVExcel"""
        try:
            self.log(f"saveresults...")
            self.update_progress(60, "saveCSVfile...")

            # DataFrame
            df = pd.DataFrame(peaks)

            # Save as CSV
            csv_file = output_path / f"{base_name}_peakdata.csv"
            df.to_csv(csv_file, index=False, encoding='utf-8-sig')
            self.log(f"CSVsave: {csv_file.name}")

            result_files = {'csv': str(csv_file)}

            # saveExcel (XLSX)
            self.update_progress(70, "saveExcelfile...")
            try:
                import openpyxl
                excel_file = output_path / f"{base_name}_peakdata.xlsx"
                df.to_excel(excel_file, index=False, sheet_name='peakdata')
                result_files['excel'] = str(excel_file)
                self.log(f"Excel(XLSX)save: {excel_file.name}")
            except ImportError:
                self.log("openpyxl, Excel(XLSX)save", "WARNING")
            except Exception as e:
                self.log(f"saveExcel(XLSX)failed: {str(e)}", "WARNING")

            # saveExcel (XLS) - data
            self.update_progress(80, "save...")
            if len(df) <= 10000: # XLS65536
                try:
                    xls_file = output_path / f"{base_name}_peakdata.xls"
                    df.to_excel(xls_file, index=False, sheet_name='peakdata')
                    result_files['xls'] = str(xls_file)
                    self.log(f"Excel(XLS)save: {xls_file.name}")
                except Exception as e:
                    self.log(f"saveExcel(XLS)failed: {str(e)}", "WARNING")

            # save
            self.update_progress(90, "...")
            stats_file = output_path / f"{base_name}_.txt"
            self.save_statistics(df, stats_file, method_used)
            result_files['stats'] = str(stats_file)

            # saveExcel
            try:
                import openpyxl
                stats_excel = output_path / f"{base_name}_.xlsx"
                self.save_statistics_excel(df, stats_excel, method_used)
                result_files['stats_excel'] = str(stats_excel)
            except:
                pass

            self.update_progress(100, "save!")
            return result_files

        except Exception as e:
            self.log(f"saveresultsfailed: {str(e)}", "ERROR")
            return {}

    def save_statistics(self, df: pd.DataFrame, stats_file: Path, method_used: str):
        """savefile"""
        try:
            with open(stats_file, 'w', encoding='utf-8') as f:
                f.write("=" * 70 + "\n")
                f.write("peak\n")
                f.write("=" * 70 + "\n\n")

                f.write(f"filename: {stats_file.stem}\n")
                f.write(f": {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"method: {method_used}\n")
                f.write(f"peak: {len(df):,}\n\n")

                if len(df) > 0:
                    f.write(":\n")
                    f.write("-" * 50 + "\n")

                    # m/z
                    f.write(f"m/z:\n")
                    f.write(f" : {df['(m/z)'].min():.6f}\n")
                    f.write(f" : {df['(m/z)'].max():.6f}\n")
                    f.write(f" : {df['(m/z)'].mean():.6f}\n")
                    f.write(f" : {df['(m/z)'].median():.6f}\n")
                    f.write(f" : {df['(m/z)'].std():.6f}\n\n")

                    # RT
                    f.write(f"RT:\n")
                    f.write(f" : {df['retention_time(RT)'].min():.2f}\n")
                    f.write(f" : {df['retention_time(RT)'].max():.2f}\n")
                    f.write(f" : {df['retention_time(RT)'].mean():.2f}\n")
                    f.write(f" : {df['retention_time(RT)'].median():.2f}\n\n")

                    # translated note
                    f.write(f":\n")
                    f.write(f" : {df[''].min():.2e}\n")
                    f.write(f" : {df[''].max():.2e}\n")
                    f.write(f" : {df[''].mean():.2e}\n")
                    f.write(f" : {df[''].median():.2e}\n")
                    f.write(f" : {df[''].sum():.2e}\n\n")

                    # peak
                    f.write(f"10peak:\n")
                    f.write("-" * 50 + "\n")
                    top_peaks = df.nlargest(10, '')
                    for i, (_, row) in enumerate(top_peaks.iterrows()):
                        f.write(
                            f"{i + 1:2d}. m/z={row['(m/z)']:9.6f}, RT={row['retention_time(RT)']:7.2f}, ={row['']:12.2e}\n")

            self.log(f"save: {stats_file.name}")

        except Exception as e:
            self.log(f"savefailed: {str(e)}", "ERROR")

    def save_statistics_excel(self, df: pd.DataFrame, excel_file: Path, method_used: str):
        """saveExcel"""
        try:
            import openpyxl

            with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
                # translated note
                basic_stats = pd.DataFrame({
                    '': ['peak', 'm/z', 'm/z', 'm/z', 'm/z',
                               'RT()', 'RT()', 'RT()', 'RT()',
                               '', '', '', ''],
                    '': [len(df),
                             df['(m/z)'].min() if len(df) > 0 else 0,
                             df['(m/z)'].max() if len(df) > 0 else 0,
                             df['(m/z)'].mean() if len(df) > 0 else 0,
                             df['(m/z)'].median() if len(df) > 0 else 0,
                             df['retention_time(RT)'].min() if len(df) > 0 else 0,
                             df['retention_time(RT)'].max() if len(df) > 0 else 0,
                             df['retention_time(RT)'].mean() if len(df) > 0 else 0,
                             df['retention_time(RT)'].median() if len(df) > 0 else 0,
                             df[''].min() if len(df) > 0 else 0,
                             df[''].max() if len(df) > 0 else 0,
                             df[''].mean() if len(df) > 0 else 0,
                             df[''].median() if len(df) > 0 else 0]
                })
                basic_stats.to_excel(writer, sheet_name='', index=False)

                # translated note
                if len(df) > 0:
                    percentiles = [1, 5, 10, 25, 50, 75, 90, 95, 99]
                    intensity_percentiles = np.percentile(df[''], percentiles)

                    percentile_df = pd.DataFrame({
                        '(%)': percentiles,
                        '': intensity_percentiles
                    })
                    percentile_df.to_excel(writer, sheet_name='', index=False)

                # peak
                if len(df) > 0:
                    top_50 = df.nlargest(50, '').copy()
                    top_50.insert(0, 'rank', range(1, len(top_50) + 1))
                    top_50.to_excel(writer, sheet_name='peakTop50', index=False)

            self.log(f"Excelsave: {excel_file.name}")

        except Exception as e:
            self.log(f"saveExcelfailed: {str(e)}", "WARNING")

    def process_file(self, input_file: str, output_dir: str, method: str = 'auto',
                     intensity_threshold: float = 1000) -> Dict[str, Any]:
        """file"""
        try:
            input_path = Path(input_file)

            # file
            self.log(f"file: {input_path.name}")
            self.update_progress(5, "file...")

            is_valid, message = self.check_mzml_file(input_path)
            if not is_valid:
                self.log(f"filefailed: {message}", "ERROR")
                return {
                    'input_file': input_file,
                    'status': 'error',
                    'message': message,
                    'peak_count': 0
                }

            self.log(f"file: {message}")

            # outputdirectory
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)

            # peakdata
            self.log(f"peakdata...")
            peaks = None
            method_used = ""
            extract_message = ""

            # method
            if method == 'auto':
                # pyOpenMS
                self.update_progress(10, "pyOpenMS...")
                peaks, msg = self.extract_peaks_with_pyopenms(input_file, intensity_threshold)
                if peaks is not None:
                    method_used = "pyOpenMS"
                    extract_message = msg
                else:
                    # pymzML
                    self.log(f"pyOpenMSfailed: {msg}", "WARNING")
                    self.update_progress(10, "pymzML...")
                    peaks, msg = self.extract_peaks_with_pymzml(input_file, intensity_threshold)
                    if peaks is not None:
                        method_used = "pymzML"
                        extract_message = msg
                    else:
                        extract_message = f"methodfailed: {msg}"

            elif method == 'pyopenms':
                self.update_progress(10, "UsepyOpenMS...")
                peaks, msg = self.extract_peaks_with_pyopenms(input_file, intensity_threshold)
                method_used = "pyOpenMS"
                extract_message = msg

            elif method == 'pymzml':
                self.update_progress(10, "UsepymzML...")
                peaks, msg = self.extract_peaks_with_pymzml(input_file, intensity_threshold)
                method_used = "pymzML"
                extract_message = msg

            # translated note
            if self.cancel_requested:
                self.log("", "WARNING")
                return {
                    'input_file': input_file,
                    'status': 'cancelled',
                    'message': 'user cancelled the operation',
                    'peak_count': 0
                }

            # saveresults
            if peaks is not None and len(peaks) > 0:
                base_name = input_path.stem
                result_files = self.save_results(peaks, output_path, base_name, method_used)

                return {
                    'input_file': input_file,
                    'output_dir': str(output_path),
                    'peak_count': len(peaks),
                    'method': method_used,
                    'status': 'success',
                    'message': extract_message,
                    'output_files': result_files
                }
            else:
                self.log(f"Extractedpeakdata: {extract_message}", "WARNING")
                return {
                    'input_file': input_file,
                    'output_dir': str(output_path),
                    'peak_count': 0,
                    'method': method_used,
                    'status': 'failed',
                    'message': extract_message
                }

        except Exception as e:
            error_msg = f"file: {str(e)}"
            self.log(error_msg, "ERROR")
            traceback.print_exc()
            return {
                'input_file': input_file,
                'status': 'error',
                'message': error_msg,
                'peak_count': 0
            }

    def process_batch(self, input_dir: str, output_dir: str, file_pattern: str = "*.mzML",
                      method: str = 'auto', intensity_threshold: float = 1000) -> List[Dict[str, Any]]:
        """file"""
        try:
            input_path = Path(input_dir)

            # inputdirectory
            if not input_path.exists():
                error_msg = f"inputdirectory: {input_dir}"
                self.log(error_msg, "ERROR")
                raise FileNotFoundError(error_msg)

            # file
            self.log(f"directoryfile: {input_dir}")
            self.update_progress(5, "file...")

            mzml_files = list(input_path.glob(file_pattern))
            if not mzml_files:
                error_msg = f" {input_dir} {file_pattern} file"
                self.log(error_msg, "ERROR")
                raise FileNotFoundError(error_msg)

            self.log(f" {len(mzml_files)} file")

            # outputdirectory
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)

            self.results_summary = []

            # file
            for i, mzml_file in enumerate(mzml_files):
                if self.cancel_requested:
                    self.log("", "WARNING")
                    break

                file_num = i + 1
                total_files = len(mzml_files)

                self.log(f"file {file_num}/{total_files}: {mzml_file.name}")

                # fileoutputdirectory
                file_output_dir = output_path / mzml_file.stem
                file_output_dir.mkdir(exist_ok=True)

                # Update progress
                progress = int(file_num / total_files * 100)
                self.update_progress(progress, f"file {file_num}/{total_files}")

                # file
                result = self.process_file(str(mzml_file), str(file_output_dir), method, intensity_threshold)

                if result:
                    self.results_summary.append(result)

            # translated note
            if self.results_summary:
                self.generate_batch_report(output_path)

            return self.results_summary

        except Exception as e:
            error_msg = f": {str(e)}"
            self.log(error_msg, "ERROR")
            traceback.print_exc()
            return []

    def generate_batch_report(self, output_path: Path):
        """"""
        try:
            self.log("...")
            self.update_progress(95, "...")

            successful = [r for r in self.results_summary if r.get('status') == 'success']
            failed = [r for r in self.results_summary if r.get('status') == 'failed']
            errors = [r for r in self.results_summary if r.get('status') == 'error']

            # save
            summary_file = output_path / ".csv"
            summary_data = []

            for result in self.results_summary:
                summary_data.append({
                    'file': Path(result['input_file']).name,
                    '': result.get('status', 'unknown'),
                    'peak': result.get('peak_count', 0),
                    'method': result.get('method', ''),
                    'message': result.get('message', '')[:100]
                })

            df_summary = pd.DataFrame(summary_data)
            df_summary.to_csv(summary_file, index=False, encoding='utf-8-sig')
            self.log(f"save: {summary_file.name}")

            # saveExcel
            try:
                import openpyxl
                excel_summary = output_path / ".xlsx"
                df_summary.to_excel(excel_summary, index=False, sheet_name='')
                self.log(f"Excelsave: {excel_summary.name}")
            except:
                pass

            # translated note
            report_file = output_path / ".txt"
            with open(report_file, 'w', encoding='utf-8') as f:
                f.write("\n")
                f.write("=" * 70 + "\n\n")
                f.write(f": {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"inputdirectory: {output_path.parent}\n")
                f.write(f"outputdirectory: {output_path}\n\n")

                f.write(f":\n")
                f.write(f" file: {len(self.results_summary)}")
                f.write(f" : {len(successful)}")
                f.write(f" failed: {len(failed)}")
                f.write(f" : {len(errors)}\n\n")

                if successful:
                    total_peaks = sum(r.get('peak_count', 0) for r in successful)
                    f.write(f"filecolumn ( {total_peaks:,} peak):\n")
                    f.write("-" * 50 + "\n")
                    for result in successful:
                        f.write(f" {Path(result['input_file']).name}: {result.get('peak_count', 0)}peak\n")
                    f.write("\n")

                if failed:
                    f.write(f"failedfilecolumn:\n")
                    f.write("-" * 50 + "\n")
                    for result in failed:
                        f.write(f" {Path(result['input_file']).name}: {result.get('message', '')}\n")
                    f.write("\n")

            self.log(f"save: {report_file.name}")
            self.update_progress(100, "!")

        except Exception as e:
            self.log(f"failed: {str(e)}", "ERROR")

    def cancel(self):
        """"""
        self.cancel_requested = True
        self.log("...", "WARNING")


class MzMLGUI:
    """mzMLfileGUI"""

    def __init__(self, root):
        """GUI"""
        self.root = root
        self.root.title("mzMLfilepeak")
        self.root.geometry("900x700")

        # translated note
        try:
            self.root.iconbitmap(default='icon.ico')
        except:
            pass

        # translated note
        self.processor = None
        self.processing_thread = None
        self.is_processing = False

        # UI
        self.setup_ui()

        # translated note
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)

    def setup_ui(self):
        """UI"""
        # translated note
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # translated note
        title_label = ttk.Label(main_frame, text="mzMLfilepeak",
                                font=("Arial", 16, "bold"))
        title_label.grid(row=0, column=0, columnspan=3, pady=10)

        # Version
        version_label = ttk.Label(main_frame, text="Version 5.1 - CSVExceloutput",
                                  font=("Arial", 10))
        version_label.grid(row=1, column=0, columnspan=3, pady=(0, 10))

        # translated note
        ttk.Separator(main_frame, orient=tk.HORIZONTAL).grid(
            row=2, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=10
        )

        # translated note
        mode_frame = ttk.LabelFrame(main_frame, text="", padding="10")
        mode_frame.grid(row=3, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5)

        self.mode_var = tk.StringVar(value="single")

        ttk.Radiobutton(mode_frame, text="file", variable=self.mode_var,
                        value="single", command=self.on_mode_change).grid(
            row=0, column=0, sticky=tk.W, padx=10)

        ttk.Radiobutton(mode_frame, text="", variable=self.mode_var,
                        value="batch", command=self.on_mode_change).grid(
            row=0, column=1, sticky=tk.W, padx=10)

        # file
        file_frame = ttk.LabelFrame(main_frame, text="file", padding="10")
        file_frame.grid(row=4, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=10)

        # file
        self.single_file_frame = ttk.Frame(file_frame)
        self.single_file_frame.grid(row=0, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5)

        ttk.Label(self.single_file_frame, text="inputfile:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5))

        self.input_file_var = tk.StringVar()
        input_file_entry = ttk.Entry(self.single_file_frame, textvariable=self.input_file_var, width=60)
        input_file_entry.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=(0, 5))

        ttk.Button(self.single_file_frame, text="...", command=self.browse_input_file).grid(
            row=0, column=2, padx=(0, 5))

        ttk.Button(self.single_file_frame, text="", command=self.clear_input_file).grid(
            row=0, column=3)

        # translated note
        self.batch_file_frame = ttk.Frame(file_frame)
        self.batch_file_frame.grid(row=1, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5)

        ttk.Label(self.batch_file_frame, text="inputdirectory:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5))

        self.input_dir_var = tk.StringVar()
        input_dir_entry = ttk.Entry(self.batch_file_frame, textvariable=self.input_dir_var, width=60)
        input_dir_entry.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=(0, 5))

        ttk.Button(self.batch_file_frame, text="...", command=self.browse_input_dir).grid(
            row=0, column=2, padx=(0, 5))

        ttk.Label(self.batch_file_frame, text="file:").grid(row=1, column=0, sticky=tk.W, padx=(0, 5), pady=(5, 0))

        self.file_pattern_var = tk.StringVar(value="*.mzML")
        file_pattern_combo = ttk.Combobox(self.batch_file_frame, textvariable=self.file_pattern_var,
                                          values=["*.mzML", "*.mzML.gz", "*.mzXML"], width=20)
        file_pattern_combo.grid(row=1, column=1, sticky=tk.W, padx=(0, 5), pady=(5, 0))

        # translated note
        self.batch_file_frame.grid_remove()

        # outputdirectory
        ttk.Label(file_frame, text="outputdirectory:").grid(row=2, column=0, sticky=tk.W, padx=(0, 5), pady=(10, 0))

        self.output_dir_var = tk.StringVar()
        output_dir_entry = ttk.Entry(file_frame, textvariable=self.output_dir_var, width=60)
        output_dir_entry.grid(row=2, column=1, sticky=(tk.W, tk.E), padx=(0, 5), pady=(10, 0))

        ttk.Button(file_frame, text="...", command=self.browse_output_dir).grid(
            row=2, column=2, pady=(10, 0))

        # (, retention)
        ttk.Label(file_frame, text=": Usefiledirectory",
                  font=("Arial", 9, "italic")).grid(
            row=3, column=0, columnspan=3, sticky=tk.W, pady=(10, 0))

        # Parameters
        param_frame = ttk.LabelFrame(main_frame, text="Parameters", padding="10")
        param_frame.grid(row=5, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=10)

        ttk.Label(param_frame, text="method:").grid(row=0, column=0, sticky=tk.W, padx=(0, 10))

        self.method_var = tk.StringVar(value="auto")
        method_combo = ttk.Combobox(param_frame, textvariable=self.method_var,
                                    values=["auto", "pyopenms", "pymzml"],
                                    state="readonly", width=15)
        method_combo.grid(row=0, column=1, sticky=tk.W)

        ttk.Label(param_frame, text=":").grid(row=0, column=2, sticky=tk.W, padx=(20, 10))

        self.threshold_var = tk.StringVar(value="1000")
        threshold_entry = ttk.Entry(param_frame, textvariable=self.threshold_var, width=15)
        threshold_entry.grid(row=0, column=3, sticky=tk.W)

        # translated note
        progress_frame = ttk.LabelFrame(main_frame, text="", padding="10")
        progress_frame.grid(row=6, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=10)

        self.progress_var = tk.StringVar(value="...")
        ttk.Label(progress_frame, textvariable=self.progress_var).grid(
            row=0, column=0, sticky=tk.W, pady=(0, 5))

        self.progress_bar = ttk.Progressbar(progress_frame, mode='determinate', length=800)
        self.progress_bar.grid(row=1, column=0, sticky=(tk.W, tk.E), pady=(0, 5))

        self.progress_percent = tk.StringVar(value="0%")
        ttk.Label(progress_frame, textvariable=self.progress_percent).grid(
            row=1, column=1, sticky=tk.E, padx=(5, 0), pady=(0, 5))

        # translated note
        log_frame = ttk.LabelFrame(main_frame, text="", padding="10")
        log_frame.grid(row=7, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), pady=10)

        # translated note
        self.log_text = scrolledtext.ScrolledText(log_frame, width=100, height=15,
                                                  wrap=tk.WORD, font=("Consolas", 10))
        self.log_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # translated note
        self.log_text.tag_config("INFO", foreground="black")
        self.log_text.tag_config("SUCCESS", foreground="green")
        self.log_text.tag_config("WARNING", foreground="orange")
        self.log_text.tag_config("ERROR", foreground="red")

        # translated note
        button_frame = ttk.Frame(main_frame)
        button_frame.grid(row=8, column=0, columnspan=3, pady=10)

        self.process_button = ttk.Button(button_frame, text="",
                                         command=self.start_processing, width=15)
        self.process_button.pack(side=tk.LEFT, padx=5)

        self.cancel_button = ttk.Button(button_frame, text="",
                                        command=self.cancel_processing, width=15,
                                        state=tk.DISABLED)
        self.cancel_button.pack(side=tk.LEFT, padx=5)

        ttk.Button(button_frame, text="", command=self.clear_log).pack(side=tk.LEFT, padx=5)

        ttk.Button(button_frame, text="outputdirectory", command=self.open_output_dir).pack(side=tk.LEFT, padx=5)

        ttk.Button(button_frame, text="", command=self.on_closing).pack(side=tk.LEFT, padx=5)

        # translated note
        self.status_var = tk.StringVar(value="")
        status_bar = ttk.Label(main_frame, textvariable=self.status_var, relief=tk.SUNKEN)
        status_bar.grid(row=9, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(10, 0))

        # translated note
        main_frame.columnconfigure(1, weight=1)
        log_frame.columnconfigure(0, weight=1)
        log_frame.rowconfigure(0, weight=1)

        # translated note
        self.on_mode_change()

    def on_mode_change(self):
        """"""
        mode = self.mode_var.get()

        if mode == "single":
            self.single_file_frame.grid()
            self.batch_file_frame.grid_remove()
        else:
            self.single_file_frame.grid_remove()
            self.batch_file_frame.grid()

    def browse_input_file(self):
        """inputfile"""
        file_types = [("mzMLfile", "*.mzML *.mzML.gz"), ("file", "*.*")]
        file_path = filedialog.askopenfilename(title="mzMLfile", filetypes=file_types)

        if file_path:
            self.input_file_var.set(file_path)

            # outputdirectory
            if not self.output_dir_var.get():
                input_path = Path(file_path)
                output_dir = input_path.parent / "output"
                self.output_dir_var.set(str(output_dir))

    def browse_input_dir(self):
        """inputdirectory"""
        dir_path = filedialog.askdirectory(title="inputdirectory")

        if dir_path:
            self.input_dir_var.set(dir_path)

            # directorymzMLfile
            file_pattern = self.file_pattern_var.get()
            input_path = Path(dir_path)
            mzml_files = list(input_path.glob(file_pattern))

            if not mzml_files:
                messagebox.showwarning("Warning", f"directory {file_pattern} file")
            else:
                messagebox.showinfo("", f" {len(mzml_files)} file")

            # outputdirectory
            if not self.output_dir_var.get():
                output_dir = input_path.parent / "batch_output"
                self.output_dir_var.set(str(output_dir))

    def browse_output_dir(self):
        """outputdirectory"""
        dir_path = filedialog.askdirectory(title="outputdirectory")

        if dir_path:
            self.output_dir_var.set(dir_path)

    def clear_input_file(self):
        """inputfile"""
        self.input_file_var.set("")

    def log_message(self, message: str, level: str = "INFO"):
        """message"""
        self.log_text.insert(tk.END, message + "\n", level)
        self.log_text.see(tk.END)
        self.root.update()

    def update_progress(self, value: int, message: str = ""):
        """Update progress"""
        self.progress_bar['value'] = value
        if message:
            self.progress_var.set(message)
        self.progress_percent.set(f"{value}%")
        self.root.update()

    def update_status(self, message: str):
        """"""
        self.status_var.set(message)
        self.root.update()

    def clear_log(self):
        """"""
        self.log_text.delete(1.0, tk.END)

    def validate_inputs(self) -> Tuple[bool, str]:
        """inputParameters"""
        mode = self.mode_var.get()

        if mode == "single":
            input_path = self.input_file_var.get().strip()
            if not input_path:
                return False, "inputfile"

            if not os.path.exists(input_path):
                return False, f"inputfile does not exist: {input_path}"

        else: # batch
            input_dir = self.input_dir_var.get().strip()
            if not input_dir:
                return False, "inputdirectory"

            if not os.path.exists(input_dir):
                return False, f"inputdirectory: {input_dir}"

            # directorymzMLfile
            file_pattern = self.file_pattern_var.get()
            input_path = Path(input_dir)
            mzml_files = list(input_path.glob(file_pattern))
            if not mzml_files:
                return False, f"directory {file_pattern} file"

        # outputdirectory
        output_dir = self.output_dir_var.get().strip()
        if not output_dir:
            return False, "outputdirectory"

        # translated note
        try:
            threshold = float(self.threshold_var.get())
            if threshold <= 0:
                return False, "0"
        except ValueError:
            return False, ""

        return True, ""

    def start_processing(self):
        """"""
        # input
        is_valid, message = self.validate_inputs()
        if not is_valid:
            messagebox.showerror("input", message)
            return

        # ,
        self.process_button.config(state=tk.DISABLED)
        self.cancel_button.config(state=tk.NORMAL)

        # translated note
        self.progress_bar['value'] = 0
        self.progress_var.set("...")
        self.progress_percent.set("0%")
        self.clear_log()

        # Parameters
        mode = self.mode_var.get()
        output_dir = self.output_dir_var.get()
        method = self.method_var.get()

        try:
            threshold = float(self.threshold_var.get())
        except ValueError:
            threshold = 1000

        # translated note
        self.processor = MzMLProcessor(
            progress_callback=self.update_progress,
            log_callback=self.log_message
        )

        # translated note
        if mode == "single":
            input_file = self.input_file_var.get()
            self.processing_thread = threading.Thread(
                target=self.process_single_file,
                args=(input_file, output_dir, method, threshold)
            )
        else:
            input_dir = self.input_dir_var.get()
            file_pattern = self.file_pattern_var.get()
            self.processing_thread = threading.Thread(
                target=self.process_batch,
                args=(input_dir, output_dir, file_pattern, method, threshold)
            )

        self.is_processing = True
        self.processing_thread.start()

        # translated note
        self.root.after(100, self.check_processing_status)

    def process_single_file(self, input_file: str, output_dir: str, method: str, threshold: float):
        """file ()"""
        try:
            result = self.processor.process_file(input_file, output_dir, method, threshold)

            # results
            self.root.after(0, self.show_processing_result, result)

        except Exception as e:
            error_msg = f": {str(e)}"
            self.root.after(0, self.log_message, error_msg, "ERROR")
            self.root.after(0, self.processing_finished)

    def process_batch(self, input_dir: str, output_dir: str, file_pattern: str, method: str, threshold: float):
        """ ()"""
        try:
            results = self.processor.process_batch(input_dir, output_dir, file_pattern, method, threshold)

            # results
            self.root.after(0, self.show_batch_result, results)

        except Exception as e:
            error_msg = f": {str(e)}"
            self.root.after(0, self.log_message, error_msg, "ERROR")
            self.root.after(0, self.processing_finished)

    def show_processing_result(self, result: Dict[str, Any]):
        """fileresults"""
        if result['status'] == 'success':
            self.log_message(f"✓ !", "SUCCESS")
            self.log_message(f" Extracted {result['peak_count']:,} peak")
            self.log_message(f" outputdirectory: {result['output_dir']}")

            if 'output_files' in result:
                self.log_message(f" file:")
                for file_type, file_path in result['output_files'].items():
                    if file_path:
                        file_name = Path(file_path).name
                        self.log_message(f"    • {file_type}: {file_name}")

            # translated note
            self.root.after(0, lambda: messagebox.showinfo(
                "",
                f" {result['peak_count']:,} peak\noutputdirectory: {result['output_dir']}"
            ))

        elif result['status'] == 'cancelled':
            self.log_message("", "WARNING")

        else:
            self.log_message(f"✗ failed: {result.get('message', '')}", "ERROR")

        self.processing_finished()

    def show_batch_result(self, results: List[Dict[str, Any]]):
        """results"""
        if not results:
            self.log_message("results", "WARNING")
            self.processing_finished()
            return

        successful = [r for r in results if r.get('status') == 'success']
        failed = [r for r in results if r.get('status') == 'failed']
        errors = [r for r in results if r.get('status') == 'error']

        self.log_message(f"!", "SUCCESS")
        self.log_message(f" file: {len(results)}")
        self.log_message(f" : {len(successful)}")
        self.log_message(f" failed: {len(failed)}")
        self.log_message(f" : {len(errors)}")

        if successful:
            total_peaks = sum(r.get('peak_count', 0) for r in successful)
            self.log_message(f" peak: {total_peaks:,}")

        # results
        result_text = f"!\n\n"
        result_text += f"file: {len(results)}\n"
        result_text += f": {len(successful)}\n"
        result_text += f"failed: {len(failed)}\n"

        if successful:
            total_peaks = sum(r.get('peak_count', 0) for r in successful)
            result_text += f"peak: {total_peaks:,}\n"

        result_text += f"\nresultsoutputdirectoryfile. "

        self.root.after(0, lambda: messagebox.showinfo("", result_text))

        self.processing_finished()

    def check_processing_status(self):
        """"""
        if self.is_processing and self.processing_thread.is_alive():
            # translated note
            self.root.after(100, self.check_processing_status)
        elif self.is_processing:
            # , results
            self.root.after(100, self.check_processing_status)

    def processing_finished(self):
        """"""
        self.is_processing = False
        self.processor = None

        # ,
        self.process_button.config(state=tk.NORMAL)
        self.cancel_button.config(state=tk.DISABLED)

        # translated note
        self.update_status("")

    def cancel_processing(self):
        """"""
        if self.processor and self.is_processing:
            self.processor.cancel()
            self.log_message("...", "WARNING")
            self.update_status("...")

    def open_output_dir(self):
        """outputdirectory"""
        output_dir = self.output_dir_var.get().strip()

        if output_dir and os.path.exists(output_dir):
            try:
                if sys.platform == 'win32':
                    os.startfile(output_dir)
                elif sys.platform == 'darwin':
                    os.system(f'open "{output_dir}"')
                else:
                    os.system(f'xdg-open "{output_dir}"')
            except Exception as e:
                messagebox.showerror("", f"directory: {str(e)}")
        else:
            messagebox.showwarning("Warning", "outputdirectory")

    def on_closing(self):
        """"""
        if self.is_processing:
            if messagebox.askyesno("", ", ？"):
                if self.processor:
                    self.processor.cancel()
                self.root.destroy()
        else:
            self.root.destroy()


def check_dependencies():
    """"""
    print("=" * 70)
    print("...")
    print("=" * 70)

    dependencies = {
        'pandas': 'data',
        'numpy': 'calculate',
        'openpyxl': 'Excel XLSX',
        'pyopenms': ' (method1)',
        'pymzml': ' (method2)'
    }

    missing = []

    for lib, desc in dependencies.items():
        try:
            __import__(lib)
            print(f"✓ {lib}: {desc}")
        except ImportError:
            print(f"✗ {lib}: {desc} - ")
            missing.append(lib)

    if missing:
        print(f"\n:")
        for lib in missing:
            print(f"  pip install {lib}")

    print("\n" + "=" * 70)
    return len(missing) == 0


def main():
    """"""

    # GUI
    if not GUI_AVAILABLE:
        print(": tkinter, GUI")
        print("tkinter:")
        print("  Ubuntu/Debian: sudo apt-get install python3-tk")
        print(" Windows/macOS: ")
        return

    # translated note
    if not check_dependencies():
        response = input("\n, ? (y/n): ")
        if response.lower() != 'y':
            return

    # translated note
    root = tk.Tk()

    # translated note
    window_width = 900
    window_height = 700
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()
    x = (screen_width - window_width) // 2
    y = (screen_height - window_height) // 2
    root.geometry(f"{window_width}x{window_height}+{x}+{y}")

    # GUI
    app = MzMLGUI(root)

    # translated note
    root.mainloop()


if __name__ == "__main__":
    print("mzMLfilepeak - GUIVersion ()")
    print("Version: 5.1")
    print(": , , CSVExceloutput")
    print(": ")
    print("=" * 70)

    main()
