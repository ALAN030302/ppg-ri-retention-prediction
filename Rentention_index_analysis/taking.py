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
Version: 7.0
Date: 2025
"""

import os
import sys
import threading
import traceback
import tempfile
import shutil
import json
import re
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Any, Union
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

    GUI_AVAILABLE = True
except ImportError:
    GUI_AVAILABLE = False
    print("Warning: tkinter is not installed, so the GUI is unavailable")


class CompoundMatcher:
    """compound matching - Excelloadcompoundmatch"""

    def __init__(self):
        """match"""
        self.compounds_df = None
        self.match_results = []
        self.match_settings = {
            'ppm_tolerance': 10, # ppm
            'rt_window': 30, # retention_time ()
            'intensity_threshold': 1000 #
        }

    def load_compounds_from_excel(self, excel_file: str) -> Tuple[bool, str, pd.DataFrame]:
        """Excelfileloadcompound"""
        try:
            file_path = Path(excel_file)

            if not file_path.exists():
                return False, "file does not exist", None

            # Try to read the Excel file
            try:
                # Try all possible worksheets
                xls = pd.ExcelFile(file_path)
                sheet_names = xls.sheet_names

                # Try the first worksheet
                df = pd.read_excel(file_path, sheet_name=sheet_names[0])

                # Find required columns
                mz_col = self._find_column(df, ['mz', 'M/Z', 'm/z', 'mass', 'Mass', 'MZ', ''])
                name_col = self._find_column(df, ['name', 'Name', 'compound', 'Compound', 'ID', 'name', 'compound_name'])
                formula_col = self._find_column(df, ['formula', 'Formula', 'molecular_formula'])
                rt_col = self._find_column(df, ['rt', 'RT', 'retention_time', 'retention_time'])

                if mz_col is None:
                    return False, "m/zcolumn, Excelfilem/z", None

                # column
                rename_dict = {}
                if mz_col:
                    rename_dict[mz_col] = 'mz'
                if name_col:
                    rename_dict[name_col] = 'name'
                if formula_col:
                    rename_dict[formula_col] = 'formula'
                if rt_col:
                    rename_dict[rt_col] = 'rt_reference'

                if rename_dict:
                    df = df.rename(columns=rename_dict)

                # column
                if 'name' not in df.columns:
                    df['name'] = [f"compound_{i + 1}" for i in range(len(df))]

                # m/zcolumn
                if 'mz' in df.columns:
                    df['mz'] = pd.to_numeric(df['mz'], errors='coerce')
                    # m/zNaN
                    df = df.dropna(subset=['mz'])

                self.compounds_df = df

                return True, f"load {len(df)} compound", df

            except Exception as e:
                return False, f"Read Excel filefailed: {str(e)}", None

        except Exception as e:
            return False, f"loadcompoundfailed: {str(e)}", None

    def _find_column(self, df: pd.DataFrame, possible_names: List[str]) -> Optional[str]:
        """DataFramecolumn"""
        for name in possible_names:
            if name in df.columns:
                return name
        return None

    def set_match_settings(self, ppm_tolerance: float = 10, rt_window: float = 30,
                           intensity_threshold: float = 1000):
        """matchParameters"""
        self.match_settings = {
            'ppm_tolerance': ppm_tolerance,
            'rt_window': rt_window,
            'intensity_threshold': intensity_threshold
        }

    def match_compounds(self, peaks_data: Union[List[Dict], pd.DataFrame, str]) -> Tuple[bool, str]:
        """matchcompoundpeakdata"""
        try:
            # compound data
            if self.compounds_df is None or len(self.compounds_df) == 0:
                return False, "loadcompound"

            # peakdataDataFrame
            if isinstance(peaks_data, str):
                # filepath
                peaks_df = pd.read_csv(peaks_data)
            elif isinstance(peaks_data, pd.DataFrame):
                peaks_df = peaks_data
            elif isinstance(peaks_data, list):
                peaks_df = pd.DataFrame(peaks_data)
            else:
                return False, "peakdata"

            # peakdatacolumn
            required_columns = ['(m/z)', 'retention_time(RT)', '']
            missing_cols = [col for col in required_columns if col not in peaks_df.columns]

            if missing_cols:
                # column
                column_mapping = {
                    '(m/z)': ['mz', 'M/Z', 'm/z'],
                    'retention_time(RT)': ['rt', 'RT', 'retention_time'],
                    '': ['intensity', 'Intensity', '']
                }

                for req_col, possible_names in column_mapping.items():
                    if req_col not in peaks_df.columns:
                        for name in possible_names:
                            if name in peaks_df.columns:
                                peaks_df = peaks_df.rename(columns={name: req_col})
                                break

            # translated note
            missing_cols = [col for col in required_columns if col not in peaks_df.columns]
            if missing_cols:
                return False, f"peakdatacolumn: {', '.join(missing_cols)}"

            # data
            peaks_df['(m/z)'] = pd.to_numeric(peaks_df['(m/z)'], errors='coerce')
            peaks_df['retention_time(RT)'] = pd.to_numeric(peaks_df['retention_time(RT)'], errors='coerce')
            peaks_df[''] = pd.to_numeric(peaks_df[''], errors='coerce')

            # peak
            peaks_df = peaks_df[peaks_df[''] >= self.match_settings['intensity_threshold']]

            self.match_results = []

            # compoundmatch
            for _, compound in self.compounds_df.iterrows():
                compound_mz = compound.get('mz')
                compound_name = compound.get('name', 'compound')
                formula = compound.get('formula', '')
                rt_reference = compound.get('rt_reference', None)

                if pd.isna(compound_mz):
                    continue

                # calculateppm
                ppm_tolerance = self.match_settings['ppm_tolerance']
                mz_tolerance = compound_mz * ppm_tolerance / 1e6

                # matchpeak
                matching_peaks = peaks_df[
                    (peaks_df['(m/z)'] >= compound_mz - mz_tolerance) &
                    (peaks_df['(m/z)'] <= compound_mz + mz_tolerance)
                    ].copy()

                # retention_time,
                if rt_reference is not None and not pd.isna(rt_reference):
                    rt_window = self.match_settings['rt_window']
                    matching_peaks = matching_peaks[
                        (matching_peaks['retention_time(RT)'] >= rt_reference - rt_window) &
                        (matching_peaks['retention_time(RT)'] <= rt_reference + rt_window)
                        ]

                if len(matching_peaks) > 0:
                    # peak
                    matching_peaks['mz_difference'] = abs(matching_peaks['(m/z)'] - compound_mz)
                    matching_peaks['mz_difference_ppm'] = matching_peaks['mz_difference'] / compound_mz * 1e6

                    # m/z
                    matching_peaks = matching_peaks.sort_values('mz_difference')

                    for i, (_, peak) in enumerate(matching_peaks.iterrows()):
                        if i >= 3: # 3match
                            break

                        result = {
                            'compound_name': compound_name,
                            'molecular_formula': formula,
                            'm/z': compound_mz,
                            'm/z': peak['(m/z)'],
                            'm/z(Da)': peak['mz_difference'],
                            'm/z(ppm)': peak['mz_difference_ppm'],
                            'retention_time(RT)': peak['retention_time(RT)'],
                            '': peak[''],
                            'match_status': 'match',
                            'match_rank': i + 1
                        }

                        # retention_time
                        if rt_reference is not None and not pd.isna(rt_reference):
                            result['RT'] = rt_reference
                            result['RT'] = abs(peak['retention_time(RT)'] - rt_reference)

                        self.match_results.append(result)
                else:
                    # match
                    result = {
                        'compound_name': compound_name,
                        'molecular_formula': formula,
                        'm/z': compound_mz,
                        'm/z': None,
                        'm/z(Da)': None,
                        'm/z(ppm)': None,
                        'retention_time(RT)': None,
                        '': None,
                        'match_status': 'unmatched',
                        'match_rank': None
                    }

                    if rt_reference is not None and not pd.isna(rt_reference):
                        result['RT'] = rt_reference

                    self.match_results.append(result)

            if len(self.match_results) == 0:
                return False, "match"

            return True, f"match: {len(self.compounds_df)} compound, {len([r for r in self.match_results if r['match_status'] == 'match'])} match"

        except Exception as e:
            return False, f"match: {str(e)}"

    def save_match_results(self, output_dir: str, base_name: str = "compound matchingresults") -> Tuple[bool, str, List[str]]:
        """savematching resultsExcelfile"""
        try:
            if not self.match_results:
                return False, "matching resultssave", []

            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)

            saved_files = []

            # DataFrame
            results_df = pd.DataFrame(self.match_results)

            # saveExcel
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

            # resultsfile
            excel_file = output_path / f"{base_name}_{timestamp}.xlsx"

            with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
                # results
                results_df.to_excel(writer, sheet_name='all matching results', index=False)

                # successful matchesresults
                successful_matches = results_df[results_df['match_status'] == 'match'].copy()
                if len(successful_matches) > 0:
                    successful_matches = successful_matches.sort_values(['compound_name', 'match_rank'])
                    successful_matches.to_excel(writer, sheet_name='successful matches', index=False)

                # compound
                compound_summary = []
                for compound_name in results_df['compound_name'].unique():
                    compound_matches = results_df[results_df['compound_name'] == compound_name]
                    successful = compound_matches[compound_matches['match_status'] == 'match']

                    summary = {
                        'compound_name': compound_name,
                        'match': len(successful),
                        'm/z(ppm)': successful['m/z(ppm)'].min() if len(successful) > 0 else None,
                        '': successful[''].max() if len(successful) > 0 else None,
                        'RT': successful['retention_time(RT)'].mean() if len(successful) > 0 else None,
                        'match_status': 'match' if len(successful) > 0 else 'unmatched'
                    }

                    compound_summary.append(summary)

                summary_df = pd.DataFrame(compound_summary)
                summary_df.to_excel(writer, sheet_name='compound', index=False)

                # translated note
                stats = {
                    '': ['compound', 'successful matches', 'unmatched', 'matchpeak',
                                 'm/z(ppm)', 'm/z(ppm)', 'm/z(ppm)'],
                    '': [
                        len(results_df['compound_name'].unique()),
                        len(results_df[results_df['match_status'] == 'match']['compound_name'].unique()),
                        len(results_df[results_df['match_status'] == 'unmatched']['compound_name'].unique()),
                        len(results_df[results_df['match_status'] == 'match']),
                        results_df['m/z(ppm)'].mean() if 'm/z(ppm)' in results_df.columns else 0,
                        results_df['m/z(ppm)'].min() if 'm/z(ppm)' in results_df.columns else 0,
                        results_df['m/z(ppm)'].max() if 'm/z(ppm)' in results_df.columns else 0
                    ]
                }
                stats_df = pd.DataFrame(stats)
                stats_df.to_excel(writer, sheet_name='', index=False)

            saved_files.append(str(excel_file))

            # Save as CSV (Use)
            csv_file = output_path / f"{base_name}_{timestamp}.csv"
            results_df.to_csv(csv_file, index=False, encoding='utf-8-sig')
            saved_files.append(str(csv_file))

            return True, f"resultssave {excel_file.name}", saved_files

        except Exception as e:
            return False, f"saveresultsfailed: {str(e)}", []


class MzMLPeakExtractor:
    """mzMLpeak - peakdatacompound matching"""

    def __init__(self, progress_callback=None, log_callback=None):
        """"""
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

    def extract_peaks_with_pymzml(self, mzml_path: str, intensity_threshold: float = 1000) -> Tuple[
        Optional[List[Dict]], str]:
        """UsepymzMLpeakdata"""
        try:
            import pymzml

            self.log(f"UsepymzMLloadfile: {Path(mzml_path).name}")
            peaks = []

            run = pymzml.run.Reader(mzml_path)
            spectrum_count = 0

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
                if spectrum_count % 100 == 0:
                    progress = min(100, spectrum_count // 10)
                    self.update_progress(progress, f" {spectrum_count}")

            if peaks:
                return peaks, f"pymzML: {len(peaks)}peak"
            else:
                return None, "pymzMLconditionpeak"

        except ImportError:
            return None, "pymzML"
        except Exception as e:
            return None, f"pymzMLfailed: {str(e)}"

    def extract_peaks_from_mzml(self, mzml_file: str, intensity_threshold: float = 1000) -> Tuple[
        bool, str, Optional[pd.DataFrame]]:
        """Extract peak data from mzML files"""
        try:
            self.log("peakdata...")
            self.update_progress(10, "peakdata...")

            # peak
            peaks, msg = self.extract_peaks_with_pymzml(mzml_file, intensity_threshold)

            if peaks is None:
                return False, msg, None

            # DataFrame
            df = pd.DataFrame(peaks)

            # columncompound matching
            column_mapping = {
                '(m/z)': '(m/z)',
                'retention_time(RT)': 'retention_time(RT)',
                '': ''
            }

            for old_col, new_col in column_mapping.items():
                if old_col in df.columns:
                    df = df.rename(columns={old_col: new_col})

            self.update_progress(100, "peak")
            return True, f" {len(df)} peak", df

        except Exception as e:
            return False, f"peakdatafailed: {str(e)}", None

    def extract_peaks_from_csv(self, csv_file: str) -> Tuple[bool, str, Optional[pd.DataFrame]]:
        """CSVfileloadpeakdata"""
        try:
            self.log("CSVfileloadpeakdata...")

            df = pd.read_csv(csv_file)

            # column
            required_columns = ['(m/z)', 'retention_time(RT)', '']

            # column
            column_mapping = {
                '(m/z)': ['mz', 'M/Z', 'm/z'],
                'retention_time(RT)': ['rt', 'RT', 'retention_time'],
                '': ['intensity', 'Intensity', '']
            }

            for req_col, possible_names in column_mapping.items():
                if req_col not in df.columns:
                    for name in possible_names:
                        if name in df.columns:
                            df = df.rename(columns={name: req_col})
                            break

            # translated note
            missing_cols = [col for col in required_columns if col not in df.columns]
            if missing_cols:
                return False, f"CSVfilecolumn: {', '.join(missing_cols)}", None

            return True, f"load {len(df)} peak", df

        except Exception as e:
            return False, f"loadCSVfilefailed: {str(e)}", None


class MzMLCompoundMatcherGUI:
    """mzMLpeakcompound matchingGUI"""

    def __init__(self, root):
        """GUI"""
        self.root = root
        self.root.title("mzMLpeakcompound matching")
        self.root.geometry("1000x800")

        # translated note
        try:
            self.root.iconbitmap(default='icon.ico')
        except:
            pass

        # translated note
        self.peak_extractor = None
        self.compound_matcher = None
        self.processing_thread = None
        self.is_processing = False

        # UI
        self.setup_ui()

        # translated note
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)

    def setup_ui(self):
        """UI"""
        # translated note
        main_frame = tk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Canvas
        canvas = tk.Canvas(main_frame)
        scrollbar = ttk.Scrollbar(main_frame, orient=tk.VERTICAL, command=canvas.yview)
        scrollable_frame = tk.Frame(canvas)

        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )

        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # translated note
        content_frame = tk.Frame(scrollable_frame, padx=20, pady=20)
        content_frame.pack(fill=tk.BOTH, expand=True)

        # translated note
        title_label = tk.Label(content_frame, text="mzMLpeakcompound matching",
                               font=("Arial", 18, "bold"))
        title_label.pack(pady=(0, 10))

        # Version
        version_label = tk.Label(content_frame, text="Version 7.0 - compound matching",
                                 font=("Arial", 10))
        version_label.pack(pady=(0, 20))

        # ==================== 1: peakdata ====================
        step1_frame = tk.LabelFrame(content_frame, text="1: peakdata", font=("Arial", 12, "bold"),
                                    padx=10, pady=10)
        step1_frame.pack(fill=tk.X, pady=(0, 20))

        # data
        self.data_source_var = tk.StringVar(value="mzml")
        ttk.Radiobutton(step1_frame, text="Extract peak data from mzML files",
                        variable=self.data_source_var, value="mzml",
                        command=self.on_data_source_change).grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)

        ttk.Radiobutton(step1_frame, text="CSVfileloadpeakdata",
                        variable=self.data_source_var, value="csv",
                        command=self.on_data_source_change).grid(row=0, column=1, sticky=tk.W, padx=5, pady=5)

        # mzMLfile
        self.mzml_frame = tk.Frame(step1_frame)
        self.mzml_frame.grid(row=1, column=0, columnspan=2, sticky=tk.W, padx=5, pady=5)

        ttk.Label(self.mzml_frame, text="mzMLfile:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5))

        self.mzml_file_var = tk.StringVar()
        mzml_file_entry = ttk.Entry(self.mzml_frame, textvariable=self.mzml_file_var, width=60)
        mzml_file_entry.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=(0, 5))

        ttk.Button(self.mzml_frame, text="...", command=self.browse_mzml_file).grid(row=0, column=2, padx=(0, 5))

        # CSVfile
        self.csv_frame = tk.Frame(step1_frame)
        self.csv_frame.grid(row=2, column=0, columnspan=2, sticky=tk.W, padx=5, pady=5)

        ttk.Label(self.csv_frame, text="CSVfile:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5))

        self.csv_file_var = tk.StringVar()
        csv_file_entry = ttk.Entry(self.csv_frame, textvariable=self.csv_file_var, width=60)
        csv_file_entry.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=(0, 5))

        ttk.Button(self.csv_frame, text="...", command=self.browse_csv_file).grid(row=0, column=2, padx=(0, 5))

        # translated note
        ttk.Label(step1_frame, text=":").grid(row=3, column=0, sticky=tk.W, padx=5, pady=5)

        self.intensity_threshold_var = tk.StringVar(value="1000")
        threshold_entry = ttk.Entry(step1_frame, textvariable=self.intensity_threshold_var, width=15)
        threshold_entry.grid(row=3, column=1, sticky=tk.W, padx=5, pady=5)

        # CSV
        self.csv_frame.grid_remove()

        # ==================== 2: compoundload ====================
        step2_frame = tk.LabelFrame(content_frame, text="2: loadcompound", font=("Arial", 12, "bold"),
                                    padx=10, pady=10)
        step2_frame.pack(fill=tk.X, pady=(0, 20))

        ttk.Label(step2_frame, text="compoundExcelfile:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)

        self.compound_file_var = tk.StringVar()
        compound_file_entry = ttk.Entry(step2_frame, textvariable=self.compound_file_var, width=60)
        compound_file_entry.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=(0, 5), pady=5)

        ttk.Button(step2_frame, text="...", command=self.browse_compound_file).grid(row=0, column=2, pady=5)

        # translated note
        ttk.Button(step2_frame, text="compound", command=self.preview_compounds).grid(row=1, column=0, sticky=tk.W,
                                                                                        padx=5, pady=5)

        # load
        self.compound_status_var = tk.StringVar(value="loadcompound")
        ttk.Label(step2_frame, textvariable=self.compound_status_var).grid(row=1, column=1, columnspan=2, sticky=tk.W,
                                                                           padx=5, pady=5)

        # ==================== 3: matchParameters ====================
        step3_frame = tk.LabelFrame(content_frame, text="3: matchParameters", font=("Arial", 12, "bold"),
                                    padx=10, pady=10)
        step3_frame.pack(fill=tk.X, pady=(0, 20))

        # m/z (ppm)
        ttk.Label(step3_frame, text="m/z (ppm):").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)

        self.ppm_tolerance_var = tk.StringVar(value="10")
        ppm_entry = ttk.Entry(step3_frame, textvariable=self.ppm_tolerance_var, width=15)
        ppm_entry.grid(row=0, column=1, sticky=tk.W, padx=5, pady=5)

        ttk.Label(step3_frame, text="ppm").grid(row=0, column=2, sticky=tk.W, padx=(0, 20), pady=5)

        # RT ()
        ttk.Label(step3_frame, text="RT ():").grid(row=0, column=3, sticky=tk.W, padx=5, pady=5)

        self.rt_window_var = tk.StringVar(value="30")
        rt_entry = ttk.Entry(step3_frame, textvariable=self.rt_window_var, width=15)
        rt_entry.grid(row=0, column=4, sticky=tk.W, padx=5, pady=5)

        ttk.Label(step3_frame, text="").grid(row=0, column=5, sticky=tk.W, padx=(0, 5), pady=5)

        # outputdirectory
        ttk.Label(step3_frame, text="outputdirectory:").grid(row=1, column=0, sticky=tk.W, padx=5, pady=(10, 5))

        self.output_dir_var = tk.StringVar()
        output_dir_entry = ttk.Entry(step3_frame, textvariable=self.output_dir_var, width=60)
        output_dir_entry.grid(row=1, column=1, columnspan=4, sticky=(tk.W, tk.E), padx=(0, 5), pady=(10, 5))

        ttk.Button(step3_frame, text="...", command=self.browse_output_dir).grid(row=1, column=5, pady=(10, 5))

        # ==================== 4: ====================
        step4_frame = tk.LabelFrame(content_frame, text="4: ", font=("Arial", 12, "bold"),
                                    padx=10, pady=10)
        step4_frame.pack(fill=tk.X, pady=(0, 20))

        self.progress_var = tk.StringVar(value="")
        ttk.Label(step4_frame, textvariable=self.progress_var).pack(anchor=tk.W, pady=(0, 5))

        self.progress_bar = ttk.Progressbar(step4_frame, mode='determinate', length=900)
        self.progress_bar.pack(fill=tk.X, pady=(0, 5))

        self.progress_percent = tk.StringVar(value="0%")
        ttk.Label(step4_frame, textvariable=self.progress_percent).pack(anchor=tk.E)

        # ==================== 5: ====================
        step5_frame = tk.LabelFrame(content_frame, text="5: ", font=("Arial", 12, "bold"),
                                    padx=10, pady=10)
        step5_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 20))

        # translated note
        self.log_text = scrolledtext.ScrolledText(step5_frame, width=120, height=15,
                                                  wrap=tk.WORD, font=("Consolas", 10))
        self.log_text.pack(fill=tk.BOTH, expand=True)

        # translated note
        self.log_text.tag_config("INFO", foreground="black")
        self.log_text.tag_config("SUCCESS", foreground="green")
        self.log_text.tag_config("WARNING", foreground="orange")
        self.log_text.tag_config("ERROR", foreground="red")

        # ==================== ====================
        button_frame = tk.Frame(content_frame)
        button_frame.pack(fill=tk.X, pady=(0, 10))

        self.process_button = ttk.Button(button_frame, text="",
                                         command=self.start_processing, width=20)
        self.process_button.pack(side=tk.LEFT, padx=5)

        self.cancel_button = ttk.Button(button_frame, text="",
                                        command=self.cancel_processing, width=20,
                                        state=tk.DISABLED)
        self.cancel_button.pack(side=tk.LEFT, padx=5)

        ttk.Button(button_frame, text="", command=self.clear_log).pack(side=tk.LEFT, padx=5)

        ttk.Button(button_frame, text="outputdirectory", command=self.open_output_dir).pack(side=tk.LEFT, padx=5)

        ttk.Button(button_frame, text="", command=self.on_closing).pack(side=tk.LEFT, padx=5)

        # translated note
        self.status_var = tk.StringVar(value="")
        status_bar = ttk.Label(content_frame, textvariable=self.status_var, relief=tk.SUNKEN)
        status_bar.pack(fill=tk.X, pady=(0, 10))

        # translated note
        self.on_data_source_change()

    def on_data_source_change(self):
        """data"""
        source = self.data_source_var.get()

        if source == "mzml":
            self.mzml_frame.grid()
            self.csv_frame.grid_remove()
        else:
            self.mzml_frame.grid_remove()
            self.csv_frame.grid()

    def browse_mzml_file(self):
        """mzMLfile"""
        file_types = [("mzMLfile", "*.mzML *.mzML.gz"), ("file", "*.*")]
        file_path = filedialog.askopenfilename(title="mzMLfile", filetypes=file_types)

        if file_path:
            self.mzml_file_var.set(file_path)

    def browse_csv_file(self):
        """CSVfile"""
        file_types = [("CSVfile", "*.csv"), ("file", "*.*")]
        file_path = filedialog.askopenfilename(title="CSVfile", filetypes=file_types)

        if file_path:
            self.csv_file_var.set(file_path)

    def browse_compound_file(self):
        """compoundExcelfile"""
        file_types = [("Excelfile", "*.xlsx *.xls"), ("file", "*.*")]
        file_path = filedialog.askopenfilename(title="compoundExcelfile", filetypes=file_types)

        if file_path:
            self.compound_file_var.set(file_path)
            # translated note
            self.preview_compounds()

    def browse_output_dir(self):
        """outputdirectory"""
        dir_path = filedialog.askdirectory(title="outputdirectory")

        if dir_path:
            self.output_dir_var.set(dir_path)

    def preview_compounds(self):
        """compound"""
        compound_file = self.compound_file_var.get()

        if not compound_file:
            messagebox.showwarning("Warning", "compoundExcelfile")
            return

        try:
            matcher = CompoundMatcher()
            success, msg, df = matcher.load_compounds_from_excel(compound_file)

            if success:
                # translated note
                preview_window = tk.Toplevel(self.root)
                preview_window.title("compound")
                preview_window.geometry("800x600")

                # Treeviewdata
                tree_frame = tk.Frame(preview_window)
                tree_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

                # translated note
                tree_scroll_y = ttk.Scrollbar(tree_frame)
                tree_scroll_y.pack(side=tk.RIGHT, fill=tk.Y)

                tree_scroll_x = ttk.Scrollbar(tree_frame, orient=tk.HORIZONTAL)
                tree_scroll_x.pack(side=tk.BOTTOM, fill=tk.X)

                # Treeview
                tree = ttk.Treeview(tree_frame, yscrollcommand=tree_scroll_y.set,
                                    xscrollcommand=tree_scroll_x.set)
                tree.pack(fill=tk.BOTH, expand=True)

                tree_scroll_y.config(command=tree.yview)
                tree_scroll_x.config(command=tree.xview)

                # column
                tree["columns"] = list(df.columns)
                tree["show"] = "headings"

                # column
                for col in df.columns:
                    tree.heading(col, text=col)
                    tree.column(col, width=100)

                # data
                for _, row in df.iterrows():
                    tree.insert("", tk.END, values=list(row))

                # translated note
                info_frame = tk.Frame(preview_window)
                info_frame.pack(fill=tk.X, padx=10, pady=(0, 10))

                tk.Label(info_frame, text=f"load {len(df)} compound").pack(side=tk.LEFT)

                if 'mz' in df.columns:
                    mz_range = f"m/z: {df['mz'].min():.4f} - {df['mz'].max():.4f}"
                    tk.Label(info_frame, text=mz_range).pack(side=tk.LEFT, padx=20)

                # translated note
                tk.Button(preview_window, text="", command=preview_window.destroy).pack(pady=10)

                self.compound_status_var.set(f"load {len(df)} compound")
                self.log_message(f"✓ load {len(df)} compound", "SUCCESS")
            else:
                messagebox.showerror("", msg)
                self.compound_status_var.set("loadfailed")
                self.log_message(f"✗ loadcompoundfailed: {msg}", "ERROR")

        except Exception as e:
            messagebox.showerror("", f"failed: {str(e)}")

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
        # data
        source = self.data_source_var.get()

        if source == "mzml":
            mzml_file = self.mzml_file_var.get().strip()
            if not mzml_file:
                return False, "mzMLfile"
            if not os.path.exists(mzml_file):
                return False, f"mzMLfile does not exist: {mzml_file}"
        else:
            csv_file = self.csv_file_var.get().strip()
            if not csv_file:
                return False, "CSVfile"
            if not os.path.exists(csv_file):
                return False, f"CSVfile does not exist: {csv_file}"

        # compoundfile
        compound_file = self.compound_file_var.get().strip()
        if not compound_file:
            return False, "compoundExcelfile"
        if not os.path.exists(compound_file):
            return False, f"compoundfile does not exist: {compound_file}"

        # outputdirectory
        output_dir = self.output_dir_var.get().strip()
        if not output_dir:
            return False, "outputdirectory"

        # Parameters
        try:
            ppm_tolerance = float(self.ppm_tolerance_var.get())
            if ppm_tolerance <= 0:
                return False, "m/z0"
        except ValueError:
            return False, "m/z"

        try:
            rt_window = float(self.rt_window_var.get())
            if rt_window < 0:
                return False, "RT"
        except ValueError:
            return False, "RT"

        try:
            intensity_threshold = float(self.intensity_threshold_var.get())
            if intensity_threshold < 0:
                return False, ""
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
        source = self.data_source_var.get()
        compound_file = self.compound_file_var.get()
        output_dir = self.output_dir_var.get()

        ppm_tolerance = float(self.ppm_tolerance_var.get())
        rt_window = float(self.rt_window_var.get())
        intensity_threshold = float(self.intensity_threshold_var.get())

        # translated note
        self.peak_extractor = MzMLPeakExtractor(
            progress_callback=self.update_progress,
            log_callback=self.log_message
        )

        self.compound_matcher = CompoundMatcher()

        # translated note
        if source == "mzml":
            mzml_file = self.mzml_file_var.get()
            self.processing_thread = threading.Thread(
                target=self.process_mzml_match,
                args=(mzml_file, compound_file, output_dir,
                      ppm_tolerance, rt_window, intensity_threshold)
            )
        else:
            csv_file = self.csv_file_var.get()
            self.processing_thread = threading.Thread(
                target=self.process_csv_match,
                args=(csv_file, compound_file, output_dir,
                      ppm_tolerance, rt_window, intensity_threshold)
            )

        self.is_processing = True
        self.processing_thread.start()

        # translated note
        self.root.after(100, self.check_processing_status)

    def process_mzml_match(self, mzml_file: str, compound_file: str, output_dir: str,
                           ppm_tolerance: float, rt_window: float, intensity_threshold: float):
        """mzMLfilecompound matching ()"""
        try:
            self.log_message("1: Extract peak data from mzML files...")
            self.update_progress(20, "peakdata...")

            # peakdata
            success, msg, peaks_df = self.peak_extractor.extract_peaks_from_mzml(
                mzml_file, intensity_threshold
            )

            if not success:
                self.root.after(0, lambda: messagebox.showerror("", msg))
                self.processing_finished()
                return

            self.log_message(f"✓ {len(peaks_df)} peak", "SUCCESS")

            # compound matching
            self.process_compound_matching(peaks_df, compound_file, output_dir,
                                           ppm_tolerance, rt_window, mzml_file)

        except Exception as e:
            error_msg = f": {str(e)}"
            self.root.after(0, lambda: self.log_message(error_msg, "ERROR"))
            self.root.after(0, self.processing_finished)

    def process_csv_match(self, csv_file: str, compound_file: str, output_dir: str,
                          ppm_tolerance: float, rt_window: float, intensity_threshold: float):
        """CSVfilecompound matching ()"""
        try:
            self.log_message("1: CSVfileloadpeakdata...")
            self.update_progress(20, "loadpeakdata...")

            # loadpeakdata
            success, msg, peaks_df = self.peak_extractor.extract_peaks_from_csv(csv_file)

            if not success:
                self.root.after(0, lambda: messagebox.showerror("", msg))
                self.processing_finished()
                return

            self.log_message(f"✓ load {len(peaks_df)} peak", "SUCCESS")

            # compound matching
            self.process_compound_matching(peaks_df, compound_file, output_dir,
                                           ppm_tolerance, rt_window, csv_file)

        except Exception as e:
            error_msg = f": {str(e)}"
            self.root.after(0, lambda: self.log_message(error_msg, "ERROR"))
            self.root.after(0, self.processing_finished)

    def process_compound_matching(self, peaks_df: pd.DataFrame, compound_file: str,
                                  output_dir: str, ppm_tolerance: float, rt_window: float,
                                  source_file: str):
        """compound matching"""
        try:
            self.log_message("2: loadcompound...")
            self.update_progress(40, "loadcompound...")

            # loadcompound
            self.compound_matcher.set_match_settings(ppm_tolerance, rt_window)

            success, msg, _ = self.compound_matcher.load_compounds_from_excel(compound_file)
            if not success:
                self.root.after(0, lambda: messagebox.showerror("", msg))
                self.processing_finished()
                return

            self.log_message(f"✓ loadcompound", "SUCCESS")

            self.log_message("3: matchcompound...")
            self.update_progress(60, "matchcompound...")

            # match
            success, msg = self.compound_matcher.match_compounds(peaks_df)
            if not success:
                self.root.after(0, lambda: messagebox.showwarning("Warning", msg))

            self.log_message(f"✓ compound matching", "SUCCESS")

            self.log_message("4: savematching results...")
            self.update_progress(80, "savematching results...")

            # saveresults
            base_name = Path(source_file).stem
            success, msg, saved_files = self.compound_matcher.save_match_results(
                output_dir, f"{base_name}_compound matching"
            )

            if success:
                self.log_message(f"✓ resultssave", "SUCCESS")
                self.log_message(f" file:", "INFO")
                for file_path in saved_files:
                    self.log_message(f"    • {Path(file_path).name}", "INFO")

                # translated note
                self.root.after(0, lambda: messagebox.showinfo(
                    "",
                    f"compound matching！\nresultssave: {output_dir}"
                ))
            else:
                self.log_message(f"✗ saveresultsfailed: {msg}", "ERROR")
                self.root.after(0, lambda: messagebox.showerror("", f"saveresultsfailed: {msg}"))

            self.update_progress(100, "")
            self.processing_finished()

        except Exception as e:
            error_msg = f"compound matching: {str(e)}"
            self.root.after(0, lambda: self.log_message(error_msg, "ERROR"))
            self.root.after(0, self.processing_finished)

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
        self.peak_extractor = None
        self.compound_matcher = None

        # ,
        self.process_button.config(state=tk.NORMAL)
        self.cancel_button.config(state=tk.DISABLED)

        # translated note
        self.update_status("")

    def cancel_processing(self):
        """"""
        if self.peak_extractor and self.is_processing:
            self.peak_extractor.cancel_requested = True
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
                if self.peak_extractor:
                    self.peak_extractor.cancel_requested = True
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
        'openpyxl': 'Excelfile',
        'pymzml': 'mzMLfile (peak)'
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
    window_width = 1000
    window_height = 800
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()
    x = (screen_width - window_width) // 2
    y = (screen_height - window_height) // 2
    root.geometry(f"{window_width}x{window_height}+{x}+{y}")

    # GUI
    app = MzMLCompoundMatcherGUI(root)

    # translated note
    root.mainloop()


if __name__ == "__main__":
    print("mzMLpeakcompound matching")
    print("Version: 7.0")
    print(":")
    print(" 1. Extract peak data from mzML files")
    print(" 2. CSVfileloadpeakdata")
    print(" 3. Excelfileloadcompound")
    print(" 4. Match peak data by compound m/z")
    print(" 5. matching resultsExcelfile")
    print("=" * 70)

    main()