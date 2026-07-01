#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PPG retention indicesfile

:
1. Load PPG standard data
2. loadPPG retention indicesfile
3. retention_indexretention_time
4. saveresults
"""

import os
import sys
import threading
import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import warnings

warnings.filterwarnings('ignore')


class IndexToRTConverter:
    """retention_indexretention_time"""

    def __init__(self):
        self.ppg_data = {} # PPG standard data
        self.standard_curves = {} # standard curveParameters
        self.index_data = None # retention_indexdata
        self.converted_data = None # data

    def load_ppg_data(self, file_path, condition="default"):
        """Load PPG standard data"""
        try:
            file_ext = Path(file_path).suffix.lower()

            if file_ext in ['.xlsx', '.xls']:
                df = pd.read_excel(file_path)
            elif file_ext == '.csv':
                df = pd.read_csv(file_path)
            else:
                return False, f"unsupported file format: {file_ext}"

            # Find required columns
            column_mapping = {
                'degree_of_polymerization': ['degree_of_polymerization', 'n', 'DP', 'degree_of_polymerizationn', 'PPG_n', 'PPG'],
                'retention_time': ['retention_time', 'RT', 'RetentionTime', 't_R', 'retention_time(RT)', 'rt']
            }

            # column
            rename_dict = {}
            for target_col, possible_names in column_mapping.items():
                for name in possible_names:
                    if name in df.columns:
                        rename_dict[name] = target_col
                        break

            if rename_dict:
                df = df.rename(columns=rename_dict)

            # column
            if 'degree_of_polymerization' not in df.columns or 'retention_time' not in df.columns:
                return False, "datafilecolumn ('degree_of_polymerization''retention_time')"

            # data
            df['degree_of_polymerization'] = pd.to_numeric(df['degree_of_polymerization'], errors='coerce')
            df['retention_time'] = pd.to_numeric(df['retention_time'], errors='coerce')
            df = df.dropna(subset=['degree_of_polymerization', 'retention_time'])

            # degree_of_polymerization
            df = df.sort_values('degree_of_polymerization')

            # data
            self.ppg_data[condition] = df

            return True, f"load {len(df)} PPG standard data"

        except Exception as e:
            return False, f"loadPPG datafailed: {str(e)}"

    def load_index_data(self, file_path):
        """loadretention_indexdata"""
        try:
            file_ext = Path(file_path).suffix.lower()

            if file_ext in ['.xlsx', '.xls']:
                self.index_data = pd.read_excel(file_path)
            elif file_ext == '.csv':
                self.index_data = pd.read_csv(file_path)
            else:
                return False, f"unsupported file format: {file_ext}"

            # PPGindexcolumn
            index_columns = []
            for col in self.index_data.columns:
                col_lower = col.lower()
                if any(keyword in col_lower for keyword in ['ppg', 'index', 'index', 'retention_index']):
                    index_columns.append(col)

            if not index_columns:
                # PPGindexcolumn,
                return True, "PPGindexcolumn, "

            return True, f"load {len(self.index_data)} data, PPGindexcolumn: {', '.join(index_columns)}"

        except Exception as e:
            return False, f"loadretention_indexdatafailed: {str(e)}"

    def fit_standard_curve(self, condition="default", model_type="linear"):
        """standard curve"""
        try:
            if condition not in self.ppg_data:
                return False, f"condition {condition} PPG data"

            df = self.ppg_data[condition]
            n_values = df['degree_of_polymerization'].values
            rt_values = df['retention_time'].values

            if model_type == "logarithmic":
                # model: RT = a + b * ln(n)
                x = np.log(n_values)
                model_name = "model (RT = a + b * ln(n))"
            elif model_type == "linear":
                # model: RT = a + b * n
                x = n_values
                model_name = "model (RT = a + b * n)"
            else:
                return False, f"model: {model_type}"

            # translated note
            slope, intercept, r_value, p_value, std_err = stats.linregress(x, rt_values)

            # store standard-curve parameters
            self.standard_curves[condition] = {
                'condition': condition,
                'model_type': model_type,
                'model_name': model_name,
                'slope': slope,
                'intercept': intercept,
                'r_squared': r_value ** 2,
                'p_value': p_value,
                'std_err': std_err
            }

            return True, f"standard curve: {model_name}, R² = {r_value ** 2:.6f}"

        except Exception as e:
            return False, f"standard curvefailed: {str(e)}"

    def convert_index_to_rt(self, index, condition="default", method="regression"):
        """retention_indexretention_time"""
        try:
            if condition not in self.ppg_data:
                return None, f"condition {condition} PPG data"

            if condition not in self.standard_curves:
                success, msg = self.fit_standard_curve(condition)
                if not success:
                    return None, f": {msg}"

            df_ppg = self.ppg_data[condition]
            ppg_rt = df_ppg['retention_time'].values
            ppg_n = df_ppg['degree_of_polymerization'].values

            # PPGindexdegree_of_polymerization
            n_value = index / 100

            if method == "interpolation":
                # translated note
                if n_value < ppg_n[0]:
                    # translated note
                    if len(ppg_n) >= 2:
                        rt_calc = ppg_rt[0] - (ppg_n[0] - n_value) / (ppg_n[1] - ppg_n[0]) * (ppg_rt[1] - ppg_rt[0])
                    else:
                        rt_calc = ppg_rt[0]
                elif n_value > ppg_n[-1]:
                    # translated note
                    if len(ppg_n) >= 2:
                        rt_calc = ppg_rt[-1] + (n_value - ppg_n[-1]) / (ppg_n[-1] - ppg_n[-2]) * (
                                ppg_rt[-1] - ppg_rt[-2])
                    else:
                        rt_calc = ppg_rt[-1]
                else:
                    # translated note
                    idx = np.searchsorted(ppg_n, n_value) - 1
                    if idx < 0:
                        idx = 0
                    elif idx >= len(ppg_n) - 1:
                        idx = len(ppg_n) - 2

                    n_i, n_j = ppg_n[idx], ppg_n[idx + 1]
                    rt_i, rt_j = ppg_rt[idx], ppg_rt[idx + 1]

                    rt_calc = rt_i + (rt_j - rt_i) * (n_value - n_i) / (n_j - n_i)

                method_used = ""

            elif method == "regression":
                # translated note
                curve = self.standard_curves[condition]

                if curve['model_type'] == "logarithmic":
                    # RT = a + b * ln(n)
                    rt_calc = curve['intercept'] + curve['slope'] * np.log(n_value)
                else:  # linear
                    # RT = a + b * n
                    rt_calc = curve['intercept'] + curve['slope'] * n_value

                method_used = ""
            else:
                return None, f"calculatemethod: {method}"

            return rt_calc, method_used

        except Exception as e:
            return None, f"failed: {str(e)}"

    def batch_convert_file(self, index_column, condition="default", method="regression"):
        """fileretention_index"""
        try:
            if self.index_data is None:
                return False, "loadretention_indexdata"

            if index_column not in self.index_data.columns:
                return False, f"column '{index_column}' data"

            # translated note
            self.converted_data = self.index_data.copy()

            # resultscolumn
            rt_column_name = f"retention_time_{condition}"
            method_column_name = f"calculatemethod_{condition}"

            # translated note
            rt_values = []
            method_values = []

            for idx, value in enumerate(self.index_data[index_column]):
                if pd.isna(value):
                    rt_values.append(np.nan)
                    method_values.append("data")
                else:
                    try:
                        rt, method_used = self.convert_index_to_rt(float(value), condition, method)
                        rt_values.append(rt)
                        method_values.append(method_used)
                    except Exception as e:
                        rt_values.append(np.nan)
                        method_values.append(f"failed: {str(e)}")

            # resultscolumn
            self.converted_data[rt_column_name] = rt_values
            self.converted_data[method_column_name] = method_values

            # translated note
            successful_conversions = sum(1 for m in method_values if m in ["", ""])
            total_conversions = len(self.index_data)

            return True, f": {successful_conversions}/{total_conversions} data"

        except Exception as e:
            return False, f"failed: {str(e)}"

    def save_converted_data(self, file_path):
        """savedata"""
        try:
            if self.converted_data is None:
                return False, "savedata"

            file_ext = Path(file_path).suffix.lower()

            if file_ext in ['.xlsx', '.xls']:
                self.converted_data.to_excel(file_path, index=False)
            elif file_ext == '.csv':
                self.converted_data.to_csv(file_path, index=False, encoding='utf-8')
            else:
                return False, f"unsupported file format: {file_ext}"

            return True, f"datasave: {file_path}"

        except Exception as e:
            return False, f"savedatafailed: {str(e)}"


class IndexFileConverterGUI:
    """retention_indexfileGUI"""

    def __init__(self, root):
        self.root = root
        self.root.title("PPG retention indicesfile")
        self.root.geometry("1200x800")

        # translated note
        self.converter = IndexToRTConverter()

        # UI
        self.setup_ui()

    def setup_ui(self):
        """UI"""
        # translated note
        main_frame = ttk.Frame(self.root, padding=10)
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Notebook ()
        self.notebook = ttk.Notebook(main_frame)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        # translated note
        self.setup_data_tab()
        self.setup_conversion_tab()
        self.setup_result_tab()

        # translated note
        self.status_var = tk.StringVar(value="")
        status_bar = ttk.Label(main_frame, textvariable=self.status_var,
                               relief=tk.SUNKEN, anchor=tk.W)
        status_bar.pack(fill=tk.X, pady=(5, 0))

    def setup_data_tab(self):
        """dataload"""
        data_tab = ttk.Frame(self.notebook)
        self.notebook.add(data_tab, text="dataload")

        # translated note
        data_frame = ttk.LabelFrame(data_tab, text="dataload", padding=15)
        data_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # ==================== PPG dataload ====================
        ppg_frame = ttk.LabelFrame(data_frame, text="1. Load PPG standard data", padding=10)
        ppg_frame.pack(fill=tk.X, pady=(0, 15))

        # file
        file_frame1 = ttk.Frame(ppg_frame)
        file_frame1.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(file_frame1, text="PPG datafile:").pack(side=tk.LEFT, padx=(0, 5))

        self.ppg_file_var = tk.StringVar()
        ppg_file_entry = ttk.Entry(file_frame1, textvariable=self.ppg_file_var, width=60)
        ppg_file_entry.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 10))

        ttk.Button(file_frame1, text="...", command=self.browse_ppg_file).pack(side=tk.LEFT)

        # condition
        condition_frame = ttk.Frame(ppg_frame)
        condition_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(condition_frame, text="conditionname:").pack(side=tk.LEFT, padx=(0, 5))

        self.ppg_condition_var = tk.StringVar(value="default")
        ppg_condition_entry = ttk.Entry(condition_frame, textvariable=self.ppg_condition_var, width=20)
        ppg_condition_entry.pack(side=tk.LEFT, padx=(0, 20))

        # load
        ttk.Button(ppg_frame, text="loadPPG data", command=self.load_ppg_data).pack(pady=10)

        # PPG data
        self.ppg_info_var = tk.StringVar(value="loadPPG data...")
        ttk.Label(ppg_frame, textvariable=self.ppg_info_var, wraplength=800).pack(anchor=tk.W)

        # ==================== retention_indexdataload ====================
        index_frame = ttk.LabelFrame(data_frame, text="2. loadretention_indexdata", padding=10)
        index_frame.pack(fill=tk.X, pady=(0, 15))

        # file
        file_frame2 = ttk.Frame(index_frame)
        file_frame2.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(file_frame2, text="retention_indexfile:").pack(side=tk.LEFT, padx=(0, 5))

        self.index_file_var = tk.StringVar()
        index_file_entry = ttk.Entry(file_frame2, textvariable=self.index_file_var, width=60)
        index_file_entry.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 10))

        ttk.Button(file_frame2, text="...", command=self.browse_index_file).pack(side=tk.LEFT)

        # load
        ttk.Button(index_frame, text="loadretention_indexdata", command=self.load_index_data).pack(pady=10)

        # retention_indexdata
        self.index_info_var = tk.StringVar(value="loadretention_indexdata...")
        ttk.Label(index_frame, textvariable=self.index_info_var, wraplength=800).pack(anchor=tk.W)

        # ==================== data ====================
        preview_frame = ttk.LabelFrame(data_frame, text="data", padding=10)
        preview_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 10))

        # translated note
        self.preview_text = scrolledtext.ScrolledText(preview_frame, width=100, height=15,
                                                      wrap=tk.WORD, font=("Consolas", 9))
        self.preview_text.pack(fill=tk.BOTH, expand=True)

        # translated note
        preview_controls = ttk.Frame(preview_frame)
        preview_controls.pack(fill=tk.X, pady=(5, 0))

        ttk.Button(preview_controls, text="PPG data",
                   command=self.preview_ppg_data).pack(side=tk.LEFT, padx=5)

        ttk.Button(preview_controls, text="retention_indexdata",
                   command=self.preview_index_data).pack(side=tk.LEFT, padx=5)

        ttk.Button(preview_controls, text="",
                   command=self.clear_preview).pack(side=tk.LEFT, padx=5)

    def setup_conversion_tab(self):
        """"""
        conversion_tab = ttk.Frame(self.notebook)
        self.notebook.add(conversion_tab, text="")

        # translated note
        conversion_frame = ttk.LabelFrame(conversion_tab, text="", padding=15)
        conversion_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # ==================== standard curve ====================
        curve_frame = ttk.LabelFrame(conversion_frame, text="1. PPGstandard curve", padding=10)
        curve_frame.pack(fill=tk.X, pady=(0, 15))

        # condition
        condition_frame = ttk.Frame(curve_frame)
        condition_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(condition_frame, text="condition:").pack(side=tk.LEFT, padx=(0, 5))

        self.curve_condition_var = tk.StringVar()
        self.curve_condition_combo = ttk.Combobox(condition_frame,
                                                  textvariable=self.curve_condition_var,
                                                  width=25, state="readonly")
        self.curve_condition_combo.pack(side=tk.LEFT, padx=(0, 20))

        # model
        ttk.Label(condition_frame, text="model:").pack(side=tk.LEFT, padx=(0, 5))

        self.curve_model_var = tk.StringVar(value="linear")
        curve_model_combo = ttk.Combobox(condition_frame, textvariable=self.curve_model_var,
                                         values=["linear", "logarithmic"],
                                         width=15, state="readonly")
        curve_model_combo.pack(side=tk.LEFT, padx=(0, 20))

        # translated note
        ttk.Button(curve_frame, text="standard curve", command=self.fit_standard_curve).pack(pady=10)

        # results
        self.curve_result_var = tk.StringVar(value="")
        ttk.Label(curve_frame, textvariable=self.curve_result_var,
                  font=("Arial", 10), wraplength=800).pack(anchor=tk.W, pady=(0, 10))

        # ==================== ====================
        settings_frame = ttk.LabelFrame(conversion_frame, text="2. ", padding=10)
        settings_frame.pack(fill=tk.X, pady=(0, 15))

        # retention_indexcolumn
        column_frame = ttk.Frame(settings_frame)
        column_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(column_frame, text="retention_indexcolumn:").pack(side=tk.LEFT, padx=(0, 5))

        self.index_column_var = tk.StringVar()
        self.index_column_combo = ttk.Combobox(column_frame,
                                               textvariable=self.index_column_var,
                                               width=30, state="readonly")
        self.index_column_combo.pack(side=tk.LEFT, padx=(0, 20))

        # column
        ttk.Button(column_frame, text="columncolumn",
                   command=self.refresh_columns).pack(side=tk.LEFT)

        # condition
        conv_condition_frame = ttk.Frame(settings_frame)
        conv_condition_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(conv_condition_frame, text="condition:").pack(side=tk.LEFT, padx=(0, 5))

        self.conv_condition_var = tk.StringVar()
        self.conv_condition_combo = ttk.Combobox(conv_condition_frame,
                                                 textvariable=self.conv_condition_var,
                                                 width=25, state="readonly")
        self.conv_condition_combo.pack(side=tk.LEFT, padx=(0, 20))

        # method
        ttk.Label(conv_condition_frame, text="calculatemethod:").pack(side=tk.LEFT, padx=(0, 5))

        self.conv_method_var = tk.StringVar(value="regression")
        conv_method_combo = ttk.Combobox(conv_condition_frame, textvariable=self.conv_method_var,
                                         values=["interpolation", "regression"],
                                         width=15, state="readonly")
        conv_method_combo.pack(side=tk.LEFT, padx=(0, 20))

        # ==================== ====================
        execute_frame = ttk.LabelFrame(conversion_frame, text="3. ", padding=10)
        execute_frame.pack(fill=tk.X, pady=(0, 15))

        # test
        test_frame = ttk.Frame(execute_frame)
        test_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(test_frame, text="test (inputPPGindex):").pack(side=tk.LEFT, padx=(0, 5))

        self.test_index_var = tk.StringVar(value="550")
        test_index_entry = ttk.Entry(test_frame, textvariable=self.test_index_var, width=15)
        test_index_entry.pack(side=tk.LEFT, padx=(0, 10))

        ttk.Button(test_frame, text="test", command=self.test_conversion).pack(side=tk.LEFT)

        # testresults
        self.test_result_var = tk.StringVar(value="")
        ttk.Label(execute_frame, textvariable=self.test_result_var,
                  font=("Arial", 10), wraplength=800).pack(anchor=tk.W, pady=(0, 10))

        # translated note
        ttk.Button(execute_frame, text="",
                   command=self.batch_convert, style="Accent.TButton").pack(pady=10)

        # translated note
        self.conversion_status_var = tk.StringVar(value="...")
        ttk.Label(execute_frame, textvariable=self.conversion_status_var,
                  font=("Arial", 10, "bold"), wraplength=800).pack(anchor=tk.W)

        # translated note
        style = ttk.Style()
        style.configure("Accent.TButton", font=("Arial", 10, "bold"), foreground="blue")

    def setup_result_tab(self):
        """results"""
        result_tab = ttk.Frame(self.notebook)
        self.notebook.add(result_tab, text="results")

        # translated note
        result_frame = ttk.LabelFrame(result_tab, text="results", padding=15)
        result_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # ==================== results ====================
        preview_frame = ttk.LabelFrame(result_frame, text="results", padding=10)
        preview_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 15))

        # translated note
        self.result_text = scrolledtext.ScrolledText(preview_frame, width=100, height=20,
                                                     wrap=tk.WORD, font=("Consolas", 9))
        self.result_text.pack(fill=tk.BOTH, expand=True)

        # ==================== results ====================
        stats_frame = ttk.LabelFrame(result_frame, text="", padding=10)
        stats_frame.pack(fill=tk.X, pady=(0, 15))

        self.stats_var = tk.StringVar(value="results")
        ttk.Label(stats_frame, textvariable=self.stats_var,
                  font=("Arial", 10), wraplength=800).pack(anchor=tk.W)

        # ==================== ====================
        export_frame = ttk.LabelFrame(result_frame, text="results", padding=10)
        export_frame.pack(fill=tk.X, pady=(0, 15))

        # file
        export_file_frame = ttk.Frame(export_frame)
        export_file_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(export_file_frame, text="outputfile:").pack(side=tk.LEFT, padx=(0, 5))

        self.export_file_var = tk.StringVar(value="")
        export_file_entry = ttk.Entry(export_file_frame, textvariable=self.export_file_var, width=60)
        export_file_entry.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 10))

        ttk.Button(export_file_frame, text="...", command=self.browse_export_file).pack(side=tk.LEFT)

        # file
        format_frame = ttk.Frame(export_frame)
        format_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(format_frame, text="file:").pack(side=tk.LEFT, padx=(0, 5))

        self.export_format_var = tk.StringVar(value="csv")
        ttk.Radiobutton(format_frame, text="CSV (.csv)",
                        variable=self.export_format_var, value="csv").pack(side=tk.LEFT, padx=(0, 20))
        ttk.Radiobutton(format_frame, text="Excel (.xlsx)",
                        variable=self.export_format_var, value="excel").pack(side=tk.LEFT)

        # translated note
        ttk.Button(export_frame, text="results",
                   command=self.export_results, style="Accent.TButton").pack(pady=10)

        # translated note
        self.export_status_var = tk.StringVar(value="...")
        ttk.Label(export_frame, textvariable=self.export_status_var,
                  font=("Arial", 10), wraplength=800).pack(anchor=tk.W)

    def browse_ppg_file(self):
        """PPGfile"""
        file_types = [("datafile", "*.csv *.xlsx *.xls"), ("CSVfile", "*.csv"),
                      ("Excelfile", "*.xlsx *.xls"), ("file", "*.*")]
        file_path = filedialog.askopenfilename(title="PPG standard datafile", filetypes=file_types)

        if file_path:
            self.ppg_file_var.set(file_path)

    def browse_index_file(self):
        """retention_indexfile"""
        file_types = [("datafile", "*.csv *.xlsx *.xls"), ("CSVfile", "*.csv"),
                      ("Excelfile", "*.xlsx *.xls"), ("file", "*.*")]
        file_path = filedialog.askopenfilename(title="retention_indexfile", filetypes=file_types)

        if file_path:
            self.index_file_var.set(file_path)

    def browse_export_file(self):
        """file"""
        default_name = "results"
        if self.export_format_var.get() == "csv":
            default_name += ".csv"
            file_types = [("CSVfile", "*.csv"), ("file", "*.*")]
        else:
            default_name += ".xlsx"
            file_types = [("Excelfile", "*.xlsx"), ("file", "*.*")]

        file_path = filedialog.asksaveasfilename(
            title="saveresults",
            filetypes=file_types,
            defaultextension=".csv" if self.export_format_var.get() == "csv" else ".xlsx",
            initialfile=default_name
        )

        if file_path:
            self.export_file_var.set(file_path)

    def load_ppg_data(self):
        """loadPPG data"""
        file_path = self.ppg_file_var.get().strip()
        condition = self.ppg_condition_var.get().strip()

        if not file_path:
            messagebox.showwarning("Warning", "PPG datafile")
            return

        if not condition:
            messagebox.showwarning("Warning", "inputconditionname")
            return

        # loaddata
        threading.Thread(target=self._load_ppg_data_thread,
                         args=(file_path, condition)).start()

    def _load_ppg_data_thread(self, file_path, condition):
        """loadPPG data"""
        self.update_status("loadPPG data...")

        success, msg = self.converter.load_ppg_data(file_path, condition)

        if success:
            # translated note
            df = self.converter.ppg_data[condition]
            info_text = f"✓ {msg}\n"
            info_text += f"degree_of_polymerization: {df['degree_of_polymerization'].min()} - {df['degree_of_polymerization'].max()}, "
            info_text += f"retention_time: {df['retention_time'].min():.2f} - {df['retention_time'].max():.2f} min"

            self.ppg_info_var.set(info_text)

            # condition
            conditions = list(self.converter.ppg_data.keys())
            self.curve_condition_combo['values'] = conditions
            self.conv_condition_combo['values'] = conditions

            if conditions:
                self.curve_condition_combo.set(conditions[0])
                self.conv_condition_combo.set(conditions[0])

            self.log_message(f"✓ PPG dataload: {condition}")
        else:
            self.ppg_info_var.set(f"✗ {msg}")
            self.log_message(f"✗ PPG dataloadfailed: {msg}")

        self.update_status("")

    def load_index_data(self):
        """loadretention_indexdata"""
        file_path = self.index_file_var.get().strip()

        if not file_path:
            messagebox.showwarning("Warning", "retention_indexfile")
            return

        # loaddata
        threading.Thread(target=self._load_index_data_thread,
                         args=(file_path,)).start()

    def _load_index_data_thread(self, file_path):
        """loadretention_indexdata"""
        self.update_status("loadretention_indexdata...")

        success, msg = self.converter.load_index_data(file_path)

        if success:
            # translated note
            self.index_info_var.set(f"✓ {msg}")

            # column
            self.refresh_columns()

            self.log_message(f"✓ retention_indexdataload")
        else:
            self.index_info_var.set(f"✗ {msg}")
            self.log_message(f"✗ retention_indexdataloadfailed: {msg}")

        self.update_status("")

    def refresh_columns(self):
        """column"""
        if self.converter.index_data is not None:
            columns = list(self.converter.index_data.columns)
            self.index_column_combo['values'] = columns

            # PPGindexcolumn
            found_column = False
            for col in columns:
                col_lower = col.lower()
                if any(keyword in col_lower for keyword in ['ppg', 'index', 'index', 'retention_index']):
                    self.index_column_combo.set(col)
                    found_column = True
                    break

            # column, Usecolumn
            if not found_column and columns:
                self.index_column_combo.set(columns[0])

    def preview_ppg_data(self):
        """PPG data"""
        condition = self.ppg_condition_var.get()

        if condition not in self.converter.ppg_data:
            messagebox.showwarning("Warning", "loadPPG data")
            return

        df = self.converter.ppg_data[condition]

        self.preview_text.delete(1.0, tk.END)
        self.preview_text.insert(tk.END, f"PPG standard data - {condition}\n")
        self.preview_text.insert(tk.END, "=" * 60 + "\n\n")
        self.preview_text.insert(tk.END, df.to_string())
        self.preview_text.insert(tk.END, f"\n\n {len(df)} data")

    def preview_index_data(self):
        """retention_indexdata"""
        if self.converter.index_data is None:
            messagebox.showwarning("Warning", "loadretention_indexdata")
            return

        df = self.converter.index_data

        self.preview_text.delete(1.0, tk.END)
        self.preview_text.insert(tk.END, "retention_indexdata\n")
        self.preview_text.insert(tk.END, "=" * 60 + "\n\n")

        # 20
        preview_df = df.head(20)
        self.preview_text.insert(tk.END, preview_df.to_string())

        if len(df) > 20:
            self.preview_text.insert(tk.END, f"\n\n... (20, {len(df)} )")
        else:
            self.preview_text.insert(tk.END, f"\n\n {len(df)} ")

        # column
        self.preview_text.insert(tk.END, "\n\ncolumn:\n")
        for i, col in enumerate(df.columns):
            self.preview_text.insert(tk.END, f"  {i + 1}. {col}\n")

    def clear_preview(self):
        """"""
        self.preview_text.delete(1.0, tk.END)

    def fit_standard_curve(self):
        """standard curve"""
        condition = self.curve_condition_var.get()
        model_type = self.curve_model_var.get()

        if not condition:
            messagebox.showwarning("Warning", "condition")
            return

        # translated note
        threading.Thread(target=self._fit_standard_curve_thread,
                         args=(condition, model_type)).start()

    def _fit_standard_curve_thread(self, condition, model_type):
        """standard curve"""
        self.update_status("standard curve...")

        success, msg = self.converter.fit_standard_curve(condition, model_type)

        if success:
            # results
            curve = self.converter.standard_curves[condition]
            result_text = f"✓ {msg}\n"
            result_text += f"model: {curve['model_name']}\n"
            result_text += f"R² = {curve['r_squared']:.6f}\n"
            result_text += f" = {curve['slope']:.4f}\n"
            result_text += f" = {curve['intercept']:.4f}"

            self.curve_result_var.set(result_text)
            self.log_message(f"✓ standard curve: {condition}")
        else:
            self.curve_result_var.set(f"✗ {msg}")
            self.log_message(f"✗ standard curvefailed: {msg}")

        self.update_status("")

    def test_conversion(self):
        """test"""
        try:
            index = float(self.test_index_var.get().strip())
            condition = self.conv_condition_var.get()
            method = self.conv_method_var.get()

            if not condition:
                messagebox.showwarning("Warning", "condition")
                return

            # translated note
            rt, method_used = self.converter.convert_index_to_rt(index, condition, method)

            if rt is not None:
                result_text = f"test!\n"
                result_text += f"PPGindex: {index:.2f} → retention_time: {rt:.4f} min\n"
                result_text += f"calculatemethod: {method_used}"

                self.test_result_var.set(result_text)
                self.log_message(f"✓ test: {index} → {rt:.4f} min")
            else:
                self.test_result_var.set(f"✗ failed: {method_used}")
                self.log_message(f"✗ testfailed")

        except ValueError:
            messagebox.showerror("", "inputPPGindex")
        except Exception as e:
            messagebox.showerror("", f"testfailed: {str(e)}")

    def batch_convert(self):
        """"""
        index_column = self.index_column_var.get()
        condition = self.conv_condition_var.get()
        method = self.conv_method_var.get()

        if not index_column:
            messagebox.showwarning("Warning", "retention_indexcolumn")
            return

        if not condition:
            messagebox.showwarning("Warning", "condition")
            return

        # loaddata
        if self.converter.index_data is None:
            messagebox.showwarning("Warning", "loadretention_indexdata")
            return

        if condition not in self.converter.ppg_data:
            messagebox.showwarning("Warning", f"loadcondition '{condition}' PPG data")
            return

        # translated note
        threading.Thread(target=self._batch_convert_thread,
                         args=(index_column, condition, method)).start()

    def _batch_convert_thread(self, index_column, condition, method):
        """"""
        self.update_status("...")
        self.conversion_status_var.set(", ...")

        success, msg = self.converter.batch_convert_file(index_column, condition, method)

        if success:
            # translated note
            self.conversion_status_var.set(f"✓ {msg}")

            # results
            self.show_conversion_results()

            # results
            self.notebook.select(2)

            self.log_message(f"✓ {msg}")
        else:
            self.conversion_status_var.set(f"✗ {msg}")
            self.log_message(f"✗ {msg}")

        self.update_status("")

    def show_conversion_results(self):
        """results"""
        if self.converter.converted_data is None:
            return

        df = self.converter.converted_data

        # results
        self.result_text.delete(1.0, tk.END)

        # 50
        preview_df = df.head(50)
        self.result_text.insert(tk.END, "results (50):\n")
        self.result_text.insert(tk.END, "=" * 80 + "\n\n")
        self.result_text.insert(tk.END, preview_df.to_string())

        if len(df) > 50:
            self.result_text.insert(tk.END, f"\n\n... (50, {len(df)} )")
        else:
            self.result_text.insert(tk.END, f"\n\n {len(df)} ")

        # translated note
        total_rows = len(df)

        # retention_timecolumn
        rt_columns = [col for col in df.columns if 'retention_time' in col]
        if rt_columns:
            rt_col = rt_columns[0]
            successful = df[rt_col].notna().sum()
            failed = total_rows - successful

            stats_text = f"\n\n:\n"
            stats_text += f"data: {total_rows}\n"
            stats_text += f": {successful} \n"
            stats_text += f"failed: {failed} \n"

            if successful > 0:
                rt_min = df[rt_col].min()
                rt_max = df[rt_col].max()
                rt_mean = df[rt_col].mean()
                stats_text += f"\nretention_time: {rt_min:.2f} - {rt_max:.2f} min\n"
                stats_text += f"retention_time: {rt_mean:.2f} min"

            self.stats_var.set(stats_text)

    def export_results(self):
        """results"""
        if self.converter.converted_data is None:
            messagebox.showwarning("Warning", "results")
            return

        file_path = self.export_file_var.get().strip()

        if not file_path:
            # Usefile
            default_dir = os.path.dirname(self.index_file_var.get()) if self.index_file_var.get() else ""
            default_name = "results"

            if self.export_format_var.get() == "csv":
                default_path = os.path.join(default_dir, f"{default_name}.csv")
            else:
                default_path = os.path.join(default_dir, f"{default_name}.xlsx")

            self.export_file_var.set(default_path)
            file_path = default_path

        # file
        if self.export_format_var.get() == "csv" and not file_path.endswith('.csv'):
            file_path += '.csv'
        elif self.export_format_var.get() == "excel" and not file_path.endswith(('.xlsx', '.xls')):
            file_path += '.xlsx'

        # savedata
        threading.Thread(target=self._export_results_thread,
                         args=(file_path,)).start()

    def _export_results_thread(self, file_path):
        """results"""
        self.update_status("results...")
        self.export_status_var.set(", ...")

        success, msg = self.converter.save_converted_data(file_path)

        if success:
            self.export_status_var.set(f"✓ {msg}")
            self.log_message(f"✓ results: {file_path}")

            # file
            if messagebox.askyesno("", f"filesave:\n{file_path}\n\nfile？"):
                try:
                    if sys.platform == 'win32':
                        os.startfile(file_path)
                    elif sys.platform == 'darwin':
                        os.system(f'open "{file_path}"')
                    else:
                        os.system(f'xdg-open "{file_path}"')
                except:
                    pass
        else:
            self.export_status_var.set(f"✗ {msg}")
            self.log_message(f"✗ resultsfailed: {msg}")

        self.update_status("")

    def log_message(self, message):
        """Record log messagesmessage"""
        print(f"[INFO] {message}")

    def update_status(self, message):
        """"""
        self.status_var.set(message)
        self.root.update()


def main():
    """"""
    root = tk.Tk()
    root.title("PPG retention indicesfile")

    # translated note
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()
    window_width = int(screen_width * 0.8)
    window_height = int(screen_height * 0.85)
    x = (screen_width - window_width) // 2
    y = (screen_height - window_height) // 2

    root.geometry(f"{window_width}x{window_height}+{x}+{y}")

    # translated note
    app = IndexFileConverterGUI(root)

    # translated note
    root.mainloop()


if __name__ == "__main__":
    print("=" * 70)
    print("PPG retention indicesfile")
    print("=" * 70)
    print(":")
    print(" 1. Load PPG standard data")
    print(" 2. loadPPG retention indicesfile")
    print(" 3. retention_indexretention_time")
    print(" 4. saveresults")
    print("=" * 70)
    print("Usemethod:")
    print(" 1. 'dataload'loadPPG dataretention_indexfile")
    print(" 2. ''Parameters")
    print(" 3. 'results'results")
    print("=" * 70)

    main()