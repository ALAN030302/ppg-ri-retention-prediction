#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PPG retention index calculation and visualization program - journal-ready complete version
:
1. Load retention-time data for PPG standards and compounds
2. Fit PPG standard curves and calculate linear relationships
3. Calculate compound PPG retention indices
4. Compare PPG indices across chromatographic conditions
5. Convert and validate PPG indices across conditions
6. comparison of multiple conversion schemes
7. visualization outputs with English labels, transparent backgrounds, adjustable font sizes, and publication-ready colors

designed for publication requirements in analytical chemistry and ES&T-style journals
"""

import os
import sys
import threading
import warnings
from datetime import datetime
from typing import Dict, List, Tuple, Optional, Union, Any
from pathlib import Path

warnings.filterwarnings('ignore')

try:
    import pandas as pd
    import numpy as np
    from scipy import stats
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    from matplotlib.figure import Figure
except ImportError as e:
    print("Error: please install the required libraries first")
    print("Installation command: pip install pandas numpy scipy matplotlib seaborn")
    sys.exit(1)

try:
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox, scrolledtext, Toplevel
    GUI_AVAILABLE = True
except ImportError:
    GUI_AVAILABLE = False
    print("Warning: tkinter is not installed, so the GUI is unavailable")

# =============================================================================
# PPGIndexCalculator core calculation class (complete version)
# =============================================================================
class PPGIndexCalculator:
    """Core class for PPG retention index calculation"""

    def __init__(self):
        self.ppg_data = {} # PPG data under different conditions
        self.compound_data = {} # compound data
        self.standard_curves = {} # standard curveParameters
        self.ppg_indices = {} # calculatePPGindex
        self.results_summary = {} # result summary
        self.conversion_results = {} # cross-condition conversion results

    def load_ppg_data(self, file_path: str, condition: str = "default") -> Tuple[bool, str]:
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
            rename_dict = {}
            for target_col, possible_names in column_mapping.items():
                for name in possible_names:
                    if name in df.columns:
                        rename_dict[name] = target_col
                        break
            if rename_dict:
                df = df.rename(columns=rename_dict)

            if 'degree_of_polymerization' not in df.columns or 'retention_time' not in df.columns:
                return False, "datafilecolumn ('degree_of_polymerization''retention_time')"

            df['degree_of_polymerization'] = pd.to_numeric(df['degree_of_polymerization'], errors='coerce')
            df['retention_time'] = pd.to_numeric(df['retention_time'], errors='coerce')
            df = df.dropna(subset=['degree_of_polymerization', 'retention_time'])
            df = df.sort_values('degree_of_polymerization')
            self.ppg_data[condition] = df
            return True, f"load {len(df)} PPG standard data (condition: {condition})"
        except Exception as e:
            return False, f"loadPPG datafailed: {str(e)}"

    def load_compound_data(self, file_path: str, category: str = "validation",
                           condition: str = "default") -> Tuple[bool, str]:
        """loadcompound data"""
        try:
            file_ext = Path(file_path).suffix.lower()
            if file_ext in ['.xlsx', '.xls']:
                df = pd.read_excel(file_path)
            elif file_ext == '.csv':
                df = pd.read_csv(file_path)
            else:
                return False, f"unsupported file format: {file_ext}"

            column_mapping = {
                'compound_name': ['compound_name', 'name', 'compound', 'Name', 'Compound', 'compound'],
                'retention_time': ['retention_time', 'RT', 'RetentionTime', 't_R', 'retention_time(RT)', 'rt'],
                'CAS': ['CAS', 'CAS', 'CAS No.', 'CAS']
            }
            rename_dict = {}
            for target_col, possible_names in column_mapping.items():
                for name in possible_names:
                    if name in df.columns:
                        rename_dict[name] = target_col
                        break
            if rename_dict:
                df = df.rename(columns=rename_dict)

            if 'compound_name' not in df.columns or 'retention_time' not in df.columns:
                return False, "datafilecolumn ('compound_name''retention_time')"

            df['retention_time'] = pd.to_numeric(df['retention_time'], errors='coerce')
            df = df.dropna(subset=['compound_name', 'retention_time'])

            key = f"{category}_{condition}"
            self.compound_data[key] = df
            return True, f"load {len(df)} compound data (: {category}, condition: {condition})"
        except Exception as e:
            return False, f"loadcompound datafailed: {str(e)}"

    def fit_standard_curve(self, condition: str = "default",
                           model_type: str = "linear") -> Tuple[bool, str]:
        """PPGstandard curve"""
        try:
            if condition not in self.ppg_data:
                return False, f"condition {condition} PPG data"
            df = self.ppg_data[condition]
            if len(df) < 3:
                return False, "PPG data, 3standard curve"

            x = df['degree_of_polymerization'].values
            y = df['retention_time'].values

            if model_type == "logarithmic":
                x_fit = np.log(x)
                model_name = "Logarithmic (RT = a + b·ln(n))"
            elif model_type == "linear":
                x_fit = x
                model_name = "Linear (RT = a + b·n)"
            else:
                return False, f"model: {model_type}"

            slope, intercept, r_value, p_value, std_err = stats.linregress(x_fit, y)
            y_pred = intercept + slope * x_fit
            residuals = y - y_pred

            self.standard_curves[condition] = {
                'condition': condition,
                'model_type': model_type,
                'slope': slope,
                'intercept': intercept,
                'r_squared': r_value ** 2,
                'p_value': p_value,
                'std_err': std_err,
                'x': x,
                'y': y,
                'y_pred': y_pred,
                'residuals': residuals,
                'model_name': model_name,
                'n_points': len(x)
            }
            return True, f"standard curve: {model_name}, R² = {r_value ** 2:.4f}"
        except Exception as e:
            return False, f"standard curvefailed: {str(e)}"

    def calculate_ppg_index(self, condition: str = "default",
                            method: str = "interpolation") -> Tuple[bool, str]:
        """calculatePPG retention indices"""
        try:
            if condition not in self.ppg_data:
                return False, f"condition {condition} PPG data"
            df_ppg = self.ppg_data[condition]

            if method == "interpolation":
                ppg_rt = df_ppg['retention_time'].values
                ppg_n = df_ppg['degree_of_polymerization'].values
                indices = {}
                for key in self.compound_data:
                    if condition in key:
                        df_comp = self.compound_data[key]
                        results = []
                        for _, row in df_comp.iterrows():
                            rt = row['retention_time']
                            compound_name = row['compound_name']
                            if rt < ppg_rt[0]:
                                n_calc = ppg_n[0] - (ppg_rt[0] - rt) / (ppg_rt[1] - ppg_rt[0]) * (ppg_n[1] - ppg_n[0])
                                if n_calc < 0:
                                    n_calc = 0
                                method_used = "Extrapolation (below min)"
                            elif rt > ppg_rt[-1]:
                                n_calc = ppg_n[-1] + (rt - ppg_rt[-1]) / (ppg_rt[-1] - ppg_rt[-2]) * (ppg_n[-1] - ppg_n[-2])
                                method_used = "Extrapolation (above max)"
                            else:
                                idx = np.searchsorted(ppg_rt, rt) - 1
                                if idx < 0:
                                    idx = 0
                                elif idx >= len(ppg_rt) - 1:
                                    idx = len(ppg_rt) - 2
                                rt_i, rt_j = ppg_rt[idx], ppg_rt[idx + 1]
                                n_i, n_j = ppg_n[idx], ppg_n[idx + 1]
                                n_calc = n_i + (n_j - n_i) * (rt - rt_i) / (rt_j - rt_i)
                                method_used = "Linear interpolation"
                            ppg_index = n_calc * 100
                            result = {'compound_name': compound_name, 'retention_time': rt, 'calculatePPGindex': ppg_index, 'calculatemethod': method_used}
                            for col in df_comp.columns:
                                if col not in result:
                                    result[col] = row[col]
                            results.append(result)
                        indices[key] = pd.DataFrame(results)

            elif method == "regression":
                if condition not in self.standard_curves:
                    success, msg = self.fit_standard_curve(condition, model_type="linear")
                    if not success:
                        return False, f"Use: {msg}"
                curve = self.standard_curves[condition]
                indices = {}
                for key in self.compound_data:
                    if condition in key:
                        df_comp = self.compound_data[key]
                        results = []
                        for _, row in df_comp.iterrows():
                            rt = row['retention_time']
                            compound_name = row['compound_name']
                            if curve['model_type'] == "logarithmic":
                                if curve['slope'] != 0:
                                    n_calc = np.exp((rt - curve['intercept']) / curve['slope'])
                                else:
                                    n_calc = np.nan
                            else:
                                if curve['slope'] != 0:
                                    n_calc = (rt - curve['intercept']) / curve['slope']
                                else:
                                    n_calc = np.nan
                            ppg_index = n_calc * 100 if not np.isnan(n_calc) else np.nan
                            result = {'compound_name': compound_name, 'retention_time': rt, 'calculatePPGindex': ppg_index, 'calculatemethod': "Regression"}
                            for col in df_comp.columns:
                                if col not in result:
                                    result[col] = row[col]
                            results.append(result)
                        indices[key] = pd.DataFrame(results)
            else:
                return False, f"calculatemethod: {method}"

            self.ppg_indices[condition] = {'method': method, 'indices': indices}
            return True, f"PPGindexcalculate (condition: {condition}, method: {method})"
        except Exception as e:
            return False, f"calculatePPGindexfailed: {str(e)}"

    def compare_conditions(self, conditions: List[str]) -> pd.DataFrame:
        """Compare PPG indices across chromatographic conditions"""
        comparison_results = []
        for key in self.compound_data:
            if any(cond in key for cond in conditions):
                category = key.split('_')[0]
                condition = key.split('_')[1] if '_' in key else "default"
                df_comp = self.compound_data[key]
                for _, row in df_comp.iterrows():
                    compound_name = row['compound_name']
                    rt = row['retention_time']
                    compound_data = {'compound_name': compound_name, 'data': category, f'{condition}_RT': rt}
                    for cond in conditions:
                        if cond in self.ppg_indices:
                            for data_key, indices_df in self.ppg_indices[cond]['indices'].items():
                                if cond in data_key:
                                    match = indices_df[indices_df['compound_name'] == compound_name]
                                    if not match.empty:
                                        compound_data[f'{cond}_PPGindex'] = match.iloc[0]['calculatePPGindex']
                                        break
                    comparison_results.append(compound_data)
        return pd.DataFrame(comparison_results)

    def calculate_conversion_error(self, from_condition: str, to_condition: str) -> pd.DataFrame:
        """calculatecondition (PPGindex)"""
        error_results = []
        if from_condition not in self.ppg_indices or to_condition not in self.ppg_indices:
            return pd.DataFrame()

        compounds_in_both = set()
        for key in self.compound_data:
            if from_condition in key:
                df = self.compound_data[key]
                compounds_in_both.update(df['compound_name'].tolist())
        for key in self.compound_data:
            if to_condition in key:
                df = self.compound_data[key]
                compounds_in_both.intersection_update(set(df['compound_name'].tolist()))

        for compound in compounds_in_both:
            from_ppg = to_ppg = None
            for key, indices_dict in self.ppg_indices.items():
                if from_condition in key:
                    for data_key, df_indices in indices_dict['indices'].items():
                        match = df_indices[df_indices['compound_name'] == compound]
                        if not match.empty:
                            from_ppg = match.iloc[0]['calculatePPGindex']
                            break
                if to_condition in key:
                    for data_key, df_indices in indices_dict['indices'].items():
                        match = df_indices[df_indices['compound_name'] == compound]
                        if not match.empty:
                            to_ppg = match.iloc[0]['calculatePPGindex']
                            break
            if from_ppg is not None and to_ppg is not None:
                absolute_error = abs(from_ppg - to_ppg)
                relative_error = (absolute_error / from_ppg * 100) if from_ppg != 0 else np.inf
                error_results.append({
                    'compound_name': compound,
                    f'{from_condition}_PPGindex': from_ppg,
                    f'{to_condition}_PPGindex': to_ppg,
                    '': absolute_error,
                    '(%)': relative_error
                })
        return pd.DataFrame(error_results)

    def convert_ppg_index_to_rt(self, from_condition: str, to_condition: str,
                                compound_names: List[str] = None) -> Tuple[pd.DataFrame, Union[Dict, str]]:
        """conditionPPGindexconditionretention_time"""
        try:
            if from_condition not in self.ppg_indices:
                return pd.DataFrame(), f"condition {from_condition} PPGindexdata"
            if to_condition not in self.standard_curves:
                success, msg = self.fit_standard_curve(to_condition, "linear")
                if not success:
                    return pd.DataFrame(), f"condition {to_condition} standard curve: {msg}"

            curve = self.standard_curves[to_condition]
            conversion_results = []

            for key, indices_df in self.ppg_indices[from_condition]['indices'].items():
                if from_condition in key:
                    for _, row in indices_df.iterrows():
                        compound_name = row['compound_name']
                        if compound_names and compound_name not in compound_names:
                            continue
                        ppg_index = row['calculatePPGindex']
                        if pd.isna(ppg_index):
                            continue
                        n_calc = ppg_index / 100
                        if curve['model_type'] == "logarithmic":
                            rt_pred = curve['intercept'] + curve['slope'] * np.log(n_calc)
                        else:
                            rt_pred = curve['intercept'] + curve['slope'] * n_calc

                        rt_actual = None
                        for comp_key, comp_df in self.compound_data.items():
                            if to_condition in comp_key:
                                match = comp_df[comp_df['compound_name'] == compound_name]
                                if not match.empty and 'retention_time' in match.columns:
                                    rt_actual = match.iloc[0]['retention_time']
                                    break

                        if rt_actual is not None:
                            absolute_error = abs(rt_pred - rt_actual)
                            relative_error = (absolute_error / rt_actual * 100) if rt_actual != 0 else np.nan
                        else:
                            absolute_error = np.nan
                            relative_error = np.nan

                        result = {
                            'compound_name': compound_name,
                            f'{from_condition}_PPGindex': ppg_index,
                            f'{from_condition}_degree_of_polymerization': n_calc,
                            f'{to_condition}_RT': rt_pred,
                            f'{to_condition}_RT': rt_actual,
                            '(min)': absolute_error,
                            '(%)': relative_error,
                            'condition': from_condition,
                            'condition': to_condition
                        }
                        conversion_results.append(result)

            conversion_df = pd.DataFrame(conversion_results)
            if not conversion_df.empty:
                valid_errors = conversion_df['(min)'].dropna()
                if len(valid_errors) > 0:
                    stats_dict = {
                        '': valid_errors.mean(),
                        '': valid_errors.std(),
                        '': valid_errors.max(),
                        '': valid_errors.min(),
                        '': valid_errors.median(),
                        '': len(valid_errors)
                    }
                    valid_rel = conversion_df['(%)'].dropna()
                    if len(valid_rel) > 0:
                        stats_dict.update({
                            '(%)': valid_rel.mean(),
                            '(%)': valid_rel.std(),
                            '(%)': valid_rel.max(),
                            '(%)': valid_rel.min()
                        })
                    return conversion_df, stats_dict
                else:
                    return conversion_df, "data"
            else:
                return conversion_df, "matchcompound data"
        except Exception as e:
            return pd.DataFrame(), f"failed: {str(e)}"

    def cross_condition_analysis(self, from_condition: str, to_condition: str,
                                 threshold: float = 0.5) -> Dict[str, Any]:
        """conditionanalysis"""
        try:
            conversion_df, stats = self.convert_ppg_index_to_rt(from_condition, to_condition)
            if conversion_df.empty:
                return {"error": "data"}

            analysis_results = {
                'condition': from_condition,
                'condition': to_condition,
                '': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                'compound': len(conversion_df),
                '': len(conversion_df['(min)'].dropna()),
                '': stats if isinstance(stats, dict) else stats,
                'data': conversion_df.to_dict('records')
            }

            if '(min)' in conversion_df.columns:
                errors = conversion_df['(min)'].dropna()
                error_bins = [0, 0.1, 0.2, 0.5, 1.0, float('inf')]
                error_labels = ['<0.1 min', '0.1-0.2 min', '0.2-0.5 min', '0.5-1.0 min', '>1.0 min']
                error_dist = {}
                for i in range(len(error_bins)-1):
                    lower, upper = error_bins[i], error_bins[i+1]
                    if i == len(error_bins)-2:
                        count = len(errors[errors >= lower])
                    else:
                        count = len(errors[(errors >= lower) & (errors < upper)])
                    error_dist[error_labels[i]] = count
                analysis_results[''] = error_dist

                passed = len(errors[errors <= threshold])
                pass_rate = (passed / len(errors) * 100) if len(errors) > 0 else 0
                analysis_results['analysis'] = {
                    '(min)': threshold,
                    '': passed,
                    '': len(errors),
                    '(%)': pass_rate
                }

            key = f"{from_condition}_to_{to_condition}"
            self.conversion_results[key] = analysis_results
            return analysis_results
        except Exception as e:
            return {"error": f"analysisfailed: {str(e)}"}

    def generate_summary_report(self) -> Dict[str, Any]:
        """analysis"""
        summary = {
            '': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'PPG datacondition': len(self.ppg_data),
            'compound data': len(self.compound_data),
            'standard curve': len(self.standard_curves),
            'PPGindexcalculateresults': len(self.ppg_indices),
            'cross-condition conversion results': len(self.conversion_results),
            'standard curve': {},
            'PPGindex': {},
            'condition': {}
        }

        for condition, curve in self.standard_curves.items():
            summary['standard curve'][condition] = {
                'model': curve['model_type'],
                'R²': curve['r_squared'],
                '': curve['slope'],
                '': curve['intercept'],
                '': curve['std_err'],
                'data': curve['n_points']
            }

        for condition, indices_data in self.ppg_indices.items():
            all_indices = []
            for df in indices_data['indices'].values():
                if 'calculatePPGindex' in df.columns:
                    all_indices.extend(df['calculatePPGindex'].dropna().tolist())
            if all_indices:
                arr = np.array(all_indices)
                summary['PPGindex'][condition] = {
                    'calculatemethod': indices_data['method'],
                    '': len(all_indices),
                    '': np.mean(arr),
                    '': np.std(arr),
                    '': np.min(arr),
                    '': np.max(arr),
                    '': np.median(arr)
                }

        for key, conv in self.conversion_results.items():
            summary['condition'][key] = {
                'condition': conv.get('condition', ''),
                'condition': conv.get('condition', ''),
                'compound': conv.get('compound', 0),
                '': conv.get('', 0),
                '': conv.get('', {}).get('', 0) if isinstance(conv.get(''), dict) else 0,
                '(%)': conv.get('analysis', {}).get('(%)', 0) if isinstance(conv.get('analysis'), dict) else 0
            }

        self.results_summary = summary
        return summary

    def save_results(self, output_dir: str) -> Tuple[bool, str, List[str]]:
        """saveresultsfile"""
        try:
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            saved_files = []

            # 1. standard curve
            if self.standard_curves:
                curves_data = []
                for condition, curve in self.standard_curves.items():
                    curves_data.append({
                        'chromatographic condition': condition,
                        'model': curve['model_type'],
                        '': curve['slope'],
                        '': curve['intercept'],
                        'R²': curve['r_squared'],
                        'p': curve['p_value'],
                        '': curve['std_err'],
                        'data': curve['n_points']
                    })
                curves_df = pd.DataFrame(curves_data)
                curves_file = output_path / f"PPG_Calibration_Curves_{timestamp}.xlsx"
                with pd.ExcelWriter(curves_file, engine='openpyxl') as writer:
                    curves_df.to_excel(writer, sheet_name='Summary', index=False)
                    for condition, curve in self.standard_curves.items():
                        detail_df = pd.DataFrame({
                            'degree_of_polymerization': curve['x'],
                            'RT': curve['y'],
                            'RT': curve['y_pred'],
                            '': curve['residuals']
                        })
                        detail_df.to_excel(writer, sheet_name=f'{condition}_details', index=False)
                saved_files.append(str(curves_file))

            # 2. PPGindex
            if self.ppg_indices:
                for condition, indices_data in self.ppg_indices.items():
                    indices_file = output_path / f"PPG_Indices_{condition}_{timestamp}.xlsx"
                    with pd.ExcelWriter(indices_file, engine='openpyxl') as writer:
                        for key, df in indices_data['indices'].items():
                            sheet_name = key.replace('_', '-')[:30]
                            df.to_excel(writer, sheet_name=sheet_name, index=False)
                    saved_files.append(str(indices_file))

            # 3. condition
            if self.conversion_results:
                conv_file = output_path / f"Conversion_Results_{timestamp}.xlsx"
                with pd.ExcelWriter(conv_file, engine='openpyxl') as writer:
                    for key, conv in self.conversion_results.items():
                        if 'data' in conv:
                            df = pd.DataFrame(conv['data'])
                            sheet_name = key[:30]
                            df.to_excel(writer, sheet_name=sheet_name, index=False)
                    # translated note
                    stats_list = []
                    for key, conv in self.conversion_results.items():
                        s = {
                            '': key,
                            'condition': conv.get('condition', ''),
                            'condition': conv.get('condition', ''),
                            'compound': conv.get('compound', 0),
                            '': conv.get('', 0)
                        }
                        if isinstance(conv.get(''), dict):
                            s.update({
                                '(min)': conv[''].get('', 0),
                                '': conv[''].get('', 0),
                                '(min)': conv[''].get('', 0)
                            })
                        if isinstance(conv.get('analysis'), dict):
                            s.update({
                                '(%)': conv['analysis'].get('(%)', 0),
                                '': conv['analysis'].get('', 0),
                                '(min)': conv['analysis'].get('(min)', 0.5)
                            })
                        stats_list.append(s)
                    if stats_list:
                        pd.DataFrame(stats_list).to_excel(writer, sheet_name='Conversion_Stats', index=False)
                saved_files.append(str(conv_file))

            # 4. analysis
            if self.results_summary:
                report_file = output_path / f"Analysis_Report_{timestamp}.txt"
                with open(report_file, 'w', encoding='utf-8') as f:
                    f.write("="*70 + "\n")
                    f.write("PPG Retention Index Analysis Report\n")
                    f.write("="*70 + "\n\n")
                    f.write(f"Generated: {self.results_summary['']}\n\n")
                    f.write(f"PPG datasets: {self.results_summary['PPG datacondition']}\n")
                    f.write(f"Compound datasets: {self.results_summary['compound data']}\n")
                    f.write(f"Calibration curves: {self.results_summary['standard curve']}\n")
                    f.write(f"PPG index calculations: {self.results_summary['PPGindexcalculateresults']}\n")
                    f.write(f"Cross‑condition conversions: {self.results_summary['cross-condition conversion results']}\n\n")
                    if self.results_summary['standard curve']:
                        f.write("Calibration curve performance:\n")
                        for cond, perf in self.results_summary['standard curve'].items():
                            f.write(f" {cond}: R²={perf['R²']:.4f}, slope={perf['']:.4f}\n")
                    if self.results_summary['PPGindex']:
                        f.write("\nPPG index statistics:\n")
                        for cond, stat in self.results_summary['PPGindex'].items():
                            f.write(f" {cond}: n={stat['']}, mean={stat['']:.2f}, sd={stat['']:.2f}\n")
                    if self.results_summary['condition']:
                        f.write("\nConversion statistics:\n")
                        for key, stat in self.results_summary['condition'].items():
                            f.write(f" {key}: valid={stat['']}, mean error={stat['']:.3f} min, pass rate={stat['(%)']:.1f}%\n")
                saved_files.append(str(report_file))

            # 5. compound data
            if self.compound_data:
                comp_file = output_path / f"Compound_Data_Summary_{timestamp}.xlsx"
                with pd.ExcelWriter(comp_file, engine='openpyxl') as writer:
                    for key, df in self.compound_data.items():
                        sheet_name = key.replace('_', '-')[:30]
                        df.to_excel(writer, sheet_name=sheet_name, index=False)
                saved_files.append(str(comp_file))

            return True, f"resultssave {output_dir}", saved_files
        except Exception as e:
            return False, f"saveresultsfailed: {str(e)}", []

# =============================================================================
# PPGVisualizer ()
# =============================================================================
class PPGVisualizer:
    """PPG data - (, , , )"""

    def __init__(self, calculator: PPGIndexCalculator):
        self.calculator = calculator
        self.figures = {}

        # translated note
        plt.rcParams['font.family'] = 'Arial'
        plt.rcParams['axes.unicode_minus'] = False
        plt.rcParams['figure.facecolor'] = 'none'
        plt.rcParams['axes.facecolor'] = 'none'
        plt.rcParams['savefig.facecolor'] = 'none'
        plt.rcParams['legend.frameon'] = False
        plt.rcParams['legend.facecolor'] = 'none'
        plt.rcParams['legend.edgecolor'] = 'none'

        # (ColorBrewer Set2)
        self.color_set = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f']

    def plot_standard_curves(self, conditions: List[str] = None, save_path: str = None,
                              fontsize: int = 12) -> plt.Figure:
        """PPGstandard curve"""
        if conditions is None:
            conditions = list(self.calculator.standard_curves.keys())
        if not conditions:
            print("Warning: No standard curve data to plot.")
            return None

        n_plots = len(conditions)
        n_cols = min(2, n_plots)
        n_rows = (n_plots + n_cols - 1) // n_cols
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows))
        if n_plots == 1:
            axes = np.array([axes])
        axes = np.atleast_2d(axes)

        for idx, condition in enumerate(conditions):
            if condition not in self.calculator.standard_curves:
                continue
            curve = self.calculator.standard_curves[condition]
            row, col = divmod(idx, n_cols)
            ax = axes[row, col]

            # data
            ax.scatter(curve['x'], curve['y'], color=self.color_set[0], s=50,
                       label='Measured', zorder=3, edgecolor='white', linewidth=0.5)
            # translated note
            if curve['model_type'] == 'logarithmic':
                x_range = np.linspace(min(curve['x']), max(curve['x']), 100)
                y_fit = curve['intercept'] + curve['slope'] * np.log(x_range)
            else:
                x_range = np.linspace(min(curve['x']), max(curve['x']), 100)
                y_fit = curve['intercept'] + curve['slope'] * x_range
            ax.plot(x_range, y_fit, color=self.color_set[1], linewidth=2, label='Fitted')

            info = (f"Model: {curve['model_name']}\n"
                    f"R² = {curve['r_squared']:.4f}\n"
                    f"Slope = {curve['slope']:.4f}\n"
                    f"Intercept = {curve['intercept']:.4f}")
            ax.text(0.05, 0.95, info, transform=ax.transAxes, va='top', fontsize=fontsize-1,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))

            ax.set_xlabel('Degree of polymerization (n)', fontsize=fontsize)
            ax.set_ylabel('Retention time (min)', fontsize=fontsize)
            ax.set_title(f'PPG calibration curve - {condition}', fontsize=fontsize+1)
            ax.legend(fontsize=fontsize-1)
            ax.tick_params(labelsize=fontsize-1)
            ax.grid(True, alpha=0.3, linestyle='--')

        # translated note
        for idx in range(len(conditions), n_rows * n_cols):
            row, col = divmod(idx, n_cols)
            axes[row, col].set_visible(False)

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)
        self.figures['standard_curves'] = fig
        return fig

    def plot_residuals(self, conditions: List[str] = None, save_path: str = None,
                       fontsize: int = 12) -> plt.Figure:
        """"""
        if conditions is None:
            conditions = list(self.calculator.standard_curves.keys())
        if not conditions:
            print("Warning: No residual data to plot.")
            return None

        n_plots = len(conditions)
        n_cols = min(2, n_plots)
        n_rows = (n_plots + n_cols - 1) // n_cols
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows))
        if n_plots == 1:
            axes = np.array([axes])
        axes = np.atleast_2d(axes)

        for idx, condition in enumerate(conditions):
            if condition not in self.calculator.standard_curves:
                continue
            curve = self.calculator.standard_curves[condition]
            row, col = divmod(idx, n_cols)
            ax = axes[row, col]

            residuals = curve['residuals']
            predicted = curve['y_pred']
            ax.scatter(predicted, residuals, color=self.color_set[2], s=50, alpha=0.7,
                       edgecolor='white', linewidth=0.5)
            ax.axhline(y=0, color='black', linestyle='--', linewidth=1)

            stats_text = (f"Residual stats:\n"
                          f"Mean = {np.mean(residuals):.4f}\n"
                          f"Std = {np.std(residuals):.4f}\n"
                          f"Max abs = {np.max(np.abs(residuals)):.4f}")
            ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, va='top', fontsize=fontsize-1,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))

            ax.set_xlabel('Predicted retention time (min)', fontsize=fontsize)
            ax.set_ylabel('Residual (min)', fontsize=fontsize)
            ax.set_title(f'Residual plot - {condition}', fontsize=fontsize+1)
            ax.tick_params(labelsize=fontsize-1)
            ax.grid(True, alpha=0.3, linestyle='--')

        for idx in range(len(conditions), n_rows * n_cols):
            row, col = divmod(idx, n_cols)
            axes[row, col].set_visible(False)

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)
        self.figures['residuals'] = fig
        return fig

    def plot_ppg_index_distribution(self, condition: str, save_path: str = None,
                                     fontsize: int = 12) -> plt.Figure:
        """PPGindex (+)"""
        if condition not in self.calculator.ppg_indices:
            print(f"Warning: No PPG index data for condition {condition}.")
            return None

        all_indices = []
        for df in self.calculator.ppg_indices[condition]['indices'].values():
            if 'calculatePPGindex' in df.columns:
                all_indices.extend(df['calculatePPGindex'].dropna().tolist())
        if not all_indices:
            print("Warning: No valid PPG index values.")
            return None

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # translated note
        ax1.hist(all_indices, bins=30, color=self.color_set[0], edgecolor='white', alpha=0.7)
        ax1.axvline(np.mean(all_indices), color='red', linestyle='--', linewidth=2,
                    label=f'Mean: {np.mean(all_indices):.1f}')
        ax1.axvline(np.median(all_indices), color='green', linestyle='--', linewidth=2,
                    label=f'Median: {np.median(all_indices):.1f}')
        ax1.set_xlabel('PPG index', fontsize=fontsize)
        ax1.set_ylabel('Frequency', fontsize=fontsize)
        ax1.set_title(f'PPG index distribution - {condition}', fontsize=fontsize+1)
        ax1.legend(fontsize=fontsize-1)
        ax1.tick_params(labelsize=fontsize-1)
        ax1.grid(True, alpha=0.3, linestyle='--')

        # translated note
        ax2.boxplot(all_indices, vert=True, patch_artist=True,
                    boxprops=dict(facecolor=self.color_set[1], color='black'),
                    medianprops=dict(color='red', linewidth=2),
                    whiskerprops=dict(color='black'),
                    capprops=dict(color='black'))
        stats_text = (f"Statistics:\n"
                      f"N = {len(all_indices)}\n"
                      f"Mean = {np.mean(all_indices):.1f}\n"
                      f"Std = {np.std(all_indices):.1f}\n"
                      f"Min = {np.min(all_indices):.1f}\n"
                      f"Max = {np.max(all_indices):.1f}\n"
                      f"Median = {np.median(all_indices):.1f}")
        ax2.text(0.05, 0.95, stats_text, transform=ax2.transAxes, va='top', fontsize=fontsize-1,
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))
        ax2.set_ylabel('PPG index', fontsize=fontsize)
        ax2.set_title(f'PPG index boxplot - {condition}', fontsize=fontsize+1)
        ax2.set_xticks([1])
        ax2.set_xticklabels([condition], fontsize=fontsize-1)
        ax2.tick_params(labelsize=fontsize-1)
        ax2.grid(True, alpha=0.3, linestyle='--')

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)
        self.figures['index_distribution'] = fig
        return fig

    def plot_condition_comparison(self, conditions: List[str], save_path: str = None,
                                   fontsize: int = 12) -> plt.Figure:
        """conditionPPGindexcompare"""
        comparison_df = self.calculator.compare_conditions(conditions)
        if comparison_df.empty:
            print("Warning: No data for comparison.")
            return None

        fig, axes = plt.subplots(2, 2, figsize=(14, 12))
        axes = axes.flatten()

        # 1. (condition)
        if len(conditions) >= 2:
            ax = axes[0]
            cond1, cond2 = conditions[0], conditions[1]
            col1 = f"{cond1}_PPGindex"
            col2 = f"{cond2}_PPGindex"
            if col1 in comparison_df.columns and col2 in comparison_df.columns:
                data = comparison_df[[col1, col2]].dropna()
                if not data.empty:
                    ax.scatter(data[col1], data[col2], color=self.color_set[0], s=50,
                               alpha=0.7, edgecolor='white', linewidth=0.5)
                    min_val = min(data.min())
                    max_val = max(data.max())
                    ax.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.5, linewidth=1.5)
                    ax.set_xlabel(f'{cond1} PPG index', fontsize=fontsize)
                    ax.set_ylabel(f'{cond2} PPG index', fontsize=fontsize)
                    ax.set_title(f'{cond1} vs {cond2}', fontsize=fontsize+1)
                    ax.tick_params(labelsize=fontsize-1)
                    ax.grid(True, alpha=0.3, linestyle='--')

        # 2. compare
        ax = axes[1]
        box_data = []
        labels = []
        for cond in conditions:
            col = f"{cond}_PPGindex"
            if col in comparison_df.columns:
                data = comparison_df[col].dropna()
                if not data.empty:
                    box_data.append(data.values)
                    labels.append(cond)
        if box_data:
            bp = ax.boxplot(box_data, labels=labels, patch_artist=True)
            for patch, color in zip(bp['boxes'], self.color_set[:len(box_data)]):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)
            ax.set_ylabel('PPG index', fontsize=fontsize)
            ax.set_title('PPG index distribution', fontsize=fontsize+1)
            ax.tick_params(labelsize=fontsize-1)
            ax.grid(True, alpha=0.3, linestyle='--')

        # 3. (condition)
        ax = axes[2]
        if len(conditions) >= 2:
            error_df = self.calculator.calculate_conversion_error(conditions[0], conditions[1])
            if not error_df.empty and '(%)' in error_df.columns:
                errors = error_df['(%)'].dropna()
                if not errors.empty:
                    ax.hist(errors, bins=20, color=self.color_set[2], edgecolor='white', alpha=0.7)
                    ax.axvline(np.mean(errors), color='red', linestyle='--', linewidth=2,
                               label=f'Mean: {np.mean(errors):.2f}%')
                    ax.set_xlabel('Relative error (%)', fontsize=fontsize)
                    ax.set_ylabel('Frequency', fontsize=fontsize)
                    ax.set_title(f'Conversion error ({conditions[0]} → {conditions[1]})', fontsize=fontsize+1)
                    ax.legend(fontsize=fontsize-1)
                    ax.tick_params(labelsize=fontsize-1)
                    ax.grid(True, alpha=0.3, linestyle='--')

        # 4.
        ax = axes[3]
        corr_data = []
        corr_labels = []
        for cond in conditions:
            col = f"{cond}_PPGindex"
            if col in comparison_df.columns:
                corr_data.append(comparison_df[col])
                corr_labels.append(cond)
        if len(corr_data) > 1:
            corr_df = pd.DataFrame(corr_data).T
            corr_df.columns = corr_labels
            corr_matrix = corr_df.corr()
            im = ax.imshow(corr_matrix, cmap='RdYlBu', vmin=0, vmax=1, aspect='auto')
            for i in range(len(corr_matrix)):
                for j in range(len(corr_matrix)):
                    ax.text(j, i, f'{corr_matrix.iloc[i, j]:.2f}',
                            ha='center', va='center', fontsize=fontsize-2,
                            color='white' if corr_matrix.iloc[i, j] > 0.5 else 'black')
            ax.set_xticks(range(len(corr_labels)))
            ax.set_yticks(range(len(corr_labels)))
            ax.set_xticklabels(corr_labels, fontsize=fontsize-1, rotation=45, ha='right')
            ax.set_yticklabels(corr_labels, fontsize=fontsize-1)
            ax.set_title('Correlation matrix', fontsize=fontsize+1)
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

        # translated note
        for i in range(4):
            if not axes[i].has_data():
                axes[i].set_visible(False)

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)
        self.figures['condition_comparison'] = fig
        return fig

    def plot_conversion_analysis(self, from_condition: str, to_condition: str,
                                  save_path: str = None, fontsize: int = 12) -> plt.Figure:
        """conditionanalysis"""
        try:
            conversion_df, _ = self.calculator.convert_ppg_index_to_rt(from_condition, to_condition)
            if conversion_df.empty or '(min)' not in conversion_df.columns:
                print("Warning: No valid conversion data.")
                return None

            valid = conversion_df.dropna(subset=['(min)',
                                                 f'{to_condition}_RT',
                                                 f'{to_condition}_RT'])
            if valid.empty:
                return None

            fig, axes = plt.subplots(2, 3, figsize=(18, 12))
            axes = axes.flatten()

            predicted = valid[f'{to_condition}_RT']
            actual = valid[f'{to_condition}_RT']
            errors = valid['(min)']
            rel_errors = valid['(%)'] if '(%)' in valid.columns else None

            # 1. vs
            ax = axes[0]
            ax.scatter(actual, predicted, color=self.color_set[0], s=50, alpha=0.7,
                       edgecolor='white', linewidth=0.5)
            lims = [min(actual.min(), predicted.min()), max(actual.max(), predicted.max())]
            ax.plot(lims, lims, 'r--', alpha=0.5, linewidth=1.5, label='y=x')
            if len(predicted) > 1:
                slope, intercept, r_val, _, _ = stats.linregress(actual, predicted)
                x_line = np.linspace(lims[0], lims[1], 50)
                ax.plot(x_line, intercept + slope * x_line, color=self.color_set[1],
                        linestyle='-', linewidth=2, label=f'Fit: R²={r_val**2:.3f}')
            ax.set_xlabel(f'Actual RT in {to_condition} (min)', fontsize=fontsize)
            ax.set_ylabel(f'Predicted RT in {to_condition} (min)', fontsize=fontsize)
            ax.set_title(f'{from_condition} → {to_condition}: Predicted vs Actual', fontsize=fontsize+1)
            ax.legend(fontsize=fontsize-1)
            ax.tick_params(labelsize=fontsize-1)
            ax.grid(True, alpha=0.3, linestyle='--')

            # 2.
            ax = axes[1]
            ax.hist(errors, bins=20, color=self.color_set[2], edgecolor='white', alpha=0.7)
            ax.axvline(errors.mean(), color='red', linestyle='--', linewidth=2,
                       label=f'Mean: {errors.mean():.3f} min')
            stats_text = (f"Error stats:\n"
                          f"N = {len(errors)}\n"
                          f"Mean = {errors.mean():.3f}\n"
                          f"Std = {errors.std():.3f}\n"
                          f"Max = {errors.max():.3f}\n"
                          f"Median = {errors.median():.3f}")
            ax.text(0.65, 0.95, stats_text, transform=ax.transAxes, va='top', fontsize=fontsize-1,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))
            ax.set_xlabel('Absolute error (min)', fontsize=fontsize)
            ax.set_ylabel('Frequency', fontsize=fontsize)
            ax.set_title('Absolute error distribution', fontsize=fontsize+1)
            ax.legend(fontsize=fontsize-1)
            ax.tick_params(labelsize=fontsize-1)
            ax.grid(True, alpha=0.3, linestyle='--')

            # 3.
            ax = axes[2]
            if rel_errors is not None and not rel_errors.dropna().empty:
                rel = rel_errors.dropna()
                ax.hist(rel, bins=20, color=self.color_set[3], edgecolor='white', alpha=0.7)
                ax.axvline(rel.mean(), color='red', linestyle='--', linewidth=2,
                           label=f'Mean: {rel.mean():.2f}%')
                ax.set_xlabel('Relative error (%)', fontsize=fontsize)
                ax.set_ylabel('Frequency', fontsize=fontsize)
                ax.set_title('Relative error distribution', fontsize=fontsize+1)
                ax.legend(fontsize=fontsize-1)
                ax.tick_params(labelsize=fontsize-1)
                ax.grid(True, alpha=0.3, linestyle='--')
            else:
                ax.set_visible(False)

            # 4. vs RT
            ax = axes[3]
            ax.scatter(actual, errors, color=self.color_set[4], s=50, alpha=0.7,
                       edgecolor='white', linewidth=0.5)
            ax.axhline(errors.mean(), color='red', linestyle='--', linewidth=1.5,
                       label=f'Mean error: {errors.mean():.3f}')
            if len(actual) > 1:
                z = np.polyfit(actual, errors, 1)
                p = np.poly1d(z)
                ax.plot(actual, p(actual), color=self.color_set[1], linestyle='-', linewidth=2,
                        label='Trend')
            ax.set_xlabel(f'Actual RT in {to_condition} (min)', fontsize=fontsize)
            ax.set_ylabel('Absolute error (min)', fontsize=fontsize)
            ax.set_title('Error vs RT', fontsize=fontsize+1)
            ax.legend(fontsize=fontsize-1)
            ax.tick_params(labelsize=fontsize-1)
            ax.grid(True, alpha=0.3, linestyle='--')

            # 5. compound
            ax = axes[4]
            top_n = min(15, len(valid))
            top_errors = valid.nlargest(top_n, '(min)')
            y_pos = np.arange(top_n)
            ax.barh(y_pos, top_errors['(min)'], color=self.color_set[5], edgecolor='white')
            names = [name[:20] + '...' if len(name) > 20 else name
                     for name in top_errors['compound_name']]
            ax.set_yticks(y_pos)
            ax.set_yticklabels(names, fontsize=fontsize-2)
            ax.invert_yaxis()
            ax.set_xlabel('Absolute error (min)', fontsize=fontsize)
            ax.set_title(f'Top {top_n} largest errors', fontsize=fontsize+1)
            ax.tick_params(labelsize=fontsize-1)
            ax.grid(True, alpha=0.3, linestyle='--', axis='x')

            # 6. Bland-Altman
            ax = axes[5]
            mean_vals = (predicted + actual) / 2
            diff = predicted - actual
            ax.scatter(mean_vals, diff, color=self.color_set[0], s=50, alpha=0.7,
                       edgecolor='white', linewidth=0.5)
            mean_diff = diff.mean()
            std_diff = diff.std()
            ax.axhline(mean_diff, color='red', linestyle='-', linewidth=2, label=f'Mean diff: {mean_diff:.3f}')
            ax.axhline(mean_diff + 1.96*std_diff, color='red', linestyle='--', linewidth=1.5,
                       label=f'+1.96SD: {mean_diff+1.96*std_diff:.3f}')
            ax.axhline(mean_diff - 1.96*std_diff, color='red', linestyle='--', linewidth=1.5,
                       label=f'-1.96SD: {mean_diff-1.96*std_diff:.3f}')
            ax.axhline(0, color='gray', linestyle='-', linewidth=0.5)
            ax.set_xlabel('Mean RT (min)', fontsize=fontsize)
            ax.set_ylabel('Predicted - Actual (min)', fontsize=fontsize)
            ax.set_title('Bland-Altman plot', fontsize=fontsize+1)
            ax.legend(fontsize=fontsize-2)
            ax.tick_params(labelsize=fontsize-1)
            ax.grid(True, alpha=0.3, linestyle='--')

            plt.tight_layout()
            if save_path:
                plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)
            self.figures['conversion_analysis'] = fig
            return fig

        except Exception as e:
            print(f"Error plotting conversion analysis: {e}")
            return None

    def plot_multiple_conversion_comparison(self, conversions: List[Tuple[str, str]],
                                             save_path: str = None, fontsize: int = 12) -> plt.Figure:
        """compare"""
        try:
            all_errors = []
            labels = []
            for from_cond, to_cond in conversions:
                df, _ = self.calculator.convert_ppg_index_to_rt(from_cond, to_cond)
                if not df.empty and '(min)' in df.columns:
                    err = df['(min)'].dropna()
                    if not err.empty:
                        all_errors.append(err)
                        labels.append(f'{from_cond}→{to_cond}')
            if not all_errors:
                print("Warning: No valid conversion data for comparison.")
                return None

            fig, axes = plt.subplots(2, 2, figsize=(14, 12))
            axes = axes.flatten()

            # 1. compare
            ax = axes[0]
            bp = ax.boxplot(all_errors, labels=labels, patch_artist=True, showfliers=False)
            for patch, color in zip(bp['boxes'], self.color_set[:len(all_errors)]):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)
            ax.set_ylabel('Absolute error (10^-1min)', fontsize=fontsize)
            ax.set_title('Error distribution', fontsize=fontsize+1)
            ax.tick_params(labelsize=fontsize-1, rotation=45)
            ax.grid(True, alpha=0.3, linestyle='--')
            for i, err in enumerate(all_errors):
                ax.text(i+1, err.max()*1.05, f'', ha='center', fontsize=fontsize-2)

            # 2.
            ax = axes[1]
            means = [e.mean() for e in all_errors]
            stds = [e.std() for e in all_errors]
            x_pos = np.arange(len(labels))
            ax.bar(x_pos, means, yerr=stds, capsize=5, color=self.color_set[0], alpha=0.7,
                   edgecolor='white', linewidth=1)
            for i, (m, s) in enumerate(zip(means, stds)):
                ax.text(i, m + 0.02, f'{m:.3f}', ha='center', fontsize=fontsize-2)
            ax.set_xticks(x_pos)
            ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=fontsize-1)
            ax.set_ylabel('Mean absolute error (10^-1min)', fontsize=fontsize)
            ax.set_title('Mean error comparison', fontsize=fontsize+1)
            ax.tick_params(labelsize=fontsize-1)
            ax.grid(True, alpha=0.3, linestyle='--', axis='y')

            # 3. compare
            ax = axes[2]
            pass_rates = [len(e[e <= 0.5]) / len(e) * 100 for e in all_errors]
            ax.bar(x_pos, pass_rates, color=self.color_set[1], alpha=0.7,
                   edgecolor='white', linewidth=1)
            for i, pr in enumerate(pass_rates):
                ax.text(i, pr + 2, f'{pr:.1f}%', ha='center', fontsize=fontsize-2)
            ax.set_xticks(x_pos)
            ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=fontsize-1)
            ax.set_ylabel('Fail rate (%)', fontsize=fontsize)
            ax.set_title('Fail rate (threshold 0.5 min)', fontsize=fontsize+1)
            ax.set_ylim(0, 105)
            ax.tick_params(labelsize=fontsize-1)
            ax.grid(True, alpha=0.3, linestyle='--', axis='y')

            # 4.
            ax = axes[3]
            for err, label in zip(all_errors, labels):
                sorted_err = np.sort(err)
                y_vals = np.arange(1, len(sorted_err)+1) / len(sorted_err) * 100
                ax.plot(sorted_err, y_vals, marker='.', label=label, linewidth=2)
            ax.set_xlabel('Absolute error (S)', fontsize=fontsize)
            ax.set_ylabel('Cumulative percentage (%)', fontsize=fontsize)
            ax.set_title('Cumulative error distribution', fontsize=fontsize+1)
            ax.legend(fontsize=fontsize-2)
            ax.tick_params(labelsize=fontsize-1)
            ax.grid(True, alpha=0.3, linestyle='--')

            plt.tight_layout()
            if save_path:
                plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)
            self.figures['multiple_conversion_comparison'] = fig
            return fig

        except Exception as e:
            print(f"Error plotting multiple conversion comparison: {e}")
            return None

# =============================================================================
# PPGIndexAnalyzerGUI ()
# =============================================================================
class PPGIndexAnalyzerGUI:
    """PPG retention indicesanalysisGUI"""

    def __init__(self, root):
        self.root = root
        self.root.title("PPG Retention Index Analyzer - Journal Version")
        self.root.geometry("1300x950")

        try:
            self.root.iconbitmap(default='icon.ico')
        except:
            pass

        self.calculator = PPGIndexCalculator()
        self.visualizer = None
        self.processing_thread = None
        self.is_processing = False

        self.loaded_conditions = set()
        self.loaded_compound_datasets = set()

        self.setup_ui()
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)

    def setup_ui(self):
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        self.notebook = ttk.Notebook(main_frame)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        self.setup_data_tab()
        self.setup_analysis_tab()
        self.setup_conversion_tab()
        self.setup_visualization_tab()
        self.setup_results_tab()
        self.setup_log_tab()

        self.status_var = tk.StringVar(value="Ready")
        status_bar = ttk.Label(main_frame, textvariable=self.status_var,
                               relief=tk.SUNKEN, anchor=tk.W)
        status_bar.pack(fill=tk.X, pady=(5, 0))

    # ---------- dataload ----------
    def setup_data_tab(self):
        data_tab = ttk.Frame(self.notebook)
        self.notebook.add(data_tab, text="Data Loading")

        data_frame = ttk.LabelFrame(data_tab, text="Data Management", padding=15)
        data_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # PPG dataload
        ppg_frame = ttk.LabelFrame(data_frame, text="PPG Standard Data", padding=10)
        ppg_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(ppg_frame, text="Condition:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        self.condition_var = tk.StringVar(value="condition1")
        ttk.Entry(ppg_frame, textvariable=self.condition_var, width=20).grid(row=0, column=1, sticky=tk.W, padx=5, pady=5)

        ttk.Label(ppg_frame, text="PPG file:").grid(row=1, column=0, sticky=tk.W, padx=5, pady=5)
        self.ppg_file_var = tk.StringVar()
        ttk.Entry(ppg_frame, textvariable=self.ppg_file_var, width=60).grid(row=1, column=1, sticky=tk.W, padx=5, pady=5)
        ttk.Button(ppg_frame, text="Browse...", command=self.browse_ppg_file).grid(row=1, column=2, padx=5, pady=5)
        ttk.Button(ppg_frame, text="Load PPG Data", command=self.load_ppg_data).grid(row=2, column=0, columnspan=3, pady=10)

        # compound dataload
        comp_frame = ttk.LabelFrame(data_frame, text="Compound Data", padding=10)
        comp_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(comp_frame, text="Category:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        self.category_var = tk.StringVar(value="validation")
        ttk.Combobox(comp_frame, textvariable=self.category_var,
                     values=["validation", "smrt", "training", "test"],
                     width=15, state="readonly").grid(row=0, column=1, sticky=tk.W, padx=5, pady=5)

        ttk.Label(comp_frame, text="Condition:").grid(row=0, column=2, sticky=tk.W, padx=5, pady=5)
        self.compound_condition_var = tk.StringVar(value="condition1")
        ttk.Entry(comp_frame, textvariable=self.compound_condition_var, width=20).grid(row=0, column=3, sticky=tk.W, padx=5, pady=5)

        ttk.Label(comp_frame, text="Compound file:").grid(row=1, column=0, sticky=tk.W, padx=5, pady=5)
        self.compound_file_var = tk.StringVar()
        ttk.Entry(comp_frame, textvariable=self.compound_file_var, width=60).grid(row=1, column=1, columnspan=3, sticky=tk.W, padx=5, pady=5)
        ttk.Button(comp_frame, text="Browse...", command=self.browse_compound_file).grid(row=1, column=4, padx=5, pady=5)
        ttk.Button(comp_frame, text="Load Compound Data", command=self.load_compound_data).grid(row=2, column=0, columnspan=5, pady=10)

        # data
        overview_frame = ttk.LabelFrame(data_frame, text="Loaded Data Overview", padding=10)
        overview_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 10))

        columns = ("Type", "Condition", "Category", "Points", "Status")
        self.data_tree = ttk.Treeview(overview_frame, columns=columns, show="headings", height=8)
        for col in columns:
            self.data_tree.heading(col, text=col)
            self.data_tree.column(col, width=100)
        scrollbar = ttk.Scrollbar(overview_frame, orient=tk.VERTICAL, command=self.data_tree.yview)
        self.data_tree.configure(yscrollcommand=scrollbar.set)
        self.data_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        btn_frame = ttk.Frame(data_frame)
        btn_frame.pack(fill=tk.X, pady=(0, 10))
        ttk.Button(btn_frame, text="Clear All", command=self.clear_all_data).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Preview PPG", command=self.preview_ppg_data).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Preview Compounds", command=self.preview_compound_data).pack(side=tk.LEFT, padx=5)

        self.data_status_var = tk.StringVar(value="Waiting for data...")
        ttk.Label(data_frame, textvariable=self.data_status_var).pack(anchor=tk.W)

    # ---------- dataanalysis ----------
    def setup_analysis_tab(self):
        analysis_tab = ttk.Frame(self.notebook)
        self.notebook.add(analysis_tab, text="Analysis")

        analysis_frame = ttk.LabelFrame(analysis_tab, text="PPG Index Analysis", padding=15)
        analysis_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # standard curve
        curve_frame = ttk.LabelFrame(analysis_frame, text="Calibration Curve Fitting", padding=10)
        curve_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(curve_frame, text="Condition:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        self.curve_condition_var = tk.StringVar()
        self.curve_condition_combo = ttk.Combobox(curve_frame, textvariable=self.curve_condition_var,
                                                  width=25, state="readonly")
        self.curve_condition_combo.grid(row=0, column=1, sticky=tk.W, padx=5, pady=5)

        ttk.Label(curve_frame, text="Model:").grid(row=0, column=2, sticky=tk.W, padx=5, pady=5)
        self.model_type_var = tk.StringVar(value="linear")
        ttk.Combobox(curve_frame, textvariable=self.model_type_var,
                     values=["linear", "logarithmic"], width=15, state="readonly").grid(row=0, column=3, sticky=tk.W, padx=5, pady=5)

        ttk.Button(curve_frame, text="Fit Curve", command=self.fit_standard_curve).grid(row=0, column=4, padx=20, pady=5)

        # PPGindexcalculate
        calc_frame = ttk.LabelFrame(analysis_frame, text="PPG Index Calculation", padding=10)
        calc_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(calc_frame, text="Condition:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        self.calc_condition_var = tk.StringVar()
        self.calc_condition_combo = ttk.Combobox(calc_frame, textvariable=self.calc_condition_var,
                                                 width=25, state="readonly")
        self.calc_condition_combo.grid(row=0, column=1, sticky=tk.W, padx=5, pady=5)

        ttk.Label(calc_frame, text="Method:").grid(row=0, column=2, sticky=tk.W, padx=5, pady=5)
        self.calc_method_var = tk.StringVar(value="interpolation")
        ttk.Combobox(calc_frame, textvariable=self.calc_method_var,
                     values=["interpolation", "regression"], width=15, state="readonly").grid(row=0, column=3, sticky=tk.W, padx=5, pady=5)

        ttk.Button(calc_frame, text="Calculate PPG Indices", command=self.calculate_ppg_index).grid(row=0, column=4, padx=20, pady=5)

        # conditioncompare
        compare_frame = ttk.LabelFrame(analysis_frame, text="Condition Comparison", padding=10)
        compare_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(compare_frame, text="Select conditions:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        self.condition_checkboxes = {}
        self.condition_checkboxes_frame = ttk.Frame(compare_frame)
        self.condition_checkboxes_frame.grid(row=1, column=0, columnspan=5, sticky=tk.W, padx=5, pady=5)

        ttk.Button(compare_frame, text="Compare Selected", command=self.compare_conditions).grid(row=2, column=0, sticky=tk.W, padx=5, pady=10)

        # results
        results_frame = ttk.LabelFrame(analysis_frame, text="Results", padding=10)
        results_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 10))

        self.analysis_text = scrolledtext.ScrolledText(results_frame, width=80, height=15,
                                                        wrap=tk.WORD, font=("Consolas", 10))
        self.analysis_text.pack(fill=tk.BOTH, expand=True)
        self.analysis_text.tag_config("INFO", foreground="black")
        self.analysis_text.tag_config("SUCCESS", foreground="green")
        self.analysis_text.tag_config("WARNING", foreground="orange")
        self.analysis_text.tag_config("ERROR", foreground="red")

        btn_frame = ttk.Frame(analysis_frame)
        btn_frame.pack(fill=tk.X, pady=(0, 10))
        ttk.Button(btn_frame, text="Generate Report", command=self.generate_report).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Clear Results", command=self.clear_analysis_text).pack(side=tk.LEFT, padx=5)

        self.analysis_status_var = tk.StringVar(value="Ready")
        ttk.Label(analysis_frame, textvariable=self.analysis_status_var).pack(anchor=tk.W)

    # ---------- condition ----------
    def setup_conversion_tab(self):
        conversion_tab = ttk.Frame(self.notebook)
        self.notebook.add(conversion_tab, text="Cross‑condition Conversion")

        conversion_frame = ttk.LabelFrame(conversion_tab, text="PPG Index Conversion & Validation", padding=15)
        conversion_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # translated note
        settings_frame = ttk.LabelFrame(conversion_frame, text="Settings", padding=10)
        settings_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(settings_frame, text="Source condition:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        self.from_condition_var = tk.StringVar()
        self.from_condition_combo = ttk.Combobox(settings_frame, textvariable=self.from_condition_var,
                                                 width=25, state="readonly")
        self.from_condition_combo.grid(row=0, column=1, sticky=tk.W, padx=5, pady=5)

        ttk.Label(settings_frame, text="Target condition:").grid(row=0, column=2, sticky=tk.W, padx=5, pady=5)
        self.to_condition_var = tk.StringVar()
        self.to_condition_combo = ttk.Combobox(settings_frame, textvariable=self.to_condition_var,
                                               width=25, state="readonly")
        self.to_condition_combo.grid(row=0, column=3, sticky=tk.W, padx=5, pady=5)

        ttk.Label(settings_frame, text="Error threshold (min):").grid(row=1, column=0, sticky=tk.W, padx=5, pady=5)
        self.threshold_var = tk.DoubleVar(value=0.5)
        ttk.Entry(settings_frame, textvariable=self.threshold_var, width=10).grid(row=1, column=1, sticky=tk.W, padx=5, pady=5)

        ttk.Button(settings_frame, text="Run Conversion", command=self.perform_conversion_analysis).grid(
            row=0, column=4, rowspan=2, padx=20, pady=5)

        # compare
        multi_frame = ttk.LabelFrame(settings_frame, text="Multiple Schemes Comparison", padding=5)
        multi_frame.grid(row=2, column=0, columnspan=5, sticky=tk.W, padx=5, pady=10)

        ttk.Label(multi_frame, text="Select conversion schemes:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=2)
        self.multi_conversion_frame = ttk.Frame(multi_frame)
        self.multi_conversion_frame.grid(row=1, column=0, sticky=tk.W, padx=5, pady=2)
        ttk.Button(multi_frame, text="Compare", command=self.compare_multiple_conversions).grid(row=1, column=1, padx=10, pady=2)

        # results
        results_frame = ttk.LabelFrame(conversion_frame, text="Conversion Results", padding=10)
        results_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 15))

        columns = ("Compound", "Source PPG", "Predicted RT", "Actual RT", "Abs Error", "Rel Error(%)")
        self.conversion_tree = ttk.Treeview(results_frame, columns=columns, show="headings", height=10)
        for col in columns:
            self.conversion_tree.heading(col, text=col)
            width = 150 if col == "Compound" else 100
            self.conversion_tree.column(col, width=width)
        scrollbar = ttk.Scrollbar(results_frame, orient=tk.VERTICAL, command=self.conversion_tree.yview)
        self.conversion_tree.configure(yscrollcommand=scrollbar.set)
        self.conversion_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        stats_frame = ttk.LabelFrame(conversion_frame, text="Statistics", padding=10)
        stats_frame.pack(fill=tk.X, pady=(0, 10))
        self.conversion_stats_text = scrolledtext.ScrolledText(stats_frame, width=80, height=6,
                                                                wrap=tk.WORD, font=("Consolas", 9))
        self.conversion_stats_text.pack(fill=tk.BOTH, expand=True)

        btn_frame = ttk.Frame(conversion_frame)
        btn_frame.pack(fill=tk.X, pady=(0, 10))
        ttk.Button(btn_frame, text="Export Results", command=self.export_conversion_results).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Clear", command=self.clear_conversion_results).pack(side=tk.LEFT, padx=5)

        self.conversion_status_var = tk.StringVar(value="Ready")
        ttk.Label(conversion_frame, textvariable=self.conversion_status_var).pack(anchor=tk.W)

    # ---------- ()----------
    def setup_visualization_tab(self):
        viz_tab = ttk.Frame(self.notebook)
        self.notebook.add(viz_tab, text="Visualization")

        viz_frame = ttk.LabelFrame(viz_tab, text="Plot Generation", padding=15)
        viz_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        options_frame = ttk.LabelFrame(viz_frame, text="Options", padding=10)
        options_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(options_frame, text="Plot type:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        self.viz_type_var = tk.StringVar(value="standard_curve")
        viz_type_combo = ttk.Combobox(options_frame, textvariable=self.viz_type_var,
                                      values=["standard_curve", "residuals", "index_distribution",
                                              "condition_comparison", "conversion_analysis",
                                              "multiple_conversion_comparison"],
                                      width=25, state="readonly")
        viz_type_combo.grid(row=0, column=1, sticky=tk.W, padx=5, pady=5)
        viz_type_combo.bind("<<ComboboxSelected>>", self.on_viz_type_change)

        ttk.Label(options_frame, text="Font size:").grid(row=0, column=2, sticky=tk.W, padx=5, pady=5)
        self.fontsize_var = tk.IntVar(value=12)
        ttk.Spinbox(options_frame, from_=8, to=24, textvariable=self.fontsize_var, width=5).grid(row=0, column=3, sticky=tk.W, padx=5, pady=5)

        ttk.Label(options_frame, text="Condition:").grid(row=0, column=4, sticky=tk.W, padx=5, pady=5)
        self.viz_condition_var = tk.StringVar()
        self.viz_condition_combo = ttk.Combobox(options_frame, textvariable=self.viz_condition_var,
                                                width=25, state="readonly")
        self.viz_condition_combo.grid(row=0, column=5, sticky=tk.W, padx=5, pady=5)

        self.viz_conditions_frame = ttk.Frame(options_frame)
        self.viz_conditions_frame.grid(row=1, column=0, columnspan=6, sticky=tk.W, padx=5, pady=5)
        self.viz_conditions_frame.grid_remove()

        self.conversion_viz_frame = ttk.Frame(options_frame)
        self.conversion_viz_frame.grid(row=1, column=0, columnspan=6, sticky=tk.W, padx=5, pady=5)
        self.conversion_viz_frame.grid_remove()

        ttk.Button(options_frame, text="Generate Plot", command=self.generate_visualization).grid(row=0, column=6, padx=20, pady=5)

        display_frame = ttk.LabelFrame(viz_frame, text="Plot Display", padding=10)
        display_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 10))

        self.viz_placeholder = ttk.Label(display_frame, text="Plot will appear here",
                                          font=("Arial", 14), foreground="gray")
        self.viz_placeholder.pack(expand=True)

        self.figure_canvas = None
        self.figure_toolbar = None
        self.current_figure = None

        btn_frame = ttk.Frame(viz_frame)
        btn_frame.pack(fill=tk.X, pady=(0, 10))
        ttk.Button(btn_frame, text="Save Plot", command=self.save_figure).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Clear Plot", command=self.clear_figure).pack(side=tk.LEFT, padx=5)

        self.viz_status_var = tk.StringVar(value="Ready")
        ttk.Label(viz_frame, textvariable=self.viz_status_var).pack(anchor=tk.W)

    # ---------- resultsoutput ----------
    def setup_results_tab(self):
        results_tab = ttk.Frame(self.notebook)
        self.notebook.add(results_tab, text="Results")

        results_frame = ttk.LabelFrame(results_tab, text="Output & Export", padding=15)
        results_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        output_frame = ttk.LabelFrame(results_frame, text="Output Settings", padding=10)
        output_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(output_frame, text="Output directory:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        self.output_dir_var = tk.StringVar(value=os.getcwd())
        ttk.Entry(output_frame, textvariable=self.output_dir_var, width=60).grid(row=0, column=1, sticky=tk.W, padx=5, pady=5)
        ttk.Button(output_frame, text="Browse...", command=self.browse_output_dir).grid(row=0, column=2, padx=5, pady=5)

        preview_frame = ttk.LabelFrame(results_frame, text="Report Preview", padding=10)
        preview_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 15))

        self.results_text = scrolledtext.ScrolledText(preview_frame, width=80, height=15,
                                                       wrap=tk.WORD, font=("Consolas", 10))
        self.results_text.pack(fill=tk.BOTH, expand=True)

        control_frame = ttk.Frame(results_frame)
        control_frame.pack(fill=tk.X, pady=(0, 10))
        ttk.Button(control_frame, text="Save All Results", command=self.save_all_results).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_frame, text="Preview Report", command=self.preview_report).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_frame, text="Open Output Folder", command=self.open_output_dir).pack(side=tk.LEFT, padx=5)

        self.output_status_var = tk.StringVar(value="Ready")
        ttk.Label(results_frame, textvariable=self.output_status_var).pack(anchor=tk.W)

    # ---------- ----------
    def setup_log_tab(self):
        log_tab = ttk.Frame(self.notebook)
        self.notebook.add(log_tab, text="Log")

        log_frame = ttk.LabelFrame(log_tab, text="Program Log", padding=15)
        log_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        self.log_text = scrolledtext.ScrolledText(log_frame, width=100, height=25,
                                                   wrap=tk.WORD, font=("Consolas", 10))
        self.log_text.pack(fill=tk.BOTH, expand=True)
        self.log_text.tag_config("INFO", foreground="black")
        self.log_text.tag_config("SUCCESS", foreground="green")
        self.log_text.tag_config("WARNING", foreground="orange")
        self.log_text.tag_config("ERROR", foreground="red")

        btn_frame = ttk.Frame(log_frame)
        btn_frame.pack(fill=tk.X, pady=(10, 0))
        ttk.Button(btn_frame, text="Clear Log", command=self.clear_log).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Save Log", command=self.save_log).pack(side=tk.LEFT, padx=5)

    # ---------- method ----------
    def browse_ppg_file(self):
        file_path = filedialog.askopenfilename(title="Select PPG data file",
                                               filetypes=[("Data files", "*.csv *.xlsx *.xls"),
                                                          ("CSV", "*.csv"), ("Excel", "*.xlsx *.xls")])
        if file_path:
            self.ppg_file_var.set(file_path)

    def browse_compound_file(self):
        file_path = filedialog.askopenfilename(title="Select compound data file",
                                               filetypes=[("Data files", "*.csv *.xlsx *.xls"),
                                                          ("CSV", "*.csv"), ("Excel", "*.xlsx *.xls")])
        if file_path:
            self.compound_file_var.set(file_path)

    def browse_output_dir(self):
        dir_path = filedialog.askdirectory(title="Select output directory")
        if dir_path:
            self.output_dir_var.set(dir_path)

    def load_ppg_data(self):
        file = self.ppg_file_var.get().strip()
        cond = self.condition_var.get().strip()
        if not file:
            messagebox.showwarning("Warning", "Please select a PPG file.")
            return
        if not cond:
            messagebox.showwarning("Warning", "Please enter a condition name.")
            return
        self.log_message(f"Loading PPG data: {file} (condition: {cond})", "INFO")
        self.update_status("Loading PPG data...")
        threading.Thread(target=self._load_ppg_data_thread, args=(file, cond), daemon=True).start()

    def _load_ppg_data_thread(self, file, cond):
        try:
            success, msg = self.calculator.load_ppg_data(file, cond)
            if success:
                self.root.after(0, lambda: self.log_message(f"✓ {msg}", "SUCCESS"))
                self.root.after(0, self.update_data_tree)
                self.root.after(0, self.update_condition_comboboxes)
                self.loaded_conditions.add(cond)
            else:
                self.root.after(0, lambda: self.log_message(f"✗ {msg}", "ERROR"))
        except Exception as e:
            self.root.after(0, lambda: self.log_message(f"✗ PPG load failed: {str(e)}", "ERROR"))
        self.root.after(0, lambda: self.update_status("Ready"))

    def load_compound_data(self):
        file = self.compound_file_var.get().strip()
        cat = self.category_var.get().strip()
        cond = self.compound_condition_var.get().strip()
        if not file:
            messagebox.showwarning("Warning", "Please select a compound file.")
            return
        if not cat:
            messagebox.showwarning("Warning", "Please select a category.")
            return
        if not cond:
            messagebox.showwarning("Warning", "Please enter a condition name.")
            return
        self.log_message(f"Loading compound data: {file} (category: {cat}, condition: {cond})", "INFO")
        self.update_status("Loading compound data...")
        threading.Thread(target=self._load_compound_data_thread, args=(file, cat, cond), daemon=True).start()

    def _load_compound_data_thread(self, file, cat, cond):
        try:
            success, msg = self.calculator.load_compound_data(file, cat, cond)
            if success:
                self.root.after(0, lambda: self.log_message(f"✓ {msg}", "SUCCESS"))
                self.root.after(0, self.update_data_tree)
                self.root.after(0, self.update_condition_comboboxes)
                self.loaded_compound_datasets.add(f"{cat}_{cond}")
            else:
                self.root.after(0, lambda: self.log_message(f"✗ {msg}", "ERROR"))
        except Exception as e:
            self.root.after(0, lambda: self.log_message(f"✗ Compound load failed: {str(e)}", "ERROR"))
        self.root.after(0, lambda: self.update_status("Ready"))

    def update_data_tree(self):
        for item in self.data_tree.get_children():
            self.data_tree.delete(item)
        for cond, df in self.calculator.ppg_data.items():
            self.data_tree.insert("", tk.END, values=("PPG", cond, "", len(df), "Loaded"))
        for key, df in self.calculator.compound_data.items():
            parts = key.split('_')
            cat, cond = parts[0], parts[1] if len(parts) > 1 else ""
            self.data_tree.insert("", tk.END, values=("Compound", cond, cat, len(df), "Loaded"))
        self.data_status_var.set(f"PPG: {len(self.calculator.ppg_data)} sets, Compounds: {len(self.calculator.compound_data)} sets")

    def update_condition_comboboxes(self):
        conditions = list(self.calculator.ppg_data.keys())
        self.curve_condition_combo['values'] = conditions
        self.calc_condition_combo['values'] = conditions
        self.viz_condition_combo['values'] = conditions
        self.from_condition_combo['values'] = conditions
        self.to_condition_combo['values'] = conditions
        if conditions:
            if not self.curve_condition_var.get():
                self.curve_condition_var.set(conditions[0])
            if not self.calc_condition_var.get():
                self.calc_condition_var.set(conditions[0])
            if not self.viz_condition_var.get():
                self.viz_condition_var.set(conditions[0])
            if not self.from_condition_var.get():
                self.from_condition_var.set(conditions[0])
            if not self.to_condition_var.get():
                self.to_condition_var.set(conditions[-1] if len(conditions) > 1 else conditions[0])
        self.update_condition_checkboxes(conditions)
        self.update_multi_conversion_checkboxes(conditions)

    def update_condition_checkboxes(self, conditions):
        for w in self.condition_checkboxes_frame.winfo_children():
            w.destroy()
        self.condition_checkboxes.clear()
        for i, cond in enumerate(conditions):
            var = tk.BooleanVar(value=(i < 2))
            cb = ttk.Checkbutton(self.condition_checkboxes_frame, text=cond, variable=var)
            cb.grid(row=i//4, column=i%4, sticky=tk.W, padx=5, pady=2)
            self.condition_checkboxes[cond] = var

    def update_multi_conversion_checkboxes(self, conditions):
        for w in self.multi_conversion_frame.winfo_children():
            w.destroy()
        self.conversion_scheme_vars = []
        if len(conditions) >= 2:
            for i in range(len(conditions)):
                for j in range(len(conditions)):
                    if i != j:
                        from_cond = conditions[i]
                        to_cond = conditions[j]
                        var = tk.BooleanVar(value=(i==0 and j==1))
                        cb = ttk.Checkbutton(self.multi_conversion_frame,
                                             text=f"{from_cond}→{to_cond}", variable=var)
                        row = (i*len(conditions)+j) // 4
                        col = (i*len(conditions)+j) % 4
                        cb.grid(row=row, column=col, sticky=tk.W, padx=2, pady=1)
                        self.conversion_scheme_vars.append((var, from_cond, to_cond))
        else:
            ttk.Label(self.multi_conversion_frame, text="Need at least 2 conditions",
                      foreground="gray").grid(row=0, column=0, sticky=tk.W)

    def clear_all_data(self):
        if messagebox.askyesno("Confirm", "Clear all loaded data?"):
            self.calculator = PPGIndexCalculator()
            self.visualizer = None
            self.loaded_conditions.clear()
            self.loaded_compound_datasets.clear()
            self.update_data_tree()
            self.update_condition_comboboxes()
            self.clear_conversion_results()
            self.log_message("All data cleared.", "INFO")

    def preview_ppg_data(self):
        if not self.calculator.ppg_data:
            messagebox.showinfo("Info", "No PPG data loaded.")
            return
        win = Toplevel(self.root)
        win.title("PPG Data Preview")
        win.geometry("800x600")
        notebook = ttk.Notebook(win)
        notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        for cond, df in self.calculator.ppg_data.items():
            frame = ttk.Frame(notebook)
            notebook.add(frame, text=cond)
            tree_frame = ttk.Frame(frame)
            tree_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
            sb_y = ttk.Scrollbar(tree_frame)
            sb_y.pack(side=tk.RIGHT, fill=tk.Y)
            sb_x = ttk.Scrollbar(tree_frame, orient=tk.HORIZONTAL)
            sb_x.pack(side=tk.BOTTOM, fill=tk.X)
            tree = ttk.Treeview(tree_frame, yscrollcommand=sb_y.set, xscrollcommand=sb_x.set)
            tree.pack(fill=tk.BOTH, expand=True)
            sb_y.config(command=tree.yview)
            sb_x.config(command=tree.xview)
            tree["columns"] = list(df.columns)
            tree["show"] = "headings"
            for col in df.columns:
                tree.heading(col, text=col)
                tree.column(col, width=100)
            for _, row in df.iterrows():
                tree.insert("", tk.END, values=list(row))
            info = ttk.Frame(frame)
            info.pack(fill=tk.X, padx=5, pady=5)
            ttk.Label(info, text=f"Points: {len(df)}").pack(side=tk.LEFT, padx=10)
            if 'degree_of_polymerization' in df.columns and 'retention_time' in df.columns:
                ttk.Label(info, text=f"n range: {df['degree_of_polymerization'].min()} - {df['degree_of_polymerization'].max()}").pack(side=tk.LEFT, padx=10)
                ttk.Label(info, text=f"RT range: {df['retention_time'].min():.2f} - {df['retention_time'].max():.2f}").pack(side=tk.LEFT, padx=10)

    def preview_compound_data(self):
        if not self.calculator.compound_data:
            messagebox.showinfo("Info", "No compound data loaded.")
            return
        win = Toplevel(self.root)
        win.title("Compound Data Preview")
        win.geometry("900x700")
        notebook = ttk.Notebook(win)
        notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        for key, df in self.calculator.compound_data.items():
            frame = ttk.Frame(notebook)
            notebook.add(frame, text=key)
            tree_frame = ttk.Frame(frame)
            tree_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
            sb_y = ttk.Scrollbar(tree_frame)
            sb_y.pack(side=tk.RIGHT, fill=tk.Y)
            sb_x = ttk.Scrollbar(tree_frame, orient=tk.HORIZONTAL)
            sb_x.pack(side=tk.BOTTOM, fill=tk.X)
            tree = ttk.Treeview(tree_frame, yscrollcommand=sb_y.set, xscrollcommand=sb_x.set)
            tree.pack(fill=tk.BOTH, expand=True)
            sb_y.config(command=tree.yview)
            sb_x.config(command=tree.xview)
            tree["columns"] = list(df.columns)
            tree["show"] = "headings"
            for col in df.columns:
                tree.heading(col, text=col)
                width = 120 if len(col) > 10 else 100
                tree.column(col, width=width)
            for _, row in df.iterrows():
                tree.insert("", tk.END, values=list(row))
            info = ttk.Frame(frame)
            info.pack(fill=tk.X, padx=5, pady=5)
            ttk.Label(info, text=f"Compounds: {len(df)}").pack(side=tk.LEFT, padx=10)
            if 'retention_time' in df.columns:
                ttk.Label(info, text=f"RT range: {df['retention_time'].min():.2f} - {df['retention_time'].max():.2f}").pack(side=tk.LEFT, padx=10)

    def fit_standard_curve(self):
        cond = self.curve_condition_var.get()
        model = self.model_type_var.get()
        if not cond:
            messagebox.showwarning("Warning", "Select a condition.")
            return
        self.analysis_message(f"Fitting curve for {cond} ({model})...", "INFO")
        self.update_status("Fitting curve...")
        threading.Thread(target=self._fit_standard_curve_thread, args=(cond, model), daemon=True).start()

    def _fit_standard_curve_thread(self, cond, model):
        try:
            success, msg = self.calculator.fit_standard_curve(cond, model)
            if success:
                curve = self.calculator.standard_curves[cond]
                self.root.after(0, lambda: self.analysis_message(f"✓ {msg}", "SUCCESS"))
                result = f"Fitting results - {cond}:\n"
                result += f"  Model: {curve['model_name']}\n"
                result += f"  R²: {curve['r_squared']:.6f}\n"
                result += f"  Slope: {curve['slope']:.6f}\n"
                result += f"  Intercept: {curve['intercept']:.6f}\n"
                result += f"  Std error: {curve['std_err']:.6f}\n"
                result += f"  Points: {curve['n_points']}\n"
                self.root.after(0, lambda: self.analysis_message(result, "INFO"))
                self.root.after(0, lambda: setattr(self, 'visualizer', PPGVisualizer(self.calculator)))
            else:
                self.root.after(0, lambda: self.analysis_message(f"✗ {msg}", "ERROR"))
        except Exception as e:
            self.root.after(0, lambda: self.analysis_message(f"✗ Fit failed: {str(e)}", "ERROR"))
        self.root.after(0, lambda: self.update_status("Ready"))

    def calculate_ppg_index(self):
        cond = self.calc_condition_var.get()
        method = self.calc_method_var.get()
        if not cond:
            messagebox.showwarning("Warning", "Select a condition.")
            return
        self.analysis_message(f"Calculating PPG indices for {cond} ({method})...", "INFO")
        self.update_status("Calculating...")
        threading.Thread(target=self._calculate_ppg_index_thread, args=(cond, method), daemon=True).start()

    def _calculate_ppg_index_thread(self, cond, method):
        try:
            success, msg = self.calculator.calculate_ppg_index(cond, method)
            if success and cond in self.calculator.ppg_indices:
                all_idx = []
                for df in self.calculator.ppg_indices[cond]['indices'].values():
                    if 'calculatePPGindex' in df.columns:
                        all_idx.extend(df['calculatePPGindex'].dropna())
                if all_idx:
                    arr = np.array(all_idx)
                    result = f"PPG index statistics - {cond}:\n"
                    result += f"  Method: {method}\n"
                    result += f"  N = {len(arr)}\n"
                    result += f"  Mean = {np.mean(arr):.2f}\n"
                    result += f"  Std = {np.std(arr):.2f}\n"
                    result += f"  Min = {np.min(arr):.2f}\n"
                    result += f"  Max = {np.max(arr):.2f}\n"
                    result += f"  Median = {np.median(arr):.2f}\n"
                    self.root.after(0, lambda: self.analysis_message(result, "INFO"))
                self.root.after(0, lambda: self.analysis_message(f"✓ {msg}", "SUCCESS"))
            else:
                self.root.after(0, lambda: self.analysis_message(f"✗ {msg}", "ERROR"))
        except Exception as e:
            self.root.after(0, lambda: self.analysis_message(f"✗ Calculation failed: {str(e)}", "ERROR"))
        self.root.after(0, lambda: self.update_status("Ready"))

    def compare_conditions(self):
        selected = [c for c, var in self.condition_checkboxes.items() if var.get()]
        if len(selected) < 2:
            messagebox.showwarning("Warning", "Select at least two conditions.")
            return
        self.analysis_message(f"Comparing: {', '.join(selected)}", "INFO")
        self.update_status("Comparing...")
        threading.Thread(target=self._compare_conditions_thread, args=(selected,), daemon=True).start()

    def _compare_conditions_thread(self, conditions):
        try:
            df = self.calculator.compare_conditions(conditions)
            if not df.empty:
                self.root.after(0, lambda: self.analysis_message(f"✓ Comparison done, {len(df)} compounds.", "SUCCESS"))
                result = f"Comparison results:\n  Conditions: {', '.join(conditions)}\n  Total compounds: {len(df)}\n\n"
                for cond in conditions:
                    col = f"{cond}_PPGindex"
                    if col in df.columns:
                        data = df[col].dropna()
                        if len(data) > 0:
                            result += f"  {cond}:\n    N={len(data)}, mean={np.mean(data):.2f}, sd={np.std(data):.2f}\n"
                if len(conditions) >= 2:
                    corr = {}
                    for i in range(len(conditions)):
                        for j in range(i+1, len(conditions)):
                            col1 = f"{conditions[i]}_PPGindex"
                            col2 = f"{conditions[j]}_PPGindex"
                            if col1 in df.columns and col2 in df.columns:
                                pair = df[[col1, col2]].dropna()
                                if not pair.empty:
                                    r = pair[col1].corr(pair[col2])
                                    corr[f"{conditions[i]} vs {conditions[j]}"] = r
                    if corr:
                        result += "\n  Correlations:\n"
                        for k, v in corr.items():
                            result += f"    {k}: {v:.4f}\n"
                self.root.after(0, lambda: self.analysis_message(result, "INFO"))
            else:
                self.root.after(0, lambda: self.analysis_message("✗ No data for comparison.", "WARNING"))
        except Exception as e:
            self.root.after(0, lambda: self.analysis_message(f"✗ Comparison failed: {str(e)}", "ERROR"))
        self.root.after(0, lambda: self.update_status("Ready"))

    def perform_conversion_analysis(self):
        from_cond = self.from_condition_var.get()
        to_cond = self.to_condition_var.get()
        thresh = self.threshold_var.get()
        if not from_cond or not to_cond or from_cond == to_cond:
            messagebox.showwarning("Warning", "Select different source and target conditions.")
            return
        self.conversion_status_var.set(f"Converting {from_cond} → {to_cond}...")
        self.update_status("Running conversion...")
        threading.Thread(target=self._perform_conversion_analysis_thread, args=(from_cond, to_cond, thresh), daemon=True).start()

    def _perform_conversion_analysis_thread(self, from_cond, to_cond, thresh):
        try:
            res = self.calculator.cross_condition_analysis(from_cond, to_cond, thresh)
            if "error" in res:
                self.root.after(0, lambda: self.conversion_status_var.set("Conversion failed"))
                self.root.after(0, lambda: self.log_message(f"✗ {res['error']}", "ERROR"))
                return
            self.root.after(0, lambda: self.display_conversion_results(res))
            self.root.after(0, lambda: self.log_message(f"✓ Conversion {from_cond}→{to_cond} complete. Valid: {res['']}, Mean error: {res[''].get('',0):.3f} min", "SUCCESS"))
            self.root.after(0, lambda: self.conversion_status_var.set("Conversion done"))
        except Exception as e:
            self.root.after(0, lambda: self.log_message(f"✗ Conversion failed: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.conversion_status_var.set("Failed"))
        self.root.after(0, lambda: self.update_status("Ready"))

    def display_conversion_results(self, res):
        for item in self.conversion_tree.get_children():
            self.conversion_tree.delete(item)
        self.conversion_stats_text.delete(1.0, tk.END)

        for d in res.get('data', []):
            vals = (
                d['compound_name'],
                f"{d.get(f'{res["condition"]}_PPGindex',''):.2f}" if isinstance(d.get(f'{res["condition"]}_PPGindex'), (int,float)) else "",
                f"{d.get(f'{res["condition"]}_RT',''):.3f}" if isinstance(d.get(f'{res["condition"]}_RT'), (int,float)) else "",
                f"{d.get(f'{res["condition"]}_RT',''):.3f}" if isinstance(d.get(f'{res["condition"]}_RT'), (int,float)) else "",
                f"{d.get('(min)',''):.3f}" if isinstance(d.get('(min)'), (int,float)) else "",
                f"{d.get('(%)',''):.2f}" if isinstance(d.get('(%)'), (int,float)) else ""
            )
            self.conversion_tree.insert("", tk.END, values=vals)

        stats = f"Conversion analysis: {res['condition']} → {res['condition']}\n"
        stats += "="*50 + "\n\n"
        stats += f"Time: {res['']}\n"
        stats += f"Total compounds: {res['compound']}\n"
        stats += f"Valid conversions: {res['']}\n\n"
        if isinstance(res[''], dict):
            stats += "Error statistics:\n"
            for k, v in res[''].items():
                stats += f"  {k}: {v:.3f}\n" if isinstance(v, float) else f"  {k}: {v}\n"
        if '' in res:
            stats += "\nError distribution:\n"
            total = res['']
            for bin_name, cnt in res[''].items():
                pct = (cnt/total*100) if total>0 else 0
                stats += f"  {bin_name}: {cnt} ({pct:.1f}%)\n"
        if 'analysis' in res:
            pa = res['analysis']
            stats += f"\nPass rate (threshold {pa['(min)']} min): {pa['']}/{pa['']} = {pa['(%)']:.1f}%\n"
        self.conversion_stats_text.insert(tk.END, stats, "INFO")

    def compare_multiple_conversions(self):
        selected = [(f, t) for var, f, t in self.conversion_scheme_vars if var.get()]
        if len(selected) < 2:
            messagebox.showwarning("Warning", "Select at least two conversion schemes.")
            return
        self.log_message(f"Comparing {len(selected)} conversion schemes...", "INFO")
        self.update_status("Comparing conversions...")
        threading.Thread(target=self._compare_multiple_conversions_thread, args=(selected,), daemon=True).start()

    def _compare_multiple_conversions_thread(self, convs):
        try:
            if self.visualizer is None:
                self.visualizer = PPGVisualizer(self.calculator)
            fig = self.visualizer.plot_multiple_conversion_comparison(convs, fontsize=self.fontsize_var.get())
            if fig:
                self.root.after(0, lambda: self.display_comparison_figure(fig, convs))
                self.root.after(0, lambda: self.log_message(f"✓ Comparison plot generated.", "SUCCESS"))
            else:
                self.root.after(0, lambda: self.log_message("✗ Could not generate comparison plot.", "ERROR"))
        except Exception as e:
            self.root.after(0, lambda: self.log_message(f"✗ Comparison failed: {str(e)}", "ERROR"))
        self.root.after(0, lambda: self.update_status("Ready"))

    def display_comparison_figure(self, fig, convs):
        self.clear_figure()
        self.viz_placeholder.pack_forget()
        canvas = FigureCanvasTkAgg(fig, master=self.viz_placeholder.master)
        canvas.draw()
        toolbar = NavigationToolbar2Tk(canvas, self.viz_placeholder.master)
        toolbar.update()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        toolbar.pack(fill=tk.X)
        self.figure_canvas = canvas
        self.figure_toolbar = toolbar
        self.current_figure = fig
        self.viz_status_var.set("Comparison plot generated")

    def export_conversion_results(self):
        if not self.calculator.conversion_results:
            messagebox.showinfo("Info", "No conversion results to export.")
            return
        file = filedialog.asksaveasfilename(defaultextension=".xlsx",
                                             filetypes=[("Excel files", "*.xlsx")])
        if not file:
            return
        try:
            with pd.ExcelWriter(file, engine='openpyxl') as writer:
                for key, conv in self.calculator.conversion_results.items():
                    if 'data' in conv:
                        df = pd.DataFrame(conv['data'])
                        df.to_excel(writer, sheet_name=key[:30], index=False)
                stats = []
                for key, conv in self.calculator.conversion_results.items():
                    s = {'Direction': key, 'Source': conv.get('condition',''), 'Target': conv.get('condition',''),
                         'Total': conv.get('compound',0), 'Valid': conv.get('',0)}
                    if isinstance(conv.get(''), dict):
                        s.update({k:v for k,v in conv[''].items() if isinstance(v,(int,float))})
                    if isinstance(conv.get('analysis'), dict):
                        s['Pass rate (%)'] = conv['analysis'].get('(%)',0)
                    stats.append(s)
                pd.DataFrame(stats).to_excel(writer, sheet_name='Statistics', index=False)
            self.log_message(f"✓ Results exported to {file}", "SUCCESS")
        except Exception as e:
            self.log_message(f"✗ Export failed: {str(e)}", "ERROR")

    def clear_conversion_results(self):
        for item in self.conversion_tree.get_children():
            self.conversion_tree.delete(item)
        self.conversion_stats_text.delete(1.0, tk.END)
        self.calculator.conversion_results.clear()
        self.conversion_status_var.set("Cleared")

    def on_viz_type_change(self, event=None):
        typ = self.viz_type_var.get()
        self.viz_conditions_frame.grid_remove()
        self.conversion_viz_frame.grid_remove()
        self.viz_condition_combo.grid()
        if typ in ["condition_comparison", "multiple_conversion_comparison"]:
            self.viz_conditions_frame.grid()
            self.viz_condition_combo.grid_remove()
            if typ == "condition_comparison":
                self.update_viz_condition_checkboxes()
            else:
                self.update_viz_conversion_checkboxes()
        elif typ == "conversion_analysis":
            self.conversion_viz_frame.grid()
            self.viz_condition_combo.grid_remove()
            self.update_conversion_viz_params()

    def update_viz_condition_checkboxes(self):
        for w in self.viz_conditions_frame.winfo_children():
            w.destroy()
        conds = list(self.calculator.ppg_data.keys())
        ttk.Label(self.viz_conditions_frame, text="Select conditions:").grid(row=0, column=0, columnspan=4, sticky=tk.W, pady=5)
        self.viz_condition_vars = {}
        for i, c in enumerate(conds):
            var = tk.BooleanVar(value=(i<2))
            cb = ttk.Checkbutton(self.viz_conditions_frame, text=c, variable=var)
            cb.grid(row=1 + i//4, column=i%4, sticky=tk.W, padx=5, pady=2)
            self.viz_condition_vars[c] = var

    def update_viz_conversion_checkboxes(self):
        for w in self.viz_conditions_frame.winfo_children():
            w.destroy()
        conds = list(self.calculator.ppg_data.keys())
        ttk.Label(self.viz_conditions_frame, text="Select conversion schemes:").grid(row=0, column=0, columnspan=4, sticky=tk.W, pady=5)
        self.viz_conversion_vars = []
        if len(conds) >= 2:
            for i in range(len(conds)):
                for j in range(len(conds)):
                    if i != j:
                        from_c, to_c = conds[i], conds[j]
                        var = tk.BooleanVar(value=(i==0 and j==1))
                        cb = ttk.Checkbutton(self.viz_conditions_frame,
                                             text=f"{from_c}→{to_c}", variable=var)
                        row = (i*len(conds)+j)//4 + 1
                        col = (i*len(conds)+j)%4
                        cb.grid(row=row, column=col, sticky=tk.W, padx=2, pady=1)
                        self.viz_conversion_vars.append((var, from_c, to_c))
        else:
            ttk.Label(self.viz_conditions_frame, text="Need at least 2 conditions",
                      foreground="gray").grid(row=1, column=0, sticky=tk.W)

    def update_conversion_viz_params(self):
        for w in self.conversion_viz_frame.winfo_children():
            w.destroy()
        conds = list(self.calculator.ppg_data.keys())
        ttk.Label(self.conversion_viz_frame, text="Source:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        self.viz_from_condition_var = tk.StringVar()
        from_combo = ttk.Combobox(self.conversion_viz_frame, textvariable=self.viz_from_condition_var,
                                   values=conds, width=15, state="readonly")
        from_combo.grid(row=0, column=1, sticky=tk.W, padx=5, pady=5)
        if conds:
            self.viz_from_condition_var.set(conds[0])
        ttk.Label(self.conversion_viz_frame, text="Target:").grid(row=0, column=2, sticky=tk.W, padx=5, pady=5)
        self.viz_to_condition_var = tk.StringVar()
        to_combo = ttk.Combobox(self.conversion_viz_frame, textvariable=self.viz_to_condition_var,
                                 values=conds, width=15, state="readonly")
        to_combo.grid(row=0, column=3, sticky=tk.W, padx=5, pady=5)
        if len(conds) > 1:
            self.viz_to_condition_var.set(conds[1])

    def generate_visualization(self):
        typ = self.viz_type_var.get()
        fs = self.fontsize_var.get()
        if typ == "condition_comparison":
            sel = [c for c, var in self.viz_condition_vars.items() if var.get()]
            if len(sel) < 2:
                messagebox.showwarning("Warning", "Select at least two conditions.")
                return
            self.create_visualization(typ, sel, fs)
        elif typ == "multiple_conversion_comparison":
            sel = [(f,t) for var,f,t in self.viz_conversion_vars if var.get()]
            if len(sel) < 2:
                messagebox.showwarning("Warning", "Select at least two conversion schemes.")
                return
            self.create_visualization(typ, sel, fs)
        elif typ == "conversion_analysis":
            f = self.viz_from_condition_var.get()
            t = self.viz_to_condition_var.get()
            if not f or not t or f == t:
                messagebox.showwarning("Warning", "Select different source and target.")
                return
            self.create_visualization(typ, (f,t), fs)
        else:
            cond = self.viz_condition_var.get()
            if not cond:
                messagebox.showwarning("Warning", "Select a condition.")
                return
            self.create_visualization(typ, cond, fs)

    def create_visualization(self, typ, params, fs):
        try:
            if self.visualizer is None:
                self.visualizer = PPGVisualizer(self.calculator)
            self.viz_status_var.set(f"Generating {typ}...")
            self.update_status("Generating plot")

            if typ == "standard_curve":
                if isinstance(params, str):
                    params = [params]
                fig = self.visualizer.plot_standard_curves(params, fontsize=fs)
            elif typ == "residuals":
                if isinstance(params, str):
                    params = [params]
                fig = self.visualizer.plot_residuals(params, fontsize=fs)
            elif typ == "index_distribution":
                fig = self.visualizer.plot_ppg_index_distribution(params, fontsize=fs)
            elif typ == "condition_comparison":
                fig = self.visualizer.plot_condition_comparison(params, fontsize=fs)
            elif typ == "conversion_analysis":
                from_c, to_c = params
                fig = self.visualizer.plot_conversion_analysis(from_c, to_c, fontsize=fs)
            elif typ == "multiple_conversion_comparison":
                fig = self.visualizer.plot_multiple_conversion_comparison(params, fontsize=fs)
            else:
                self.viz_status_var.set("Unsupported type")
                return

            if fig is None:
                self.viz_status_var.set("Failed to generate. Check data.")
                return

            self.clear_figure()
            self.viz_placeholder.pack_forget()

            canvas = FigureCanvasTkAgg(fig, master=self.viz_placeholder.master)
            canvas.draw()
            toolbar = NavigationToolbar2Tk(canvas, self.viz_placeholder.master)
            toolbar.update()
            canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
            toolbar.pack(fill=tk.X)

            self.figure_canvas = canvas
            self.figure_toolbar = toolbar
            self.current_figure = fig

            self.viz_status_var.set("Plot generated")
            self.update_status("Ready")
        except Exception as e:
            self.viz_status_var.set(f"Error: {str(e)}")
            self.log_message(f"✗ Plot failed: {str(e)}", "ERROR")

    def save_figure(self):
        if self.current_figure is None:
            messagebox.showwarning("Warning", "No figure to save.")
            return
        file = filedialog.asksaveasfilename(defaultextension=".png",
                                             filetypes=[("PNG", "*.png"), ("PDF", "*.pdf"), ("SVG", "*.svg")])
        if file:
            try:
                self.current_figure.savefig(file, dpi=300, bbox_inches='tight', transparent=True)
                self.log_message(f"✓ Figure saved to {file}", "SUCCESS")
            except Exception as e:
                self.log_message(f"✗ Save failed: {str(e)}", "ERROR")

    def clear_figure(self):
        if self.figure_canvas:
            self.figure_canvas.get_tk_widget().destroy()
            self.figure_toolbar.destroy()
            self.figure_canvas = None
            self.figure_toolbar = None
        if self.current_figure:
            plt.close(self.current_figure)
            self.current_figure = None
        self.viz_placeholder.pack(expand=True)

    def generate_report(self):
        try:
            summary = self.calculator.generate_summary_report()
            self.analysis_message("="*60, "INFO")
            self.analysis_message("PPG Retention Index Analysis Report", "INFO")
            self.analysis_message("="*60, "INFO")
            self.analysis_message(f"Generated: {summary['']}", "INFO")
            self.analysis_message("", "INFO")
            self.analysis_message(f"PPG datasets: {summary['PPG datacondition']}", "INFO")
            self.analysis_message(f"Compound datasets: {summary['compound data']}", "INFO")
            self.analysis_message(f"Calibration curves: {summary['standard curve']}", "INFO")
            self.analysis_message(f"PPG index calculations: {summary['PPGindexcalculateresults']}", "INFO")
            self.analysis_message(f"Cross‑condition conversions: {summary['cross-condition conversion results']}", "INFO")
            self.analysis_message("", "INFO")
            if summary['standard curve']:
                self.analysis_message("Calibration curve performance:", "INFO")
                for cond, perf in summary['standard curve'].items():
                    self.analysis_message(f" {cond}: R²={perf['R²']:.4f}, slope={perf['']:.4f}", "INFO")
            if summary['PPGindex']:
                self.analysis_message("\nPPG index statistics:", "INFO")
                for cond, stat in summary['PPGindex'].items():
                    self.analysis_message(f" {cond}: n={stat['']}, mean={stat['']:.2f}, sd={stat['']:.2f}", "INFO")
            if summary['condition']:
                self.analysis_message("\nConversion statistics:", "INFO")
                for key, stat in summary['condition'].items():
                    self.analysis_message(f" {key}: valid={stat['']}, mean error={stat['']:.3f} min, pass rate={stat['(%)']:.1f}%", "INFO")
            self.analysis_message("\nReport generation complete.", "SUCCESS")
            self.results_text.delete(1.0, tk.END)
            self.results_text.insert(tk.END, "PPG Retention Index Analysis Report\n")
            self.results_text.insert(tk.END, "="*60 + "\n\n")
            self.results_text.insert(tk.END, f"Generated: {summary['']}\n\n")
            self.results_text.insert(tk.END, f"PPG datasets: {summary['PPG datacondition']}\n")
            self.results_text.insert(tk.END, f"Compound datasets: {summary['compound data']}\n")
            self.results_text.insert(tk.END, f"Calibration curves: {summary['standard curve']}\n")
            self.results_text.insert(tk.END, f"PPG index calculations: {summary['PPGindexcalculateresults']}\n")
            self.results_text.insert(tk.END, f"Cross‑condition conversions: {summary['cross-condition conversion results']}\n\n")
            self.analysis_status_var.set("Report generated")
        except Exception as e:
            self.analysis_message(f"✗ Report generation failed: {str(e)}", "ERROR")

    def save_all_results(self):
        outdir = self.output_dir_var.get().strip()
        if not outdir:
            messagebox.showwarning("Warning", "Select an output directory.")
            return
        self.output_status_var.set("Saving results...")
        self.update_status("Saving...")
        threading.Thread(target=self._save_all_results_thread, args=(outdir,), daemon=True).start()

    def _save_all_results_thread(self, outdir):
        try:
            success, msg, files = self.calculator.save_results(outdir)
            if success:
                self.root.after(0, lambda: self.output_status_var.set(f"Saved {len(files)} files"))
                self.root.after(0, lambda: self.log_message(f"✓ {msg}", "SUCCESS"))
                self.root.after(0, lambda: self.results_text.delete(1.0, tk.END))
                self.root.after(0, lambda: self.results_text.insert(tk.END, "Saved files:\n" + "\n".join([Path(f).name for f in files])))
            else:
                self.root.after(0, lambda: self.output_status_var.set("Save failed"))
                self.root.after(0, lambda: self.log_message(f"✗ {msg}", "ERROR"))
        except Exception as e:
            self.root.after(0, lambda: self.log_message(f"✗ Save failed: {str(e)}", "ERROR"))
        self.root.after(0, lambda: self.update_status("Ready"))

    def preview_report(self):
        self.generate_report()
        self.notebook.select(4)

    def open_output_dir(self):
        d = self.output_dir_var.get().strip()
        if d and os.path.exists(d):
            try:
                if sys.platform == 'win32':
                    os.startfile(d)
                elif sys.platform == 'darwin':
                    os.system(f'open "{d}"')
                else:
                    os.system(f'xdg-open "{d}"')
            except Exception as e:
                messagebox.showerror("Error", f"Cannot open directory: {str(e)}")
        else:
            messagebox.showwarning("Warning", "Output directory does not exist or not set.")

    def log_message(self, msg, level="INFO"):
        ts = datetime.now().strftime("%H:%M:%S")
        self.log_text.insert(tk.END, f"[{ts}] {msg}\n", level)
        self.log_text.see(tk.END)
        self.root.update()

    def analysis_message(self, msg, level="INFO"):
        self.analysis_text.insert(tk.END, msg + "\n", level)
        self.analysis_text.see(tk.END)
        self.root.update()

    def clear_log(self):
        self.log_text.delete(1.0, tk.END)

    def clear_analysis_text(self):
        self.analysis_text.delete(1.0, tk.END)

    def save_log(self):
        file = filedialog.asksaveasfilename(defaultextension=".txt",
                                             filetypes=[("Text files", "*.txt")])
        if file:
            try:
                with open(file, 'w', encoding='utf-8') as f:
                    f.write(self.log_text.get(1.0, tk.END))
                self.log_message(f"✓ Log saved to {file}", "SUCCESS")
            except Exception as e:
                self.log_message(f"✗ Save log failed: {str(e)}", "ERROR")

    def update_status(self, msg):
        self.status_var.set(msg)
        self.root.update()

    def on_closing(self):
        if messagebox.askyesno("Confirm", "Exit program?"):
            if self.current_figure:
                plt.close(self.current_figure)
            self.root.destroy()

# =============================================================================
# translated note
# =============================================================================
def main():
    print("PPG Retention Index Analyzer - Journal Version")
    print("Features:")
    print("  1. Load PPG standard and compound data")
    print("  2. Fit calibration curves and calculate PPG indices")
    print("  3. Compare PPG indices under different conditions")
    print("  4. Cross-condition conversion and validation")
    print("  5. Publication-quality plots (English, transparent, adjustable font)")
    print("="*70)

    if not GUI_AVAILABLE:
        print("Error: tkinter not available. Cannot start GUI.")
        return

    root = tk.Tk()
    root.title("PPG Retention Index Analyzer - Journal Version")
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()
    window_width = int(screen_width * 0.8)
    window_height = int(screen_height * 0.8)
    x = (screen_width - window_width) // 2
    y = (screen_height - window_height) // 2
    root.geometry(f"{window_width}x{window_height}+{x}+{y}")

    app = PPGIndexAnalyzerGUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()
