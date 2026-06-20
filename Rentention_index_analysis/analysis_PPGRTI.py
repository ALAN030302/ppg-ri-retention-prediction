
# !/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PPG retention index calculation and visualization program - enhanced version with cross-condition conversion
:
1. Load retention-time data for PPG standards and compounds
2. Fit PPG standard curves and calculate linear relationships
3. Calculate compound PPG retention indices
4. Compare PPG indices across chromatographic conditions
5. Convert and validate PPG indices across conditions
6. Visualize analysis results and generate reports

Implement a unified LC retention-index framework based on the experimental design
"""

import os
import sys
import threading
import traceback
import warnings
from datetime import datetime
from typing import Dict, List, Tuple, Optional, Union, Any
from pathlib import Path

# Suppress warnings
warnings.filterwarnings('ignore')

try:
    import pandas as pd
    import numpy as np
    from scipy import stats
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    from matplotlib.figure import Figure
except ImportError:
    print("Error: please install the required libraries first")
    print("Installation command: pip install pandas numpy scipy matplotlib seaborn")
    sys.exit(1)

# Try to import GUI libraries
try:
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox, scrolledtext, Toplevel

    GUI_AVAILABLE = True
except ImportError:
    GUI_AVAILABLE = False
    print("Warning: tkinter is not installed, so the GUI is unavailable")


class PPGIndexCalculator:
    """Core class for PPG retention index calculation"""

    def __init__(self):
        """Initialize the calculator"""
        self.ppg_data = {} # store PPG data under different conditions
        self.compound_data = {} # store compound data
        self.standard_curves = {} # store standard-curve parameters
        self.ppg_indices = {} # store calculated PPG indices
        self.results_summary = {} # store result summaries
        self.conversion_results = {} # store cross-condition conversion results

    def load_ppg_data(self, file_path: str, condition: str = "default") -> Tuple[bool, str]:
        """
        Load PPG standard data

        Parameters:
            file_path: data file path (supports Excel and CSV)
            condition: chromatographic condition identifier

        Returns:
            (success flag, message)
        """
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

            return True, f"load {len(df)} PPG standard data (condition: {condition})"

        except Exception as e:
            return False, f"loadPPG datafailed: {str(e)}"

    def load_compound_data(self, file_path: str, category: str = "validation",
                           condition: str = "default") -> Tuple[bool, str]:
        """
        loadcompound data

        Parameters:
            file_path: data file path
            category: data ('validation'validation set, 'smrt'training)
            condition: chromatographic condition identifier

        Returns:
            (success flag, message)
        """
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
                'compound_name': ['compound_name', 'name', 'compound', 'Name', 'Compound', 'compound'],
                'retention_time': ['retention_time', 'RT', 'RetentionTime', 't_R', 'retention_time(RT)', 'rt'],
                'CAS': ['CAS', 'CAS', 'CAS No.', 'CAS']
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
            if 'compound_name' not in df.columns or 'retention_time' not in df.columns:
                return False, "datafilecolumn ('compound_name''retention_time')"

            # data
            df['retention_time'] = pd.to_numeric(df['retention_time'], errors='coerce')
            df = df.dropna(subset=['compound_name', 'retention_time'])

            # data
            key = f"{category}_{condition}"
            self.compound_data[key] = df

            return True, f"load {len(df)} compound data (: {category}, condition: {condition})"

        except Exception as e:
            return False, f"loadcompound datafailed: {str(e)}"

    def fit_standard_curve(self, condition: str = "default",
                           model_type: str = "logarithmic") -> Tuple[bool, str]:
        """
        PPGstandard curve

        Parameters:
            condition: chromatographic condition identifier
            model_type: model ('linear', 'logarithmic')

        Returns:
            (success flag, message)
        """
        try:
            if condition not in self.ppg_data:
                return False, f"condition {condition} PPG data"

            df = self.ppg_data[condition]

            if len(df) < 3:
                return False, "PPG data, 3standard curve"

            x = df['degree_of_polymerization'].values
            y = df['retention_time'].values

            if model_type == "logarithmic":
                # model: RT = a + b * ln(n)
                x_fit = np.log(x)
                model_name = "model (RT = a + b * ln(n))"
            elif model_type == "linear":
                # model: RT = a + b * n
                x_fit = x
                model_name = "model (RT = a + b * n)"
            else:
                return False, f"model: {model_type}"

            # translated note
            slope, intercept, r_value, p_value, std_err = stats.linregress(x_fit, y)

            # calculate
            y_pred = intercept + slope * x_fit
            residuals = y - y_pred

            # store standard-curve parameters
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
        """
        calculatePPG retention indices

        Parameters:
            condition: chromatographic condition identifier
            method: calculatemethod ('interpolation', 'regression')

        Returns:
            (success flag, message)
        """
        try:
            if condition not in self.ppg_data:
                return False, f"condition {condition} PPG data"

            df_ppg = self.ppg_data[condition]

            if method == "interpolation":
                # translated note
                ppg_rt = df_ppg['retention_time'].values
                ppg_n = df_ppg['degree_of_polymerization'].values

                # compound datacalculatePPGindex
                indices = {}

                for key in self.compound_data:
                    if condition in key:
                        df_comp = self.compound_data[key]
                        results = []

                        for _, row in df_comp.iterrows():
                            rt = row['retention_time']
                            compound_name = row['compound_name']

                            # translated note
                            if rt < ppg_rt[0]:
                                # translated note
                                n_calc = ppg_n[0] - (ppg_rt[0] - rt) / (ppg_rt[1] - ppg_rt[0]) * (ppg_n[1] - ppg_n[0])
                                if n_calc < 0:
                                    n_calc = 0
                                method_used = " (PPG RT)"
                            elif rt > ppg_rt[-1]:
                                # translated note
                                n_calc = ppg_n[-1] + (rt - ppg_rt[-1]) / (ppg_rt[-1] - ppg_rt[-2]) * (
                                        ppg_n[-1] - ppg_n[-2])
                                method_used = " (PPG RT)"
                            else:
                                # translated note
                                idx = np.searchsorted(ppg_rt, rt) - 1
                                if idx < 0:
                                    idx = 0
                                elif idx >= len(ppg_rt) - 1:
                                    idx = len(ppg_rt) - 2

                                rt_i, rt_j = ppg_rt[idx], ppg_rt[idx + 1]
                                n_i, n_j = ppg_n[idx], ppg_n[idx + 1]

                                n_calc = n_i + (n_j - n_i) * (rt - rt_i) / (rt_j - rt_i)
                                method_used = ""

                            # retention_index (100)
                            ppg_index = n_calc * 100

                            result = {
                                'compound_name': compound_name,
                                'retention_time': rt,
                                'calculatePPGindex': ppg_index,
                                'calculatemethod': method_used,
                            }

                            # column
                            for col in df_comp.columns:
                                if col not in result:
                                    result[col] = row[col]

                            results.append(result)

                        indices[key] = pd.DataFrame(results)

            elif method == "regression":
                # (Usestandard curve)
                if condition not in self.standard_curves:
                    success, msg = self.fit_standard_curve(condition, model_type="logarithmic")
                    if not success:
                        return False, f"Use: {msg}"

                curve = self.standard_curves[condition]

                # compoundcalculatePPGindex
                indices = {}

                for key in self.compound_data:
                    if condition in key:
                        df_comp = self.compound_data[key]
                        results = []

                        for _, row in df_comp.iterrows():
                            rt = row['retention_time']
                            compound_name = row['compound_name']

                            # RTn
                            if curve['model_type'] == "logarithmic":
                                # RT = a + b * ln(n) => n = exp((RT - a) / b)
                                if curve['slope'] != 0:
                                    n_calc = np.exp((rt - curve['intercept']) / curve['slope'])
                                else:
                                    n_calc = np.nan
                            else:  # linear
                                # RT = a + b * n => n = (RT - a) / b
                                if curve['slope'] != 0:
                                    n_calc = (rt - curve['intercept']) / curve['slope']
                                else:
                                    n_calc = np.nan

                            ppg_index = n_calc * 100 if not np.isnan(n_calc) else np.nan

                            result = {
                                'compound_name': compound_name,
                                'retention_time': rt,
                                'calculatePPGindex': ppg_index,
                                'calculatemethod': "",
                            }

                            # column
                            for col in df_comp.columns:
                                if col not in result:
                                    result[col] = row[col]

                            results.append(result)

                        indices[key] = pd.DataFrame(results)
            else:
                return False, f"calculatemethod: {method}"

            # results
            self.ppg_indices[condition] = {
                'method': method,
                'indices': indices
            }

            return True, f"PPGindexcalculate (condition: {condition}, method: {method})"

        except Exception as e:
            return False, f"calculatePPGindexfailed: {str(e)}"

    def compare_conditions(self, conditions: List[str]) -> pd.DataFrame:
        """
        Compare PPG indices across chromatographic conditions

        Parameters:
            conditions: comparechromatographic conditioncolumn

        Returns:
            compareresultsDataFrame
        """
        comparison_results = []

        for key in self.compound_data:
            # comparecondition
            if any(cond in key for cond in conditions):
                category = key.split('_')[0]
                condition = key.split('_')[1] if '_' in key else "default"

                # conditioncompound data
                df_comp = self.compound_data[key]

                for _, row in df_comp.iterrows():
                    compound_name = row['compound_name']
                    rt = row['retention_time']

                    # compoundconditiondata
                    compound_data = {
                        'compound_name': compound_name,
                        'data': category,
                        f'{condition}_RT': rt
                    }

                    # PPGindex (calculate)
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
        """
        calculatecondition

        Parameters:
            from_condition: chromatographic condition
            to_condition: chromatographic condition

        Returns:
            analysis results
        """
        error_results = []

        # conditionPPGindex
        if from_condition not in self.ppg_indices or to_condition not in self.ppg_indices:
            return pd.DataFrame()

        # conditioncompound
        compounds_in_both = set()

        for key in self.compound_data:
            if from_condition in key:
                df = self.compound_data[key]
                compounds_in_both.update(df['compound_name'].tolist())

        for key in self.compound_data:
            if to_condition in key:
                df = self.compound_data[key]
                compounds_in_both.intersection_update(set(df['compound_name'].tolist()))

        # calculate
        for compound in compounds_in_both:
            # conditionPPGindex
            from_ppg = None
            to_ppg = None

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
        """
        conditionPPGindexconditionretention_time

        Parameters:
            from_condition: chromatographic condition (PPGindex)
            to_condition: chromatographic condition (condition)
            compound_names: compound_namecolumn (Nonecompound)

        Returns:
            (resultsDataFrame, message)
        """
        try:
            # conditionPPGindex
            if from_condition not in self.ppg_indices:
                return pd.DataFrame(), f"condition {from_condition} PPGindexdata"

            # conditionstandard curve
            if to_condition not in self.standard_curves:
                # standard curve
                success, msg = self.fit_standard_curve(to_condition, "logarithmic")
                if not success:
                    return pd.DataFrame(), f"condition {to_condition} standard curve: {msg}"

            # conditionstandard curve
            curve = self.standard_curves[to_condition]

            # results
            conversion_results = []

            # conditionPPGindex
            for key, indices_df in self.ppg_indices[from_condition]['indices'].items():
                if from_condition in key:
                    for _, row in indices_df.iterrows():
                        compound_name = row['compound_name']

                        # compoundcolumn, compound
                        if compound_names and compound_name not in compound_names:
                            continue

                        # conditionPPGindex
                        ppg_index = row['calculatePPGindex']
                        if pd.isna(ppg_index):
                            continue

                        # PPGindexdegree_of_polymerizationn (100)
                        n_calc = ppg_index / 100

                        # Useconditionstandard curvecalculateretention_time
                        if curve['model_type'] == "logarithmic":
                            # model: RT = a + b * ln(n)
                            rt_pred = curve['intercept'] + curve['slope'] * np.log(n_calc)
                        else:  # linear
                            # model: RT = a + b * n
                            rt_pred = curve['intercept'] + curve['slope'] * n_calc

                        # conditionretention_time
                        rt_actual = None
                        for comp_key, comp_df in self.compound_data.items():
                            if to_condition in comp_key:
                                match = comp_df[comp_df['compound_name'] == compound_name]
                                if not match.empty and 'retention_time' in match.columns:
                                    rt_actual = match.iloc[0]['retention_time']
                                    break

                        # calculate
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

            # calculate
            if not conversion_df.empty:
                valid_errors = conversion_df['(min)'].dropna()
                if len(valid_errors) > 0:
                    stats = {
                        '': valid_errors.mean(),
                        '': valid_errors.std(),
                        '': valid_errors.max(),
                        '': valid_errors.min(),
                        '': valid_errors.median(),
                        '': len(valid_errors)
                    }

                    # calculate
                    valid_rel_errors = conversion_df['(%)'].dropna()
                    if len(valid_rel_errors) > 0:
                        stats.update({
                            '(%)': valid_rel_errors.mean(),
                            '(%)': valid_rel_errors.std(),
                            '(%)': valid_rel_errors.max(),
                            '(%)': valid_rel_errors.min()
                        })

                    return conversion_df, stats
                else:
                    return conversion_df, "data"
            else:
                return conversion_df, "matchcompound data"

        except Exception as e:
            return pd.DataFrame(), f"failed: {str(e)}"

    def cross_condition_analysis(self, from_condition: str, to_condition: str,
                                 threshold: float = 0.5) -> Dict[str, Any]:
        """
        conditionanalysis

        Parameters:
            from_condition: chromatographic condition
            to_condition: chromatographic condition
            threshold: ()

        Returns:
            analysis results
        """
        try:
            # translated note
            conversion_df, stats = self.convert_ppg_index_to_rt(from_condition, to_condition)

            if conversion_df.empty:
                return {"error": "data"}

            # translated note
            analysis_results = {
                'condition': from_condition,
                'condition': to_condition,
                '': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                'compound': len(conversion_df),
                '': len(conversion_df['(min)'].dropna()),
                '': stats if isinstance(stats, dict) else stats,
                'data': conversion_df.to_dict('records')
            }

            # translated note
            if '(min)' in conversion_df.columns:
                errors = conversion_df['(min)'].dropna()

                # translated note
                error_bins = [0, 0.1, 0.2, 0.5, 1.0, float('inf')]
                error_labels = ['<0.1min', '0.1-0.2min', '0.2-0.5min', '0.5-1.0min', '>1.0min']

                error_distribution = {}
                for i in range(len(error_bins) - 1):
                    lower = error_bins[i]
                    upper = error_bins[i + 1]
                    if i == len(error_bins) - 2:
                        count = len(errors[errors >= lower])
                    else:
                        count = len(errors[(errors >= lower) & (errors < upper)])
                    error_distribution[error_labels[i]] = count

                analysis_results[''] = error_distribution

                # analysis ()
                passed = len(errors[errors <= threshold])
                pass_rate = (passed / len(errors) * 100) if len(errors) > 0 else 0

                analysis_results['analysis'] = {
                    '(min)': threshold,
                    '': passed,
                    '': len(errors),
                    '(%)': pass_rate
                }

            # results
            key = f"{from_condition}_to_{to_condition}"
            self.conversion_results[key] = analysis_results

            return analysis_results

        except Exception as e:
            return {"error": f"analysisfailed: {str(e)}"}

    def generate_summary_report(self) -> Dict[str, Any]:
        """
        analysis

        Returns:
            translated note
        """
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

        # standard curve
        for condition, curve in self.standard_curves.items():
            summary['standard curve'][condition] = {
                'model': curve['model_type'],
                'R²': curve['r_squared'],
                '': curve['slope'],
                '': curve['intercept'],
                '': curve['std_err'],
                'data': curve['n_points']
            }

        # PPGindex
        for condition, indices_data in self.ppg_indices.items():
            all_indices = []
            for key, df in indices_data['indices'].items():
                if 'calculatePPGindex' in df.columns:
                    valid_indices = df['calculatePPGindex'].dropna()
                    all_indices.extend(valid_indices.tolist())

            if all_indices:
                indices_array = np.array(all_indices)
                summary['PPGindex'][condition] = {
                    'calculatemethod': indices_data['method'],
                    '': len(all_indices),
                    '': np.mean(indices_array),
                    '': np.std(indices_array),
                    '': np.min(indices_array),
                    '': np.max(indices_array),
                    '': np.median(indices_array)
                }

        # condition
        for key, conversion in self.conversion_results.items():
            summary['condition'][key] = {
                'condition': conversion.get('condition', ''),
                'condition': conversion.get('condition', ''),
                'compound': conversion.get('compound', 0),
                '': conversion.get('', 0),
                '': conversion.get('', {}).get('', 0) if isinstance(
                    conversion.get(''), dict) else 0,
                '(%)': conversion.get('analysis', {}).get('(%)', 0) if isinstance(
                    conversion.get('analysis'), dict) else 0
            }

        self.results_summary = summary
        return summary

    def save_results(self, output_dir: str) -> Tuple[bool, str, List[str]]:
        """
        saveresultsfile

        Parameters:
            output_dir: outputdirectory

        Returns:
            (success flag, message, savefilecolumn)
        """
        try:
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)

            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            saved_files = []

            # 1. savePPGstandard curvedata
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
                curves_file = output_path / f"PPGstandard curve_{timestamp}.xlsx"
                with pd.ExcelWriter(curves_file, engine='openpyxl') as writer:
                    curves_df.to_excel(writer, sheet_name='standard curve', index=False)

                    # conditiondata
                    for condition, curve in self.standard_curves.items():
                        detail_df = pd.DataFrame({
                            'degree_of_polymerization': curve['x'],
                            'RT': curve['y'],
                            'RT': curve['y_pred'],
                            '': curve['residuals']
                        })
                        detail_df.to_excel(writer, sheet_name=f'{condition}_data', index=False)

                saved_files.append(str(curves_file))

            # 2. savePPGindexcalculateresults
            if self.ppg_indices:
                for condition, indices_data in self.ppg_indices.items():
                    indices_file = output_path / f"PPGindex_{condition}_{timestamp}.xlsx"

                    with pd.ExcelWriter(indices_file, engine='openpyxl') as writer:
                        for key, df in indices_data['indices'].items():
                            # sheetname
                            sheet_name = key.replace('_', '-')[:30]
                            if sheet_name in writer.sheets:
                                sheet_name = f"{sheet_name[:25]}_{hash(key) % 1000:03d}"

                            df.to_excel(writer, sheet_name=sheet_name, index=False)

                    saved_files.append(str(indices_file))

            # 3. savecross-condition conversion results
            if self.conversion_results:
                conversion_file = output_path / f"cross-condition conversion results_{timestamp}.xlsx"
                with pd.ExcelWriter(conversion_file, engine='openpyxl') as writer:
                    for key, conversion in self.conversion_results.items():
                        if 'data' in conversion:
                            df = pd.DataFrame(conversion['data'])
                            # sheetname
                            sheet_name = key[:30]
                            if sheet_name in writer.sheets:
                                sheet_name = f"{sheet_name[:25]}_{hash(key) % 1000:03d}"
                            df.to_excel(writer, sheet_name=sheet_name, index=False)

                    # save
                    conversion_stats = []
                    for key, conversion in self.conversion_results.items():
                        stats = {
                            '': key,
                            'condition': conversion.get('condition', ''),
                            'condition': conversion.get('condition', ''),
                            'compound': conversion.get('compound', 0),
                            '': conversion.get('', 0)
                        }

                        if isinstance(conversion.get(''), dict):
                            stats.update({
                                '(min)': conversion[''].get('', 0),
                                '': conversion[''].get('', 0),
                                '(min)': conversion[''].get('', 0)
                            })

                        if isinstance(conversion.get('analysis'), dict):
                            stats.update({
                                '(%)': conversion['analysis'].get('(%)', 0),
                                '': conversion['analysis'].get('', 0),
                                '(min)': conversion['analysis'].get('(min)', 0.5)
                            })

                        conversion_stats.append(stats)

                    if conversion_stats:
                        stats_df = pd.DataFrame(conversion_stats)
                        stats_df.to_excel(writer, sheet_name='', index=False)

                saved_files.append(str(conversion_file))

            # 4. saveanalysis
            if self.results_summary:
                report_file = output_path / f"analysis_{timestamp}.txt"
                with open(report_file, 'w', encoding='utf-8') as f:
                    f.write("=" * 70 + "\n")
                    f.write("PPG retention indicesanalysis\n")
                    f.write("=" * 70 + "\n\n")

                    f.write(f": {self.results_summary['']}\n\n")

                    f.write("data:\n")
                    f.write(f" - PPG datacondition: {self.results_summary['PPG datacondition']}\n")
                    f.write(f" - compound data: {self.results_summary['compound data']}\n")
                    f.write(f" - standard curve: {self.results_summary['standard curve']}\n")
                    f.write(f" - PPGindexcalculateresults: {self.results_summary['PPGindexcalculateresults']}\n")
                    f.write(f" - cross-condition conversion results: {self.results_summary['cross-condition conversion results']}\n\n")

                    if self.results_summary['standard curve']:
                        f.write("standard curve:\n")
                        for condition, perf in self.results_summary['standard curve'].items():
                            f.write(f"  {condition}:\n")
                            f.write(f" - model: {perf['model']}\n")
                            f.write(f"    - R²: {perf['R²']:.4f}\n")
                            f.write(f" - : {perf['']:.4f}\n")
                            f.write(f" - : {perf['']:.4f}\n")
                            f.write(f" - : {perf['']:.4f}\n")
                            f.write(f" - data: {perf['data']}\n")
                        f.write("\n")

                    if self.results_summary['PPGindex']:
                        f.write("PPGindex:\n")
                        for condition, stats in self.results_summary['PPGindex'].items():
                            f.write(f"  {condition}:\n")
                            f.write(f" - calculatemethod: {stats['calculatemethod']}\n")
                            f.write(f" - : {stats['']}\n")
                            f.write(f" - : {stats['']:.2f}\n")
                            f.write(f" - : {stats['']:.2f}\n")
                            f.write(f" - : {stats['']:.2f} - {stats['']:.2f}\n")
                            f.write(f" - : {stats['']:.2f}\n")
                        f.write("\n")

                    if self.results_summary['condition']:
                        f.write("condition:\n")
                        for key, stats in self.results_summary['condition'].items():
                            f.write(f"  {key}:\n")
                            f.write(f" - condition: {stats['condition']}\n")
                            f.write(f" - condition: {stats['condition']}\n")
                            f.write(f" - compound: {stats['compound']}\n")
                            f.write(f" - : {stats['']}\n")
                            f.write(f" - : {stats['']:.3f} min\n")
                            f.write(f" - : {stats['(%)']:.1f}%\n")
                        f.write("\n")

                saved_files.append(str(report_file))

            # 5. savecompound data
            if self.compound_data:
                compounds_file = output_path / f"compound data_{timestamp}.xlsx"
                with pd.ExcelWriter(compounds_file, engine='openpyxl') as writer:
                    for key, df in self.compound_data.items():
                        sheet_name = key.replace('_', '-')[:30]
                        if sheet_name in writer.sheets:
                            sheet_name = f"{sheet_name[:25]}_{hash(key) % 1000:03d}"
                        df.to_excel(writer, sheet_name=sheet_name, index=False)

                saved_files.append(str(compounds_file))

            return True, f"resultssave {output_dir}", saved_files

        except Exception as e:
            return False, f"saveresultsfailed: {str(e)}", []


class PPGVisualizer:
    """PPG data"""

    def __init__(self, calculator: PPGIndexCalculator):
        """"""
        self.calculator = calculator
        self.figures = {}

        # Set fonts
        plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
        plt.rcParams['axes.unicode_minus'] = False

    def plot_standard_curves(self, conditions: List[str] = None,
                             save_path: str = None) -> plt.Figure:
        """
        PPGstandard curve

        Parameters:
            conditions: chromatographic conditioncolumn (None)
            save_path: savepath ()

        Returns:
            matplotlib Figure
        """
        if conditions is None:
            conditions = list(self.calculator.standard_curves.keys())

        if not conditions:
            print("Warning: standard curvedata")
            return None

        n_plots = len(conditions)
        n_cols = min(2, n_plots)
        n_rows = (n_plots + n_cols - 1) // n_cols

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows))
        if n_plots == 1:
            axes = np.array([axes])
        if axes.ndim == 1:
            axes = axes.reshape(-1, n_cols)

        for idx, condition in enumerate(conditions):
            if condition not in self.calculator.standard_curves:
                continue

            curve = self.calculator.standard_curves[condition]
            row = idx // n_cols
            col = idx % n_cols

            ax = axes[row, col]

            # raw data
            ax.scatter(curve['x'], curve['y'], color='blue', s=50, label='data', zorder=3)

            # translated note
            if curve['model_type'] == 'logarithmic':
                x_fit = np.log(curve['x'])
                x_range = np.linspace(min(curve['x']), max(curve['x']), 100)
                x_fit_range = np.log(x_range)
            else:
                x_fit = curve['x']
                x_range = np.linspace(min(curve['x']), max(curve['x']), 100)
                x_fit_range = x_range

            y_fit_range = curve['intercept'] + curve['slope'] * x_fit_range
            ax.plot(x_range, y_fit_range, 'r-', label='', linewidth=2)

            # translated note
            info_text = f"model: {curve['model_name']}\n"
            info_text += f"R² = {curve['r_squared']:.4f}\n"
            info_text += f" = {curve['slope']:.4f}\n"
            info_text += f" = {curve['intercept']:.4f}"

            ax.text(0.05, 0.95, info_text, transform=ax.transAxes,
                    verticalalignment='top', fontsize=10,
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

            ax.set_xlabel('degree_of_polymerization (n)')
            ax.set_ylabel('retention_time (min)')
            ax.set_title(f'PPGstandard curve - {condition}')
            ax.legend()
            ax.grid(True, alpha=0.3)

        # translated note
        for idx in range(len(conditions), n_rows * n_cols):
            row = idx // n_cols
            col = idx % n_cols
            axes[row, col].set_visible(False)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"standard curvesave: {save_path}")

        self.figures['standard_curves'] = fig
        return fig

    def plot_residuals(self, conditions: List[str] = None,
                       save_path: str = None) -> plt.Figure:
        """
        translated note

        Parameters:
            conditions: chromatographic conditioncolumn (None)
            save_path: savepath ()

        Returns:
            matplotlib Figure
        """
        if conditions is None:
            conditions = list(self.calculator.standard_curves.keys())

        if not conditions:
            print("Warning: data")
            return None

        n_plots = len(conditions)
        n_cols = min(2, n_plots)
        n_rows = (n_plots + n_cols - 1) // n_cols

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows))
        if n_plots == 1:
            axes = np.array([axes])
        if axes.ndim == 1:
            axes = axes.reshape(-1, n_cols)

        for idx, condition in enumerate(conditions):
            if condition not in self.calculator.standard_curves:
                continue

            curve = self.calculator.standard_curves[condition]
            row = idx // n_cols
            col = idx % n_cols

            ax = axes[row, col]

            # translated note
            residuals = curve['residuals']
            predicted = curve['y_pred']

            ax.scatter(predicted, residuals, color='blue', s=50, alpha=0.7)
            ax.axhline(y=0, color='red', linestyle='--', linewidth=1)

            # translated note
            mean_residual = np.mean(residuals)
            std_residual = np.std(residuals)

            info_text = f":\n"
            info_text += f" = {mean_residual:.4f}\n"
            info_text += f" = {std_residual:.4f}\n"
            info_text += f" = {max(abs(residuals)):.4f}"

            ax.text(0.05, 0.95, info_text, transform=ax.transAxes,
                    verticalalignment='top', fontsize=9,
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

            ax.set_xlabel('retention_time (min)')
            ax.set_ylabel(' (min)')
            ax.set_title(f' - {condition}')
            ax.grid(True, alpha=0.3)

        # translated note
        for idx in range(len(conditions), n_rows * n_cols):
            row = idx // n_cols
            col = idx % n_cols
            axes[row, col].set_visible(False)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"save: {save_path}")

        self.figures['residuals'] = fig
        return fig

    def plot_ppg_index_distribution(self, condition: str,
                                    save_path: str = None) -> plt.Figure:
        """
        PPGindex

        Parameters:
            condition: chromatographic condition
            save_path: savepath ()

        Returns:
            matplotlib Figure
        """
        if condition not in self.calculator.ppg_indices:
            print(f"Warning: condition {condition} PPGindexdata")
            return None

        # PPGindex
        all_indices = []
        for key, df in self.calculator.ppg_indices[condition]['indices'].items():
            if 'calculatePPGindex' in df.columns:
                indices = df['calculatePPGindex'].dropna().tolist()
                all_indices.extend(indices)

        if not all_indices:
            print(f"Warning: condition {condition} PPGindex")
            return None

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # translated note
        ax1.hist(all_indices, bins=30, color='steelblue', edgecolor='black', alpha=0.7)
        ax1.axvline(x=np.mean(all_indices), color='red', linestyle='--',
                    linewidth=2, label=f': {np.mean(all_indices):.1f}')
        ax1.axvline(x=np.median(all_indices), color='green', linestyle='--',
                    linewidth=2, label=f': {np.median(all_indices):.1f}')

        ax1.set_xlabel('PPG retention indices')
        ax1.set_ylabel('')
        ax1.set_title(f'PPGindex - {condition}')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # translated note
        ax2.boxplot(all_indices, vert=True, patch_artist=True,
                    boxprops=dict(facecolor='lightblue', color='blue'),
                    medianprops=dict(color='red', linewidth=2))

        # translated note
        stats_text = f":\n"
        stats_text += f": {len(all_indices)}\n"
        stats_text += f": {np.mean(all_indices):.1f}\n"
        stats_text += f": {np.std(all_indices):.1f}\n"
        stats_text += f": {np.min(all_indices):.1f}\n"
        stats_text += f": {np.max(all_indices):.1f}\n"
        stats_text += f": {np.median(all_indices):.1f}"

        ax2.text(0.05, 0.95, stats_text, transform=ax2.transAxes,
                 verticalalignment='top', fontsize=10,
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

        ax2.set_ylabel('PPG retention indices')
        ax2.set_title(f'PPGindex - {condition}')
        ax2.set_xticks([1])
        ax2.set_xticklabels([condition])
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"PPGindexsave: {save_path}")

        self.figures['index_distribution'] = fig
        return fig

    def plot_condition_comparison(self, conditions: List[str],
                                  save_path: str = None) -> plt.Figure:
        """
        conditionPPGindexcompare

        Parameters:
            conditions: comparechromatographic conditioncolumn
            save_path: savepath ()

        Returns:
            matplotlib Figure
        """
        # comparedata
        comparison_df = self.calculator.compare_conditions(conditions)

        if comparison_df.empty:
            print("Warning: comparedata")
            return None

        # conditionPPGindex
        fig, axes = plt.subplots(2, 2, figsize=(14, 12))
        axes = axes.flatten()

        # 1. compare (compare)
        if len(conditions) >= 2:
            ax = axes[0]
            cond1, cond2 = conditions[0], conditions[1]

            # conditioncompound
            compounds = set()
            for cond in [cond1, cond2]:
                col_name = f"{cond}_PPGindex"
                if col_name in comparison_df.columns:
                    valid_data = comparison_df[comparison_df[col_name].notna()]
                    compounds.update(valid_data['compound_name'].tolist())

            # Extract data
            data_points = []
            for compound in compounds:
                row = comparison_df[comparison_df['compound_name'] == compound]
                if not row.empty:
                    val1 = row[f"{cond1}_PPGindex"].iloc[0] if f"{cond1}_PPGindex" in row.columns else np.nan
                    val2 = row[f"{cond2}_PPGindex"].iloc[0] if f"{cond2}_PPGindex" in row.columns else np.nan
                    if not (np.isnan(val1) or np.isnan(val2)):
                        data_points.append((val1, val2))

            if data_points:
                val1_vals, val2_vals = zip(*data_points)
                ax.scatter(val1_vals, val2_vals, alpha=0.6, s=50)

                # translated note
                min_val = min(min(val1_vals), min(val2_vals))
                max_val = max(max(val1_vals), max(val2_vals))
                ax.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.5)

                ax.set_xlabel(f'{cond1} PPGindex')
                ax.set_ylabel(f'{cond2} PPGindex')
                ax.set_title(f'{cond1} vs {cond2} PPGindexcompare')
                ax.grid(True, alpha=0.3)

        # 2. conditioncompare
        ax = axes[1]
        box_data = []
        labels = []

        for cond in conditions:
            col_name = f"{cond}_PPGindex"
            if col_name in comparison_df.columns:
                data = comparison_df[col_name].dropna().tolist()
                if data:
                    box_data.append(data)
                    labels.append(cond)

        if box_data:
            bp = ax.boxplot(box_data, labels=labels, patch_artist=True)

            # translated note
            colors = ['lightblue', 'lightgreen', 'lightcoral', 'lightsalmon']
            for patch, color in zip(bp['boxes'], colors[:len(box_data)]):
                patch.set_facecolor(color)

            ax.set_ylabel('PPG retention indices')
            ax.set_title('conditionPPGindexcompare')
            ax.grid(True, alpha=0.3)

        # 3. (data)
        ax = axes[2]
        if len(conditions) >= 2:
            error_df = self.calculator.calculate_conversion_error(conditions[0], conditions[1])
            if not error_df.empty and '(%)' in error_df.columns:
                errors = error_df['(%)'].dropna()
                if len(errors) > 0:
                    ax.hist(errors, bins=20, color='coral', edgecolor='black', alpha=0.7)
                    ax.axvline(x=np.mean(errors), color='red', linestyle='--',
                               linewidth=2, label=f': {np.mean(errors):.2f}%')

                    ax.set_xlabel(' (%)')
                    ax.set_ylabel('')
                    ax.set_title(f'{conditions[0]} → {conditions[1]} ')
                    ax.legend()
                    ax.grid(True, alpha=0.3)

        # 4.
        ax = axes[3]
        if len(conditions) >= 2:
            # translated note
            corr_data = []
            for cond in conditions:
                col_name = f"{cond}_PPGindex"
                if col_name in comparison_df.columns:
                    corr_data.append(comparison_df[col_name])

            if corr_data and len(corr_data) > 1:
                corr_df = pd.DataFrame(corr_data).T
                corr_df.columns = conditions[:len(corr_data)]
                corr_matrix = corr_df.corr()

                # translated note
                im = ax.imshow(corr_matrix, cmap='RdYlBu', vmin=0, vmax=1)

                # translated note
                for i in range(len(corr_matrix)):
                    for j in range(len(corr_matrix)):
                        ax.text(j, i, f'{corr_matrix.iloc[i, j]:.3f}',
                                ha='center', va='center',
                                color='black' if abs(corr_matrix.iloc[i, j]) < 0.7 else 'white')

                ax.set_xticks(range(len(corr_matrix)))
                ax.set_yticks(range(len(corr_matrix)))
                ax.set_xticklabels(corr_matrix.columns)
                ax.set_yticklabels(corr_matrix.columns)
                ax.set_title('PPGindex')

                # translated note
                plt.colorbar(im, ax=ax)

        # Use
        for i in range(len(conditions), 4):
            axes[i].set_visible(False)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"conditioncomparesave: {save_path}")

        self.figures['condition_comparison'] = fig
        return fig

    def plot_conversion_analysis(self, from_condition: str, to_condition: str,
                                 save_path: str = None) -> plt.Figure:
        """
        conditionanalysis

        Parameters:
            from_condition: chromatographic condition
            to_condition: chromatographic condition
            save_path: savepath ()

        Returns:
            matplotlib Figure
        """
        try:
            # analysis
            conversion_df, stats = self.calculator.convert_ppg_index_to_rt(from_condition, to_condition)

            if conversion_df.empty or '(min)' not in conversion_df.columns:
                print("Warning: data")
                return None

            # data
            valid_df = conversion_df.dropna(
                subset=['(min)', f'{to_condition}_RT', f'{to_condition}_RT'])

            if valid_df.empty:
                print("Warning: data")
                return None

            # translated note
            fig, axes = plt.subplots(2, 3, figsize=(18, 12))
            axes = axes.flatten()

            # 1. RT vs RT
            ax = axes[0]
            predicted = valid_df[f'{to_condition}_RT']
            actual = valid_df[f'{to_condition}_RT']

            ax.scatter(actual, predicted, alpha=0.6, s=50, color='steelblue')

            # translated note
            min_val = min(predicted.min(), actual.min())
            max_val = max(predicted.max(), actual.max())
            ax.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.5, label='y=x')

            # translated note
            if len(predicted) > 1:
                slope, intercept, r_value, p_value, std_err = stats.linregress(actual, predicted)
                x_range = np.linspace(min_val, max_val, 100)
                y_pred = intercept + slope * x_range
                ax.plot(x_range, y_pred, 'g-', alpha=0.7,
                        label=f': R²={r_value ** 2:.3f}')

            ax.set_xlabel(f'{to_condition} retention_time (min)')
            ax.set_ylabel(f'{to_condition} retention_time (min)')
            ax.set_title(f'{from_condition}→{to_condition}: vs retention_time')
            ax.legend()
            ax.grid(True, alpha=0.3)

            # 2.
            ax = axes[1]
            errors = valid_df['(min)']

            ax.hist(errors, bins=20, color='coral', edgecolor='black', alpha=0.7)
            ax.axvline(x=errors.mean(), color='red', linestyle='--',
                       linewidth=2, label=f': {errors.mean():.3f}min')

            # translated note
            stats_text = f":\n"
            stats_text += f": {len(errors)}\n"
            stats_text += f": {errors.mean():.3f}min\n"
            stats_text += f": {errors.std():.3f}min\n"
            stats_text += f": {errors.max():.3f}min\n"
            stats_text += f": {errors.median():.3f}min"

            ax.text(0.65, 0.95, stats_text, transform=ax.transAxes,
                    verticalalignment='top', fontsize=10,
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

            ax.set_xlabel(' (min)')
            ax.set_ylabel('')
            ax.set_title(f'{from_condition}→{to_condition}: ')
            ax.legend()
            ax.grid(True, alpha=0.3)

            # 3.
            ax = axes[2]
            if '(%)' in valid_df.columns:
                rel_errors = valid_df['(%)'].dropna()
                if len(rel_errors) > 0:
                    ax.hist(rel_errors, bins=20, color='lightgreen', edgecolor='black', alpha=0.7)
                    ax.axvline(x=rel_errors.mean(), color='green', linestyle='--',
                               linewidth=2, label=f': {rel_errors.mean():.2f}%')

                    ax.set_xlabel(' (%)')
                    ax.set_ylabel('')
                    ax.set_title(f'{from_condition}→{to_condition}: ')
                    ax.legend()
                    ax.grid(True, alpha=0.3)

            # 4. vs RT
            ax = axes[3]
            ax.scatter(actual, errors, alpha=0.6, s=50, color='purple')
            ax.axhline(y=errors.mean(), color='red', linestyle='--',
                       linewidth=1, label=f': {errors.mean():.3f}min')

            # translated note
            if len(actual) > 1:
                z = np.polyfit(actual, errors, 1)
                p = np.poly1d(z)
                ax.plot(actual, p(actual), "b--", alpha=0.5, label='')

            ax.set_xlabel(f'{to_condition} retention_time (min)')
            ax.set_ylabel(' (min)')
            ax.set_title(f'{from_condition}→{to_condition}: vs retention_time')
            ax.legend()
            ax.grid(True, alpha=0.3)

            # 5. rank (compound)
            ax = axes[4]
            top_n = min(15, len(valid_df))
            top_errors = valid_df.nlargest(top_n, '(min)')

            y_pos = np.arange(top_n)
            ax.barh(y_pos, top_errors['(min)'], color='tomato', alpha=0.7)

            # compound_name
            compound_names = []
            for name in top_errors['compound_name']:
                if len(name) > 20:
                    compound_names.append(name[:17] + '...')
                else:
                    compound_names.append(name)

            ax.set_yticks(y_pos)
            ax.set_yticklabels(compound_names)
            ax.invert_yaxis() #
            ax.set_xlabel(' (min)')
            ax.set_title(f'{from_condition}→{to_condition}: {top_n}compound')
            ax.grid(True, alpha=0.3, axis='x')

            # 6. Bland-Altman (analysis)
            ax = axes[5]
            mean_values = (predicted + actual) / 2
            differences = predicted - actual

            ax.scatter(mean_values, differences, alpha=0.6, s=50, color='orange')
            ax.axhline(y=differences.mean(), color='red', linestyle='-',
                       linewidth=2, label=f': {differences.mean():.3f}')

            # 95%
            mean_diff = differences.mean()
            std_diff = differences.std()
            upper_limit = mean_diff + 1.96 * std_diff
            lower_limit = mean_diff - 1.96 * std_diff

            ax.axhline(y=upper_limit, color='red', linestyle='--',
                       linewidth=1, label=f'+1.96SD: {upper_limit:.3f}')
            ax.axhline(y=lower_limit, color='red', linestyle='--',
                       linewidth=1, label=f'-1.96SD: {lower_limit:.3f}')

            ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.5, alpha=0.5)

            ax.set_xlabel('retention_time (min)')
            ax.set_ylabel(' - (min)')
            ax.set_title(f'{from_condition}→{to_condition}: Bland-Altman')
            ax.legend(fontsize=9)
            ax.grid(True, alpha=0.3)

            plt.tight_layout()

            if save_path:
                plt.savefig(save_path, dpi=300, bbox_inches='tight')
                print(f"analysissave: {save_path}")

            self.figures['conversion_analysis'] = fig
            return fig

        except Exception as e:
            print(f"analysisfailed: {str(e)}")
            return None

    def plot_multiple_conversion_comparison(self, conversions: List[Tuple[str, str]],
                                            save_path: str = None) -> plt.Figure:
        """
        compare

        Parameters:
            conversions: column, (from_condition, to_condition)
            save_path: savepath ()

        Returns:
            matplotlib Figure
        """
        try:
            # data
            all_errors = []
            labels = []

            for from_cond, to_cond in conversions:
                conversion_df, _ = self.calculator.convert_ppg_index_to_rt(from_cond, to_cond)

                if not conversion_df.empty and '(min)' in conversion_df.columns:
                    errors = conversion_df['(min)'].dropna()
                    if len(errors) > 0:
                        all_errors.append(errors)
                        labels.append(f'{from_cond}→{to_cond}')

            if not all_errors:
                print("Warning: datacompare")
                return None

            # compare
            fig, axes = plt.subplots(2, 2, figsize=(14, 12))
            axes = axes.flatten()

            # 1. compare
            ax = axes[0]
            bp = ax.boxplot(all_errors, labels=labels, patch_artist=True, showfliers=False)

            # translated note
            colors = ['lightblue', 'lightgreen', 'lightcoral', 'lightsalmon', 'lightyellow']
            for patch, color in zip(bp['boxes'], colors[:len(all_errors)]):
                patch.set_facecolor(color)

            ax.set_ylabel(' (min)')
            ax.set_title('compare')
            ax.grid(True, alpha=0.3)

            # data
            for i, errors in enumerate(all_errors):
                ax.text(i + 1, errors.max() + 0.05, f'n={len(errors)}',
                        ha='center', va='bottom', fontsize=9)

            # 2.
            ax = axes[1]
            means = [err.mean() for err in all_errors]
            stds = [err.std() for err in all_errors]

            x_pos = np.arange(len(labels))
            bars = ax.bar(x_pos, means, yerr=stds, capsize=5, alpha=0.7, color='steelblue')

            # translated note
            for bar, mean in zip(bars, means):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width() / 2., height + 0.01,
                        f'{mean:.3f}', ha='center', va='bottom', fontsize=9)

            ax.set_xlabel('')
            ax.set_ylabel(' (min)')
            ax.set_title('compare')
            ax.set_xticks(x_pos)
            ax.set_xticklabels(labels, rotation=45, ha='right')
            ax.grid(True, alpha=0.3, axis='y')

            # 3. compare (0.5min)
            ax = axes[2]
            pass_rates = []
            for errors in all_errors:
                passed = len(errors[errors <= 0.5])
                pass_rate = (passed / len(errors) * 100) if len(errors) > 0 else 0
                pass_rates.append(pass_rate)

            bars = ax.bar(x_pos, pass_rates, alpha=0.7, color='lightgreen')

            # translated note
            for bar, rate in zip(bars, pass_rates):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width() / 2., height + 1,
                        f'{rate:.1f}%', ha='center', va='bottom', fontsize=9)

            ax.set_xlabel('')
            ax.set_ylabel(' (%)')
            ax.set_title('compare (: 0.5min)')
            ax.set_xticks(x_pos)
            ax.set_xticklabels(labels, rotation=45, ha='right')
            ax.set_ylim([0, 105])
            ax.grid(True, alpha=0.3, axis='y')

            # 4.
            ax = axes[3]
            for errors, label in zip(all_errors, labels):
                sorted_errors = np.sort(errors)
                y_vals = np.arange(1, len(sorted_errors) + 1) / len(sorted_errors) * 100

                ax.plot(sorted_errors, y_vals, marker='.', label=label, linewidth=2)

            ax.set_xlabel(' (min)')
            ax.set_ylabel(' (%)')
            ax.set_title('')
            ax.legend()
            ax.grid(True, alpha=0.3)

            plt.tight_layout()

            if save_path:
                plt.savefig(save_path, dpi=300, bbox_inches='tight')
                print(f"comparesave: {save_path}")

            self.figures['multiple_conversion_comparison'] = fig
            return fig

        except Exception as e:
            print(f"comparefailed: {str(e)}")
            return None


class PPGIndexAnalyzerGUI:
    """PPG retention indicesanalysisGUI"""

    def __init__(self, root):
        """GUI"""
        self.root = root
        self.root.title("PPG retention indicescalculateanalysis - ")
        self.root.geometry("1300x950")

        # translated note
        try:
            self.root.iconbitmap(default='icon.ico')
        except:
            pass

        # translated note
        self.calculator = PPGIndexCalculator()
        self.visualizer = None
        self.processing_thread = None
        self.is_processing = False

        # loadcondition
        self.loaded_conditions = set()
        self.loaded_compound_datasets = set()

        # UI
        self.setup_ui()

        # translated note
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)

    def setup_ui(self):
        """UI"""
        # translated note
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Notebook ()
        self.notebook = ttk.Notebook(main_frame)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        # translated note
        self.setup_data_tab()
        self.setup_analysis_tab()
        self.setup_conversion_tab()
        self.setup_visualization_tab()
        self.setup_results_tab()
        self.setup_log_tab()

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
        data_frame = ttk.LabelFrame(data_tab, text="Data loading and management", padding=15)
        data_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # ==================== PPG dataload ====================
        ppg_frame = ttk.LabelFrame(data_frame, text="PPG standard dataload", padding=10)
        ppg_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(ppg_frame, text="conditionname:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.condition_var = tk.StringVar(value="condition1")
        condition_entry = ttk.Entry(ppg_frame, textvariable=self.condition_var, width=20)
        condition_entry.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Label(ppg_frame, text="PPG datafile:").grid(row=1, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.ppg_file_var = tk.StringVar()
        ppg_file_entry = ttk.Entry(ppg_frame, textvariable=self.ppg_file_var, width=60)
        ppg_file_entry.grid(row=1, column=1, sticky=(tk.W, tk.E), padx=(0, 10), pady=5)

        ttk.Button(ppg_frame, text="...", command=self.browse_ppg_file).grid(row=1, column=2, pady=5)

        ttk.Button(ppg_frame, text="loadPPG data", command=self.load_ppg_data).grid(row=2, column=0, columnspan=3,
                                                                                   pady=10)

        # ==================== compound dataload ====================
        compound_frame = ttk.LabelFrame(data_frame, text="compound dataload", padding=10)
        compound_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(compound_frame, text="data:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.category_var = tk.StringVar(value="validation")
        category_combo = ttk.Combobox(compound_frame, textvariable=self.category_var,
                                      values=["validation", "smrt", "training", "test"],
                                      width=15, state="readonly")
        category_combo.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Label(compound_frame, text="conditionname:").grid(row=0, column=2, sticky=tk.W, padx=(0, 5), pady=5)
        self.compound_condition_var = tk.StringVar(value="condition1")
        compound_condition_entry = ttk.Entry(compound_frame, textvariable=self.compound_condition_var, width=20)
        compound_condition_entry.grid(row=0, column=3, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Label(compound_frame, text="compound datafile:").grid(row=1, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.compound_file_var = tk.StringVar()
        compound_file_entry = ttk.Entry(compound_frame, textvariable=self.compound_file_var, width=60)
        compound_file_entry.grid(row=1, column=1, columnspan=3, sticky=(tk.W, tk.E), padx=(0, 10), pady=5)

        ttk.Button(compound_frame, text="...", command=self.browse_compound_file).grid(row=1, column=4, pady=5)

        ttk.Button(compound_frame, text="loadcompound data", command=self.load_compound_data).grid(row=2, column=0,
                                                                                                columnspan=5, pady=10)

        # ==================== loaddata ====================
        overview_frame = ttk.LabelFrame(data_frame, text="loaddata", padding=10)
        overview_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 10))

        # Treeviewloaddata
        columns = ("", "condition", "data", "data", "")
        self.data_tree = ttk.Treeview(overview_frame, columns=columns, show="headings", height=8)

        # column
        for col in columns:
            self.data_tree.heading(col, text=col)
            self.data_tree.column(col, width=100)

        # translated note
        scrollbar = ttk.Scrollbar(overview_frame, orient=tk.VERTICAL, command=self.data_tree.yview)
        self.data_tree.configure(yscrollcommand=scrollbar.set)

        self.data_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # data
        button_frame = ttk.Frame(data_frame)
        button_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Button(button_frame, text="data", command=self.clear_all_data).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="PPG data", command=self.preview_ppg_data).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="compound data", command=self.preview_compound_data).pack(side=tk.LEFT, padx=5)

        # translated note
        self.data_status_var = tk.StringVar(value="loaddata...")
        ttk.Label(data_frame, textvariable=self.data_status_var).pack(anchor=tk.W)

    def setup_analysis_tab(self):
        """analysis"""
        analysis_tab = ttk.Frame(self.notebook)
        self.notebook.add(analysis_tab, text="dataanalysis")

        # translated note
        analysis_frame = ttk.LabelFrame(analysis_tab, text="PPG retention indicesanalysis", padding=15)
        analysis_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # ==================== standard curve ====================
        curve_frame = ttk.LabelFrame(analysis_frame, text="PPGstandard curve", padding=10)
        curve_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(curve_frame, text="condition:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.curve_condition_var = tk.StringVar()
        self.curve_condition_combo = ttk.Combobox(curve_frame, textvariable=self.curve_condition_var,
                                                  width=25, state="readonly")
        self.curve_condition_combo.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Label(curve_frame, text="model:").grid(row=0, column=2, sticky=tk.W, padx=(0, 5), pady=5)
        self.model_type_var = tk.StringVar(value="logarithmic")
        model_combo = ttk.Combobox(curve_frame, textvariable=self.model_type_var,
                                   values=["logarithmic", "linear"],
                                   width=15, state="readonly")
        model_combo.grid(row=0, column=3, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Button(curve_frame, text="standard curve", command=self.fit_standard_curve).grid(row=0, column=4,
                                                                                           padx=(20, 0), pady=5)

        # ==================== PPGindexcalculate ====================
        calc_frame = ttk.LabelFrame(analysis_frame, text="PPGindexcalculate", padding=10)
        calc_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(calc_frame, text="condition:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.calc_condition_var = tk.StringVar()
        self.calc_condition_combo = ttk.Combobox(calc_frame, textvariable=self.calc_condition_var,
                                                 width=25, state="readonly")
        self.calc_condition_combo.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Label(calc_frame, text="calculatemethod:").grid(row=0, column=2, sticky=tk.W, padx=(0, 5), pady=5)
        self.calc_method_var = tk.StringVar(value="interpolation")
        method_combo = ttk.Combobox(calc_frame, textvariable=self.calc_method_var,
                                    values=["interpolation", "regression"],
                                    width=15, state="readonly")
        method_combo.grid(row=0, column=3, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Button(calc_frame, text="calculatePPGindex", command=self.calculate_ppg_index).grid(row=0, column=4, padx=(20, 0),
                                                                                          pady=5)

        # ==================== conditioncompare ====================
        compare_frame = ttk.LabelFrame(analysis_frame, text="conditioncompareanalysis", padding=10)
        compare_frame.pack(fill=tk.X, pady=(0, 15))

        # condition
        ttk.Label(compare_frame, text="comparecondition:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)

        # condition
        self.condition_checkboxes = {}
        self.condition_checkboxes_frame = ttk.Frame(compare_frame)
        self.condition_checkboxes_frame.grid(row=1, column=0, columnspan=5, sticky=tk.W, padx=5, pady=5)

        # compare
        ttk.Button(compare_frame, text="comparecondition", command=self.compare_conditions).grid(row=2, column=0,
                                                                                             sticky=tk.W, padx=5,
                                                                                             pady=10)

        # ==================== analysis results ====================
        results_frame = ttk.LabelFrame(analysis_frame, text="analysis results", padding=10)
        results_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 10))

        # translated note
        self.analysis_text = scrolledtext.ScrolledText(results_frame, width=80, height=15,
                                                       wrap=tk.WORD, font=("Consolas", 10))
        self.analysis_text.pack(fill=tk.BOTH, expand=True)

        # translated note
        self.analysis_text.tag_config("INFO", foreground="black")
        self.analysis_text.tag_config("SUCCESS", foreground="green")
        self.analysis_text.tag_config("WARNING", foreground="orange")
        self.analysis_text.tag_config("ERROR", foreground="red")

        # analysis
        button_frame = ttk.Frame(analysis_frame)
        button_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Button(button_frame, text="analysis", command=self.generate_report).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="results", command=self.clear_analysis_text).pack(side=tk.LEFT, padx=5)

        # analysis
        self.analysis_status_var = tk.StringVar(value="analysis...")
        ttk.Label(analysis_frame, textvariable=self.analysis_status_var).pack(anchor=tk.W)

    def setup_conversion_tab(self):
        """condition"""
        conversion_tab = ttk.Frame(self.notebook)
        self.notebook.add(conversion_tab, text="condition")

        # translated note
        conversion_frame = ttk.LabelFrame(conversion_tab, text="conditionPPGindexanalysis", padding=15)
        conversion_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # ==================== ====================
        settings_frame = ttk.LabelFrame(conversion_frame, text="", padding=10)
        settings_frame.pack(fill=tk.X, pady=(0, 15))

        # conditioncondition
        ttk.Label(settings_frame, text="condition (PPGindex):").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.from_condition_var = tk.StringVar()
        self.from_condition_combo = ttk.Combobox(settings_frame, textvariable=self.from_condition_var,
                                                 width=25, state="readonly")
        self.from_condition_combo.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Label(settings_frame, text="condition ():").grid(row=0, column=2, sticky=tk.W, padx=(0, 5), pady=5)
        self.to_condition_var = tk.StringVar()
        self.to_condition_combo = ttk.Combobox(settings_frame, textvariable=self.to_condition_var,
                                               width=25, state="readonly")
        self.to_condition_combo.grid(row=0, column=3, sticky=tk.W, padx=(0, 10), pady=5)

        # translated note
        ttk.Label(settings_frame, text="(min):").grid(row=1, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.threshold_var = tk.DoubleVar(value=0.5)
        threshold_entry = ttk.Entry(settings_frame, textvariable=self.threshold_var, width=10)
        threshold_entry.grid(row=1, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        # translated note
        ttk.Button(settings_frame, text="analysis", command=self.perform_conversion_analysis).grid(
            row=0, column=4, rowspan=2, padx=(20, 0), pady=5, sticky=tk.NS)

        # compare
        multi_compare_frame = ttk.LabelFrame(settings_frame, text="compare", padding=5)
        multi_compare_frame.grid(row=2, column=0, columnspan=5, sticky=tk.W, padx=5, pady=10)

        ttk.Label(multi_compare_frame, text=":").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=2)

        # translated note
        self.multi_conversion_frame = ttk.Frame(multi_compare_frame)
        self.multi_conversion_frame.grid(row=1, column=0, columnspan=3, sticky=tk.W, padx=5, pady=2)

        ttk.Button(multi_compare_frame, text="compare", command=self.compare_multiple_conversions).grid(
            row=1, column=3, padx=(10, 0), pady=2)

        # ==================== results ====================
        results_frame = ttk.LabelFrame(conversion_frame, text="results", padding=10)
        results_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 15))

        # Treeviewresults
        columns = ("compound_name", "PPGindex", "RT", "RT", "", "(%)")
        self.conversion_tree = ttk.Treeview(results_frame, columns=columns, show="headings", height=10)

        # column
        for col in columns:
            self.conversion_tree.heading(col, text=col)
            if col == "compound_name":
                self.conversion_tree.column(col, width=150)
            else:
                self.conversion_tree.column(col, width=100)

        # translated note
        scrollbar = ttk.Scrollbar(results_frame, orient=tk.VERTICAL, command=self.conversion_tree.yview)
        self.conversion_tree.configure(yscrollcommand=scrollbar.set)

        self.conversion_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # translated note
        stats_frame = ttk.LabelFrame(conversion_frame, text="", padding=10)
        stats_frame.pack(fill=tk.X, pady=(0, 10))

        self.conversion_stats_text = scrolledtext.ScrolledText(stats_frame, width=80, height=6,
                                                               wrap=tk.WORD, font=("Consolas", 9))
        self.conversion_stats_text.pack(fill=tk.BOTH, expand=True)

        # translated note
        self.conversion_stats_text.tag_config("INFO", foreground="black")
        self.conversion_stats_text.tag_config("SUCCESS", foreground="green")
        self.conversion_stats_text.tag_config("WARNING", foreground="orange")
        self.conversion_stats_text.tag_config("ERROR", foreground="red")

        # translated note
        button_frame = ttk.Frame(conversion_frame)
        button_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Button(button_frame, text="results", command=self.export_conversion_results).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="results", command=self.clear_conversion_results).pack(side=tk.LEFT, padx=5)

        # translated note
        self.conversion_status_var = tk.StringVar(value="analysis...")
        ttk.Label(conversion_frame, textvariable=self.conversion_status_var).pack(anchor=tk.W)

    def setup_visualization_tab(self):
        """"""
        viz_tab = ttk.Frame(self.notebook)
        self.notebook.add(viz_tab, text="data")

        # translated note
        viz_frame = ttk.LabelFrame(viz_tab, text="data", padding=15)
        viz_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # ==================== ====================
        options_frame = ttk.LabelFrame(viz_frame, text="", padding=10)
        options_frame.pack(fill=tk.X, pady=(0, 15))

        # translated note
        ttk.Label(options_frame, text=":").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.viz_type_var = tk.StringVar(value="standard_curve")
        viz_type_combo = ttk.Combobox(options_frame, textvariable=self.viz_type_var,
                                      values=["standard_curve", "residuals", "index_distribution",
                                              "condition_comparison", "conversion_analysis",
                                              "multiple_conversion_comparison"],
                                      width=25, state="readonly")
        viz_type_combo.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)
        viz_type_combo.bind("<<ComboboxSelected>>", self.on_viz_type_change)

        # condition ()
        ttk.Label(options_frame, text="condition:").grid(row=0, column=2, sticky=tk.W, padx=(0, 5), pady=5)
        self.viz_condition_var = tk.StringVar()
        self.viz_condition_combo = ttk.Combobox(options_frame, textvariable=self.viz_condition_var,
                                                width=25, state="readonly")
        self.viz_condition_combo.grid(row=0, column=3, sticky=tk.W, padx=(0, 10), pady=5)

        # condition (conditioncompareanalysis)
        self.viz_conditions_frame = ttk.Frame(options_frame)
        self.viz_conditions_frame.grid(row=1, column=0, columnspan=4, sticky=tk.W, padx=5, pady=5)
        self.viz_conditions_frame.grid_remove() #

        # analysisParameters
        self.conversion_viz_frame = ttk.Frame(options_frame)
        self.conversion_viz_frame.grid(row=1, column=0, columnspan=4, sticky=tk.W, padx=5, pady=5)
        self.conversion_viz_frame.grid_remove() #

        # chart
        ttk.Button(options_frame, text="chart", command=self.generate_visualization).grid(row=0, column=4,
                                                                                             padx=(20, 0), pady=5)

        # ==================== chart ====================
        display_frame = ttk.LabelFrame(viz_frame, text="chart", padding=10)
        display_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 10))

        # translated note
        self.figure_canvas = None
        self.figure_toolbar = None
        self.current_figure = None

        # translated note
        self.viz_placeholder = ttk.Label(display_frame, text="chart",
                                         font=("Arial", 14), foreground="gray")
        self.viz_placeholder.pack(expand=True)

        # translated note
        button_frame = ttk.Frame(viz_frame)
        button_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Button(button_frame, text="savechart", command=self.save_figure).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="chart", command=self.clear_figure).pack(side=tk.LEFT, padx=5)

        # translated note
        self.viz_status_var = tk.StringVar(value="chart...")
        ttk.Label(viz_frame, textvariable=self.viz_status_var).pack(anchor=tk.W)

    def setup_results_tab(self):
        """results"""
        results_tab = ttk.Frame(self.notebook)
        self.notebook.add(results_tab, text="resultsoutput")

        # translated note
        results_frame = ttk.LabelFrame(results_tab, text="resultsoutput", padding=15)
        results_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # ==================== output ====================
        output_frame = ttk.LabelFrame(results_frame, text="output", padding=10)
        output_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(output_frame, text="outputdirectory:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.output_dir_var = tk.StringVar(value=os.getcwd())
        output_dir_entry = ttk.Entry(output_frame, textvariable=self.output_dir_var, width=60)
        output_dir_entry.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=(0, 10), pady=5)

        ttk.Button(output_frame, text="...", command=self.browse_output_dir).grid(row=0, column=2, pady=5)

        ttk.Label(output_frame, text="outputfile:").grid(row=1, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.output_prefix_var = tk.StringVar(value="PPG_Analysis")
        ttk.Entry(output_frame, textvariable=self.output_prefix_var, width=30).grid(row=1, column=1, sticky=tk.W,
                                                                                    padx=(0, 10), pady=5)

        # ==================== results ====================
        preview_frame = ttk.LabelFrame(results_frame, text="results", padding=10)
        preview_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 15))

        # translated note
        self.results_text = scrolledtext.ScrolledText(preview_frame, width=80, height=15,
                                                      wrap=tk.WORD, font=("Consolas", 10))
        self.results_text.pack(fill=tk.BOTH, expand=True)

        # ==================== output ====================
        control_frame = ttk.Frame(results_frame)
        control_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Button(control_frame, text="saveresults", command=self.save_all_results).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_frame, text="", command=self.preview_report).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_frame, text="outputdirectory", command=self.open_output_dir).pack(side=tk.LEFT, padx=5)

        # output
        self.output_status_var = tk.StringVar(value="outputresults...")
        ttk.Label(results_frame, textvariable=self.output_status_var).pack(anchor=tk.W)

    def setup_log_tab(self):
        """"""
        log_tab = ttk.Frame(self.notebook)
        self.notebook.add(log_tab, text="")

        # translated note
        log_frame = ttk.LabelFrame(log_tab, text="", padding=15)
        log_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # translated note
        self.log_text = scrolledtext.ScrolledText(log_frame, width=100, height=25,
                                                  wrap=tk.WORD, font=("Consolas", 10))
        self.log_text.pack(fill=tk.BOTH, expand=True)

        # translated note
        self.log_text.tag_config("INFO", foreground="black")
        self.log_text.tag_config("SUCCESS", foreground="green")
        self.log_text.tag_config("WARNING", foreground="orange")
        self.log_text.tag_config("ERROR", foreground="red")

        # translated note
        button_frame = ttk.Frame(log_frame)
        button_frame.pack(fill=tk.X, pady=(10, 0))

        ttk.Button(button_frame, text="", command=self.clear_log).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="save", command=self.save_log).pack(side=tk.LEFT, padx=5)

    def browse_ppg_file(self):
        """PPG datafile"""
        file_types = [("datafile", "*.csv *.xlsx *.xls"), ("CSVfile", "*.csv"),
                      ("Excelfile", "*.xlsx *.xls"), ("file", "*.*")]
        file_path = filedialog.askopenfilename(title="PPG datafile", filetypes=file_types)

        if file_path:
            self.ppg_file_var.set(file_path)

    def browse_compound_file(self):
        """compound datafile"""
        file_types = [("datafile", "*.csv *.xlsx *.xls"), ("CSVfile", "*.csv"),
                      ("Excelfile", "*.xlsx *.xls"), ("file", "*.*")]
        file_path = filedialog.askopenfilename(title="compound datafile", filetypes=file_types)

        if file_path:
            self.compound_file_var.set(file_path)

    def browse_output_dir(self):
        """outputdirectory"""
        dir_path = filedialog.askdirectory(title="outputdirectory")

        if dir_path:
            self.output_dir_var.set(dir_path)

    def load_ppg_data(self):
        """loadPPG data"""
        ppg_file = self.ppg_file_var.get().strip()
        condition = self.condition_var.get().strip()

        if not ppg_file:
            messagebox.showwarning("Warning", "PPG datafile")
            return

        if not condition:
            messagebox.showwarning("Warning", "inputconditionname")
            return

        # load
        self.log_message(f"loadPPG data: {ppg_file} (condition: {condition})", "INFO")
        self.update_status(f"loadPPG data: {Path(ppg_file).name}")

        # loaddata
        self.processing_thread = threading.Thread(
            target=self._load_ppg_data_thread,
            args=(ppg_file, condition)
        )
        self.processing_thread.start()

    def _load_ppg_data_thread(self, ppg_file, condition):
        """loadPPG data"""
        try:
            success, msg = self.calculator.load_ppg_data(ppg_file, condition)

            if success:
                self.root.after(0, lambda: self.log_message(f"✓ {msg}", "SUCCESS"))
                self.root.after(0, self.update_data_tree)
                self.root.after(0, lambda: self.update_condition_comboboxes())
                self.loaded_conditions.add(condition)
            else:
                self.root.after(0, lambda: self.log_message(f"✗ {msg}", "ERROR"))

            self.root.after(0, lambda: self.update_status(""))

        except Exception as e:
            self.root.after(0, lambda: self.log_message(f"✗ loadPPG datafailed: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("loadfailed"))

    def load_compound_data(self):
        """loadcompound data"""
        compound_file = self.compound_file_var.get().strip()
        category = self.category_var.get().strip()
        condition = self.compound_condition_var.get().strip()

        if not compound_file:
            messagebox.showwarning("Warning", "compound datafile")
            return

        if not category:
            messagebox.showwarning("Warning", "data")
            return

        if not condition:
            messagebox.showwarning("Warning", "inputconditionname")
            return

        # load
        self.log_message(f"loadcompound data: {compound_file} (: {category}, condition: {condition})", "INFO")
        self.update_status(f"loadcompound data: {Path(compound_file).name}")

        # loaddata
        self.processing_thread = threading.Thread(
            target=self._load_compound_data_thread,
            args=(compound_file, category, condition)
        )
        self.processing_thread.start()

    def _load_compound_data_thread(self, compound_file, category, condition):
        """loadcompound data"""
        try:
            success, msg = self.calculator.load_compound_data(compound_file, category, condition)

            if success:
                self.root.after(0, lambda: self.log_message(f"✓ {msg}", "SUCCESS"))
                self.root.after(0, self.update_data_tree)
                self.root.after(0, lambda: self.update_condition_comboboxes())
                self.loaded_compound_datasets.add(f"{category}_{condition}")
            else:
                self.root.after(0, lambda: self.log_message(f"✗ {msg}", "ERROR"))

            self.root.after(0, lambda: self.update_status(""))

        except Exception as e:
            self.root.after(0, lambda: self.log_message(f"✗ loadcompound datafailed: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("loadfailed"))

    def update_data_tree(self):
        """data"""
        # data
        for item in self.data_tree.get_children():
            self.data_tree.delete(item)

        # PPG data
        for condition, df in self.calculator.ppg_data.items():
            self.data_tree.insert("", tk.END, values=("PPGstandard", condition, "", len(df), "load"))

        # compound data
        for key, df in self.calculator.compound_data.items():
            parts = key.split('_')
            if len(parts) >= 2:
                category, condition = parts[0], parts[1]
                self.data_tree.insert("", tk.END, values=("compound", condition, category, len(df), "load"))

        # translated note
        self.data_status_var.set(
            f"load {len(self.calculator.ppg_data)} PPG data, {len(self.calculator.compound_data)} compound data")

    def update_condition_comboboxes(self):
        """condition"""
        conditions = list(self.calculator.ppg_data.keys())

        # condition
        self.curve_condition_combo['values'] = conditions
        if conditions and not self.curve_condition_var.get():
            self.curve_condition_var.set(conditions[0])

        # calculatecondition
        self.calc_condition_combo['values'] = conditions
        if conditions and not self.calc_condition_var.get():
            self.calc_condition_var.set(conditions[0])

        # condition
        self.viz_condition_combo['values'] = conditions
        if conditions and not self.viz_condition_var.get():
            self.viz_condition_var.set(conditions[0])

        # condition
        self.from_condition_combo['values'] = conditions
        self.to_condition_combo['values'] = conditions
        if conditions:
            if not self.from_condition_var.get():
                self.from_condition_var.set(conditions[0])
            if not self.to_condition_var.get():
                self.to_condition_var.set(conditions[-1] if len(conditions) > 1 else conditions[0])

        # conditioncompare
        self.update_condition_checkboxes(conditions)

        # translated note
        self.update_multi_conversion_checkboxes(conditions)

    def update_condition_checkboxes(self, conditions):
        """conditioncompare"""
        # translated note
        for widget in self.condition_checkboxes_frame.winfo_children():
            widget.destroy()
        self.condition_checkboxes.clear()

        # translated note
        for i, condition in enumerate(conditions):
            var = tk.BooleanVar(value=(i < 2)) # condition
            cb = ttk.Checkbutton(self.condition_checkboxes_frame, text=condition, variable=var)
            cb.grid(row=i // 4, column=i % 4, sticky=tk.W, padx=5, pady=2)
            self.condition_checkboxes[condition] = var

    def update_multi_conversion_checkboxes(self, conditions):
        """"""
        # translated note
        for widget in self.multi_conversion_frame.winfo_children():
            widget.destroy()

        # translated note
        self.conversion_scheme_vars = []

        if len(conditions) >= 2:
            for i in range(len(conditions)):
                for j in range(len(conditions)):
                    if i != j:
                        from_cond = conditions[i]
                        to_cond = conditions[j]
                        var = tk.BooleanVar(value=(i == 0 and j == 1)) #
                        cb = ttk.Checkbutton(self.multi_conversion_frame,
                                             text=f"{from_cond}→{to_cond}",
                                             variable=var)
                        row = (i * len(conditions) + j) // 4
                        col = (i * len(conditions) + j) % 4
                        cb.grid(row=row, column=col, sticky=tk.W, padx=2, pady=1)
                        self.conversion_scheme_vars.append((var, from_cond, to_cond))
        else:
            ttk.Label(self.multi_conversion_frame, text="2conditioncompare",
                      foreground="gray").grid(row=0, column=0, sticky=tk.W)

    def clear_all_data(self):
        """data"""
        if messagebox.askyesno("", "loaddata？"):
            self.calculator = PPGIndexCalculator()
            self.visualizer = None
            self.loaded_conditions.clear()
            self.loaded_compound_datasets.clear()
            self.update_data_tree()
            self.update_condition_comboboxes()
            self.clear_conversion_results()
            self.log_message("data", "INFO")

    def preview_ppg_data(self):
        """PPG data"""
        conditions = list(self.calculator.ppg_data.keys())
        if not conditions:
            messagebox.showinfo("", "loadPPG data")
            return

        # translated note
        preview_window = Toplevel(self.root)
        preview_window.title("PPG data")
        preview_window.geometry("800x600")

        # Notebook
        notebook = ttk.Notebook(preview_window)
        notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # condition
        for condition in conditions:
            df = self.calculator.ppg_data[condition]

            # translated note
            frame = ttk.Frame(notebook)
            notebook.add(frame, text=condition)

            # Treeview
            tree_frame = ttk.Frame(frame)
            tree_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

            # translated note
            scrollbar_y = ttk.Scrollbar(tree_frame)
            scrollbar_y.pack(side=tk.RIGHT, fill=tk.Y)

            scrollbar_x = ttk.Scrollbar(tree_frame, orient=tk.HORIZONTAL)
            scrollbar_x.pack(side=tk.BOTTOM, fill=tk.X)

            # Treeview
            tree = ttk.Treeview(tree_frame, yscrollcommand=scrollbar_y.set,
                                xscrollcommand=scrollbar_x.set)
            tree.pack(fill=tk.BOTH, expand=True)

            scrollbar_y.config(command=tree.yview)
            scrollbar_x.config(command=tree.xview)

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
            info_frame = ttk.Frame(frame)
            info_frame.pack(fill=tk.X, padx=5, pady=(0, 5))

            ttk.Label(info_frame, text=f"data: {len(df)}").pack(side=tk.LEFT, padx=10)

            if 'degree_of_polymerization' in df.columns and 'retention_time' in df.columns:
                n_range = f"degree_of_polymerization: {df['degree_of_polymerization'].min()} - {df['degree_of_polymerization'].max()}"
                rt_range = f"retention_time: {df['retention_time'].min():.2f} - {df['retention_time'].max():.2f}"
                ttk.Label(info_frame, text=n_range).pack(side=tk.LEFT, padx=10)
                ttk.Label(info_frame, text=rt_range).pack(side=tk.LEFT, padx=10)

    def preview_compound_data(self):
        """compound data"""
        if not self.calculator.compound_data:
            messagebox.showinfo("", "loadcompound data")
            return

        # translated note
        preview_window = Toplevel(self.root)
        preview_window.title("compound data")
        preview_window.geometry("900x700")

        # Notebook
        notebook = ttk.Notebook(preview_window)
        notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # data
        for key, df in self.calculator.compound_data.items():
            # translated note
            frame = ttk.Frame(notebook)
            notebook.add(frame, text=key)

            # Treeview
            tree_frame = ttk.Frame(frame)
            tree_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

            # translated note
            scrollbar_y = ttk.Scrollbar(tree_frame)
            scrollbar_y.pack(side=tk.RIGHT, fill=tk.Y)

            scrollbar_x = ttk.Scrollbar(tree_frame, orient=tk.HORIZONTAL)
            scrollbar_x.pack(side=tk.BOTTOM, fill=tk.X)

            # Treeview
            tree = ttk.Treeview(tree_frame, yscrollcommand=scrollbar_y.set,
                                xscrollcommand=scrollbar_x.set)
            tree.pack(fill=tk.BOTH, expand=True)

            scrollbar_y.config(command=tree.yview)
            scrollbar_x.config(command=tree.xview)

            # column
            tree["columns"] = list(df.columns)
            tree["show"] = "headings"

            # column
            for col in df.columns:
                tree.heading(col, text=col)
                tree.column(col, width=120 if len(col) > 10 else 100)

            # data
            for _, row in df.iterrows():
                tree.insert("", tk.END, values=list(row))

            # translated note
            info_frame = ttk.Frame(frame)
            info_frame.pack(fill=tk.X, padx=5, pady=(0, 5))

            ttk.Label(info_frame, text=f"compound: {len(df)}").pack(side=tk.LEFT, padx=10)

            if 'retention_time' in df.columns:
                rt_range = f"retention_time: {df['retention_time'].min():.2f} - {df['retention_time'].max():.2f}"
                ttk.Label(info_frame, text=rt_range).pack(side=tk.LEFT, padx=10)

    def fit_standard_curve(self):
        """standard curve"""
        condition = self.curve_condition_var.get()
        model_type = self.model_type_var.get()

        if not condition:
            messagebox.showwarning("Warning", "condition")
            return

        # analysis
        self.analysis_message(f"standard curve (condition: {condition}, model: {model_type})", "INFO")
        self.update_status(f"standard curve: {condition}")

        # translated note
        self.processing_thread = threading.Thread(
            target=self._fit_standard_curve_thread,
            args=(condition, model_type)
        )
        self.processing_thread.start()

    def _fit_standard_curve_thread(self, condition, model_type):
        """standard curve"""
        try:
            success, msg = self.calculator.fit_standard_curve(condition, model_type)

            if success:
                self.root.after(0, lambda: self.analysis_message(f"✓ {msg}", "SUCCESS"))

                # results
                curve = self.calculator.standard_curves[condition]
                result_text = f"standard curveresults - {condition}:\n"
                result_text += f" model: {curve['model_name']}\n"
                result_text += f"  R²: {curve['r_squared']:.6f}\n"
                result_text += f" : {curve['slope']:.6f}\n"
                result_text += f" : {curve['intercept']:.6f}\n"
                result_text += f" : {curve['std_err']:.6f}\n"
                result_text += f" p: {curve['p_value']:.6f}\n"
                result_text += f" data: {curve['n_points']}\n"

                self.root.after(0, lambda: self.analysis_message(result_text, "INFO"))

                # translated note
                self.visualizer = PPGVisualizer(self.calculator)
            else:
                self.root.after(0, lambda: self.analysis_message(f"✗ {msg}", "ERROR"))

            self.root.after(0, lambda: self.update_status(""))

        except Exception as e:
            self.root.after(0, lambda: self.analysis_message(f"✗ standard curvefailed: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("failed"))

    def calculate_ppg_index(self):
        """calculatePPGindex"""
        condition = self.calc_condition_var.get()
        method = self.calc_method_var.get()

        if not condition:
            messagebox.showwarning("Warning", "condition")
            return

        # analysis
        self.analysis_message(f"calculatePPGindex (condition: {condition}, method: {method})", "INFO")
        self.update_status(f"calculatePPGindex: {condition}")

        # calculate
        self.processing_thread = threading.Thread(
            target=self._calculate_ppg_index_thread,
            args=(condition, method)
        )
        self.processing_thread.start()

    def _calculate_ppg_index_thread(self, condition, method):
        """calculatePPGindex"""
        try:
            success, msg = self.calculator.calculate_ppg_index(condition, method)

            if success:
                self.root.after(0, lambda: self.analysis_message(f"✓ {msg}", "SUCCESS"))

                # results
                if condition in self.calculator.ppg_indices:
                    indices_data = self.calculator.ppg_indices[condition]
                    all_indices = []

                    for key, df in indices_data['indices'].items():
                        if 'calculatePPGindex' in df.columns:
                            indices = df['calculatePPGindex'].dropna().tolist()
                            all_indices.extend(indices)

                    if all_indices:
                        indices_array = np.array(all_indices)
                        result_text = f"PPGindexresults - {condition}:\n"
                        result_text += f" calculatemethod: {method}\n"
                        result_text += f" compound: {len(all_indices)}\n"
                        result_text += f" : {np.mean(indices_array):.2f}\n"
                        result_text += f" : {np.std(indices_array):.2f}\n"
                        result_text += f" : {np.min(indices_array):.2f}\n"
                        result_text += f" : {np.max(indices_array):.2f}\n"
                        result_text += f" : {np.median(indices_array):.2f}\n"

                        self.root.after(0, lambda: self.analysis_message(result_text, "INFO"))
            else:
                self.root.after(0, lambda: self.analysis_message(f"✗ {msg}", "ERROR"))

            self.root.after(0, lambda: self.update_status(""))

        except Exception as e:
            self.root.after(0, lambda: self.analysis_message(f"✗ calculatePPGindexfailed: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("calculatefailed"))

    def compare_conditions(self):
        """comparecondition"""
        # condition
        selected_conditions = []
        for condition, var in self.condition_checkboxes.items():
            if var.get():
                selected_conditions.append(condition)

        if len(selected_conditions) < 2:
            messagebox.showwarning("Warning", "conditioncompare")
            return

        # analysis
        self.analysis_message(f"comparecondition: {', '.join(selected_conditions)}", "INFO")
        self.update_status(f"comparecondition")

        # compare
        self.processing_thread = threading.Thread(
            target=self._compare_conditions_thread,
            args=(selected_conditions,)
        )
        self.processing_thread.start()

    def _compare_conditions_thread(self, conditions):
        """comparecondition"""
        try:
            comparison_df = self.calculator.compare_conditions(conditions)

            if not comparison_df.empty:
                self.root.after(0, lambda: self.analysis_message(f"✓ conditioncompare, compare {len(comparison_df)} compound",
                                                                 "SUCCESS"))

                # compareresults
                result_text = f"conditioncompareresults:\n"
                result_text += f" comparecondition: {', '.join(conditions)}\n"
                result_text += f" compound: {len(comparison_df)}\n\n"

                # calculateconditionPPGindex
                for condition in conditions:
                    col_name = f"{condition}_PPGindex"
                    if col_name in comparison_df.columns:
                        data = comparison_df[col_name].dropna()
                        if len(data) > 0:
                            result_text += f"  {condition}:\n"
                            result_text += f" data: {len(data)}\n"
                            result_text += f" : {np.mean(data):.2f}\n"
                            result_text += f" : {np.std(data):.2f}\n"

                # calculatecondition
                if len(conditions) >= 2:
                    result_text += f"\n condition:\n"

                    # Extract data
                    corr_data = {}
                    for condition in conditions:
                        col_name = f"{condition}_PPGindex"
                        if col_name in comparison_df.columns:
                            corr_data[condition] = comparison_df[col_name]

                    if len(corr_data) >= 2:
                        corr_df = pd.DataFrame(corr_data)
                        corr_matrix = corr_df.corr()

                        for i, cond1 in enumerate(conditions):
                            for j, cond2 in enumerate(conditions):
                                if i < j and cond1 in corr_matrix.columns and cond2 in corr_matrix.columns:
                                    corr_value = corr_matrix.loc[cond1, cond2]
                                    result_text += f"    {cond1} vs {cond2}: {corr_value:.4f}\n"

                self.root.after(0, lambda: self.analysis_message(result_text, "INFO"))
            else:
                self.root.after(0, lambda: self.analysis_message("✗ comparedata", "WARNING"))

            self.root.after(0, lambda: self.update_status(""))

        except Exception as e:
            self.root.after(0, lambda: self.analysis_message(f"✗ conditioncomparefailed: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("comparefailed"))

    def perform_conversion_analysis(self):
        """conditionanalysis"""
        from_condition = self.from_condition_var.get()
        to_condition = self.to_condition_var.get()
        threshold = self.threshold_var.get()

        if not from_condition or not to_condition:
            messagebox.showwarning("Warning", "conditioncondition")
            return

        if from_condition == to_condition:
            messagebox.showwarning("Warning", "conditioncondition")
            return

        # translated note
        self.conversion_status_var.set(f": {from_condition} → {to_condition}")
        self.update_status(f"analysis")

        # translated note
        self.processing_thread = threading.Thread(
            target=self._perform_conversion_analysis_thread,
            args=(from_condition, to_condition, threshold)
        )
        self.processing_thread.start()

    def _perform_conversion_analysis_thread(self, from_condition, to_condition, threshold):
        """analysis"""
        try:
            # analysis
            analysis_results = self.calculator.cross_condition_analysis(from_condition, to_condition, threshold)

            if "error" in analysis_results:
                self.root.after(0, lambda: self.conversion_status_var.set("failed"))
                self.root.after(0, lambda: self.log_message(f"✗ analysisfailed: {analysis_results['error']}", "ERROR"))
                return

            # results
            self.root.after(0, lambda: self.display_conversion_results(analysis_results))

            # translated note
            self.root.after(0, lambda: self.log_message(
                f"✓ analysis: {from_condition} → {to_condition}, " +
                f": {analysis_results['']}, " +
                f": {analysis_results[''].get('', 0):.3f}min",
                "SUCCESS"))

            self.root.after(0, lambda: self.conversion_status_var.set("analysis"))
            self.root.after(0, lambda: self.update_status(""))

        except Exception as e:
            self.root.after(0, lambda: self.conversion_status_var.set("failed"))
            self.root.after(0, lambda: self.log_message(f"✗ analysisfailed: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("failed"))

    def display_conversion_results(self, analysis_results):
        """results"""
        # results
        for item in self.conversion_tree.get_children():
            self.conversion_tree.delete(item)

        # translated note
        self.conversion_stats_text.delete(1.0, tk.END)

        # data
        if 'data' in analysis_results:
            for item in analysis_results['data']:
                compound_name = item['compound_name']
                source_ppg = item.get(f"{analysis_results['condition']}_PPGindex", "")
                pred_rt = item.get(f"{analysis_results['condition']}_RT", "")
                actual_rt = item.get(f"{analysis_results['condition']}_RT", "")
                abs_error = item.get('(min)', "")
                rel_error = item.get('(%)', "")

                # translated note
                values = (
                    compound_name,
                    f"{source_ppg:.2f}" if isinstance(source_ppg, (int, float)) else source_ppg,
                    f"{pred_rt:.3f}" if isinstance(pred_rt, (int, float)) else pred_rt,
                    f"{actual_rt:.3f}" if isinstance(actual_rt, (int, float)) else actual_rt,
                    f"{abs_error:.3f}" if isinstance(abs_error, (int, float)) else abs_error,
                    f"{rel_error:.2f}" if isinstance(rel_error, (int, float)) else rel_error
                )

                self.conversion_tree.insert("", tk.END, values=values)

        # translated note
        stats_text = f"analysis results: {analysis_results['condition']} → {analysis_results['condition']}\n"
        stats_text += "=" * 50 + "\n\n"

        stats_text += f": {analysis_results['']}\n"
        stats_text += f"compound: {analysis_results['compound']}\n"
        stats_text += f": {analysis_results['']}\n\n"

        if isinstance(analysis_results[''], dict):
            stats_text += ":\n"
            for key, value in analysis_results[''].items():
                if isinstance(value, float):
                    if '' in key or '' in key:
                        stats_text += f"  {key}: {value:.3f}\n"
                    else:
                        stats_text += f"  {key}: {value:.2f}\n"
                else:
                    stats_text += f"  {key}: {value}\n"

        if '' in analysis_results:
            stats_text += "\n:\n"
            for bin_name, count in analysis_results[''].items():
                percentage = (count / analysis_results[''] * 100) if analysis_results[''] > 0 else 0
                stats_text += f"  {bin_name}: {count} ({percentage:.1f}%)\n"

        if 'analysis' in analysis_results:
            pass_analysis = analysis_results['analysis']
            stats_text += f"\nanalysis (: {pass_analysis['(min)']}min):\n"
            stats_text += f" : {pass_analysis['']}\n"
            stats_text += f" : {pass_analysis['']}\n"
            stats_text += f" : {pass_analysis['(%)']:.1f}%\n"

        self.conversion_stats_text.insert(tk.END, stats_text, "INFO")

    def compare_multiple_conversions(self):
        """compare"""
        # translated note
        selected_conversions = []
        for var, from_cond, to_cond in self.conversion_scheme_vars:
            if var.get():
                selected_conversions.append((from_cond, to_cond))

        if len(selected_conversions) < 2:
            messagebox.showwarning("Warning", "compare")
            return

        # analysis
        self.log_message(f"compare: {len(selected_conversions)} ", "INFO")
        self.update_status(f"compare")

        # compare
        self.processing_thread = threading.Thread(
            target=self._compare_multiple_conversions_thread,
            args=(selected_conversions,)
        )
        self.processing_thread.start()

    def _compare_multiple_conversions_thread(self, conversions):
        """compare"""
        try:
            # translated note
            if self.visualizer is None:
                self.visualizer = PPGVisualizer(self.calculator)

            # comparechart
            fig = self.visualizer.plot_multiple_conversion_comparison(conversions)

            if fig is not None:
                # chart
                self.root.after(0, lambda: self.display_comparison_figure(fig, conversions))

                self.root.after(0, lambda: self.log_message(
                    f"✓ compare, compare {len(conversions)} ", "SUCCESS"))
            else:
                self.root.after(0, lambda: self.log_message("✗ comparechart", "ERROR"))

            self.root.after(0, lambda: self.update_status(""))

        except Exception as e:
            self.root.after(0, lambda: self.log_message(f"✗ comparefailed: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("comparefailed"))

    def display_comparison_figure(self, fig, conversions):
        """comparechart"""
        # chart
        self.clear_figure()

        # translated note
        self.viz_placeholder.pack_forget()

        # translated note
        canvas = FigureCanvasTkAgg(fig, master=self.viz_placeholder.master)
        canvas.draw()

        # translated note
        toolbar = NavigationToolbar2Tk(canvas, self.viz_placeholder.master)
        toolbar.update()

        # translated note
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        toolbar.pack(fill=tk.X)

        # save
        self.figure_canvas = canvas
        self.figure_toolbar = toolbar
        self.current_figure = fig

        # translated note
        self.viz_status_var.set(f"compare")

    def export_conversion_results(self):
        """results"""
        if not self.calculator.conversion_results:
            messagebox.showinfo("", "results")
            return

        file_path = filedialog.asksaveasfilename(
            title="results",
            defaultextension=".xlsx",
            filetypes=[("Excelfile", "*.xlsx"), ("file", "*.*")]
        )

        if file_path:
            try:
                with pd.ExcelWriter(file_path, engine='openpyxl') as writer:
                    # saveresults
                    for key, conversion in self.calculator.conversion_results.items():
                        if 'data' in conversion:
                            df = pd.DataFrame(conversion['data'])
                            sheet_name = key[:30] # Excel sheetname
                            df.to_excel(writer, sheet_name=sheet_name, index=False)

                    # save
                    conversion_stats = []
                    for key, conversion in self.calculator.conversion_results.items():
                        stats = {
                            '': key,
                            'condition': conversion.get('condition', ''),
                            'condition': conversion.get('condition', ''),
                            'compound': conversion.get('compound', 0),
                            '': conversion.get('', 0)
                        }

                        if isinstance(conversion.get(''), dict):
                            stats.update({
                                '(min)': conversion[''].get('', 0),
                                '': conversion[''].get('', 0),
                                '(min)': conversion[''].get('', 0)
                            })

                        if isinstance(conversion.get('analysis'), dict):
                            stats.update({
                                '(%)': conversion['analysis'].get('(%)', 0),
                                '': conversion['analysis'].get('', 0),
                                '(min)': conversion['analysis'].get('(min)', 0.5)
                            })

                        conversion_stats.append(stats)

                    if conversion_stats:
                        stats_df = pd.DataFrame(conversion_stats)
                        stats_df.to_excel(writer, sheet_name='', index=False)

                self.log_message(f"✓ results: {file_path}", "SUCCESS")
                self.conversion_status_var.set(f"results: {Path(file_path).name}")

            except Exception as e:
                self.log_message(f"✗ resultsfailed: {str(e)}", "ERROR")

    def clear_conversion_results(self):
        """results"""
        # Treeview
        for item in self.conversion_tree.get_children():
            self.conversion_tree.delete(item)

        # translated note
        self.conversion_stats_text.delete(1.0, tk.END)

        # calculateresults
        self.calculator.conversion_results.clear()

        self.conversion_status_var.set("results")
        self.log_message("results", "INFO")

    def generate_report(self):
        """analysis"""
        try:
            summary = self.calculator.generate_summary_report()

            # analysis results
            self.analysis_message("=" * 60, "INFO")
            self.analysis_message("PPG retention indicesanalysis", "INFO")
            self.analysis_message("=" * 60, "INFO")
            self.analysis_message(f": {summary['']}", "INFO")
            self.analysis_message("", "INFO")

            self.analysis_message("data:", "INFO")
            self.analysis_message(f" - PPG datacondition: {summary['PPG datacondition']}", "INFO")
            self.analysis_message(f" - compound data: {summary['compound data']}", "INFO")
            self.analysis_message(f" - standard curve: {summary['standard curve']}", "INFO")
            self.analysis_message(f" - PPGindexcalculateresults: {summary['PPGindexcalculateresults']}", "INFO")
            self.analysis_message(f" - cross-condition conversion results: {summary['cross-condition conversion results']}", "INFO")
            self.analysis_message("", "INFO")

            if summary['standard curve']:
                self.analysis_message("standard curve:", "INFO")
                for condition, perf in summary['standard curve'].items():
                    self.analysis_message(f"  {condition}:", "INFO")
                    self.analysis_message(f" - model: {perf['model']}", "INFO")
                    self.analysis_message(f"    - R²: {perf['R²']:.4f}", "INFO")
                    self.analysis_message(f" - : {perf['']:.4f}", "INFO")
                    self.analysis_message(f" - : {perf['']:.4f}", "INFO")
                    self.analysis_message(f" - : {perf['']:.4f}", "INFO")
                self.analysis_message("", "INFO")

            if summary['PPGindex']:
                self.analysis_message("PPGindex:", "INFO")
                for condition, stats in summary['PPGindex'].items():
                    self.analysis_message(f"  {condition}:", "INFO")
                    self.analysis_message(f" - calculatemethod: {stats['calculatemethod']}", "INFO")
                    self.analysis_message(f" - : {stats['']}", "INFO")
                    self.analysis_message(f" - : {stats['']:.2f}", "INFO")
                    self.analysis_message(f" - : {stats['']:.2f}", "INFO")
                    self.analysis_message(f" - : {stats['']:.2f} - {stats['']:.2f}", "INFO")
                self.analysis_message("", "INFO")

            if summary['condition']:
                self.analysis_message("condition:", "INFO")
                for key, stats in summary['condition'].items():
                    self.analysis_message(f"  {key}:", "INFO")
                    self.analysis_message(f" - condition: {stats['condition']}", "INFO")
                    self.analysis_message(f" - condition: {stats['condition']}", "INFO")
                    self.analysis_message(f" - compound: {stats['compound']}", "INFO")
                    self.analysis_message(f" - : {stats['']}", "INFO")
                    self.analysis_message(f" - : {stats['']:.3f} min", "INFO")
                    self.analysis_message(f" - : {stats['(%)']:.1f}%", "INFO")
                self.analysis_message("", "INFO")

            self.analysis_message("", "SUCCESS")

            # results
            self.results_text.delete(1.0, tk.END)
            self.results_text.insert(tk.END, "PPG retention indicesanalysis\n")
            self.results_text.insert(tk.END, "=" * 60 + "\n\n")
            self.results_text.insert(tk.END, f": {summary['']}\n\n")

            self.results_text.insert(tk.END, "data:\n")
            self.results_text.insert(tk.END, f" - PPG datacondition: {summary['PPG datacondition']}\n")
            self.results_text.insert(tk.END, f" - compound data: {summary['compound data']}\n")
            self.results_text.insert(tk.END, f" - standard curve: {summary['standard curve']}\n")
            self.results_text.insert(tk.END, f" - PPGindexcalculateresults: {summary['PPGindexcalculateresults']}\n")
            self.results_text.insert(tk.END, f" - cross-condition conversion results: {summary['cross-condition conversion results']}\n\n")

            self.analysis_status_var.set("")

        except Exception as e:
            self.analysis_message(f"✗ failed: {str(e)}", "ERROR")

    def on_viz_type_change(self, event=None):
        """"""
        viz_type = self.viz_type_var.get()

        # Parameters
        self.viz_conditions_frame.grid_remove()
        self.conversion_viz_frame.grid_remove()
        self.viz_condition_combo.grid()

        if viz_type in ["condition_comparison", "multiple_conversion_comparison"]:
            # condition
            self.viz_conditions_frame.grid()
            self.viz_condition_combo.grid_remove()

            # condition
            if viz_type == "condition_comparison":
                self.update_viz_condition_checkboxes()
            else:
                self.update_viz_conversion_checkboxes()

        elif viz_type == "conversion_analysis":
            # analysisParameters
            self.conversion_viz_frame.grid()
            self.viz_condition_combo.grid_remove()

            # analysisParameters
            self.update_conversion_viz_params()

    def update_viz_condition_checkboxes(self):
        """condition"""
        # translated note
        for widget in self.viz_conditions_frame.winfo_children():
            widget.destroy()

        # condition
        conditions = list(self.calculator.ppg_data.keys())

        # translated note
        ttk.Label(self.viz_conditions_frame, text="comparecondition:").grid(row=0, column=0, columnspan=4, sticky=tk.W,
                                                                            pady=(0, 5))

        self.viz_condition_vars = {}
        for i, condition in enumerate(conditions):
            var = tk.BooleanVar(value=(i < 2)) # condition
            cb = ttk.Checkbutton(self.viz_conditions_frame, text=condition, variable=var)
            cb.grid(row=1 + i // 4, column=i % 4, sticky=tk.W, padx=5, pady=2)
            self.viz_condition_vars[condition] = var

    def update_viz_conversion_checkboxes(self):
        """"""
        # translated note
        for widget in self.viz_conditions_frame.winfo_children():
            widget.destroy()

        # condition
        conditions = list(self.calculator.ppg_data.keys())

        # translated note
        ttk.Label(self.viz_conditions_frame, text="compare:").grid(row=0, column=0, columnspan=4,
                                                                                sticky=tk.W,
                                                                                pady=(0, 5))

        self.viz_conversion_vars = []

        if len(conditions) >= 2:
            for i in range(len(conditions)):
                for j in range(len(conditions)):
                    if i != j:
                        from_cond = conditions[i]
                        to_cond = conditions[j]
                        var = tk.BooleanVar(value=(i == 0 and j == 1)) #
                        cb = ttk.Checkbutton(self.viz_conditions_frame,
                                             text=f"{from_cond}→{to_cond}",
                                             variable=var)
                        row = (i * len(conditions) + j) // 4 + 1
                        col = (i * len(conditions) + j) % 4
                        cb.grid(row=row, column=col, sticky=tk.W, padx=2, pady=1)
                        self.viz_conversion_vars.append((var, from_cond, to_cond))
        else:
            ttk.Label(self.viz_conditions_frame, text="2conditioncompare",
                      foreground="gray").grid(row=1, column=0, sticky=tk.W)

    def update_conversion_viz_params(self):
        """Parameters"""
        # translated note
        for widget in self.conversion_viz_frame.winfo_children():
            widget.destroy()

        # condition
        conditions = list(self.calculator.ppg_data.keys())

        # conditioncondition
        ttk.Label(self.conversion_viz_frame, text="condition:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.viz_from_condition_var = tk.StringVar()
        viz_from_combo = ttk.Combobox(self.conversion_viz_frame, textvariable=self.viz_from_condition_var,
                                      values=conditions, width=15, state="readonly")
        viz_from_combo.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)
        if conditions:
            self.viz_from_condition_var.set(conditions[0])

        ttk.Label(self.conversion_viz_frame, text="condition:").grid(row=0, column=2, sticky=tk.W, padx=(0, 5), pady=5)
        self.viz_to_condition_var = tk.StringVar()
        viz_to_combo = ttk.Combobox(self.conversion_viz_frame, textvariable=self.viz_to_condition_var,
                                    values=conditions, width=15, state="readonly")
        viz_to_combo.grid(row=0, column=3, sticky=tk.W, padx=(0, 10), pady=5)
        if len(conditions) > 1:
            self.viz_to_condition_var.set(conditions[1])

    def generate_visualization(self):
        """visualization charts"""
        viz_type = self.viz_type_var.get()

        if viz_type == "condition_comparison":
            # condition
            selected_conditions = []
            if hasattr(self, 'viz_condition_vars'):
                for condition, var in self.viz_condition_vars.items():
                    if var.get():
                        selected_conditions.append(condition)

            if len(selected_conditions) < 2:
                messagebox.showwarning("Warning", "conditioncompare")
                return

            self.create_visualization(viz_type, selected_conditions)

        elif viz_type == "multiple_conversion_comparison":
            # translated note
            selected_conversions = []
            if hasattr(self, 'viz_conversion_vars'):
                for var, from_cond, to_cond in self.viz_conversion_vars:
                    if var.get():
                        selected_conversions.append((from_cond, to_cond))

            if len(selected_conversions) < 2:
                messagebox.showwarning("Warning", "compare")
                return

            self.create_visualization(viz_type, selected_conversions)

        elif viz_type == "conversion_analysis":
            # Parameters
            from_condition = self.viz_from_condition_var.get()
            to_condition = self.viz_to_condition_var.get()

            if not from_condition or not to_condition:
                messagebox.showwarning("Warning", "conditioncondition")
                return

            if from_condition == to_condition:
                messagebox.showwarning("Warning", "conditioncondition")
                return

            self.create_visualization(viz_type, (from_condition, to_condition))

        else:
            condition = self.viz_condition_var.get()
            if not condition:
                messagebox.showwarning("Warning", "condition")
                return

            self.create_visualization(viz_type, condition)

    def create_visualization(self, viz_type, conditions):
        """visualization charts"""
        try:
            # translated note
            if self.visualizer is None and self.calculator.standard_curves:
                self.visualizer = PPGVisualizer(self.calculator)
            elif self.visualizer is None:
                messagebox.showwarning("Warning", "standard curve")
                return

            # translated note
            self.viz_status_var.set(f" {viz_type} chart...")
            self.update_status(f"chart: {viz_type}")

            # chart
            if viz_type == "standard_curve":
                if isinstance(conditions, str):
                    conditions = [conditions]
                fig = self.visualizer.plot_standard_curves(conditions)
            elif viz_type == "residuals":
                if isinstance(conditions, str):
                    conditions = [conditions]
                fig = self.visualizer.plot_residuals(conditions)
            elif viz_type == "index_distribution":
                fig = self.visualizer.plot_ppg_index_distribution(conditions)
            elif viz_type == "condition_comparison":
                fig = self.visualizer.plot_condition_comparison(conditions)
            elif viz_type == "conversion_analysis":
                from_condition, to_condition = conditions
                fig = self.visualizer.plot_conversion_analysis(from_condition, to_condition)
            elif viz_type == "multiple_conversion_comparison":
                fig = self.visualizer.plot_multiple_conversion_comparison(conditions)
            else:
                self.viz_status_var.set("")
                return

            if fig is None:
                self.viz_status_var.set("chart, data")
                return

            # chart
            self.clear_figure()

            # translated note
            self.viz_placeholder.pack_forget()

            # translated note
            canvas = FigureCanvasTkAgg(fig, master=self.viz_placeholder.master)
            canvas.draw()

            # translated note
            toolbar = NavigationToolbar2Tk(canvas, self.viz_placeholder.master)
            toolbar.update()

            # translated note
            canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
            toolbar.pack(fill=tk.X)

            # save
            self.figure_canvas = canvas
            self.figure_toolbar = toolbar
            self.current_figure = fig

            self.viz_status_var.set(f"{viz_type} chart")
            self.update_status("")

        except Exception as e:
            self.viz_status_var.set(f"chartfailed: {str(e)}")
            self.log_message(f"✗ chartfailed: {str(e)}", "ERROR")

    def save_figure(self):
        """savechart"""
        if self.current_figure is None:
            messagebox.showwarning("Warning", "savechart")
            return

        file_types = [("PNGfile", "*.png"), ("PDFfile", "*.pdf"),
                      ("SVGfile", "*.svg"), ("file", "*.*")]

        file_path = filedialog.asksaveasfilename(
            title="savechart",
            filetypes=file_types,
            defaultextension=".png"
        )

        if file_path:
            try:
                self.current_figure.savefig(file_path, dpi=300, bbox_inches='tight')
                self.log_message(f"✓ chartsave: {file_path}", "SUCCESS")
                self.viz_status_var.set(f"chartsave: {Path(file_path).name}")
            except Exception as e:
                self.log_message(f"✗ savechartfailed: {str(e)}", "ERROR")

    def clear_figure(self):
        """chart"""
        if self.figure_canvas:
            self.figure_canvas.get_tk_widget().destroy()
            self.figure_toolbar.destroy()
            self.figure_canvas = None
            self.figure_toolbar = None

        if self.current_figure:
            plt.close(self.current_figure)
            self.current_figure = None

        # translated note
        self.viz_placeholder.pack(expand=True)

    def save_all_results(self):
        """saveresults"""
        output_dir = self.output_dir_var.get().strip()

        if not output_dir:
            messagebox.showwarning("Warning", "outputdirectory")
            return

        # save
        self.output_status_var.set("saveresults...")
        self.update_status("saveresults")

        # save
        self.processing_thread = threading.Thread(
            target=self._save_all_results_thread,
            args=(output_dir,)
        )
        self.processing_thread.start()

    def _save_all_results_thread(self, output_dir):
        """saveresults"""
        try:
            success, msg, saved_files = self.calculator.save_results(output_dir)

            if success:
                self.root.after(0, lambda: self.output_status_var.set(f"resultssave: {len(saved_files)} file"))
                self.root.after(0, lambda: self.log_message(f"✓ {msg}", "SUCCESS"))

                # resultssavefile
                self.results_text.delete(1.0, tk.END)
                self.results_text.insert(tk.END, "savefilecolumn:\n")
                self.results_text.insert(tk.END, "=" * 60 + "\n\n")

                for file_path in saved_files:
                    self.results_text.insert(tk.END, f"• {Path(file_path).name}\n")

                self.results_text.insert(tk.END, f"\nfilesave: {output_dir}")
            else:
                self.root.after(0, lambda: self.output_status_var.set("savefailed"))
                self.root.after(0, lambda: self.log_message(f"✗ {msg}", "ERROR"))

            self.root.after(0, lambda: self.update_status(""))

        except Exception as e:
            self.root.after(0, lambda: self.output_status_var.set("savefailed"))
            self.root.after(0, lambda: self.log_message(f"✗ saveresultsfailed: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("savefailed"))

    def preview_report(self):
        """"""
        # translated note
        self.generate_report()

        # results
        self.notebook.select(4) # results5 (4)

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

    def log_message(self, message: str, level: str = "INFO"):
        """message"""
        timestamp = datetime.now().strftime("%H:%M:%S")
        log_message = f"[{timestamp}] {message}"

        self.log_text.insert(tk.END, log_message + "\n", level)
        self.log_text.see(tk.END)
        self.root.update()

    def analysis_message(self, message: str, level: str = "INFO"):
        """analysismessage"""
        self.analysis_text.insert(tk.END, message + "\n", level)
        self.analysis_text.see(tk.END)
        self.root.update()

    def clear_log(self):
        """"""
        self.log_text.delete(1.0, tk.END)

    def clear_analysis_text(self):
        """analysis"""
        self.analysis_text.delete(1.0, tk.END)

    def save_log(self):
        """save"""
        file_path = filedialog.asksaveasfilename(
            title="save",
            defaultextension=".txt",
            filetypes=[("file", "*.txt"), ("file", "*.*")]
        )

        if file_path:
            try:
                with open(file_path, 'w', encoding='utf-8') as f:
                    f.write(self.log_text.get(1.0, tk.END))
                self.log_message(f"✓ save: {file_path}", "SUCCESS")
            except Exception as e:
                self.log_message(f"✗ savefailed: {str(e)}", "ERROR")

    def update_status(self, message: str):
        """"""
        self.status_var.set(message)
        self.root.update()

    def on_closing(self):
        """"""
        if messagebox.askyesno("", "？"):
            # translated note
            if self.current_figure:
                plt.close(self.current_figure)
            self.root.destroy()


def check_dependencies():
    """"""
    print("=" * 70)
    print("...")
    print("=" * 70)

    dependencies = {
        'pandas': 'data',
        'numpy': 'calculate',
        'scipy': 'calculate',
        'matplotlib': '',
        'seaborn': '',
        'openpyxl': 'Excelfile',
        'tkinter': 'GUI'
    }

    missing = []

    for lib, desc in dependencies.items():
        try:
            if lib == 'tkinter':
                import tkinter
            else:
                __import__(lib)
            print(f"✓ {lib}: {desc}")
        except ImportError:
            print(f"✗ {lib}: {desc} - ")
            missing.append(lib)

    if missing:
        print(f"\n:")
        for lib in missing:
            if lib == 'tkinter':
                print(" : tkinterPython, :")
                print("    Ubuntu/Debian: sudo apt-get install python3-tk")
                print(" Windows/macOS: Pythontkinter")
            else:
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
        print(" Windows/macOS: Pythontkinter")
        return

    # translated note
    if not check_dependencies():
        response = input("\n, ? (y/n): ")
        if response.lower() != 'y':
            return

    # translated note
    root = tk.Tk()

    # translated note
    root.title("PPG retention indicescalculateanalysis - ")

    # translated note
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()

    # (80%)
    window_width = int(screen_width * 0.8)
    window_height = int(screen_height * 0.8)

    # calculate ()
    x = (screen_width - window_width) // 2
    y = (screen_height - window_height) // 2

    root.geometry(f"{window_width}x{window_height}+{x}+{y}")

    # GUI
    app = PPGIndexAnalyzerGUI(root)

    # translated note
    root.mainloop()


if __name__ == "__main__":
    print("PPG retention indicescalculateanalysis - ")
    print("Version: 3.0 (condition)")
    print(":")
    print(" 1. loadPPGstandardcompoundretention_timedata")
    print(" 2. Fit PPG standard curves and calculate linear relationships")
    print(" 3. Calculate compound PPG retention indices")
    print(" 4. Compare PPG indices across chromatographic conditions")
    print(" 5. Convert and validate PPG indices across conditions")
    print(" 6. comparison of multiple conversion schemes")
    print(" 7. Visualize analysis results and generate reports")
    print("=" * 70)

    main()
