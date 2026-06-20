#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PPG retention index analysis system - complete GUI version

:
1. Data loading and management (PPG data, SMRT dataset, validation set)
2. PPG retention index calculation and conversion
3. analysis and comparison of five calibration methods
4. model performance evaluation
5. visualization generation
6. experimental report export
"""

import os
import sys
import threading
import traceback
import warnings
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any

# Suppress warnings
warnings.filterwarnings('ignore')

# Import required libraries
try:
    import pandas as pd
    import numpy as np
    from scipy import stats, optimize
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    from matplotlib.figure import Figure
    import seaborn as sns
    from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
    from sklearn.model_selection import train_test_split
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LinearRegression
    import matplotlib

    matplotlib.use('TkAgg') # Set the matplotlib backend
except ImportError as e:
    print(f"Error: please install the required libraries first: {e}")
    print("Installation command: pip install pandas numpy scipy matplotlib seaborn scikit-learn")
    sys.exit(1)

# Try to import GUI libraries
try:
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox, scrolledtext, Toplevel, StringVar, BooleanVar, IntVar

    GUI_AVAILABLE = True
except ImportError:
    GUI_AVAILABLE = False
    print("Warning: tkinter is not installed, so the GUI is unavailable")

# Set fonts
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
sns.set_style("whitegrid")


class PPGExperimentAnalyzer:
    """Core class for the PPG experimental-design analyzer"""

    def __init__(self):
        self.ppg_data = {} # PPG data under different conditions
        self.smrt_data = {} # SMRT dataset
        self.validation_data = {} # validation setdata
        self.calibration_data = {} # calibration compound data
        self.ppg_indices = {} # PPG retention indices
        self.models = {} # trained models
        self.results = {} # analysis results
        self.calibration_methods = {} # calibration methods
        self.visualizations = {} # visualization charts

        # experimental parameters
        self.experiment_params = {
            'ppg_model_type': 'logarithmic', # PPG model type
            'calibration_method': 'linear', # calibration methods
            'validation_split': 0.2, # validation split
            'random_seed': 42 #
        }

    def load_ppg_data(self, file_path: str, condition: str) -> Tuple[bool, str]:
        """loadPPG data"""
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

    def load_compound_data(self, file_path: str, data_type: str, condition: str = "default") -> Tuple[bool, str]:
        """loadcompound data"""
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
                'CAS': ['CAS', 'CAS', 'CAS No.', 'CAS'],
                'molecule': ['molecule', 'MW', 'MolecularWeight', 'molecule(Da)'],
                'logP': ['logP', 'LogP', 'log P', '']
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
            key = f"{data_type}_{condition}"

            if data_type == "smrt":
                self.smrt_data[key] = df
            elif data_type == "validation":
                self.validation_data[key] = df
            elif data_type == "calibration":
                self.calibration_data[key] = df
            else:
                return False, f"data: {data_type}"

            return True, f"load {len(df)} {data_type}compound data (condition: {condition})"

        except Exception as e:
            return False, f"loadcompound datafailed: {str(e)}"

    def calculate_ppg_index(self, rt: float, condition: str) -> float:
        """calculatePPG retention indices"""
        if condition not in self.ppg_data:
            raise ValueError(f"condition {condition} PPG data")

        df_ppg = self.ppg_data[condition]
        ppg_rt = df_ppg['retention_time'].values
        ppg_n = df_ppg['degree_of_polymerization'].values

        # translated note
        if rt < ppg_rt[0]:
            # translated note
            if len(ppg_rt) >= 2:
                n_calc = ppg_n[0] - (ppg_rt[0] - rt) / (ppg_rt[1] - ppg_rt[0]) * (ppg_n[1] - ppg_n[0])
            else:
                n_calc = ppg_n[0]
        elif rt > ppg_rt[-1]:
            # translated note
            if len(ppg_rt) >= 2:
                n_calc = ppg_n[-1] + (rt - ppg_rt[-1]) / (ppg_rt[-1] - ppg_rt[-2]) * (ppg_n[-1] - ppg_n[-2])
            else:
                n_calc = ppg_n[-1]
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

        # retention_index (100)
        return n_calc * 100

    def calculate_rt_from_index(self, index: float, condition: str, method: str = "interpolation") -> float:
        """PPG retention indicescalculateretention_time"""
        if condition not in self.ppg_data:
            raise ValueError(f"condition {condition} PPG data")

        df_ppg = self.ppg_data[condition]
        ppg_rt = df_ppg['retention_time'].values
        ppg_n = df_ppg['degree_of_polymerization'].values
        n_target = index / 100 # indexdegree_of_polymerization

        if method == "interpolation":
            # translated note
            if n_target < ppg_n[0]:
                # translated note
                if len(ppg_n) >= 2:
                    rt_calc = ppg_rt[0] - (ppg_n[0] - n_target) / (ppg_n[1] - ppg_n[0]) * (ppg_rt[1] - ppg_rt[0])
                else:
                    rt_calc = ppg_rt[0]
            elif n_target > ppg_n[-1]:
                # translated note
                if len(ppg_n) >= 2:
                    rt_calc = ppg_rt[-1] + (n_target - ppg_n[-1]) / (ppg_n[-1] - ppg_n[-2]) * (ppg_rt[-1] - ppg_rt[-2])
                else:
                    rt_calc = ppg_rt[-1]
            else:
                # translated note
                idx = np.searchsorted(ppg_n, n_target) - 1
                if idx < 0:
                    idx = 0
                elif idx >= len(ppg_n) - 1:
                    idx = len(ppg_n) - 2

                n_i, n_j = ppg_n[idx], ppg_n[idx + 1]
                rt_i, rt_j = ppg_rt[idx], ppg_rt[idx + 1]

                rt_calc = rt_i + (rt_j - rt_i) * (n_target - n_i) / (n_j - n_i)

        elif method == "regression":
            # - UsePPGstandard curve
            # model: RT = a + b * ln(n)
            log_n = np.log(ppg_n)

            # translated note
            slope, intercept, r_value, p_value, std_err = stats.linregress(log_n, ppg_rt)

            # retention_time
            rt_calc = intercept + slope * np.log(n_target)

            # Parameters
            if 'regression_params' not in self.calibration_methods:
                self.calibration_methods['regression_params'] = {}
            self.calibration_methods['regression_params'][condition] = {
                'slope': slope,
                'intercept': intercept,
                'r_squared': r_value ** 2,
                'std_err': std_err
            }

        else:
            raise ValueError(f"calculatemethod: {method}")

        return rt_calc

    def fit_ppg_standard_curve(self, condition: str) -> Tuple[bool, str, Dict]:
        """PPGstandard curve"""
        if condition not in self.ppg_data:
            return False, f"condition {condition} PPG data", {}

        df_ppg = self.ppg_data[condition]
        n_values = df_ppg['degree_of_polymerization'].values
        rt_values = df_ppg['retention_time'].values

        # model
        models = {}

        # 1. model: RT = a + b * ln(n)
        log_n = np.log(n_values)
        try:
            slope_log, intercept_log, r_log, p_log, std_err_log = stats.linregress(log_n, rt_values)
            y_pred_log = intercept_log + slope_log * log_n
            ss_res_log = np.sum((rt_values - y_pred_log) ** 2)
            ss_tot_log = np.sum((rt_values - np.mean(rt_values)) ** 2)
            r_squared_log = 1 - (ss_res_log / ss_tot_log) if ss_tot_log != 0 else 0

            models['logarithmic'] = {
                'params': (intercept_log, slope_log),
                'r_squared': r_squared_log,
                'std_err': std_err_log,
                'p_value': p_log,
                'func': lambda x, a, b: a + b * np.log(x)
            }
        except:
            models['logarithmic'] = {'r_squared': -1}

        # 2. model: RT = a + b * n
        try:
            slope_lin, intercept_lin, r_lin, p_lin, std_err_lin = stats.linregress(n_values, rt_values)
            y_pred_lin = intercept_lin + slope_lin * n_values
            ss_res_lin = np.sum((rt_values - y_pred_lin) ** 2)
            ss_tot_lin = np.sum((rt_values - np.mean(rt_values)) ** 2)
            r_squared_lin = 1 - (ss_res_lin / ss_tot_lin) if ss_tot_lin != 0 else 0

            models['linear'] = {
                'params': (intercept_lin, slope_lin),
                'r_squared': r_squared_lin,
                'std_err': std_err_lin,
                'p_value': p_lin,
                'func': lambda x, a, b: a + b * x
            }
        except:
            models['linear'] = {'r_squared': -1}

        # model
        best_model_name = None
        best_r2 = -1

        for model_name, model_data in models.items():
            if 'r_squared' in model_data and model_data['r_squared'] > best_r2:
                best_r2 = model_data['r_squared']
                best_model_name = model_name

        if best_model_name is None:
            return False, "model", {}

        best_model = models[best_model_name]

        # standard curve
        if 'standard_curves' not in self.calibration_methods:
            self.calibration_methods['standard_curves'] = {}

        self.calibration_methods['standard_curves'][condition] = {
            'model_type': best_model_name,
            'params': best_model['params'],
            'r_squared': best_model['r_squared'],
            'n_values': n_values,
            'rt_values': rt_values
        }

        return True, f"PPGstandard curve: {best_model_name}model, R²={best_model['r_squared']:.6f}", best_model

    def apply_calibration_methods(self, source_condition: str, target_condition: str) -> Tuple[bool, str]:
        """calibration"""
        print(f"calibrationanalysis ({source_condition} → {target_condition}):")

        if len(self.smrt_data) == 0:
            return False, "SMRT dataset"

        # SMRTdata
        smrt_key = f"smrt_{source_condition}"
        if smrt_key not in self.smrt_data:
            # SMRTdata
            smrt_keys = [k for k in self.smrt_data.keys() if k.startswith('smrt_')]
            if not smrt_keys:
                return False, "SMRTdata"
            smrt_key = smrt_keys[0]

        smrt_df = self.smrt_data[smrt_key]

        # data
        validation_key = f"validation_{target_condition}"
        if validation_key not in self.validation_data:
            # data
            validation_keys = [k for k in self.validation_data.keys() if k.startswith('validation_')]
            if not validation_keys:
                # data, UseSMRTdata
                validation_df = smrt_df
            else:
                validation_key = validation_keys[0]
                validation_df = self.validation_data[validation_key]
        else:
            validation_df = self.validation_data[validation_key]

        # datacolumn
        if 'retention_time' not in smrt_df.columns:
            return False, "SMRTdata'retention_time'column"

        # method1: Use, calibration
        print("1. method1: UseSMRT dataset, calibration")
        # conditioncondition, dataRT
        if 'RT_target' in smrt_df.columns and 'RT_source' in smrt_df.columns:
            rt_source = smrt_df['RT_source'].values
            rt_target = smrt_df['RT_target'].values

            mae = mean_absolute_error(rt_target, rt_source)
            rmse = np.sqrt(mean_squared_error(rt_target, rt_source))
            r2 = r2_score(rt_target, rt_source)

            self.results['method1'] = {
                'MAE': mae,
                'RMSE': rmse,
                'R2': r2,
                'description': 'UseSMRT dataset, calibration'
            }
            print(f"   MAE: {mae:.4f}, RMSE: {rmse:.4f}, R²: {r2:.4f}")
        else:
            self.results['method1'] = {
                'MAE': np.nan,
                'RMSE': np.nan,
                'R2': np.nan,
                'description': 'UseSMRT dataset, calibration',
                'error': 'retention_timedata'
            }

        # method2: Use28compoundcalibration
        print("2. method2: Use28compoundcalibration")
        # calibration compound data
        calib_key = f"calibration_{source_condition}"
        if calib_key in self.calibration_data:
            calib_df = self.calibration_data[calib_key]

            if 'RT_source' in calib_df.columns and 'RT_target' in calib_df.columns:
                rt_source_calib = calib_df['RT_source'].values
                rt_target_calib = calib_df['RT_target'].values

                # calibrationmodel
                try:
                    slope, intercept, r_value, p_value, std_err = stats.linregress(rt_source_calib, rt_target_calib)

                    # calibrationSMRTcompound
                    rt_source_all = smrt_df['retention_time'].values
                    rt_corrected = intercept + slope * rt_source_all

                    # data,
                    if 'RT_target' in smrt_df.columns:
                        rt_actual = smrt_df['RT_target'].values
                        mae = mean_absolute_error(rt_actual, rt_corrected)
                        rmse = np.sqrt(mean_squared_error(rt_actual, rt_corrected))
                        r2 = r2_score(rt_actual, rt_corrected)
                    else:
                        mae = rmse = r2 = np.nan

                    self.results['method2'] = {
                        'MAE': mae,
                        'RMSE': rmse,
                        'R2': r2,
                        'correction_params': {'slope': slope, 'intercept': intercept, 'r_squared': r_value ** 2},
                        'description': 'Usecompoundcalibration'
                    }
                    print(f"   MAE: {mae:.4f}, RMSE: {rmse:.4f}, R²: {r2:.4f}")
                except Exception as e:
                    self.results['method2'] = {
                        'MAE': np.nan,
                        'RMSE': np.nan,
                        'R2': np.nan,
                        'description': 'Usecompoundcalibration',
                        'error': f'calibrationfailed: {str(e)}'
                    }
            else:
                self.results['method2'] = {
                    'MAE': np.nan,
                    'RMSE': np.nan,
                    'R2': np.nan,
                    'description': 'Usecompoundcalibration',
                    'error': 'calibrationdataRT_sourceRT_targetcolumn'
                }
        else:
            self.results['method2'] = {
                'MAE': np.nan,
                'RMSE': np.nan,
                'R2': np.nan,
                'description': 'Usecompoundcalibration',
                'error': 'calibration compound data'
            }

        # method3: UsePPGcalibration
        print("3. method3: UsePPGcalibration")
        if source_condition in self.ppg_data and target_condition in self.ppg_data:
            # datacalculatePPGindex
            ppg_indices = []
            rt_target_actual = []
            rt_target_predicted = []

            for _, row in validation_df.iterrows():
                if 'retention_time' in row:
                    rt_source = row['retention_time']

                    try:
                        # calculatePPGindex
                        ppg_index = self.calculate_ppg_index(rt_source, source_condition)

                        # indexcalculateRT
                        rt_pred = self.calculate_rt_from_index(ppg_index, target_condition, method='regression')

                        ppg_indices.append(ppg_index)
                        rt_target_predicted.append(rt_pred)

                        # RT
                        if 'RT_target' in row:
                            rt_target_actual.append(row['RT_target'])
                    except Exception as e:
                        print(f" calculatefailed: {e}")
                        continue

            if rt_target_predicted:
                rt_target_predicted = np.array(rt_target_predicted)

                if rt_target_actual:
                    rt_target_actual = np.array(rt_target_actual)

                    mae = mean_absolute_error(rt_target_actual, rt_target_predicted)
                    rmse = np.sqrt(mean_squared_error(rt_target_actual, rt_target_predicted))
                    r2 = r2_score(rt_target_actual, rt_target_predicted)
                else:
                    mae = rmse = r2 = np.nan

                self.results['method3'] = {
                    'MAE': mae,
                    'RMSE': rmse,
                    'R2': r2,
                    'ppg_indices': ppg_indices,
                    'rt_predicted': rt_target_predicted.tolist(),
                    'rt_actual': rt_target_actual.tolist() if rt_target_actual else [],
                    'description': 'UsePPGcalibration'
                }
                print(f"   MAE: {mae:.4f}, RMSE: {rmse:.4f}, R²: {r2:.4f}")
            else:
                self.results['method3'] = {
                    'MAE': np.nan,
                    'RMSE': np.nan,
                    'R2': np.nan,
                    'description': 'UsePPGcalibration',
                    'error': 'calculatePPGindexRT'
                }
        else:
            self.results['method3'] = {
                'MAE': np.nan,
                'RMSE': np.nan,
                'R2': np.nan,
                'description': 'UsePPGcalibration',
                'error': f'PPG data: condition={source_condition}, condition={target_condition}'
            }

        # method4: UsePPGcalibrationmodelretention_time
        print("4. method4: UsePPGcalibrationmodelretention_time")
        # PPGstandard curve
        success, msg, ppg_model = self.fit_ppg_standard_curve(target_condition)

        if success and ppg_model and 'func' in ppg_model:
            rt_target_predicted = []
            rt_target_actual = []

            for _, row in validation_df.iterrows():
                if 'retention_time' in row:
                    rt_source = row['retention_time']

                    try:
                        # calculatePPGindex
                        ppg_index = self.calculate_ppg_index(rt_source, source_condition)
                        n_value = ppg_index / 100

                        # UsePPGmodelRT
                        if 'params' in ppg_model:
                            rt_pred = ppg_model['func'](n_value, *ppg_model['params'])
                        else:
                            rt_pred = ppg_model['func'](n_value)

                        rt_target_predicted.append(rt_pred)

                        # RT
                        if 'RT_target' in row:
                            rt_target_actual.append(row['RT_target'])
                    except Exception as e:
                        print(f" failed: {e}")
                        continue

            if rt_target_predicted:
                rt_target_predicted = np.array(rt_target_predicted)

                if rt_target_actual:
                    rt_target_actual = np.array(rt_target_actual)

                    mae = mean_absolute_error(rt_target_actual, rt_target_predicted)
                    rmse = np.sqrt(mean_squared_error(rt_target_actual, rt_target_predicted))
                    r2 = r2_score(rt_target_actual, rt_target_predicted)
                else:
                    mae = rmse = r2 = np.nan

                self.results['method4'] = {
                    'MAE': mae,
                    'RMSE': rmse,
                    'R2': r2,
                    'description': 'UsePPGcalibrationmodelretention_time',
                    'model_type': ppg_model.get('model_type', 'unknown')
                }
                print(f"   MAE: {mae:.4f}, RMSE: {rmse:.4f}, R²: {r2:.4f}")
            else:
                self.results['method4'] = {
                    'MAE': np.nan,
                    'RMSE': np.nan,
                    'R2': np.nan,
                    'description': 'UsePPGcalibrationmodelretention_time',
                    'error': 'retention_time'
                }
        else:
            self.results['method4'] = {
                'MAE': np.nan,
                'RMSE': np.nan,
                'R2': np.nan,
                'description': 'UsePPGcalibrationmodelretention_time',
                'error': f'PPGstandard curve: {msg}'
            }

        # method5: PPGcalibrationcompoundcalibration
        print("5. method5: PPGcalibrationcompoundcalibration")
        # method2method3results
        comparison_results = {}

        if 'method2' in self.results and 'method3' in self.results:
            mae_diff = abs(self.results['method2'].get('MAE', np.nan) - self.results['method3'].get('MAE', np.nan))
            r2_diff = abs(self.results['method2'].get('R2', np.nan) - self.results['method3'].get('R2', np.nan))

            comparison_results = {
                'MAE_difference': mae_diff,
                'R2_difference': r2_diff,
                'recommendation': 'PPGcalibration' if self.results['method3'].get('R2', 0) > self.results['method2'].get('R2',
                                                                                                                  0) else 'compoundcalibration'
            }

        self.results['method5'] = {
            'description': 'PPGcalibrationcompoundcalibration',
            'comparison': comparison_results
        }

        print("calibrationanalysis!")
        return True, "calibrationanalysis"

    def evaluate_model_performance(self, model_type: str = 'linear') -> Tuple[bool, str]:
        """model"""
        print(f"{model_type}model...")

        # model
        # model

        # data
        all_data = []

        # SMRTdata
        for key, df in self.smrt_data.items():
            for _, row in df.iterrows():
                if 'retention_time' in row and 'logP' in row:
                    all_data.append({
                        'rt': row['retention_time'],
                        'logp': row['logP'] if not pd.isna(row['logP']) else 0
                    })

        # data
        for key, df in self.validation_data.items():
            for _, row in df.iterrows():
                if 'retention_time' in row and 'logP' in row:
                    all_data.append({
                        'rt': row['retention_time'],
                        'logp': row['logP'] if not pd.isna(row['logP']) else 0
                    })

        if len(all_data) < 10:
            return False, f"data, {len(all_data)} "

        # translated note
        features = np.array([[d['logp']] for d in all_data])
        targets = np.array([d['rt'] for d in all_data])

        # trainingtest
        X_train, X_test, y_train, y_test = train_test_split(
            features, targets, test_size=0.2, random_state=42
        )

        # translated note
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)

        # model
        model = LinearRegression()
        model.fit(X_train_scaled, y_train)

        # translated note
        y_pred = model.predict(X_test_scaled)

        # translated note
        mae = mean_absolute_error(y_test, y_pred)
        rmse = np.sqrt(mean_squared_error(y_test, y_pred))
        r2 = r2_score(y_test, y_pred)

        print(f" test:")
        print(f"    MAE: {mae:.4f}")
        print(f"    RMSE: {rmse:.4f}")
        print(f"    R²: {r2:.4f}")

        self.models[model_type] = {
            'model': model,
            'scaler': scaler,
            'performance': {'MAE': mae, 'RMSE': rmse, 'R2': r2},
            'n_samples': len(all_data)
        }

        return True, f"model: MAE={mae:.4f}, RMSE={rmse:.4f}, R²={r2:.4f}"

    def generate_visualizations(self, output_dir: str = None) -> Tuple[bool, str, Dict]:
        """visualization charts"""
        print("visualization charts...")

        visualizations = {}

        # 1. PPGstandard curve
        if self.ppg_data:
            try:
                fig_ppg = self._create_ppg_curves_plot()
                visualizations['ppg_curves'] = fig_ppg
            except Exception as e:
                print(f" PPGstandard curvefailed: {e}")

        # 2. calibration methodscompare
        if self.results:
            try:
                fig_comparison = self._create_calibration_comparison_plot()
                visualizations['calibration_comparison'] = fig_comparison
            except Exception as e:
                print(f" calibration methodscomparefailed: {e}")

        # 3. PPGindex
        if 'method3' in self.results and 'ppg_indices' in self.results['method3']:
            try:
                fig_distribution = self._create_ppg_index_distribution_plot()
                visualizations['index_distribution'] = fig_distribution
            except Exception as e:
                print(f" PPGindexfailed: {e}")

        # 4. analysis
        try:
            fig_error = self._create_error_analysis_plot()
            visualizations['error_analysis'] = fig_error
        except Exception as e:
            print(f" analysisfailed: {e}")

        # 5. vs
        if 'method3' in self.results and 'rt_predicted' in self.results['method3'] and 'rt_actual' in self.results[
            'method3']:
            try:
                fig_pred_vs_actual = self._create_prediction_vs_actual_plot()
                visualizations['prediction_vs_actual'] = fig_pred_vs_actual
            except Exception as e:
                print(f" vsfailed: {e}")

        # savefile
        if output_dir:
            try:
                Path(output_dir).mkdir(parents=True, exist_ok=True)

                for name, fig in visualizations.items():
                    file_path = os.path.join(output_dir, f"{name}.png")
                    fig.savefig(file_path, dpi=300, bbox_inches='tight')
                    print(f" savechart: {file_path}")

                # savePDF
                pdf_path = os.path.join(output_dir, "all_visualizations.pdf")
                with PdfPages(pdf_path) as pdf:
                    for name, fig in visualizations.items():
                        pdf.savefig(fig)
                    print(f" savePDF: {pdf_path}")

            except Exception as e:
                print(f" savechartfailed: {e}")

        self.visualizations = visualizations
        return True, f" {len(visualizations)} visualization charts", visualizations

    def _create_ppg_curves_plot(self) -> Figure:
        """PPGstandard curve"""
        n_conditions = len(self.ppg_data)
        n_cols = min(3, n_conditions)
        n_rows = (n_conditions + n_cols - 1) // n_cols

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows))
        if n_conditions == 1:
            axes = np.array([axes])
        if axes.ndim == 1:
            axes = axes.reshape(-1, n_cols)

        conditions = list(self.ppg_data.keys())

        for idx, condition in enumerate(conditions[:n_rows * n_cols]):
            row = idx // n_cols
            col = idx % n_cols
            ax = axes[row, col] if n_rows > 1 or n_cols > 1 else axes[col]

            df_ppg = self.ppg_data[condition]
            n_values = df_ppg['degree_of_polymerization'].values
            rt_values = df_ppg['retention_time'].values

            # translated note
            ax.scatter(n_values, rt_values, alpha=0.7, label='data', s=50)

            # translated note
            try:
                log_n = np.log(n_values)
                slope, intercept, r_value, _, _ = stats.linregress(log_n, rt_values)
                n_fit = np.linspace(n_values.min(), n_values.max(), 100)
                rt_fit = intercept + slope * np.log(n_fit)

                ax.plot(n_fit, rt_fit, 'r-', label=f' (R²={r_value ** 2:.4f})', linewidth=2)
            except:
                pass

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
        return fig

    def _create_calibration_comparison_plot(self) -> Figure:
        """calibration methodscompare"""
        methods = []
        mae_values = []
        rmse_values = []
        r2_values = []

        for method_key in ['method1', 'method2', 'method3', 'method4']:
            if method_key in self.results:
                result = self.results[method_key]
                methods.append(result.get('description', method_key))
                mae_values.append(result.get('MAE', np.nan))
                rmse_values.append(result.get('RMSE', np.nan))
                r2_values.append(result.get('R2', np.nan))

        if not methods:
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.text(0.5, 0.5, 'calibration methodsresults', ha='center', va='center', fontsize=14)
            ax.set_axis_off()
            return fig

        fig, axes = plt.subplots(1, 3, figsize=(15, 5))

        # MAEcompare
        bars1 = axes[0].bar(range(len(methods)), mae_values, color='skyblue', alpha=0.8)
        axes[0].set_title('calibration methodsMAEcompare')
        axes[0].set_ylabel('MAE (min)')
        axes[0].set_xticks(range(len(methods)))
        axes[0].set_xticklabels([m[:15] + '...' if len(m) > 15 else m for m in methods], rotation=45, ha='right')

        # translated note
        for bar, val in zip(bars1, mae_values):
            if not np.isnan(val):
                height = bar.get_height()
                axes[0].text(bar.get_x() + bar.get_width() / 2., height + 0.01,
                             f'{val:.3f}', ha='center', va='bottom', fontsize=9)

        # RMSEcompare
        bars2 = axes[1].bar(range(len(methods)), rmse_values, color='lightcoral', alpha=0.8)
        axes[1].set_title('calibration methodsRMSEcompare')
        axes[1].set_ylabel('RMSE (min)')
        axes[1].set_xticks(range(len(methods)))
        axes[1].set_xticklabels([m[:15] + '...' if len(m) > 15 else m for m in methods], rotation=45, ha='right')

        # translated note
        for bar, val in zip(bars2, rmse_values):
            if not np.isnan(val):
                height = bar.get_height()
                axes[1].text(bar.get_x() + bar.get_width() / 2., height + 0.01,
                             f'{val:.3f}', ha='center', va='bottom', fontsize=9)

        # R²compare
        bars3 = axes[2].bar(range(len(methods)), r2_values, color='lightgreen', alpha=0.8)
        axes[2].set_title('calibration methodsR²compare')
        axes[2].set_ylabel('R²')
        axes[2].set_xticks(range(len(methods)))
        axes[2].set_xticklabels([m[:15] + '...' if len(m) > 15 else m for m in methods], rotation=45, ha='right')
        axes[2].axhline(y=0.9, color='r', linestyle='--', alpha=0.5, label='R²=0.9')
        axes[2].legend()

        # translated note
        for bar, val in zip(bars3, r2_values):
            if not np.isnan(val):
                height = bar.get_height()
                axes[2].text(bar.get_x() + bar.get_width() / 2., height + 0.01,
                             f'{val:.3f}', ha='center', va='bottom', fontsize=9)

        plt.tight_layout()
        return fig

    def _create_ppg_index_distribution_plot(self) -> Figure:
        """PPGindex"""
        if 'method3' not in self.results or 'ppg_indices' not in self.results['method3']:
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.text(0.5, 0.5, 'PPGindexdata', ha='center', va='center', fontsize=14)
            ax.set_axis_off()
            return fig

        indices = self.results['method3']['ppg_indices']

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        # translated note
        ax1.hist(indices, bins=20, edgecolor='black', alpha=0.7, color='steelblue')
        ax1.axvline(x=np.mean(indices), color='r', linestyle='--', linewidth=2, label=f': {np.mean(indices):.1f}')
        ax1.axvline(x=np.median(indices), color='g', linestyle='--', linewidth=2,
                    label=f': {np.median(indices):.1f}')
        ax1.set_xlabel('PPG retention indices')
        ax1.set_ylabel('')
        ax1.set_title('PPGindex')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # translated note
        bp = ax2.boxplot(indices, vert=True, patch_artist=True,
                         boxprops=dict(facecolor='lightblue', color='blue'),
                         medianprops=dict(color='red', linewidth=2))

        # translated note
        stats_text = f":\n"
        stats_text += f": {len(indices)}\n"
        stats_text += f": {np.mean(indices):.1f}\n"
        stats_text += f": {np.std(indices):.1f}\n"
        stats_text += f": {np.min(indices):.1f}\n"
        stats_text += f": {np.max(indices):.1f}\n"
        stats_text += f": {np.median(indices):.1f}"

        ax2.text(0.05, 0.95, stats_text, transform=ax2.transAxes,
                 verticalalignment='top', fontsize=10,
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

        ax2.set_ylabel('PPG retention indices')
        ax2.set_title('PPGindex')
        ax2.set_xticks([1])
        ax2.set_xticklabels(['PPGindex'])
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        return fig

    def _create_error_analysis_plot(self) -> Figure:
        """analysis"""
        methods = []
        errors = []

        for method_key in ['method1', 'method2', 'method3', 'method4']:
            if method_key in self.results and 'MAE' in self.results[method_key]:
                result = self.results[method_key]
                methods.append(result.get('description', method_key))
                errors.append(result['MAE'])

        if not methods:
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.text(0.5, 0.5, 'data', ha='center', va='center', fontsize=14)
            ax.set_axis_off()
            return fig

        fig, ax = plt.subplots(figsize=(10, 6))

        colors = ['skyblue', 'lightcoral', 'lightgreen', 'gold']
        bars = ax.bar(range(len(methods)), errors, color=colors[:len(methods)], alpha=0.8)

        ax.set_xlabel('calibration methods')
        ax.set_ylabel('MAE (min)')
        ax.set_title('calibration methods')
        ax.set_xticks(range(len(methods)))
        ax.set_xticklabels([m[:20] + '...' if len(m) > 20 else m for m in methods], rotation=45, ha='right')
        ax.grid(True, alpha=0.3, axis='y')

        # translated note
        for bar, error in zip(bars, errors):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2., height + 0.001,
                    f'{error:.4f}', ha='center', va='bottom')

        plt.tight_layout()
        return fig

    def _create_prediction_vs_actual_plot(self) -> Figure:
        """vs"""
        if 'method3' not in self.results:
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.text(0.5, 0.5, 'data', ha='center', va='center', fontsize=14)
            ax.set_axis_off()
            return fig

        result = self.results['method3']

        if 'rt_predicted' not in result or 'rt_actual' not in result:
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.text(0.5, 0.5, 'data', ha='center', va='center', fontsize=14)
            ax.set_axis_off()
            return fig

        rt_predicted = result['rt_predicted']
        rt_actual = result['rt_actual']

        if not rt_predicted or not rt_actual:
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.text(0.5, 0.5, 'data', ha='center', va='center', fontsize=14)
            ax.set_axis_off()
            return fig

        # translated note
        min_len = min(len(rt_predicted), len(rt_actual))
        rt_predicted = rt_predicted[:min_len]
        rt_actual = rt_actual[:min_len]

        fig, ax = plt.subplots(figsize=(8, 6))

        # translated note
        ax.scatter(rt_actual, rt_predicted, alpha=0.6, s=50)

        # translated note
        min_val = min(min(rt_actual), min(rt_predicted))
        max_val = max(max(rt_actual), max(rt_predicted))
        ax.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.5, label='')

        # translated note
        if len(rt_actual) >= 2:
            slope, intercept, r_value, _, _ = stats.linregress(rt_actual, rt_predicted)
            x_line = np.array([min_val, max_val])
            y_line = intercept + slope * x_line
            ax.plot(x_line, y_line, 'g-', alpha=0.7, label=f' (R²={r_value ** 2:.4f})')

        ax.set_xlabel('retention_time (min)')
        ax.set_ylabel('retention_time (min)')
        ax.set_title('vsretention_time')
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        return fig

    def generate_report(self, output_dir: str = None) -> Tuple[bool, str, str]:
        """analysis"""
        print("analysis...")

        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        report = "=" * 70 + "\n"
        report += "PPG retention indicesanalysis\n"
        report += "=" * 70 + "\n\n"

        report += f": {timestamp}\n\n"

        report += "1. data\n"
        report += "-" * 40 + "\n"
        report += f"PPG data: {len(self.ppg_data)}\n"

        smrt_count = sum(len(df) for df in self.smrt_data.values())
        report += f"SMRT dataset: {smrt_count}\n"

        validation_count = sum(len(df) for df in self.validation_data.values())
        report += f"validation set: {validation_count}\n"

        calibration_count = sum(len(df) for df in self.calibration_data.values())
        report += f"calibrationcompound: {calibration_count}\n\n"

        report += "2. PPGstandard curve\n"
        report += "-" * 40 + "\n"
        if 'standard_curves' in self.calibration_methods:
            for condition, curve in self.calibration_methods['standard_curves'].items():
                report += f"{condition}:\n"
                report += f" model: {curve.get('model_type', '')}\n"
                report += f"  R²: {curve.get('r_squared', 0):.6f}\n"
                report += f" data: {len(curve.get('n_values', []))}\n"
        else:
            report += "PPGstandard curve\n"
        report += "\n"

        report += "3. calibration methodscompare\n"
        report += "-" * 40 + "\n"

        for method_key in ['method1', 'method2', 'method3', 'method4']:
            if method_key in self.results:
                result = self.results[method_key]
                report += f"{result.get('description', method_key)}:\n"

                if 'MAE' in result and not np.isnan(result['MAE']):
                    report += f"  MAE: {result['MAE']:.4f}\n"
                if 'RMSE' in result and not np.isnan(result['RMSE']):
                    report += f"  RMSE: {result['RMSE']:.4f}\n"
                if 'R2' in result and not np.isnan(result['R2']):
                    report += f"  R²: {result['R2']:.4f}\n"

                if 'error' in result:
                    report += f" : {result['error']}\n"

                report += "\n"

        report += "4. model performance evaluation\n"
        report += "-" * 40 + "\n"
        if self.models:
            for model_name, model_data in self.models.items():
                report += f"{model_name}model:\n"
                perf = model_data.get('performance', {})
                report += f" : {model_data.get('n_samples', 0)}\n"
                report += f"  MAE: {perf.get('MAE', 0):.4f}\n"
                report += f"  RMSE: {perf.get('RMSE', 0):.4f}\n"
                report += f"  R²: {perf.get('R2', 0):.4f}\n"
        else:
            report += "trainingmodel\n"
        report += "\n"

        report += "5. \n"
        report += "-" * 40 + "\n"

        # method
        best_method = None
        best_r2 = -1

        for method_key in ['method1', 'method2', 'method3', 'method4']:
            if method_key in self.results:
                result = self.results[method_key]
                r2 = result.get('R2', -1)
                if not np.isnan(r2) and r2 > best_r2:
                    best_r2 = r2
                    best_method = result.get('description', method_key)

        if best_method:
            report += f"calibration methods: {best_method} (R² = {best_r2:.4f})\n"

        # translated note
        report += "\n:\n"
        report += "1. UsePPGcalibration methodsretention_time\n"
        report += "2. PPGstandard curve\n"
        report += "3. chromatographic condition, UsePPGindex\n"
        report += "4. model\n"

        report += "\n" + "=" * 70 + "\n"
        report += "\n"
        report += "=" * 70 + "\n"

        # savefile
        if output_dir:
            try:
                Path(output_dir).mkdir(parents=True, exist_ok=True)
                report_path = os.path.join(output_dir,
                                           f"PPG_analysis_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt")

                with open(report_path, 'w', encoding='utf-8') as f:
                    f.write(report)

                return True, "", report_path
            except Exception as e:
                return False, f"savefailed: {str(e)}", report

        return True, "", report

    def save_all_results(self, output_dir: str) -> Tuple[bool, str, List[str]]:
        """saveresults"""
        try:
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)

            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            saved_files = []

            # 1. savePPG data
            if self.ppg_data:
                for condition, df in self.ppg_data.items():
                    file_name = f"PPG_{condition}_{timestamp}.csv"
                    file_path = output_path / file_name
                    df.to_csv(file_path, index=False, encoding='utf-8')
                    saved_files.append(str(file_path))

            # 2. savecompound data
            for data_type, data_dict in [('SMRT', self.smrt_data), ('Validation', self.validation_data),
                                         ('Calibration', self.calibration_data)]:
                if data_dict:
                    for key, df in data_dict.items():
                        file_name = f"{data_type}_{key}_{timestamp}.csv"
                        file_path = output_path / file_name
                        df.to_csv(file_path, index=False, encoding='utf-8')
                        saved_files.append(str(file_path))

            # 3. saveanalysis results
            if self.results:
                results_df = pd.DataFrame(self.results).T
                file_name = f"Analysis_Results_{timestamp}.csv"
                file_path = output_path / file_name
                results_df.to_csv(file_path, encoding='utf-8')
                saved_files.append(str(file_path))

            # 4. save
            success, msg, report_content = self.generate_report(output_dir)
            if success and isinstance(report_content, str) and os.path.exists(report_content):
                saved_files.append(report_content)

            # 5. savevisualization charts
            if self.visualizations:
                for name, fig in self.visualizations.items():
                    file_name = f"Visualization_{name}_{timestamp}.png"
                    file_path = output_path / file_name
                    fig.savefig(file_path, dpi=300, bbox_inches='tight')
                    saved_files.append(str(file_path))

            return True, f"resultssave {output_dir}", saved_files

        except Exception as e:
            return False, f"saveresultsfailed: {str(e)}", []


class PPGExperimentGUI:
    """PPGanalysisGUI"""

    def __init__(self, root):
        """GUI"""
        self.root = root
        self.root.title("PPG retention indicesanalysis - ")
        self.root.geometry("1400x900")

        # translated note
        try:
            self.root.iconbitmap(default='icon.ico')
        except:
            pass

        # analysis
        self.analyzer = PPGExperimentAnalyzer()

        # translated note
        self.processing_thread = None
        self.is_processing = False

        # loadcondition
        self.loaded_conditions = set()
        self.loaded_datasets = set()

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
        self.setup_visualization_tab()
        self.setup_report_tab()
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
        self.ppg_condition_var = tk.StringVar(value="C18_gradient1")
        ppg_condition_entry = ttk.Entry(ppg_frame, textvariable=self.ppg_condition_var, width=20)
        ppg_condition_entry.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Label(ppg_frame, text="PPG datafile:").grid(row=1, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.ppg_file_var = tk.StringVar()
        ppg_file_entry = ttk.Entry(ppg_frame, textvariable=self.ppg_file_var, width=60)
        ppg_file_entry.grid(row=1, column=1, sticky=(tk.W, tk.E), padx=(0, 10), pady=5)

        ttk.Button(ppg_frame, text="...", command=self.browse_ppg_file).grid(row=1, column=2, pady=5)

        ttk.Button(ppg_frame, text="loadPPG data", command=self.load_ppg_data).grid(row=2, column=0, columnspan=3,
                                                                                   pady=10)

        # ==================== SMRTdataload ====================
        smrt_frame = ttk.LabelFrame(data_frame, text="SMRT datasetload", padding=10)
        smrt_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(smrt_frame, text="conditionname:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.smrt_condition_var = tk.StringVar(value="default")
        smrt_condition_entry = ttk.Entry(smrt_frame, textvariable=self.smrt_condition_var, width=20)
        smrt_condition_entry.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Label(smrt_frame, text="SMRTdatafile:").grid(row=1, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.smrt_file_var = tk.StringVar()
        smrt_file_entry = ttk.Entry(smrt_frame, textvariable=self.smrt_file_var, width=60)
        smrt_file_entry.grid(row=1, column=1, sticky=(tk.W, tk.E), padx=(0, 10), pady=5)

        ttk.Button(smrt_frame, text="...", command=self.browse_smrt_file).grid(row=1, column=2, pady=5)

        ttk.Button(smrt_frame, text="loadSMRTdata", command=self.load_smrt_data).grid(row=2, column=0, columnspan=3,
                                                                                      pady=10)

        # ==================== validation setdataload ====================
        validation_frame = ttk.LabelFrame(data_frame, text="validation setdataload", padding=10)
        validation_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(validation_frame, text="conditionname:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.validation_condition_var = tk.StringVar(value="default")
        validation_condition_entry = ttk.Entry(validation_frame, textvariable=self.validation_condition_var, width=20)
        validation_condition_entry.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Label(validation_frame, text="validation setfile:").grid(row=1, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.validation_file_var = tk.StringVar()
        validation_file_entry = ttk.Entry(validation_frame, textvariable=self.validation_file_var, width=60)
        validation_file_entry.grid(row=1, column=1, sticky=(tk.W, tk.E), padx=(0, 10), pady=5)

        ttk.Button(validation_frame, text="...", command=self.browse_validation_file).grid(row=1, column=2, pady=5)

        ttk.Button(validation_frame, text="loadvalidation setdata", command=self.load_validation_data).grid(row=2, column=0,
                                                                                                    columnspan=3,
                                                                                                    pady=10)

        # ==================== calibration compound dataload ====================
        calibration_frame = ttk.LabelFrame(data_frame, text="calibration compound dataload", padding=10)
        calibration_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(calibration_frame, text="conditionname:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.calibration_condition_var = tk.StringVar(value="default")
        calibration_condition_entry = ttk.Entry(calibration_frame, textvariable=self.calibration_condition_var,
                                                width=20)
        calibration_condition_entry.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Label(calibration_frame, text="calibrationcompoundfile:").grid(row=1, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.calibration_file_var = tk.StringVar()
        calibration_file_entry = ttk.Entry(calibration_frame, textvariable=self.calibration_file_var, width=60)
        calibration_file_entry.grid(row=1, column=1, sticky=(tk.W, tk.E), padx=(0, 10), pady=5)

        ttk.Button(calibration_frame, text="...", command=self.browse_calibration_file).grid(row=1, column=2,
                                                                                                 pady=5)

        ttk.Button(calibration_frame, text="loadcalibration compound data", command=self.load_calibration_data).grid(row=2,
                                                                                                          column=0,
                                                                                                          columnspan=3,
                                                                                                          pady=10)

        # ==================== loaddata ====================
        overview_frame = ttk.LabelFrame(data_frame, text="loaddata", padding=10)
        overview_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 10))

        # Treeviewloaddata
        columns = ("data", "condition", "data", "")
        self.data_tree = ttk.Treeview(overview_frame, columns=columns, show="headings", height=8)

        # column
        for col in columns:
            self.data_tree.heading(col, text=col)
            self.data_tree.column(col, width=120)

        # translated note
        scrollbar = ttk.Scrollbar(overview_frame, orient=tk.VERTICAL, command=self.data_tree.yview)
        self.data_tree.configure(yscrollcommand=scrollbar.set)

        self.data_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # data
        button_frame = ttk.Frame(data_frame)
        button_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Button(button_frame, text="data", command=self.clear_all_data).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="data", command=self.generate_example_data).pack(side=tk.LEFT, padx=5)

        # data
        self.data_status_var = tk.StringVar(value="loaddata...")
        ttk.Label(data_frame, textvariable=self.data_status_var).pack(anchor=tk.W)

    def setup_analysis_tab(self):
        """analysis"""
        analysis_tab = ttk.Frame(self.notebook)
        self.notebook.add(analysis_tab, text="dataanalysis")

        # translated note
        analysis_frame = ttk.LabelFrame(analysis_tab, text="PPGanalysis", padding=15)
        analysis_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # ==================== PPGstandard curve ====================
        curve_frame = ttk.LabelFrame(analysis_frame, text="PPGstandard curve", padding=10)
        curve_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(curve_frame, text="condition:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.curve_condition_var = tk.StringVar()
        self.curve_condition_combo = ttk.Combobox(curve_frame, textvariable=self.curve_condition_var,
                                                  width=25, state="readonly")
        self.curve_condition_combo.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Label(curve_frame, text="model:").grid(row=0, column=2, sticky=tk.W, padx=(0, 5), pady=5)
        self.curve_model_var = tk.StringVar(value="logarithmic")
        curve_model_combo = ttk.Combobox(curve_frame, textvariable=self.curve_model_var,
                                         values=["logarithmic", "linear"],
                                         width=15, state="readonly")
        curve_model_combo.grid(row=0, column=3, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Button(curve_frame, text="PPGstandard curve", command=self.fit_ppg_curve).grid(row=0, column=4, padx=(20, 0),
                                                                                         pady=5)

        # ==================== PPGindexcalculate ====================
        calc_frame = ttk.LabelFrame(analysis_frame, text="PPGindexcalculate", padding=10)
        calc_frame.pack(fill=tk.X, pady=(0, 15))

        # calculatePPGindex
        ttk.Label(calc_frame, text="retention_time(min):").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.calc_rt_var = tk.StringVar(value="5.5")
        calc_rt_entry = ttk.Entry(calc_frame, textvariable=self.calc_rt_var, width=15)
        calc_rt_entry.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Label(calc_frame, text="condition:").grid(row=0, column=2, sticky=tk.W, padx=(0, 5), pady=5)
        self.calc_source_var = tk.StringVar()
        self.calc_source_combo = ttk.Combobox(calc_frame, textvariable=self.calc_source_var,
                                              width=20, state="readonly")
        self.calc_source_combo.grid(row=0, column=3, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Button(calc_frame, text="calculatePPGindex", command=self.calculate_ppg_index).grid(row=0, column=4, padx=(10, 0),
                                                                                          pady=5)

        # PPGindex
        self.ppg_index_var = tk.StringVar(value="")
        ttk.Label(calc_frame, text="PPGindex:").grid(row=1, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        ttk.Label(calc_frame, textvariable=self.ppg_index_var, font=("Arial", 10, "bold")).grid(row=1, column=1,
                                                                                                sticky=tk.W,
                                                                                                padx=(0, 10), pady=5)

        # PPGindexcalculateretention_time
        ttk.Label(calc_frame, text="PPGindex:").grid(row=2, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.rt_from_index_var = tk.StringVar(value="550")
        rt_from_index_entry = ttk.Entry(calc_frame, textvariable=self.rt_from_index_var, width=15)
        rt_from_index_entry.grid(row=2, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Label(calc_frame, text="condition:").grid(row=2, column=2, sticky=tk.W, padx=(0, 5), pady=5)
        self.calc_target_var = tk.StringVar()
        self.calc_target_combo = ttk.Combobox(calc_frame, textvariable=self.calc_target_var,
                                              width=20, state="readonly")
        self.calc_target_combo.grid(row=2, column=3, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Label(calc_frame, text="calculatemethod:").grid(row=2, column=4, sticky=tk.W, padx=(10, 5), pady=5)
        self.calc_method_var = tk.StringVar(value="regression")
        calc_method_combo = ttk.Combobox(calc_frame, textvariable=self.calc_method_var,
                                         values=["interpolation", "regression"],
                                         width=12, state="readonly")
        calc_method_combo.grid(row=2, column=5, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Button(calc_frame, text="calculateretention_time", command=self.calculate_rt_from_index).grid(row=2, column=6,
                                                                                               padx=(10, 0), pady=5)

        # calculateresults
        self.calc_rt_result_var = tk.StringVar(value="")
        ttk.Label(calc_frame, text="retention_time:").grid(row=3, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        ttk.Label(calc_frame, textvariable=self.calc_rt_result_var, font=("Arial", 10, "bold")).grid(row=3, column=1,
                                                                                                     sticky=tk.W,
                                                                                                     padx=(0, 10),
                                                                                                     pady=5)

        # ==================== calibration methodsanalysis ====================
        calibration_frame = ttk.LabelFrame(analysis_frame, text="calibration methodsanalysis", padding=10)
        calibration_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(calibration_frame, text="condition:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.calibration_source_var = tk.StringVar()
        self.calibration_source_combo = ttk.Combobox(calibration_frame, textvariable=self.calibration_source_var,
                                                     width=20, state="readonly")
        self.calibration_source_combo.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Label(calibration_frame, text="condition:").grid(row=0, column=2, sticky=tk.W, padx=(0, 5), pady=5)
        self.calibration_target_var = tk.StringVar()
        self.calibration_target_combo = ttk.Combobox(calibration_frame, textvariable=self.calibration_target_var,
                                                     width=20, state="readonly")
        self.calibration_target_combo.grid(row=0, column=3, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Button(calibration_frame, text="calibration", command=self.apply_calibration_methods).grid(row=0,
                                                                                                            column=4,
                                                                                                            padx=(
                                                                                                            20, 0),
                                                                                                            pady=5)

        # ==================== model performance evaluation ====================
        model_frame = ttk.LabelFrame(analysis_frame, text="model performance evaluation", padding=10)
        model_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(model_frame, text="model:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.model_type_var = tk.StringVar(value="linear")
        model_type_combo = ttk.Combobox(model_frame, textvariable=self.model_type_var,
                                        values=["linear", "ridge", "lasso"],
                                        width=15, state="readonly")
        model_type_combo.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Button(model_frame, text="model", command=self.evaluate_model_performance).grid(row=0, column=2,
                                                                                                   padx=(20, 0), pady=5)

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

        ttk.Button(button_frame, text="results", command=self.clear_analysis_text).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="analysis results", command=self.export_analysis_results).pack(side=tk.LEFT, padx=5)

        # analysis
        self.analysis_status_var = tk.StringVar(value="analysis...")
        ttk.Label(analysis_frame, textvariable=self.analysis_status_var).pack(anchor=tk.W)

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
        self.viz_type_var = tk.StringVar(value="ppg_curves")
        viz_type_combo = ttk.Combobox(options_frame, textvariable=self.viz_type_var,
                                      values=["ppg_curves", "calibration_comparison",
                                              "index_distribution", "error_analysis",
                                              "prediction_vs_actual", "all"],
                                      width=20, state="readonly")
        viz_type_combo.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        # chart
        ttk.Button(options_frame, text="chart", command=self.generate_visualization).grid(row=0, column=2,
                                                                                             padx=(20, 0), pady=5)

        # savechart
        ttk.Button(options_frame, text="savechart", command=self.save_visualization).grid(row=0, column=3, padx=(10, 0),
                                                                                         pady=5)

        # ==================== chart ====================
        display_frame = ttk.LabelFrame(viz_frame, text="chart", padding=10)
        display_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 10))

        # translated note
        self.figure_canvas = None
        self.figure_toolbar = None
        self.current_figure = None
        self.current_viz_type = None

        # translated note
        self.viz_placeholder = ttk.Label(display_frame, text="chart",
                                         font=("Arial", 14), foreground="gray")
        self.viz_placeholder.pack(expand=True)

        # translated note
        button_frame = ttk.Frame(viz_frame)
        button_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Button(button_frame, text="chart", command=self.clear_visualization).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="savechart", command=self.save_all_visualizations).pack(side=tk.LEFT, padx=5)

        # translated note
        self.viz_status_var = tk.StringVar(value="chart...")
        ttk.Label(viz_frame, textvariable=self.viz_status_var).pack(anchor=tk.W)

    def setup_report_tab(self):
        """"""
        report_tab = ttk.Frame(self.notebook)
        self.notebook.add(report_tab, text="output")

        # translated note
        report_frame = ttk.LabelFrame(report_tab, text="output", padding=15)
        report_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # ==================== output ====================
        output_frame = ttk.LabelFrame(report_frame, text="output", padding=10)
        output_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(output_frame, text="outputdirectory:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.output_dir_var = tk.StringVar(value=os.path.join(os.getcwd(), "PPG_Results"))
        output_dir_entry = ttk.Entry(output_frame, textvariable=self.output_dir_var, width=60)
        output_dir_entry.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=(0, 10), pady=5)

        ttk.Button(output_frame, text="...", command=self.browse_output_dir).grid(row=0, column=2, pady=5)

        # ==================== ====================
        preview_frame = ttk.LabelFrame(report_frame, text="", padding=10)
        preview_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 15))

        # translated note
        self.report_text = scrolledtext.ScrolledText(preview_frame, width=80, height=15,
                                                     wrap=tk.WORD, font=("Consolas", 10))
        self.report_text.pack(fill=tk.BOTH, expand=True)

        # ==================== output ====================
        control_frame = ttk.Frame(report_frame)
        control_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Button(control_frame, text="", command=self.generate_report).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_frame, text="saveresults", command=self.save_all_results).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_frame, text="outputdirectory", command=self.open_output_dir).pack(side=tk.LEFT, padx=5)

        # output
        self.output_status_var = tk.StringVar(value="outputresults...")
        ttk.Label(report_frame, textvariable=self.output_status_var).pack(anchor=tk.W)

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

    def browse_smrt_file(self):
        """SMRTdatafile"""
        file_types = [("datafile", "*.csv *.xlsx *.xls"), ("CSVfile", "*.csv"),
                      ("Excelfile", "*.xlsx *.xls"), ("file", "*.*")]
        file_path = filedialog.askopenfilename(title="SMRTdatafile", filetypes=file_types)

        if file_path:
            self.smrt_file_var.set(file_path)

    def browse_validation_file(self):
        """validation setfile"""
        file_types = [("datafile", "*.csv *.xlsx *.xls"), ("CSVfile", "*.csv"),
                      ("Excelfile", "*.xlsx *.xls"), ("file", "*.*")]
        file_path = filedialog.askopenfilename(title="validation setfile", filetypes=file_types)

        if file_path:
            self.validation_file_var.set(file_path)

    def browse_calibration_file(self):
        """calibrationcompoundfile"""
        file_types = [("datafile", "*.csv *.xlsx *.xls"), ("CSVfile", "*.csv"),
                      ("Excelfile", "*.xlsx *.xls"), ("file", "*.*")]
        file_path = filedialog.askopenfilename(title="calibrationcompoundfile", filetypes=file_types)

        if file_path:
            self.calibration_file_var.set(file_path)

    def browse_output_dir(self):
        """outputdirectory"""
        dir_path = filedialog.askdirectory(title="outputdirectory")

        if dir_path:
            self.output_dir_var.set(dir_path)

    def load_ppg_data(self):
        """loadPPG data"""
        ppg_file = self.ppg_file_var.get().strip()
        condition = self.ppg_condition_var.get().strip()

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
            success, msg = self.analyzer.load_ppg_data(ppg_file, condition)

            if success:
                self.root.after(0, lambda: self.log_message(f"✓ {msg}", "SUCCESS"))
                self.root.after(0, self.update_data_tree)
                self.root.after(0, self.update_condition_comboboxes)
                self.loaded_conditions.add(condition)
            else:
                self.root.after(0, lambda: self.log_message(f"✗ {msg}", "ERROR"))

            self.root.after(0, lambda: self.update_status(""))

        except Exception as e:
            self.root.after(0, lambda: self.log_message(f"✗ loadPPG datafailed: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("loadfailed"))

    def load_smrt_data(self):
        """loadSMRTdata"""
        smrt_file = self.smrt_file_var.get().strip()
        condition = self.smrt_condition_var.get().strip()

        if not smrt_file:
            messagebox.showwarning("Warning", "SMRTdatafile")
            return

        if not condition:
            messagebox.showwarning("Warning", "inputconditionname")
            return

        # load
        self.log_message(f"loadSMRTdata: {smrt_file} (condition: {condition})", "INFO")
        self.update_status(f"loadSMRTdata: {Path(smrt_file).name}")

        # loaddata
        self.processing_thread = threading.Thread(
            target=self._load_smrt_data_thread,
            args=(smrt_file, condition)
        )
        self.processing_thread.start()

    def _load_smrt_data_thread(self, smrt_file, condition):
        """loadSMRTdata"""
        try:
            success, msg = self.analyzer.load_compound_data(smrt_file, "smrt", condition)

            if success:
                self.root.after(0, lambda: self.log_message(f"✓ {msg}", "SUCCESS"))
                self.root.after(0, self.update_data_tree)
                self.root.after(0, self.update_condition_comboboxes)
                self.loaded_datasets.add(f"smrt_{condition}")
            else:
                self.root.after(0, lambda: self.log_message(f"✗ {msg}", "ERROR"))

            self.root.after(0, lambda: self.update_status(""))

        except Exception as e:
            self.root.after(0, lambda: self.log_message(f"✗ loadSMRTdatafailed: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("loadfailed"))

    def load_validation_data(self):
        """loadvalidation setdata"""
        validation_file = self.validation_file_var.get().strip()
        condition = self.validation_condition_var.get().strip()

        if not validation_file:
            messagebox.showwarning("Warning", "validation setfile")
            return

        if not condition:
            messagebox.showwarning("Warning", "inputconditionname")
            return

        # load
        self.log_message(f"loadvalidation setdata: {validation_file} (condition: {condition})", "INFO")
        self.update_status(f"loadvalidation setdata: {Path(validation_file).name}")

        # loaddata
        self.processing_thread = threading.Thread(
            target=self._load_validation_data_thread,
            args=(validation_file, condition)
        )
        self.processing_thread.start()

    def _load_validation_data_thread(self, validation_file, condition):
        """loadvalidation setdata"""
        try:
            success, msg = self.analyzer.load_compound_data(validation_file, "validation", condition)

            if success:
                self.root.after(0, lambda: self.log_message(f"✓ {msg}", "SUCCESS"))
                self.root.after(0, self.update_data_tree)
                self.root.after(0, self.update_condition_comboboxes)
                self.loaded_datasets.add(f"validation_{condition}")
            else:
                self.root.after(0, lambda: self.log_message(f"✗ {msg}", "ERROR"))

            self.root.after(0, lambda: self.update_status(""))

        except Exception as e:
            self.root.after(0, lambda: self.log_message(f"✗ loadvalidation setdatafailed: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("loadfailed"))

    def load_calibration_data(self):
        """loadcalibration compound data"""
        calibration_file = self.calibration_file_var.get().strip()
        condition = self.calibration_condition_var.get().strip()

        if not calibration_file:
            messagebox.showwarning("Warning", "calibrationcompoundfile")
            return

        if not condition:
            messagebox.showwarning("Warning", "inputconditionname")
            return

        # load
        self.log_message(f"loadcalibration compound data: {calibration_file} (condition: {condition})", "INFO")
        self.update_status(f"loadcalibration compound data: {Path(calibration_file).name}")

        # loaddata
        self.processing_thread = threading.Thread(
            target=self._load_calibration_data_thread,
            args=(calibration_file, condition)
        )
        self.processing_thread.start()

    def _load_calibration_data_thread(self, calibration_file, condition):
        """loadcalibration compound data"""
        try:
            success, msg = self.analyzer.load_compound_data(calibration_file, "calibration", condition)

            if success:
                self.root.after(0, lambda: self.log_message(f"✓ {msg}", "SUCCESS"))
                self.root.after(0, self.update_data_tree)
                self.root.after(0, self.update_condition_comboboxes)
                self.loaded_datasets.add(f"calibration_{condition}")
            else:
                self.root.after(0, lambda: self.log_message(f"✗ {msg}", "ERROR"))

            self.root.after(0, lambda: self.update_status(""))

        except Exception as e:
            self.root.after(0, lambda: self.log_message(f"✗ loadcalibration compound datafailed: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("loadfailed"))

    def update_data_tree(self):
        """data"""
        # data
        for item in self.data_tree.get_children():
            self.data_tree.delete(item)

        # PPG data
        for condition, df in self.analyzer.ppg_data.items():
            self.data_tree.insert("", tk.END, values=("PPGstandard", condition, len(df), "load"))

        # SMRTdata
        for key, df in self.analyzer.smrt_data.items():
            condition = key.replace('smrt_', '')
            self.data_tree.insert("", tk.END, values=("SMRT dataset", condition, len(df), "load"))

        # data
        for key, df in self.analyzer.validation_data.items():
            condition = key.replace('validation_', '')
            self.data_tree.insert("", tk.END, values=("validation set", condition, len(df), "load"))

        # calibration compound data
        for key, df in self.analyzer.calibration_data.items():
            condition = key.replace('calibration_', '')
            self.data_tree.insert("", tk.END, values=("calibrationcompound", condition, len(df), "load"))

        # translated note
        total_data = (len(self.analyzer.ppg_data) + len(self.analyzer.smrt_data) +
                      len(self.analyzer.validation_data) + len(self.analyzer.calibration_data))
        self.data_status_var.set(f"load {total_data} data")

    def update_condition_comboboxes(self):
        """condition"""
        conditions = list(self.analyzer.ppg_data.keys())

        # condition
        self.curve_condition_combo['values'] = conditions
        if conditions and not self.curve_condition_var.get():
            self.curve_condition_var.set(conditions[0])

        # calculatecondition
        self.calc_source_combo['values'] = conditions
        if conditions and not self.calc_source_var.get():
            self.calc_source_var.set(conditions[0])

        # calculatecondition
        self.calc_target_combo['values'] = conditions
        if conditions and not self.calc_target_var.get():
            self.calc_target_var.set(
                conditions[0] if len(conditions) == 1 else conditions[1] if len(conditions) > 1 else "")

        # calibrationcondition
        self.calibration_source_combo['values'] = conditions
        if conditions and not self.calibration_source_var.get():
            self.calibration_source_var.set(conditions[0])

        # calibrationcondition
        self.calibration_target_combo['values'] = conditions
        if conditions and not self.calibration_target_var.get():
            self.calibration_target_var.set(
                conditions[0] if len(conditions) == 1 else conditions[1] if len(conditions) > 1 else "")

    def clear_all_data(self):
        """data"""
        if messagebox.askyesno("", "loaddata？"):
            self.analyzer = PPGExperimentAnalyzer()
            self.loaded_conditions.clear()
            self.loaded_datasets.clear()
            self.update_data_tree()
            self.update_condition_comboboxes()
            self.log_message("data", "INFO")

    def generate_example_data(self):
        """data"""
        if messagebox.askyesno("", "data, data？"):
            self.clear_all_data()
            self.log_message("data...", "INFO")

            # PPG data
            for i, condition in enumerate(['C18_gradient1', 'C18_gradient2', 'C18_gradient3']):
                n_values = np.arange(2, 31)  # PPG2-30
                if i == 0:
                    rt_values = 2.0 + 0.5 * np.log(n_values) + np.random.normal(0, 0.1, len(n_values))
                elif i == 1:
                    rt_values = 2.5 + 0.6 * np.log(n_values) + np.random.normal(0, 0.1, len(n_values))
                else:
                    rt_values = 3.0 + 0.7 * np.log(n_values) + np.random.normal(0, 0.1, len(n_values))

                df_ppg = pd.DataFrame({'degree_of_polymerization': n_values, 'retention_time': rt_values})
                self.analyzer.ppg_data[condition] = df_ppg

            # SMRTdata
            n_smrt = 28
            compound_names = [f'Compound_{i + 1}' for i in range(n_smrt)]
            rt_source = np.random.uniform(2, 15, n_smrt)
            rt_target = rt_source * 1.1 + np.random.normal(0, 0.2, n_smrt)
            logp_values = np.random.uniform(-2, 5, n_smrt)

            self.analyzer.smrt_data['smrt_default'] = pd.DataFrame({
                'compound_name': compound_names,
                'retention_time': rt_source,
                'RT_target': rt_target,
                'logP': logp_values,
                'molecule': np.random.uniform(100, 500, n_smrt)
            })

            # data
            n_validation = 50
            validation_names = [f'Val_Compound_{i + 1}' for i in range(n_validation)]
            rt_val_source = np.random.uniform(2, 18, n_validation)
            rt_val_target = rt_val_source * 1.05 + np.random.normal(0, 0.15, n_validation)
            logp_val = np.random.uniform(-2, 5, n_validation)

            self.analyzer.validation_data['validation_default'] = pd.DataFrame({
                'compound_name': validation_names,
                'retention_time': rt_val_source,
                'RT_target': rt_val_target,
                'logP': logp_val,
                'molecule': np.random.uniform(100, 500, n_validation)
            })

            # calibration compound data
            n_calib = 28
            calib_names = [f'Calib_Compound_{i + 1}' for i in range(n_calib)]
            rt_calib_source = np.random.uniform(2, 15, n_calib)
            rt_calib_target = rt_calib_source * 1.08 + np.random.normal(0, 0.1, n_calib)

            self.analyzer.calibration_data['calibration_default'] = pd.DataFrame({
                'compound_name': calib_names,
                'RT_source': rt_calib_source,
                'RT_target': rt_calib_target
            })

            self.update_data_tree()
            self.update_condition_comboboxes()
            self.log_message("data", "SUCCESS")

    def fit_ppg_curve(self):
        """PPGstandard curve"""
        condition = self.curve_condition_var.get()

        if not condition:
            messagebox.showwarning("Warning", "condition")
            return

        # analysis
        self.analysis_message(f"PPGstandard curve (condition: {condition})", "INFO")
        self.update_status(f"PPGstandard curve: {condition}")

        # translated note
        self.processing_thread = threading.Thread(
            target=self._fit_ppg_curve_thread,
            args=(condition,)
        )
        self.processing_thread.start()

    def _fit_ppg_curve_thread(self, condition):
        """PPGstandard curve"""
        try:
            success, msg, curve_data = self.analyzer.fit_ppg_standard_curve(condition)

            if success:
                self.root.after(0, lambda: self.analysis_message(f"✓ {msg}", "SUCCESS"))

                # results
                result_text = f"PPGstandard curveresults - {condition}:\n"
                result_text += f" model: {curve_data.get('model_type', '')}\n"
                result_text += f"  R²: {curve_data.get('r_squared', 0):.6f}\n"

                if 'params' in curve_data:
                    params = curve_data['params']
                    if curve_data.get('model_type') == 'logarithmic':
                        result_text += f" : {params[0]:.6f}\n"
                        result_text += f" : {params[1]:.6f}\n"
                        result_text += f" : RT = {params[0]:.4f} + {params[1]:.4f} * ln(n)\n"
                    else:
                        result_text += f" : {params[0]:.6f}\n"
                        result_text += f" : {params[1]:.6f}\n"
                        result_text += f" : RT = {params[0]:.4f} + {params[1]:.4f} * n\n"

                self.root.after(0, lambda: self.analysis_message(result_text, "INFO"))
            else:
                self.root.after(0, lambda: self.analysis_message(f"✗ {msg}", "ERROR"))

            self.root.after(0, lambda: self.update_status(""))

        except Exception as e:
            self.root.after(0, lambda: self.analysis_message(f"✗ PPGstandard curvefailed: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("failed"))

    def calculate_ppg_index(self):
        """calculatePPGindex"""
        try:
            rt = float(self.calc_rt_var.get().strip())
            condition = self.calc_source_var.get()

            if not condition:
                messagebox.showwarning("Warning", "condition")
                return

            # calculatePPGindex
            ppg_index = self.analyzer.calculate_ppg_index(rt, condition)

            # results
            self.ppg_index_var.set(f"{ppg_index:.1f}")
            self.analysis_message(f"retention_time {rt} min → PPGindex: {ppg_index:.1f}", "SUCCESS")

        except ValueError:
            messagebox.showerror("", "inputretention_time")
        except Exception as e:
            messagebox.showerror("", f"calculatePPGindexfailed: {str(e)}")

    def calculate_rt_from_index(self):
        """PPGindexcalculateretention_time"""
        try:
            index = float(self.rt_from_index_var.get().strip())
            condition = self.calc_target_var.get()
            method = self.calc_method_var.get()

            if not condition:
                messagebox.showwarning("Warning", "condition")
                return

            # calculateretention_time
            rt = self.analyzer.calculate_rt_from_index(index, condition, method)

            # results
            self.calc_rt_result_var.set(f"{rt:.2f} min")
            self.analysis_message(f"PPGindex {index:.1f} → retention_time: {rt:.2f} min (method: {method})", "SUCCESS")

        except ValueError:
            messagebox.showerror("", "inputPPGindex")
        except Exception as e:
            messagebox.showerror("", f"calculateretention_timefailed: {str(e)}")

    def apply_calibration_methods(self):
        """calibration methods"""
        source_condition = self.calibration_source_var.get()
        target_condition = self.calibration_target_var.get()

        if not source_condition:
            messagebox.showwarning("Warning", "condition")
            return

        if not target_condition:
            messagebox.showwarning("Warning", "condition")
            return

        if source_condition == target_condition:
            if not messagebox.askyesno("", "conditioncondition, ？"):
                return

        # analysis
        self.analysis_message(f"calibration ({source_condition} → {target_condition})...", "INFO")
        self.update_status(f"calibration methods: {source_condition}→{target_condition}")

        # calibration
        self.processing_thread = threading.Thread(
            target=self._apply_calibration_methods_thread,
            args=(source_condition, target_condition)
        )
        self.processing_thread.start()

    def _apply_calibration_methods_thread(self, source_condition, target_condition):
        """calibration methods"""
        try:
            success, msg = self.analyzer.apply_calibration_methods(source_condition, target_condition)

            if success:
                self.root.after(0, lambda: self.analysis_message(f"✓ {msg}", "SUCCESS"))

                # results
                result_text = "calibration methodsanalysis results:\n"
                result_text += "-" * 40 + "\n"

                for method_key in ['method1', 'method2', 'method3', 'method4']:
                    if method_key in self.analyzer.results:
                        result = self.analyzer.results[method_key]
                        result_text += f"{result.get('description', method_key)}:\n"

                        if 'MAE' in result and not np.isnan(result['MAE']):
                            result_text += f"  MAE: {result['MAE']:.4f}\n"
                        if 'RMSE' in result and not np.isnan(result['RMSE']):
                            result_text += f"  RMSE: {result['RMSE']:.4f}\n"
                        if 'R2' in result and not np.isnan(result['R2']):
                            result_text += f"  R²: {result['R2']:.4f}\n"

                        if 'error' in result:
                            result_text += f" : {result['error']}\n"

                        result_text += "\n"

                self.root.after(0, lambda: self.analysis_message(result_text, "INFO"))
            else:
                self.root.after(0, lambda: self.analysis_message(f"✗ {msg}", "ERROR"))

            self.root.after(0, lambda: self.update_status(""))

        except Exception as e:
            self.root.after(0, lambda: self.analysis_message(f"✗ calibration methodsfailed: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("calibrationfailed"))

    def evaluate_model_performance(self):
        """model"""
        model_type = self.model_type_var.get()

        # analysis
        self.analysis_message(f"{model_type}model...", "INFO")
        self.update_status(f"model: {model_type}")

        # model
        self.processing_thread = threading.Thread(
            target=self._evaluate_model_performance_thread,
            args=(model_type,)
        )
        self.processing_thread.start()

    def _evaluate_model_performance_thread(self, model_type):
        """model"""
        try:
            success, msg = self.analyzer.evaluate_model_performance(model_type)

            if success:
                self.root.after(0, lambda: self.analysis_message(f"✓ {msg}", "SUCCESS"))

                # model
                if model_type in self.analyzer.models:
                    model_data = self.analyzer.models[model_type]
                    perf = model_data.get('performance', {})

                    result_text = f"{model_type}model:\n"
                    result_text += f" : {model_data.get('n_samples', 0)}\n"
                    result_text += f"  MAE: {perf.get('MAE', 0):.4f}\n"
                    result_text += f"  RMSE: {perf.get('RMSE', 0):.4f}\n"
                    result_text += f"  R²: {perf.get('R2', 0):.4f}\n"

                    self.root.after(0, lambda: self.analysis_message(result_text, "INFO"))
            else:
                self.root.after(0, lambda: self.analysis_message(f"✗ {msg}", "ERROR"))

            self.root.after(0, lambda: self.update_status(""))

        except Exception as e:
            self.root.after(0, lambda: self.analysis_message(f"✗ modelfailed: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("failed"))

    def generate_visualization(self):
        """visualization charts"""
        viz_type = self.viz_type_var.get()

        if viz_type == "all":
            # chart
            viz_types = ["ppg_curves", "calibration_comparison", "index_distribution",
                         "error_analysis", "prediction_vs_actual"]
        else:
            viz_types = [viz_type]

        # translated note
        self.viz_status_var.set(f" {viz_type} chart...")
        self.update_status(f"chart: {viz_type}")

        # chart
        self.processing_thread = threading.Thread(
            target=self._generate_visualization_thread,
            args=(viz_types,)
        )
        self.processing_thread.start()

    def _generate_visualization_thread(self, viz_types):
        """visualization charts"""
        try:
            # chart
            success, msg, visualizations = self.analyzer.generate_visualizations()

            if success and visualizations:
                # chart
                if viz_types:
                    first_type = viz_types[0]
                    if first_type in visualizations:
                        self.root.after(0, lambda: self.display_figure(visualizations[first_type], first_type))
                    elif "ppg_curves" in visualizations:
                        self.root.after(0, lambda: self.display_figure(visualizations["ppg_curves"], "ppg_curves"))

                self.root.after(0, lambda: self.viz_status_var.set(f" {len(visualizations)} chart"))
                self.root.after(0, lambda: self.log_message(f"✓ {msg}", "SUCCESS"))
            else:
                self.root.after(0, lambda: self.viz_status_var.set("chartfailed"))
                self.root.after(0, lambda: self.log_message(f"✗ {msg}", "ERROR"))

            self.root.after(0, lambda: self.update_status(""))

        except Exception as e:
            self.root.after(0, lambda: self.viz_status_var.set(f"chartfailed: {str(e)}"))
            self.root.after(0, lambda: self.log_message(f"✗ chartfailed: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("failed"))

    def display_figure(self, fig, viz_type):
        """chart"""
        # chart
        self.clear_visualization()

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
        self.current_viz_type = viz_type

        # translated note
        self.viz_status_var.set(f"chart: {viz_type}")

    def save_visualization(self):
        """savechart"""
        if self.current_figure is None:
            messagebox.showwarning("Warning", "savechart")
            return

        file_types = [("PNGfile", "*.png"), ("PDFfile", "*.pdf"),
                      ("SVGfile", "*.svg"), ("file", "*.*")]

        default_name = f"PPG_{self.current_viz_type if self.current_viz_type else 'chart'}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

        file_path = filedialog.asksaveasfilename(
            title="savechart",
            filetypes=file_types,
            defaultextension=".png",
            initialfile=default_name
        )

        if file_path:
            try:
                self.current_figure.savefig(file_path, dpi=300, bbox_inches='tight')
                self.log_message(f"✓ chartsave: {file_path}", "SUCCESS")
                self.viz_status_var.set(f"chartsave: {Path(file_path).name}")
            except Exception as e:
                self.log_message(f"✗ savechartfailed: {str(e)}", "ERROR")

    def save_all_visualizations(self):
        """savechart"""
        if not self.analyzer.visualizations:
            messagebox.showwarning("Warning", "savechart")
            return

        output_dir = filedialog.askdirectory(title="savedirectory")

        if output_dir:
            try:
                # chart
                success, msg, visualizations = self.analyzer.generate_visualizations(output_dir)

                if success:
                    self.log_message(f"✓ chartsave: {output_dir}", "SUCCESS")
                    self.viz_status_var.set(f"save {len(visualizations)} chart")
                else:
                    self.log_message(f"✗ {msg}", "ERROR")

            except Exception as e:
                self.log_message(f"✗ savechartfailed: {str(e)}", "ERROR")

    def clear_visualization(self):
        """chart"""
        if self.figure_canvas:
            self.figure_canvas.get_tk_widget().destroy()
            self.figure_canvas = None

        if self.figure_toolbar:
            self.figure_toolbar.destroy()
            self.figure_toolbar = None

        if self.current_figure:
            plt.close(self.current_figure)
            self.current_figure = None
            self.current_viz_type = None

        # translated note
        self.viz_placeholder.pack(expand=True)

    def generate_report(self):
        """"""
        # translated note
        self.update_status("...")

        # translated note
        self.processing_thread = threading.Thread(
            target=self._generate_report_thread
        )
        self.processing_thread.start()

    def _generate_report_thread(self):
        """"""
        try:
            # translated note
            success, msg, report_content = self.analyzer.generate_report()

            if success:
                # GUI
                if isinstance(report_content, str) and os.path.exists(report_content):
                    # file
                    with open(report_content, 'r', encoding='utf-8') as f:
                        report_text = f.read()
                else:
                    # Use
                    report_text = report_content

                self.root.after(0, lambda: self.report_text.delete(1.0, tk.END))
                self.root.after(0, lambda: self.report_text.insert(tk.END, report_text))

                self.root.after(0, lambda: self.output_status_var.set(""))
                self.root.after(0, lambda: self.log_message(f"✓ {msg}", "SUCCESS"))
            else:
                self.root.after(0, lambda: self.output_status_var.set("failed"))
                self.root.after(0, lambda: self.log_message(f"✗ {msg}", "ERROR"))

            self.root.after(0, lambda: self.update_status(""))

        except Exception as e:
            self.root.after(0, lambda: self.output_status_var.set("failed"))
            self.root.after(0, lambda: self.log_message(f"✗ failed: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("failed"))

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
            success, msg, saved_files = self.analyzer.save_all_results(output_dir)

            if success:
                self.root.after(0, lambda: self.output_status_var.set(f"resultssave: {len(saved_files)} file"))
                self.root.after(0, lambda: self.log_message(f"✓ {msg}", "SUCCESS"))

                # savefile
                self.root.after(0, lambda: self.report_text.delete(1.0, tk.END))
                self.root.after(0, lambda: self.report_text.insert(tk.END, "savefilecolumn:\n"))
                self.root.after(0, lambda: self.report_text.insert(tk.END, "=" * 60 + "\n\n"))

                for file_path in saved_files:
                    self.root.after(0, lambda fp=file_path: self.report_text.insert(tk.END, f"• {Path(fp).name}\n"))

                self.root.after(0, lambda: self.report_text.insert(tk.END, f"\nfilesave: {output_dir}"))
            else:
                self.root.after(0, lambda: self.output_status_var.set("savefailed"))
                self.root.after(0, lambda: self.log_message(f"✗ {msg}", "ERROR"))

            self.root.after(0, lambda: self.update_status(""))

        except Exception as e:
            self.root.after(0, lambda: self.output_status_var.set("savefailed"))
            self.root.after(0, lambda: self.log_message(f"✗ saveresultsfailed: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("savefailed"))

    def export_analysis_results(self):
        """analysis results"""
        if not self.analyzer.results:
            messagebox.showwarning("Warning", "analysis results")
            return

        file_path = filedialog.asksaveasfilename(
            title="analysis results",
            defaultextension=".csv",
            filetypes=[("CSVfile", "*.csv"), ("Excelfile", "*.xlsx"), ("file", "*.*")]
        )

        if file_path:
            try:
                # resultsDataFrame
                results_list = []
                for method_key, result in self.analyzer.results.items():
                    if method_key.startswith('method'):
                        row = {
                            'method': result.get('description', method_key),
                            'MAE': result.get('MAE', np.nan),
                            'RMSE': result.get('RMSE', np.nan),
                            'R2': result.get('R2', np.nan)
                        }
                        results_list.append(row)

                if results_list:
                    results_df = pd.DataFrame(results_list)

                    if file_path.endswith('.xlsx'):
                        results_df.to_excel(file_path, index=False)
                    else:
                        results_df.to_csv(file_path, index=False, encoding='utf-8')

                    self.log_message(f"✓ analysis results: {file_path}", "SUCCESS")

            except Exception as e:
                self.log_message(f"✗ analysis resultsfailed: {str(e)}", "ERROR")

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

    def clear_analysis_text(self):
        """analysis"""
        self.analysis_text.delete(1.0, tk.END)

    def clear_log(self):
        """"""
        self.log_text.delete(1.0, tk.END)

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
        'scikit-learn': '',
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
    root.title("PPG retention indicesanalysis - ")

    # translated note
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()

    # (85%)
    window_width = int(screen_width * 0.85)
    window_height = int(screen_height * 0.85)

    # calculate ()
    x = (screen_width - window_width) // 2
    y = (screen_height - window_height) // 2

    root.geometry(f"{window_width}x{window_height}+{x}+{y}")

    # GUI
    app = PPGExperimentGUI(root)

    # translated note
    root.mainloop()


if __name__ == "__main__":
    print("PPG retention indicesanalysis - ")
    print("=" * 70)
    print(":")
    print(" 1. Data loading and management (PPG, SMRT, validation set, calibrationcompound)")
    print(" 2. PPG retention index calculation and conversion")
    print(" 3. analysis and comparison of five calibration methods")
    print(" 4. model performance evaluation")
    print(" 5. visualization generation")
    print(" 6. experimental report export")
    print("=" * 70)

    main()