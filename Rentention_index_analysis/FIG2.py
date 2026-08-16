#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PPG保留指数分析系统 - 完整GUI版本

功能：
1. 数据加载与管理（PPG数据、SMRT数据集、验证集）
2. PPG保留指数计算与转换
3. 五种校正方法分析比较
4. 模型性能评估
5. 可视化图表生成
6. 实验报告导出
"""

import os
import sys
import threading
import traceback
import warnings
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any

# 抑制警告
warnings.filterwarnings('ignore')

# 导入必要的库
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

    matplotlib.use('TkAgg')  # 设置matplotlib后端
except ImportError as e:
    print(f"错误: 请先安装必要的库: {e}")
    print("安装命令: pip install pandas numpy scipy matplotlib seaborn scikit-learn")
    sys.exit(1)

# 尝试导入GUI库
try:
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox, scrolledtext, Toplevel, StringVar, BooleanVar, IntVar

    GUI_AVAILABLE = True
except ImportError:
    GUI_AVAILABLE = False
    print("警告: tkinter未安装，GUI不可用")

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
sns.set_style("whitegrid")


class PPGExperimentAnalyzer:
    """PPG试验方案分析器核心类"""

    def __init__(self):
        self.ppg_data = {}  # 不同条件下的PPG数据
        self.smrt_data = {}  # SMRT数据集
        self.validation_data = {}  # 验证集数据
        self.calibration_data = {}  # 校正化合物数据
        self.ppg_indices = {}  # PPG保留指数
        self.models = {}  # 训练的模型
        self.results = {}  # 分析结果
        self.calibration_methods = {}  # 校正方法
        self.visualizations = {}  # 可视化图表

        # 实验参数
        self.experiment_params = {
            'ppg_model_type': 'logarithmic',  # PPG模型类型
            'calibration_method': 'linear',  # 校正方法
            'validation_split': 0.2,  # 验证集比例
            'random_seed': 42  # 随机种子
        }

    def load_ppg_data(self, file_path: str, condition: str) -> Tuple[bool, str]:
        """加载PPG数据"""
        try:
            file_ext = Path(file_path).suffix.lower()

            if file_ext in ['.xlsx', '.xls']:
                df = pd.read_excel(file_path)
            elif file_ext == '.csv':
                df = pd.read_csv(file_path)
            else:
                return False, f"不支持的文件格式: {file_ext}"

            # 查找必要的列
            column_mapping = {
                '聚合度': ['聚合度', 'n', 'DP', '聚合度n', 'PPG_n', 'PPG'],
                '保留时间': ['保留时间', 'RT', 'RetentionTime', 't_R', '保留时间(RT)', 'rt']
            }

            # 重命名列
            rename_dict = {}
            for target_col, possible_names in column_mapping.items():
                for name in possible_names:
                    if name in df.columns:
                        rename_dict[name] = target_col
                        break

            if rename_dict:
                df = df.rename(columns=rename_dict)

            # 检查必要列
            if '聚合度' not in df.columns or '保留时间' not in df.columns:
                return False, "数据文件中缺少必要的列（需要'聚合度'和'保留时间'）"

            # 确保数据类型正确
            df['聚合度'] = pd.to_numeric(df['聚合度'], errors='coerce')
            df['保留时间'] = pd.to_numeric(df['保留时间'], errors='coerce')
            df = df.dropna(subset=['聚合度', '保留时间'])

            # 按聚合度排序
            df = df.sort_values('聚合度')

            # 存储数据
            self.ppg_data[condition] = df

            return True, f"成功加载 {len(df)} 个PPG标准品数据（条件: {condition}）"

        except Exception as e:
            return False, f"加载PPG数据失败: {str(e)}"

    def load_compound_data(self, file_path: str, data_type: str, condition: str = "default") -> Tuple[bool, str]:
        """加载化合物数据"""
        try:
            file_ext = Path(file_path).suffix.lower()

            if file_ext in ['.xlsx', '.xls']:
                df = pd.read_excel(file_path)
            elif file_ext == '.csv':
                df = pd.read_csv(file_path)
            else:
                return False, f"不支持的文件格式: {file_ext}"

            # 查找必要的列
            column_mapping = {
                '化合物名称': ['化合物名称', '名称', '化合物', 'Name', 'Compound', '化合物名'],
                '保留时间': ['保留时间', 'RT', 'RetentionTime', 't_R', '保留时间(RT)', 'rt'],
                'CAS': ['CAS', 'CAS号', 'CAS No.', 'CAS号'],
                '分子量': ['分子量', 'MW', 'MolecularWeight', '分子量(Da)'],
                'logP': ['logP', 'LogP', 'log P', '疏水性']
            }

            # 重命名列
            rename_dict = {}
            for target_col, possible_names in column_mapping.items():
                for name in possible_names:
                    if name in df.columns:
                        rename_dict[name] = target_col
                        break

            if rename_dict:
                df = df.rename(columns=rename_dict)

            # 检查必要列
            if '化合物名称' not in df.columns or '保留时间' not in df.columns:
                return False, "数据文件中缺少必要的列（需要'化合物名称'和'保留时间'）"

            # 确保数据类型正确
            df['保留时间'] = pd.to_numeric(df['保留时间'], errors='coerce')
            df = df.dropna(subset=['化合物名称', '保留时间'])

            # 存储数据
            key = f"{data_type}_{condition}"

            if data_type == "smrt":
                self.smrt_data[key] = df
            elif data_type == "validation":
                self.validation_data[key] = df
            elif data_type == "calibration":
                self.calibration_data[key] = df
            else:
                return False, f"不支持的数据类型: {data_type}"

            return True, f"成功加载 {len(df)} 个{data_type}化合物数据（条件: {condition}）"

        except Exception as e:
            return False, f"加载化合物数据失败: {str(e)}"

    def calculate_ppg_index(self, rt: float, condition: str) -> float:
        """计算PPG保留指数"""
        if condition not in self.ppg_data:
            raise ValueError(f"未找到条件 {condition} 的PPG数据")

        df_ppg = self.ppg_data[condition]
        ppg_rt = df_ppg['保留时间'].values
        ppg_n = df_ppg['聚合度'].values

        # 边界处理
        if rt < ppg_rt[0]:
            # 外推
            if len(ppg_rt) >= 2:
                n_calc = ppg_n[0] - (ppg_rt[0] - rt) / (ppg_rt[1] - ppg_rt[0]) * (ppg_n[1] - ppg_n[0])
            else:
                n_calc = ppg_n[0]
        elif rt > ppg_rt[-1]:
            # 外推
            if len(ppg_rt) >= 2:
                n_calc = ppg_n[-1] + (rt - ppg_rt[-1]) / (ppg_rt[-1] - ppg_rt[-2]) * (ppg_n[-1] - ppg_n[-2])
            else:
                n_calc = ppg_n[-1]
        else:
            # 线性插值
            idx = np.searchsorted(ppg_rt, rt) - 1
            if idx < 0:
                idx = 0
            elif idx >= len(ppg_rt) - 1:
                idx = len(ppg_rt) - 2

            rt_i, rt_j = ppg_rt[idx], ppg_rt[idx + 1]
            n_i, n_j = ppg_n[idx], ppg_n[idx + 1]

            n_calc = n_i + (n_j - n_i) * (rt - rt_i) / (rt_j - rt_i)

        # 标准化为保留指数（乘以100）
        return n_calc * 100

    def calculate_rt_from_index(self, index: float, condition: str, method: str = "interpolation") -> float:
        """从PPG保留指数计算回保留时间"""
        if condition not in self.ppg_data:
            raise ValueError(f"未找到条件 {condition} 的PPG数据")

        df_ppg = self.ppg_data[condition]
        ppg_rt = df_ppg['保留时间'].values
        ppg_n = df_ppg['聚合度'].values
        n_target = index / 100  # 将指数转换回聚合度

        if method == "interpolation":
            # 线性插值法
            if n_target < ppg_n[0]:
                # 外推
                if len(ppg_n) >= 2:
                    rt_calc = ppg_rt[0] - (ppg_n[0] - n_target) / (ppg_n[1] - ppg_n[0]) * (ppg_rt[1] - ppg_rt[0])
                else:
                    rt_calc = ppg_rt[0]
            elif n_target > ppg_n[-1]:
                # 外推
                if len(ppg_n) >= 2:
                    rt_calc = ppg_rt[-1] + (n_target - ppg_n[-1]) / (ppg_n[-1] - ppg_n[-2]) * (ppg_rt[-1] - ppg_rt[-2])
                else:
                    rt_calc = ppg_rt[-1]
            else:
                # 线性插值
                idx = np.searchsorted(ppg_n, n_target) - 1
                if idx < 0:
                    idx = 0
                elif idx >= len(ppg_n) - 1:
                    idx = len(ppg_n) - 2

                n_i, n_j = ppg_n[idx], ppg_n[idx + 1]
                rt_i, rt_j = ppg_rt[idx], ppg_rt[idx + 1]

                rt_calc = rt_i + (rt_j - rt_i) * (n_target - n_i) / (n_j - n_i)

        elif method == "regression":
            # 回归法 - 使用PPG标准曲线
            # 拟合对数模型: RT = a + b * ln(n)
            log_n = np.log(ppg_n)

            # 线性回归
            slope, intercept, r_value, p_value, std_err = stats.linregress(log_n, ppg_rt)

            # 预测保留时间
            rt_calc = intercept + slope * np.log(n_target)

            # 存储回归参数
            if 'regression_params' not in self.calibration_methods:
                self.calibration_methods['regression_params'] = {}
            self.calibration_methods['regression_params'][condition] = {
                'slope': slope,
                'intercept': intercept,
                'r_squared': r_value ** 2,
                'std_err': std_err
            }

        else:
            raise ValueError(f"不支持的计算方法: {method}")

        return rt_calc

    def fit_ppg_standard_curve(self, condition: str) -> Tuple[bool, str, Dict]:
        """拟合PPG标准曲线"""
        if condition not in self.ppg_data:
            return False, f"未找到条件 {condition} 的PPG数据", {}

        df_ppg = self.ppg_data[condition]
        n_values = df_ppg['聚合度'].values
        rt_values = df_ppg['保留时间'].values

        # 尝试多种模型拟合
        models = {}

        # 1. 对数模型: RT = a + b * ln(n)
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

        # 2. 线性模型: RT = a + b * n
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

        # 选择最佳模型
        best_model_name = None
        best_r2 = -1

        for model_name, model_data in models.items():
            if 'r_squared' in model_data and model_data['r_squared'] > best_r2:
                best_r2 = model_data['r_squared']
                best_model_name = model_name

        if best_model_name is None:
            return False, "无法拟合任何模型", {}

        best_model = models[best_model_name]

        # 存储标准曲线
        if 'standard_curves' not in self.calibration_methods:
            self.calibration_methods['standard_curves'] = {}

        self.calibration_methods['standard_curves'][condition] = {
            'model_type': best_model_name,
            'params': best_model['params'],
            'r_squared': best_model['r_squared'],
            'n_values': n_values,
            'rt_values': rt_values
        }

        return True, f"成功拟合PPG标准曲线: {best_model_name}模型, R²={best_model['r_squared']:.6f}", best_model

    def apply_calibration_methods(self, source_condition: str, target_condition: str) -> Tuple[bool, str]:
        """应用五种校正方式"""
        print(f"应用校正方式分析 ({source_condition} → {target_condition}):")

        if len(self.smrt_data) == 0:
            return False, "没有SMRT数据集"

        # 获取SMRT数据
        smrt_key = f"smrt_{source_condition}"
        if smrt_key not in self.smrt_data:
            # 尝试查找任何SMRT数据
            smrt_keys = [k for k in self.smrt_data.keys() if k.startswith('smrt_')]
            if not smrt_keys:
                return False, "没有找到SMRT数据"
            smrt_key = smrt_keys[0]

        smrt_df = self.smrt_data[smrt_key]

        # 获取验证数据
        validation_key = f"validation_{target_condition}"
        if validation_key not in self.validation_data:
            # 尝试查找任何验证数据
            validation_keys = [k for k in self.validation_data.keys() if k.startswith('validation_')]
            if not validation_keys:
                # 如果没有验证数据，使用SMRT数据
                validation_df = smrt_df
            else:
                validation_key = validation_keys[0]
                validation_df = self.validation_data[validation_key]
        else:
            validation_df = self.validation_data[validation_key]

        # 确保数据包含必要的列
        if '保留时间' not in smrt_df.columns:
            return False, "SMRT数据缺少'保留时间'列"

        # 方法1: 直接使用，不校正
        print("1. 方法1: 直接使用SMRT数据集，不校正")
        # 假设源条件和目标条件相同，或者数据已经包含目标RT
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
                'description': '直接使用SMRT数据集，不校正'
            }
            print(f"   MAE: {mae:.4f}, RMSE: {rmse:.4f}, R²: {r2:.4f}")
        else:
            self.results['method1'] = {
                'MAE': np.nan,
                'RMSE': np.nan,
                'R2': np.nan,
                'description': '直接使用SMRT数据集，不校正',
                'error': '缺少目标保留时间数据'
            }

        # 方法2: 使用28个化合物校正
        print("2. 方法2: 使用28个化合物校正")
        # 查找校正化合物数据
        calib_key = f"calibration_{source_condition}"
        if calib_key in self.calibration_data:
            calib_df = self.calibration_data[calib_key]

            if 'RT_source' in calib_df.columns and 'RT_target' in calib_df.columns:
                rt_source_calib = calib_df['RT_source'].values
                rt_target_calib = calib_df['RT_target'].values

                # 线性校正模型
                try:
                    slope, intercept, r_value, p_value, std_err = stats.linregress(rt_source_calib, rt_target_calib)

                    # 应用校正到所有SMRT化合物
                    rt_source_all = smrt_df['保留时间'].values
                    rt_corrected = intercept + slope * rt_source_all

                    # 如果有目标数据，评估性能
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
                        'description': '使用化合物校正'
                    }
                    print(f"   MAE: {mae:.4f}, RMSE: {rmse:.4f}, R²: {r2:.4f}")
                except Exception as e:
                    self.results['method2'] = {
                        'MAE': np.nan,
                        'RMSE': np.nan,
                        'R2': np.nan,
                        'description': '使用化合物校正',
                        'error': f'校正失败: {str(e)}'
                    }
            else:
                self.results['method2'] = {
                    'MAE': np.nan,
                    'RMSE': np.nan,
                    'R2': np.nan,
                    'description': '使用化合物校正',
                    'error': '校正数据缺少RT_source或RT_target列'
                }
        else:
            self.results['method2'] = {
                'MAE': np.nan,
                'RMSE': np.nan,
                'R2': np.nan,
                'description': '使用化合物校正',
                'error': '没有校正化合物数据'
            }

        # 方法3: 使用PPG校正
        print("3. 方法3: 使用PPG校正")
        if source_condition in self.ppg_data and target_condition in self.ppg_data:
            # 对验证数据计算PPG指数
            ppg_indices = []
            rt_target_actual = []
            rt_target_predicted = []

            for _, row in validation_df.iterrows():
                if '保留时间' in row:
                    rt_source = row['保留时间']

                    try:
                        # 计算PPG指数
                        ppg_index = self.calculate_ppg_index(rt_source, source_condition)

                        # 从指数计算目标RT
                        rt_pred = self.calculate_rt_from_index(ppg_index, target_condition, method='regression')

                        ppg_indices.append(ppg_index)
                        rt_target_predicted.append(rt_pred)

                        # 如果有实际目标RT
                        if 'RT_target' in row:
                            rt_target_actual.append(row['RT_target'])
                    except Exception as e:
                        print(f"   计算失败: {e}")
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
                    'description': '使用PPG校正'
                }
                print(f"   MAE: {mae:.4f}, RMSE: {rmse:.4f}, R²: {r2:.4f}")
            else:
                self.results['method3'] = {
                    'MAE': np.nan,
                    'RMSE': np.nan,
                    'R2': np.nan,
                    'description': '使用PPG校正',
                    'error': '无法计算PPG指数或预测RT'
                }
        else:
            self.results['method3'] = {
                'MAE': np.nan,
                'RMSE': np.nan,
                'R2': np.nan,
                'description': '使用PPG校正',
                'error': f'缺少PPG数据: 源条件={source_condition}, 目标条件={target_condition}'
            }

        # 方法4: 使用PPG校正模型预测保留时间
        print("4. 方法4: 使用PPG校正模型预测保留时间")
        # 首先拟合PPG标准曲线
        success, msg, ppg_model = self.fit_ppg_standard_curve(target_condition)

        if success and ppg_model and 'func' in ppg_model:
            rt_target_predicted = []
            rt_target_actual = []

            for _, row in validation_df.iterrows():
                if '保留时间' in row:
                    rt_source = row['保留时间']

                    try:
                        # 计算PPG指数
                        ppg_index = self.calculate_ppg_index(rt_source, source_condition)
                        n_value = ppg_index / 100

                        # 使用PPG模型预测RT
                        if 'params' in ppg_model:
                            rt_pred = ppg_model['func'](n_value, *ppg_model['params'])
                        else:
                            rt_pred = ppg_model['func'](n_value)

                        rt_target_predicted.append(rt_pred)

                        # 如果有实际目标RT
                        if 'RT_target' in row:
                            rt_target_actual.append(row['RT_target'])
                    except Exception as e:
                        print(f"   预测失败: {e}")
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
                    'description': '使用PPG校正模型预测保留时间',
                    'model_type': ppg_model.get('model_type', 'unknown')
                }
                print(f"   MAE: {mae:.4f}, RMSE: {rmse:.4f}, R²: {r2:.4f}")
            else:
                self.results['method4'] = {
                    'MAE': np.nan,
                    'RMSE': np.nan,
                    'R2': np.nan,
                    'description': '使用PPG校正模型预测保留时间',
                    'error': '无法预测保留时间'
                }
        else:
            self.results['method4'] = {
                'MAE': np.nan,
                'RMSE': np.nan,
                'R2': np.nan,
                'description': '使用PPG校正模型预测保留时间',
                'error': f'无法拟合PPG标准曲线: {msg}'
            }

        # 方法5: 对比PPG校正和化合物校正的区别
        print("5. 方法5: 对比PPG校正和化合物校正的区别")
        # 这里可以对比方法2和方法3的结果
        comparison_results = {}

        if 'method2' in self.results and 'method3' in self.results:
            mae_diff = abs(self.results['method2'].get('MAE', np.nan) - self.results['method3'].get('MAE', np.nan))
            r2_diff = abs(self.results['method2'].get('R2', np.nan) - self.results['method3'].get('R2', np.nan))

            comparison_results = {
                'MAE_difference': mae_diff,
                'R2_difference': r2_diff,
                'recommendation': 'PPG校正' if self.results['method3'].get('R2', 0) > self.results['method2'].get('R2',
                                                                                                                  0) else '化合物校正'
            }

        self.results['method5'] = {
            'description': '对比PPG校正和化合物校正的区别',
            'comparison': comparison_results
        }

        print("校正方式分析完成!")
        return True, "校正方式分析完成"

    def evaluate_model_performance(self, model_type: str = 'linear') -> Tuple[bool, str]:
        """评估模型性能"""
        print(f"评估{model_type}模型性能...")

        # 这里可以实现不同模型的性能评估
        # 暂时用简单回归模型演示

        # 收集所有数据
        all_data = []

        # 从SMRT数据收集
        for key, df in self.smrt_data.items():
            for _, row in df.iterrows():
                if '保留时间' in row and 'logP' in row:
                    all_data.append({
                        'rt': row['保留时间'],
                        'logp': row['logP'] if not pd.isna(row['logP']) else 0
                    })

        # 从验证数据收集
        for key, df in self.validation_data.items():
            for _, row in df.iterrows():
                if '保留时间' in row and 'logP' in row:
                    all_data.append({
                        'rt': row['保留时间'],
                        'logp': row['logP'] if not pd.isna(row['logP']) else 0
                    })

        if len(all_data) < 10:
            return False, f"数据不足，只有 {len(all_data)} 个有效样本"

        # 准备特征和目标
        features = np.array([[d['logp']] for d in all_data])
        targets = np.array([d['rt'] for d in all_data])

        # 划分训练集和测试集
        X_train, X_test, y_train, y_test = train_test_split(
            features, targets, test_size=0.2, random_state=42
        )

        # 标准化
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)

        # 线性回归模型
        model = LinearRegression()
        model.fit(X_train_scaled, y_train)

        # 预测
        y_pred = model.predict(X_test_scaled)

        # 评估
        mae = mean_absolute_error(y_test, y_pred)
        rmse = np.sqrt(mean_squared_error(y_test, y_pred))
        r2 = r2_score(y_test, y_pred)

        print(f"  测试集性能:")
        print(f"    MAE: {mae:.4f}")
        print(f"    RMSE: {rmse:.4f}")
        print(f"    R²: {r2:.4f}")

        self.models[model_type] = {
            'model': model,
            'scaler': scaler,
            'performance': {'MAE': mae, 'RMSE': rmse, 'R2': r2},
            'n_samples': len(all_data)
        }

        return True, f"模型评估完成: MAE={mae:.4f}, RMSE={rmse:.4f}, R²={r2:.4f}"

    def generate_visualizations(self, output_dir: str = None) -> Tuple[bool, str, Dict]:
        """生成可视化图表"""
        print("生成可视化图表...")

        visualizations = {}

        # 1. PPG标准曲线图
        if self.ppg_data:
            try:
                fig_ppg = self._create_ppg_curves_plot()
                visualizations['ppg_curves'] = fig_ppg
            except Exception as e:
                print(f"  生成PPG标准曲线图失败: {e}")

        # 2. 校正方法比较图
        if self.results:
            try:
                fig_comparison = self._create_calibration_comparison_plot()
                visualizations['calibration_comparison'] = fig_comparison
            except Exception as e:
                print(f"  生成校正方法比较图失败: {e}")

        # 3. PPG指数分布图
        if 'method3' in self.results and 'ppg_indices' in self.results['method3']:
            try:
                fig_distribution = self._create_ppg_index_distribution_plot()
                visualizations['index_distribution'] = fig_distribution
            except Exception as e:
                print(f"  生成PPG指数分布图失败: {e}")

        # 4. 误差分析图
        try:
            fig_error = self._create_error_analysis_plot()
            visualizations['error_analysis'] = fig_error
        except Exception as e:
            print(f"  生成误差分析图失败: {e}")

        # 5. 预测vs实际图
        if 'method3' in self.results and 'rt_predicted' in self.results['method3'] and 'rt_actual' in self.results[
            'method3']:
            try:
                fig_pred_vs_actual = self._create_prediction_vs_actual_plot()
                visualizations['prediction_vs_actual'] = fig_pred_vs_actual
            except Exception as e:
                print(f"  生成预测vs实际图失败: {e}")

        # 保存到文件
        if output_dir:
            try:
                Path(output_dir).mkdir(parents=True, exist_ok=True)

                for name, fig in visualizations.items():
                    file_path = os.path.join(output_dir, f"{name}.png")
                    fig.savefig(file_path, dpi=300, bbox_inches='tight')
                    print(f"  保存图表: {file_path}")

                # 保存为PDF
                pdf_path = os.path.join(output_dir, "all_visualizations.pdf")
                with PdfPages(pdf_path) as pdf:
                    for name, fig in visualizations.items():
                        pdf.savefig(fig)
                    print(f"  保存PDF: {pdf_path}")

            except Exception as e:
                print(f"  保存图表失败: {e}")

        self.visualizations = visualizations
        return True, f"生成 {len(visualizations)} 个可视化图表", visualizations

    def _create_ppg_curves_plot(self) -> Figure:
        """创建PPG标准曲线图"""
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
            n_values = df_ppg['聚合度'].values
            rt_values = df_ppg['保留时间'].values

            # 绘制散点
            ax.scatter(n_values, rt_values, alpha=0.7, label='实测数据', s=50)

            # 尝试拟合对数曲线
            try:
                log_n = np.log(n_values)
                slope, intercept, r_value, _, _ = stats.linregress(log_n, rt_values)
                n_fit = np.linspace(n_values.min(), n_values.max(), 100)
                rt_fit = intercept + slope * np.log(n_fit)

                ax.plot(n_fit, rt_fit, 'r-', label=f'对数拟合 (R²={r_value ** 2:.4f})', linewidth=2)
            except:
                pass

            ax.set_xlabel('聚合度 (n)')
            ax.set_ylabel('保留时间 (min)')
            ax.set_title(f'PPG标准曲线 - {condition}')
            ax.legend()
            ax.grid(True, alpha=0.3)

        # 隐藏多余的子图
        for idx in range(len(conditions), n_rows * n_cols):
            row = idx // n_cols
            col = idx % n_cols
            axes[row, col].set_visible(False)

        plt.tight_layout()
        return fig

    def _create_calibration_comparison_plot(self) -> Figure:
        """创建校正方法比较图"""
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
            ax.text(0.5, 0.5, '无校正方法结果', ha='center', va='center', fontsize=14)
            ax.set_axis_off()
            return fig

        fig, axes = plt.subplots(1, 3, figsize=(15, 5))

        # MAE比较
        bars1 = axes[0].bar(range(len(methods)), mae_values, color='skyblue', alpha=0.8)
        axes[0].set_title('校正方法MAE比较')
        axes[0].set_ylabel('MAE (min)')
        axes[0].set_xticks(range(len(methods)))
        axes[0].set_xticklabels([m[:15] + '...' if len(m) > 15 else m for m in methods], rotation=45, ha='right')

        # 添加数值标签
        for bar, val in zip(bars1, mae_values):
            if not np.isnan(val):
                height = bar.get_height()
                axes[0].text(bar.get_x() + bar.get_width() / 2., height + 0.01,
                             f'{val:.3f}', ha='center', va='bottom', fontsize=9)

        # RMSE比较
        bars2 = axes[1].bar(range(len(methods)), rmse_values, color='lightcoral', alpha=0.8)
        axes[1].set_title('校正方法RMSE比较')
        axes[1].set_ylabel('RMSE (min)')
        axes[1].set_xticks(range(len(methods)))
        axes[1].set_xticklabels([m[:15] + '...' if len(m) > 15 else m for m in methods], rotation=45, ha='right')

        # 添加数值标签
        for bar, val in zip(bars2, rmse_values):
            if not np.isnan(val):
                height = bar.get_height()
                axes[1].text(bar.get_x() + bar.get_width() / 2., height + 0.01,
                             f'{val:.3f}', ha='center', va='bottom', fontsize=9)

        # R²比较
        bars3 = axes[2].bar(range(len(methods)), r2_values, color='lightgreen', alpha=0.8)
        axes[2].set_title('校正方法R²比较')
        axes[2].set_ylabel('R²')
        axes[2].set_xticks(range(len(methods)))
        axes[2].set_xticklabels([m[:15] + '...' if len(m) > 15 else m for m in methods], rotation=45, ha='right')
        axes[2].axhline(y=0.9, color='r', linestyle='--', alpha=0.5, label='R²=0.9')
        axes[2].legend()

        # 添加数值标签
        for bar, val in zip(bars3, r2_values):
            if not np.isnan(val):
                height = bar.get_height()
                axes[2].text(bar.get_x() + bar.get_width() / 2., height + 0.01,
                             f'{val:.3f}', ha='center', va='bottom', fontsize=9)

        plt.tight_layout()
        return fig

    def _create_ppg_index_distribution_plot(self) -> Figure:
        """创建PPG指数分布图"""
        if 'method3' not in self.results or 'ppg_indices' not in self.results['method3']:
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.text(0.5, 0.5, '无PPG指数数据', ha='center', va='center', fontsize=14)
            ax.set_axis_off()
            return fig

        indices = self.results['method3']['ppg_indices']

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        # 直方图
        ax1.hist(indices, bins=20, edgecolor='black', alpha=0.7, color='steelblue')
        ax1.axvline(x=np.mean(indices), color='r', linestyle='--', linewidth=2, label=f'均值: {np.mean(indices):.1f}')
        ax1.axvline(x=np.median(indices), color='g', linestyle='--', linewidth=2,
                    label=f'中位数: {np.median(indices):.1f}')
        ax1.set_xlabel('PPG保留指数')
        ax1.set_ylabel('频数')
        ax1.set_title('PPG指数分布')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # 箱线图
        bp = ax2.boxplot(indices, vert=True, patch_artist=True,
                         boxprops=dict(facecolor='lightblue', color='blue'),
                         medianprops=dict(color='red', linewidth=2))

        # 添加统计信息
        stats_text = f"统计信息:\n"
        stats_text += f"样本数: {len(indices)}\n"
        stats_text += f"均值: {np.mean(indices):.1f}\n"
        stats_text += f"标准差: {np.std(indices):.1f}\n"
        stats_text += f"最小值: {np.min(indices):.1f}\n"
        stats_text += f"最大值: {np.max(indices):.1f}\n"
        stats_text += f"中位数: {np.median(indices):.1f}"

        ax2.text(0.05, 0.95, stats_text, transform=ax2.transAxes,
                 verticalalignment='top', fontsize=10,
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

        ax2.set_ylabel('PPG保留指数')
        ax2.set_title('PPG指数箱线图')
        ax2.set_xticks([1])
        ax2.set_xticklabels(['PPG指数'])
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        return fig

    def _create_error_analysis_plot(self) -> Figure:
        """创建误差分析图"""
        methods = []
        errors = []

        for method_key in ['method1', 'method2', 'method3', 'method4']:
            if method_key in self.results and 'MAE' in self.results[method_key]:
                result = self.results[method_key]
                methods.append(result.get('description', method_key))
                errors.append(result['MAE'])

        if not methods:
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.text(0.5, 0.5, '无误差数据', ha='center', va='center', fontsize=14)
            ax.set_axis_off()
            return fig

        fig, ax = plt.subplots(figsize=(10, 6))

        colors = ['skyblue', 'lightcoral', 'lightgreen', 'gold']
        bars = ax.bar(range(len(methods)), errors, color=colors[:len(methods)], alpha=0.8)

        ax.set_xlabel('校正方法')
        ax.set_ylabel('MAE (min)')
        ax.set_title('不同校正方法的预测误差')
        ax.set_xticks(range(len(methods)))
        ax.set_xticklabels([m[:20] + '...' if len(m) > 20 else m for m in methods], rotation=45, ha='right')
        ax.grid(True, alpha=0.3, axis='y')

        # 添加数值标签
        for bar, error in zip(bars, errors):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2., height + 0.001,
                    f'{error:.4f}', ha='center', va='bottom')

        plt.tight_layout()
        return fig

    def _create_prediction_vs_actual_plot(self) -> Figure:
        """创建预测vs实际图"""
        if 'method3' not in self.results:
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.text(0.5, 0.5, '无预测数据', ha='center', va='center', fontsize=14)
            ax.set_axis_off()
            return fig

        result = self.results['method3']

        if 'rt_predicted' not in result or 'rt_actual' not in result:
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.text(0.5, 0.5, '无预测或实际数据', ha='center', va='center', fontsize=14)
            ax.set_axis_off()
            return fig

        rt_predicted = result['rt_predicted']
        rt_actual = result['rt_actual']

        if not rt_predicted or not rt_actual:
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.text(0.5, 0.5, '数据为空', ha='center', va='center', fontsize=14)
            ax.set_axis_off()
            return fig

        # 确保长度一致
        min_len = min(len(rt_predicted), len(rt_actual))
        rt_predicted = rt_predicted[:min_len]
        rt_actual = rt_actual[:min_len]

        fig, ax = plt.subplots(figsize=(8, 6))

        # 散点图
        ax.scatter(rt_actual, rt_predicted, alpha=0.6, s=50)

        # 添加对角线
        min_val = min(min(rt_actual), min(rt_predicted))
        max_val = max(max(rt_actual), max(rt_predicted))
        ax.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.5, label='理想预测线')

        # 添加回归线
        if len(rt_actual) >= 2:
            slope, intercept, r_value, _, _ = stats.linregress(rt_actual, rt_predicted)
            x_line = np.array([min_val, max_val])
            y_line = intercept + slope * x_line
            ax.plot(x_line, y_line, 'g-', alpha=0.7, label=f'回归线 (R²={r_value ** 2:.4f})')

        ax.set_xlabel('实际保留时间 (min)')
        ax.set_ylabel('预测保留时间 (min)')
        ax.set_title('预测vs实际保留时间')
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        return fig

    def generate_report(self, output_dir: str = None) -> Tuple[bool, str, str]:
        """生成分析报告"""
        print("生成分析报告...")

        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        report = "=" * 70 + "\n"
        report += "PPG保留指数综合分析报告\n"
        report += "=" * 70 + "\n\n"

        report += f"报告生成时间: {timestamp}\n\n"

        report += "1. 数据概览\n"
        report += "-" * 40 + "\n"
        report += f"PPG数据集数量: {len(self.ppg_data)}\n"

        smrt_count = sum(len(df) for df in self.smrt_data.values())
        report += f"SMRT数据集大小: {smrt_count}\n"

        validation_count = sum(len(df) for df in self.validation_data.values())
        report += f"验证集大小: {validation_count}\n"

        calibration_count = sum(len(df) for df in self.calibration_data.values())
        report += f"校正化合物数量: {calibration_count}\n\n"

        report += "2. PPG标准曲线\n"
        report += "-" * 40 + "\n"
        if 'standard_curves' in self.calibration_methods:
            for condition, curve in self.calibration_methods['standard_curves'].items():
                report += f"{condition}:\n"
                report += f"  模型类型: {curve.get('model_type', '未知')}\n"
                report += f"  R²: {curve.get('r_squared', 0):.6f}\n"
                report += f"  数据点数: {len(curve.get('n_values', []))}\n"
        else:
            report += "未拟合PPG标准曲线\n"
        report += "\n"

        report += "3. 校正方法性能比较\n"
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
                    report += f"  错误: {result['error']}\n"

                report += "\n"

        report += "4. 模型性能评估\n"
        report += "-" * 40 + "\n"
        if self.models:
            for model_name, model_data in self.models.items():
                report += f"{model_name}模型:\n"
                perf = model_data.get('performance', {})
                report += f"  样本数: {model_data.get('n_samples', 0)}\n"
                report += f"  MAE: {perf.get('MAE', 0):.4f}\n"
                report += f"  RMSE: {perf.get('RMSE', 0):.4f}\n"
                report += f"  R²: {perf.get('R2', 0):.4f}\n"
        else:
            report += "未训练模型\n"
        report += "\n"

        report += "5. 结论与建议\n"
        report += "-" * 40 + "\n"

        # 找出最佳方法
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
            report += f"最佳校正方法: {best_method} (R² = {best_r2:.4f})\n"

        # 建议
        report += "\n建议:\n"
        report += "1. 推荐使用PPG校正方法进行保留时间标准化\n"
        report += "2. 建议定期验证PPG标准曲线的线性度\n"
        report += "3. 对于新色谱条件，使用PPG指数系统进行迁移\n"
        report += "4. 结合机器学习模型进一步提高预测精度\n"

        report += "\n" + "=" * 70 + "\n"
        report += "报告生成完成\n"
        report += "=" * 70 + "\n"

        # 保存到文件
        if output_dir:
            try:
                Path(output_dir).mkdir(parents=True, exist_ok=True)
                report_path = os.path.join(output_dir,
                                           f"PPG_analysis_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt")

                with open(report_path, 'w', encoding='utf-8') as f:
                    f.write(report)

                return True, "报告生成成功", report_path
            except Exception as e:
                return False, f"保存报告失败: {str(e)}", report

        return True, "报告生成成功", report

    def save_all_results(self, output_dir: str) -> Tuple[bool, str, List[str]]:
        """保存所有结果"""
        try:
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)

            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            saved_files = []

            # 1. 保存PPG数据
            if self.ppg_data:
                for condition, df in self.ppg_data.items():
                    file_name = f"PPG_{condition}_{timestamp}.csv"
                    file_path = output_path / file_name
                    df.to_csv(file_path, index=False, encoding='utf-8')
                    saved_files.append(str(file_path))

            # 2. 保存化合物数据
            for data_type, data_dict in [('SMRT', self.smrt_data), ('Validation', self.validation_data),
                                         ('Calibration', self.calibration_data)]:
                if data_dict:
                    for key, df in data_dict.items():
                        file_name = f"{data_type}_{key}_{timestamp}.csv"
                        file_path = output_path / file_name
                        df.to_csv(file_path, index=False, encoding='utf-8')
                        saved_files.append(str(file_path))

            # 3. 保存分析结果
            if self.results:
                results_df = pd.DataFrame(self.results).T
                file_name = f"Analysis_Results_{timestamp}.csv"
                file_path = output_path / file_name
                results_df.to_csv(file_path, encoding='utf-8')
                saved_files.append(str(file_path))

            # 4. 生成并保存报告
            success, msg, report_content = self.generate_report(output_dir)
            if success and isinstance(report_content, str) and os.path.exists(report_content):
                saved_files.append(report_content)

            # 5. 保存可视化图表
            if self.visualizations:
                for name, fig in self.visualizations.items():
                    file_name = f"Visualization_{name}_{timestamp}.png"
                    file_path = output_path / file_name
                    fig.savefig(file_path, dpi=300, bbox_inches='tight')
                    saved_files.append(str(file_path))

            return True, f"结果已保存到 {output_dir}", saved_files

        except Exception as e:
            return False, f"保存结果失败: {str(e)}", []


class PPGExperimentGUI:
    """PPG试验方案分析GUI"""

    def __init__(self, root):
        """初始化GUI"""
        self.root = root
        self.root.title("PPG保留指数分析系统 - 试验方案实现")
        self.root.geometry("1400x900")

        # 设置图标
        try:
            self.root.iconbitmap(default='icon.ico')
        except:
            pass

        # 初始化分析器
        self.analyzer = PPGExperimentAnalyzer()

        # 处理线程相关
        self.processing_thread = None
        self.is_processing = False

        # 存储加载的条件信息
        self.loaded_conditions = set()
        self.loaded_datasets = set()

        # 创建UI
        self.setup_ui()

        # 绑定关闭事件
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)

    def setup_ui(self):
        """设置UI界面"""
        # 创建主框架
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # 创建Notebook（标签页）
        self.notebook = ttk.Notebook(main_frame)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        # 创建各个标签页
        self.setup_data_tab()
        self.setup_analysis_tab()
        self.setup_visualization_tab()
        self.setup_report_tab()
        self.setup_log_tab()

        # 状态栏
        self.status_var = tk.StringVar(value="就绪")
        status_bar = ttk.Label(main_frame, textvariable=self.status_var,
                               relief=tk.SUNKEN, anchor=tk.W)
        status_bar.pack(fill=tk.X, pady=(5, 0))

    def setup_data_tab(self):
        """设置数据加载标签页"""
        data_tab = ttk.Frame(self.notebook)
        self.notebook.add(data_tab, text="数据加载")

        # 主框架
        data_frame = ttk.LabelFrame(data_tab, text="数据加载与管理", padding=15)
        data_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # ==================== PPG数据加载 ====================
        ppg_frame = ttk.LabelFrame(data_frame, text="PPG标准品数据加载", padding=10)
        ppg_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(ppg_frame, text="条件名称:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.ppg_condition_var = tk.StringVar(value="C18_gradient1")
        ppg_condition_entry = ttk.Entry(ppg_frame, textvariable=self.ppg_condition_var, width=20)
        ppg_condition_entry.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Label(ppg_frame, text="PPG数据文件:").grid(row=1, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.ppg_file_var = tk.StringVar()
        ppg_file_entry = ttk.Entry(ppg_frame, textvariable=self.ppg_file_var, width=60)
        ppg_file_entry.grid(row=1, column=1, sticky=(tk.W, tk.E), padx=(0, 10), pady=5)

        ttk.Button(ppg_frame, text="浏览...", command=self.browse_ppg_file).grid(row=1, column=2, pady=5)

        ttk.Button(ppg_frame, text="加载PPG数据", command=self.load_ppg_data).grid(row=2, column=0, columnspan=3,
                                                                                   pady=10)

        # ==================== SMRT数据加载 ====================
        smrt_frame = ttk.LabelFrame(data_frame, text="SMRT数据集加载", padding=10)
        smrt_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(smrt_frame, text="条件名称:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.smrt_condition_var = tk.StringVar(value="default")
        smrt_condition_entry = ttk.Entry(smrt_frame, textvariable=self.smrt_condition_var, width=20)
        smrt_condition_entry.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Label(smrt_frame, text="SMRT数据文件:").grid(row=1, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.smrt_file_var = tk.StringVar()
        smrt_file_entry = ttk.Entry(smrt_frame, textvariable=self.smrt_file_var, width=60)
        smrt_file_entry.grid(row=1, column=1, sticky=(tk.W, tk.E), padx=(0, 10), pady=5)

        ttk.Button(smrt_frame, text="浏览...", command=self.browse_smrt_file).grid(row=1, column=2, pady=5)

        ttk.Button(smrt_frame, text="加载SMRT数据", command=self.load_smrt_data).grid(row=2, column=0, columnspan=3,
                                                                                      pady=10)

        # ==================== 验证集数据加载 ====================
        validation_frame = ttk.LabelFrame(data_frame, text="验证集数据加载", padding=10)
        validation_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(validation_frame, text="条件名称:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.validation_condition_var = tk.StringVar(value="default")
        validation_condition_entry = ttk.Entry(validation_frame, textvariable=self.validation_condition_var, width=20)
        validation_condition_entry.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Label(validation_frame, text="验证集文件:").grid(row=1, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.validation_file_var = tk.StringVar()
        validation_file_entry = ttk.Entry(validation_frame, textvariable=self.validation_file_var, width=60)
        validation_file_entry.grid(row=1, column=1, sticky=(tk.W, tk.E), padx=(0, 10), pady=5)

        ttk.Button(validation_frame, text="浏览...", command=self.browse_validation_file).grid(row=1, column=2, pady=5)

        ttk.Button(validation_frame, text="加载验证集数据", command=self.load_validation_data).grid(row=2, column=0,
                                                                                                    columnspan=3,
                                                                                                    pady=10)

        # ==================== 校正化合物数据加载 ====================
        calibration_frame = ttk.LabelFrame(data_frame, text="校正化合物数据加载", padding=10)
        calibration_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(calibration_frame, text="条件名称:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.calibration_condition_var = tk.StringVar(value="default")
        calibration_condition_entry = ttk.Entry(calibration_frame, textvariable=self.calibration_condition_var,
                                                width=20)
        calibration_condition_entry.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Label(calibration_frame, text="校正化合物文件:").grid(row=1, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.calibration_file_var = tk.StringVar()
        calibration_file_entry = ttk.Entry(calibration_frame, textvariable=self.calibration_file_var, width=60)
        calibration_file_entry.grid(row=1, column=1, sticky=(tk.W, tk.E), padx=(0, 10), pady=5)

        ttk.Button(calibration_frame, text="浏览...", command=self.browse_calibration_file).grid(row=1, column=2,
                                                                                                 pady=5)

        ttk.Button(calibration_frame, text="加载校正化合物数据", command=self.load_calibration_data).grid(row=2,
                                                                                                          column=0,
                                                                                                          columnspan=3,
                                                                                                          pady=10)

        # ==================== 加载的数据概览 ====================
        overview_frame = ttk.LabelFrame(data_frame, text="已加载数据概览", padding=10)
        overview_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 10))

        # 创建Treeview显示加载的数据
        columns = ("数据类型", "条件", "数据点数", "状态")
        self.data_tree = ttk.Treeview(overview_frame, columns=columns, show="headings", height=8)

        # 设置列标题
        for col in columns:
            self.data_tree.heading(col, text=col)
            self.data_tree.column(col, width=120)

        # 添加滚动条
        scrollbar = ttk.Scrollbar(overview_frame, orient=tk.VERTICAL, command=self.data_tree.yview)
        self.data_tree.configure(yscrollcommand=scrollbar.set)

        self.data_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # 数据管理按钮
        button_frame = ttk.Frame(data_frame)
        button_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Button(button_frame, text="清空所有数据", command=self.clear_all_data).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="生成示例数据", command=self.generate_example_data).pack(side=tk.LEFT, padx=5)

        # 数据状态
        self.data_status_var = tk.StringVar(value="等待加载数据...")
        ttk.Label(data_frame, textvariable=self.data_status_var).pack(anchor=tk.W)

    def setup_analysis_tab(self):
        """设置分析标签页"""
        analysis_tab = ttk.Frame(self.notebook)
        self.notebook.add(analysis_tab, text="数据分析")

        # 主框架
        analysis_frame = ttk.LabelFrame(analysis_tab, text="PPG试验方案分析", padding=15)
        analysis_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # ==================== PPG标准曲线拟合 ====================
        curve_frame = ttk.LabelFrame(analysis_frame, text="PPG标准曲线拟合", padding=10)
        curve_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(curve_frame, text="选择条件:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.curve_condition_var = tk.StringVar()
        self.curve_condition_combo = ttk.Combobox(curve_frame, textvariable=self.curve_condition_var,
                                                  width=25, state="readonly")
        self.curve_condition_combo.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Label(curve_frame, text="模型类型:").grid(row=0, column=2, sticky=tk.W, padx=(0, 5), pady=5)
        self.curve_model_var = tk.StringVar(value="logarithmic")
        curve_model_combo = ttk.Combobox(curve_frame, textvariable=self.curve_model_var,
                                         values=["logarithmic", "linear"],
                                         width=15, state="readonly")
        curve_model_combo.grid(row=0, column=3, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Button(curve_frame, text="拟合PPG标准曲线", command=self.fit_ppg_curve).grid(row=0, column=4, padx=(20, 0),
                                                                                         pady=5)

        # ==================== PPG指数计算 ====================
        calc_frame = ttk.LabelFrame(analysis_frame, text="PPG指数计算与转换", padding=10)
        calc_frame.pack(fill=tk.X, pady=(0, 15))

        # 计算PPG指数
        ttk.Label(calc_frame, text="保留时间(min):").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.calc_rt_var = tk.StringVar(value="5.5")
        calc_rt_entry = ttk.Entry(calc_frame, textvariable=self.calc_rt_var, width=15)
        calc_rt_entry.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Label(calc_frame, text="源条件:").grid(row=0, column=2, sticky=tk.W, padx=(0, 5), pady=5)
        self.calc_source_var = tk.StringVar()
        self.calc_source_combo = ttk.Combobox(calc_frame, textvariable=self.calc_source_var,
                                              width=20, state="readonly")
        self.calc_source_combo.grid(row=0, column=3, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Button(calc_frame, text="计算PPG指数", command=self.calculate_ppg_index).grid(row=0, column=4, padx=(10, 0),
                                                                                          pady=5)

        # PPG指数显示
        self.ppg_index_var = tk.StringVar(value="")
        ttk.Label(calc_frame, text="PPG指数:").grid(row=1, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        ttk.Label(calc_frame, textvariable=self.ppg_index_var, font=("Arial", 10, "bold")).grid(row=1, column=1,
                                                                                                sticky=tk.W,
                                                                                                padx=(0, 10), pady=5)

        # 从PPG指数计算保留时间
        ttk.Label(calc_frame, text="PPG指数:").grid(row=2, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.rt_from_index_var = tk.StringVar(value="550")
        rt_from_index_entry = ttk.Entry(calc_frame, textvariable=self.rt_from_index_var, width=15)
        rt_from_index_entry.grid(row=2, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Label(calc_frame, text="目标条件:").grid(row=2, column=2, sticky=tk.W, padx=(0, 5), pady=5)
        self.calc_target_var = tk.StringVar()
        self.calc_target_combo = ttk.Combobox(calc_frame, textvariable=self.calc_target_var,
                                              width=20, state="readonly")
        self.calc_target_combo.grid(row=2, column=3, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Label(calc_frame, text="计算方法:").grid(row=2, column=4, sticky=tk.W, padx=(10, 5), pady=5)
        self.calc_method_var = tk.StringVar(value="regression")
        calc_method_combo = ttk.Combobox(calc_frame, textvariable=self.calc_method_var,
                                         values=["interpolation", "regression"],
                                         width=12, state="readonly")
        calc_method_combo.grid(row=2, column=5, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Button(calc_frame, text="计算保留时间", command=self.calculate_rt_from_index).grid(row=2, column=6,
                                                                                               padx=(10, 0), pady=5)

        # 计算结果显示
        self.calc_rt_result_var = tk.StringVar(value="")
        ttk.Label(calc_frame, text="保留时间:").grid(row=3, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        ttk.Label(calc_frame, textvariable=self.calc_rt_result_var, font=("Arial", 10, "bold")).grid(row=3, column=1,
                                                                                                     sticky=tk.W,
                                                                                                     padx=(0, 10),
                                                                                                     pady=5)

        # ==================== 校正方法分析 ====================
        calibration_frame = ttk.LabelFrame(analysis_frame, text="校正方法分析", padding=10)
        calibration_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(calibration_frame, text="源条件:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.calibration_source_var = tk.StringVar()
        self.calibration_source_combo = ttk.Combobox(calibration_frame, textvariable=self.calibration_source_var,
                                                     width=20, state="readonly")
        self.calibration_source_combo.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Label(calibration_frame, text="目标条件:").grid(row=0, column=2, sticky=tk.W, padx=(0, 5), pady=5)
        self.calibration_target_var = tk.StringVar()
        self.calibration_target_combo = ttk.Combobox(calibration_frame, textvariable=self.calibration_target_var,
                                                     width=20, state="readonly")
        self.calibration_target_combo.grid(row=0, column=3, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Button(calibration_frame, text="应用五种校正方式", command=self.apply_calibration_methods).grid(row=0,
                                                                                                            column=4,
                                                                                                            padx=(
                                                                                                            20, 0),
                                                                                                            pady=5)

        # ==================== 模型性能评估 ====================
        model_frame = ttk.LabelFrame(analysis_frame, text="模型性能评估", padding=10)
        model_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(model_frame, text="模型类型:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.model_type_var = tk.StringVar(value="linear")
        model_type_combo = ttk.Combobox(model_frame, textvariable=self.model_type_var,
                                        values=["linear", "ridge", "lasso"],
                                        width=15, state="readonly")
        model_type_combo.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Button(model_frame, text="评估模型性能", command=self.evaluate_model_performance).grid(row=0, column=2,
                                                                                                   padx=(20, 0), pady=5)

        # ==================== 分析结果显示 ====================
        results_frame = ttk.LabelFrame(analysis_frame, text="分析结果", padding=10)
        results_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 10))

        # 创建文本显示区域
        self.analysis_text = scrolledtext.ScrolledText(results_frame, width=80, height=15,
                                                       wrap=tk.WORD, font=("Consolas", 10))
        self.analysis_text.pack(fill=tk.BOTH, expand=True)

        # 配置标签颜色
        self.analysis_text.tag_config("INFO", foreground="black")
        self.analysis_text.tag_config("SUCCESS", foreground="green")
        self.analysis_text.tag_config("WARNING", foreground="orange")
        self.analysis_text.tag_config("ERROR", foreground="red")

        # 分析按钮
        button_frame = ttk.Frame(analysis_frame)
        button_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Button(button_frame, text="清空结果显示", command=self.clear_analysis_text).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="导出分析结果", command=self.export_analysis_results).pack(side=tk.LEFT, padx=5)

        # 分析状态
        self.analysis_status_var = tk.StringVar(value="等待分析...")
        ttk.Label(analysis_frame, textvariable=self.analysis_status_var).pack(anchor=tk.W)

    def setup_visualization_tab(self):
        """设置可视化标签页"""
        viz_tab = ttk.Frame(self.notebook)
        self.notebook.add(viz_tab, text="数据可视化")

        # 主框架
        viz_frame = ttk.LabelFrame(viz_tab, text="数据可视化", padding=15)
        viz_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # ==================== 可视化选项 ====================
        options_frame = ttk.LabelFrame(viz_frame, text="可视化选项", padding=10)
        options_frame.pack(fill=tk.X, pady=(0, 15))

        # 可视化类型选择
        ttk.Label(options_frame, text="可视化类型:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.viz_type_var = tk.StringVar(value="ppg_curves")
        viz_type_combo = ttk.Combobox(options_frame, textvariable=self.viz_type_var,
                                      values=["ppg_curves", "calibration_comparison",
                                              "index_distribution", "error_analysis",
                                              "prediction_vs_actual", "all"],
                                      width=20, state="readonly")
        viz_type_combo.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        # 生成图表按钮
        ttk.Button(options_frame, text="生成图表", command=self.generate_visualization).grid(row=0, column=2,
                                                                                             padx=(20, 0), pady=5)

        # 保存图表按钮
        ttk.Button(options_frame, text="保存图表", command=self.save_visualization).grid(row=0, column=3, padx=(10, 0),
                                                                                         pady=5)

        # ==================== 图表显示区域 ====================
        display_frame = ttk.LabelFrame(viz_frame, text="图表显示", padding=10)
        display_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 10))

        # 创建画布框架
        self.figure_canvas = None
        self.figure_toolbar = None
        self.current_figure = None
        self.current_viz_type = None

        # 创建占位标签
        self.viz_placeholder = ttk.Label(display_frame, text="图表将在此处显示",
                                         font=("Arial", 14), foreground="gray")
        self.viz_placeholder.pack(expand=True)

        # 可视化控制按钮
        button_frame = ttk.Frame(viz_frame)
        button_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Button(button_frame, text="清除图表", command=self.clear_visualization).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="保存所有图表", command=self.save_all_visualizations).pack(side=tk.LEFT, padx=5)

        # 可视化状态
        self.viz_status_var = tk.StringVar(value="等待生成图表...")
        ttk.Label(viz_frame, textvariable=self.viz_status_var).pack(anchor=tk.W)

    def setup_report_tab(self):
        """设置报告标签页"""
        report_tab = ttk.Frame(self.notebook)
        self.notebook.add(report_tab, text="报告输出")

        # 主框架
        report_frame = ttk.LabelFrame(report_tab, text="报告输出与管理", padding=15)
        report_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # ==================== 输出设置 ====================
        output_frame = ttk.LabelFrame(report_frame, text="输出设置", padding=10)
        output_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(output_frame, text="输出目录:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.output_dir_var = tk.StringVar(value=os.path.join(os.getcwd(), "PPG_Results"))
        output_dir_entry = ttk.Entry(output_frame, textvariable=self.output_dir_var, width=60)
        output_dir_entry.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=(0, 10), pady=5)

        ttk.Button(output_frame, text="浏览...", command=self.browse_output_dir).grid(row=0, column=2, pady=5)

        # ==================== 报告预览 ====================
        preview_frame = ttk.LabelFrame(report_frame, text="报告预览", padding=10)
        preview_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 15))

        # 创建文本显示区域
        self.report_text = scrolledtext.ScrolledText(preview_frame, width=80, height=15,
                                                     wrap=tk.WORD, font=("Consolas", 10))
        self.report_text.pack(fill=tk.BOTH, expand=True)

        # ==================== 输出控制 ====================
        control_frame = ttk.Frame(report_frame)
        control_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Button(control_frame, text="生成报告", command=self.generate_report).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_frame, text="保存所有结果", command=self.save_all_results).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_frame, text="打开输出目录", command=self.open_output_dir).pack(side=tk.LEFT, padx=5)

        # 输出状态
        self.output_status_var = tk.StringVar(value="等待输出结果...")
        ttk.Label(report_frame, textvariable=self.output_status_var).pack(anchor=tk.W)

    def setup_log_tab(self):
        """设置日志标签页"""
        log_tab = ttk.Frame(self.notebook)
        self.notebook.add(log_tab, text="运行日志")

        # 主框架
        log_frame = ttk.LabelFrame(log_tab, text="程序运行日志", padding=15)
        log_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # 日志显示区域
        self.log_text = scrolledtext.ScrolledText(log_frame, width=100, height=25,
                                                  wrap=tk.WORD, font=("Consolas", 10))
        self.log_text.pack(fill=tk.BOTH, expand=True)

        # 配置标签颜色
        self.log_text.tag_config("INFO", foreground="black")
        self.log_text.tag_config("SUCCESS", foreground="green")
        self.log_text.tag_config("WARNING", foreground="orange")
        self.log_text.tag_config("ERROR", foreground="red")

        # 日志控制按钮
        button_frame = ttk.Frame(log_frame)
        button_frame.pack(fill=tk.X, pady=(10, 0))

        ttk.Button(button_frame, text="清空日志", command=self.clear_log).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="保存日志", command=self.save_log).pack(side=tk.LEFT, padx=5)

    def browse_ppg_file(self):
        """浏览PPG数据文件"""
        file_types = [("数据文件", "*.csv *.xlsx *.xls"), ("CSV文件", "*.csv"),
                      ("Excel文件", "*.xlsx *.xls"), ("所有文件", "*.*")]
        file_path = filedialog.askopenfilename(title="选择PPG数据文件", filetypes=file_types)

        if file_path:
            self.ppg_file_var.set(file_path)

    def browse_smrt_file(self):
        """浏览SMRT数据文件"""
        file_types = [("数据文件", "*.csv *.xlsx *.xls"), ("CSV文件", "*.csv"),
                      ("Excel文件", "*.xlsx *.xls"), ("所有文件", "*.*")]
        file_path = filedialog.askopenfilename(title="选择SMRT数据文件", filetypes=file_types)

        if file_path:
            self.smrt_file_var.set(file_path)

    def browse_validation_file(self):
        """浏览验证集文件"""
        file_types = [("数据文件", "*.csv *.xlsx *.xls"), ("CSV文件", "*.csv"),
                      ("Excel文件", "*.xlsx *.xls"), ("所有文件", "*.*")]
        file_path = filedialog.askopenfilename(title="选择验证集文件", filetypes=file_types)

        if file_path:
            self.validation_file_var.set(file_path)

    def browse_calibration_file(self):
        """浏览校正化合物文件"""
        file_types = [("数据文件", "*.csv *.xlsx *.xls"), ("CSV文件", "*.csv"),
                      ("Excel文件", "*.xlsx *.xls"), ("所有文件", "*.*")]
        file_path = filedialog.askopenfilename(title="选择校正化合物文件", filetypes=file_types)

        if file_path:
            self.calibration_file_var.set(file_path)

    def browse_output_dir(self):
        """浏览输出目录"""
        dir_path = filedialog.askdirectory(title="选择输出目录")

        if dir_path:
            self.output_dir_var.set(dir_path)

    def load_ppg_data(self):
        """加载PPG数据"""
        ppg_file = self.ppg_file_var.get().strip()
        condition = self.ppg_condition_var.get().strip()

        if not ppg_file:
            messagebox.showwarning("警告", "请选择PPG数据文件")
            return

        if not condition:
            messagebox.showwarning("警告", "请输入条件名称")
            return

        # 显示加载状态
        self.log_message(f"正在加载PPG数据: {ppg_file} (条件: {condition})", "INFO")
        self.update_status(f"加载PPG数据: {Path(ppg_file).name}")

        # 在线程中加载数据
        self.processing_thread = threading.Thread(
            target=self._load_ppg_data_thread,
            args=(ppg_file, condition)
        )
        self.processing_thread.start()

    def _load_ppg_data_thread(self, ppg_file, condition):
        """加载PPG数据的线程函数"""
        try:
            success, msg = self.analyzer.load_ppg_data(ppg_file, condition)

            if success:
                self.root.after(0, lambda: self.log_message(f"✓ {msg}", "SUCCESS"))
                self.root.after(0, self.update_data_tree)
                self.root.after(0, self.update_condition_comboboxes)
                self.loaded_conditions.add(condition)
            else:
                self.root.after(0, lambda: self.log_message(f"✗ {msg}", "ERROR"))

            self.root.after(0, lambda: self.update_status("就绪"))

        except Exception as e:
            self.root.after(0, lambda: self.log_message(f"✗ 加载PPG数据失败: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("加载失败"))

    def load_smrt_data(self):
        """加载SMRT数据"""
        smrt_file = self.smrt_file_var.get().strip()
        condition = self.smrt_condition_var.get().strip()

        if not smrt_file:
            messagebox.showwarning("警告", "请选择SMRT数据文件")
            return

        if not condition:
            messagebox.showwarning("警告", "请输入条件名称")
            return

        # 显示加载状态
        self.log_message(f"正在加载SMRT数据: {smrt_file} (条件: {condition})", "INFO")
        self.update_status(f"加载SMRT数据: {Path(smrt_file).name}")

        # 在线程中加载数据
        self.processing_thread = threading.Thread(
            target=self._load_smrt_data_thread,
            args=(smrt_file, condition)
        )
        self.processing_thread.start()

    def _load_smrt_data_thread(self, smrt_file, condition):
        """加载SMRT数据的线程函数"""
        try:
            success, msg = self.analyzer.load_compound_data(smrt_file, "smrt", condition)

            if success:
                self.root.after(0, lambda: self.log_message(f"✓ {msg}", "SUCCESS"))
                self.root.after(0, self.update_data_tree)
                self.root.after(0, self.update_condition_comboboxes)
                self.loaded_datasets.add(f"smrt_{condition}")
            else:
                self.root.after(0, lambda: self.log_message(f"✗ {msg}", "ERROR"))

            self.root.after(0, lambda: self.update_status("就绪"))

        except Exception as e:
            self.root.after(0, lambda: self.log_message(f"✗ 加载SMRT数据失败: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("加载失败"))

    def load_validation_data(self):
        """加载验证集数据"""
        validation_file = self.validation_file_var.get().strip()
        condition = self.validation_condition_var.get().strip()

        if not validation_file:
            messagebox.showwarning("警告", "请选择验证集文件")
            return

        if not condition:
            messagebox.showwarning("警告", "请输入条件名称")
            return

        # 显示加载状态
        self.log_message(f"正在加载验证集数据: {validation_file} (条件: {condition})", "INFO")
        self.update_status(f"加载验证集数据: {Path(validation_file).name}")

        # 在线程中加载数据
        self.processing_thread = threading.Thread(
            target=self._load_validation_data_thread,
            args=(validation_file, condition)
        )
        self.processing_thread.start()

    def _load_validation_data_thread(self, validation_file, condition):
        """加载验证集数据的线程函数"""
        try:
            success, msg = self.analyzer.load_compound_data(validation_file, "validation", condition)

            if success:
                self.root.after(0, lambda: self.log_message(f"✓ {msg}", "SUCCESS"))
                self.root.after(0, self.update_data_tree)
                self.root.after(0, self.update_condition_comboboxes)
                self.loaded_datasets.add(f"validation_{condition}")
            else:
                self.root.after(0, lambda: self.log_message(f"✗ {msg}", "ERROR"))

            self.root.after(0, lambda: self.update_status("就绪"))

        except Exception as e:
            self.root.after(0, lambda: self.log_message(f"✗ 加载验证集数据失败: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("加载失败"))

    def load_calibration_data(self):
        """加载校正化合物数据"""
        calibration_file = self.calibration_file_var.get().strip()
        condition = self.calibration_condition_var.get().strip()

        if not calibration_file:
            messagebox.showwarning("警告", "请选择校正化合物文件")
            return

        if not condition:
            messagebox.showwarning("警告", "请输入条件名称")
            return

        # 显示加载状态
        self.log_message(f"正在加载校正化合物数据: {calibration_file} (条件: {condition})", "INFO")
        self.update_status(f"加载校正化合物数据: {Path(calibration_file).name}")

        # 在线程中加载数据
        self.processing_thread = threading.Thread(
            target=self._load_calibration_data_thread,
            args=(calibration_file, condition)
        )
        self.processing_thread.start()

    def _load_calibration_data_thread(self, calibration_file, condition):
        """加载校正化合物数据的线程函数"""
        try:
            success, msg = self.analyzer.load_compound_data(calibration_file, "calibration", condition)

            if success:
                self.root.after(0, lambda: self.log_message(f"✓ {msg}", "SUCCESS"))
                self.root.after(0, self.update_data_tree)
                self.root.after(0, self.update_condition_comboboxes)
                self.loaded_datasets.add(f"calibration_{condition}")
            else:
                self.root.after(0, lambda: self.log_message(f"✗ {msg}", "ERROR"))

            self.root.after(0, lambda: self.update_status("就绪"))

        except Exception as e:
            self.root.after(0, lambda: self.log_message(f"✗ 加载校正化合物数据失败: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("加载失败"))

    def update_data_tree(self):
        """更新数据树显示"""
        # 清空现有数据
        for item in self.data_tree.get_children():
            self.data_tree.delete(item)

        # 添加PPG数据
        for condition, df in self.analyzer.ppg_data.items():
            self.data_tree.insert("", tk.END, values=("PPG标准品", condition, len(df), "已加载"))

        # 添加SMRT数据
        for key, df in self.analyzer.smrt_data.items():
            condition = key.replace('smrt_', '')
            self.data_tree.insert("", tk.END, values=("SMRT数据集", condition, len(df), "已加载"))

        # 添加验证数据
        for key, df in self.analyzer.validation_data.items():
            condition = key.replace('validation_', '')
            self.data_tree.insert("", tk.END, values=("验证集", condition, len(df), "已加载"))

        # 添加校正化合物数据
        for key, df in self.analyzer.calibration_data.items():
            condition = key.replace('calibration_', '')
            self.data_tree.insert("", tk.END, values=("校正化合物", condition, len(df), "已加载"))

        # 更新状态
        total_data = (len(self.analyzer.ppg_data) + len(self.analyzer.smrt_data) +
                      len(self.analyzer.validation_data) + len(self.analyzer.calibration_data))
        self.data_status_var.set(f"已加载 {total_data} 个数据集")

    def update_condition_comboboxes(self):
        """更新条件选择下拉框"""
        conditions = list(self.analyzer.ppg_data.keys())

        # 更新曲线拟合条件
        self.curve_condition_combo['values'] = conditions
        if conditions and not self.curve_condition_var.get():
            self.curve_condition_var.set(conditions[0])

        # 更新计算源条件
        self.calc_source_combo['values'] = conditions
        if conditions and not self.calc_source_var.get():
            self.calc_source_var.set(conditions[0])

        # 更新计算目标条件
        self.calc_target_combo['values'] = conditions
        if conditions and not self.calc_target_var.get():
            self.calc_target_var.set(
                conditions[0] if len(conditions) == 1 else conditions[1] if len(conditions) > 1 else "")

        # 更新校正源条件
        self.calibration_source_combo['values'] = conditions
        if conditions and not self.calibration_source_var.get():
            self.calibration_source_var.set(conditions[0])

        # 更新校正目标条件
        self.calibration_target_combo['values'] = conditions
        if conditions and not self.calibration_target_var.get():
            self.calibration_target_var.set(
                conditions[0] if len(conditions) == 1 else conditions[1] if len(conditions) > 1 else "")

    def clear_all_data(self):
        """清空所有数据"""
        if messagebox.askyesno("确认", "确定要清空所有已加载的数据吗？"):
            self.analyzer = PPGExperimentAnalyzer()
            self.loaded_conditions.clear()
            self.loaded_datasets.clear()
            self.update_data_tree()
            self.update_condition_comboboxes()
            self.log_message("已清空所有数据", "INFO")

    def generate_example_data(self):
        """生成示例数据"""
        if messagebox.askyesno("确认", "将生成示例数据用于演示，这会覆盖当前数据吗？"):
            self.clear_all_data()
            self.log_message("生成示例数据...", "INFO")

            # 生成示例PPG数据
            for i, condition in enumerate(['C18_gradient1', 'C18_gradient2', 'C18_gradient3']):
                n_values = np.arange(2, 31)  # PPG2-30
                if i == 0:
                    rt_values = 2.0 + 0.5 * np.log(n_values) + np.random.normal(0, 0.1, len(n_values))
                elif i == 1:
                    rt_values = 2.5 + 0.6 * np.log(n_values) + np.random.normal(0, 0.1, len(n_values))
                else:
                    rt_values = 3.0 + 0.7 * np.log(n_values) + np.random.normal(0, 0.1, len(n_values))

                df_ppg = pd.DataFrame({'聚合度': n_values, '保留时间': rt_values})
                self.analyzer.ppg_data[condition] = df_ppg

            # 生成示例SMRT数据
            n_smrt = 28
            compound_names = [f'Compound_{i + 1}' for i in range(n_smrt)]
            rt_source = np.random.uniform(2, 15, n_smrt)
            rt_target = rt_source * 1.1 + np.random.normal(0, 0.2, n_smrt)
            logp_values = np.random.uniform(-2, 5, n_smrt)

            self.analyzer.smrt_data['smrt_default'] = pd.DataFrame({
                '化合物名称': compound_names,
                '保留时间': rt_source,
                'RT_target': rt_target,
                'logP': logp_values,
                '分子量': np.random.uniform(100, 500, n_smrt)
            })

            # 生成示例验证数据
            n_validation = 50
            validation_names = [f'Val_Compound_{i + 1}' for i in range(n_validation)]
            rt_val_source = np.random.uniform(2, 18, n_validation)
            rt_val_target = rt_val_source * 1.05 + np.random.normal(0, 0.15, n_validation)
            logp_val = np.random.uniform(-2, 5, n_validation)

            self.analyzer.validation_data['validation_default'] = pd.DataFrame({
                '化合物名称': validation_names,
                '保留时间': rt_val_source,
                'RT_target': rt_val_target,
                'logP': logp_val,
                '分子量': np.random.uniform(100, 500, n_validation)
            })

            # 生成示例校正化合物数据
            n_calib = 28
            calib_names = [f'Calib_Compound_{i + 1}' for i in range(n_calib)]
            rt_calib_source = np.random.uniform(2, 15, n_calib)
            rt_calib_target = rt_calib_source * 1.08 + np.random.normal(0, 0.1, n_calib)

            self.analyzer.calibration_data['calibration_default'] = pd.DataFrame({
                '化合物名称': calib_names,
                'RT_source': rt_calib_source,
                'RT_target': rt_calib_target
            })

            self.update_data_tree()
            self.update_condition_comboboxes()
            self.log_message("示例数据生成完成", "SUCCESS")

    def fit_ppg_curve(self):
        """拟合PPG标准曲线"""
        condition = self.curve_condition_var.get()

        if not condition:
            messagebox.showwarning("警告", "请选择条件")
            return

        # 显示分析状态
        self.analysis_message(f"正在拟合PPG标准曲线 (条件: {condition})", "INFO")
        self.update_status(f"拟合PPG标准曲线: {condition}")

        # 在线程中拟合曲线
        self.processing_thread = threading.Thread(
            target=self._fit_ppg_curve_thread,
            args=(condition,)
        )
        self.processing_thread.start()

    def _fit_ppg_curve_thread(self, condition):
        """拟合PPG标准曲线的线程函数"""
        try:
            success, msg, curve_data = self.analyzer.fit_ppg_standard_curve(condition)

            if success:
                self.root.after(0, lambda: self.analysis_message(f"✓ {msg}", "SUCCESS"))

                # 显示拟合结果
                result_text = f"PPG标准曲线拟合结果 - {condition}:\n"
                result_text += f"  模型类型: {curve_data.get('model_type', '未知')}\n"
                result_text += f"  R²: {curve_data.get('r_squared', 0):.6f}\n"

                if 'params' in curve_data:
                    params = curve_data['params']
                    if curve_data.get('model_type') == 'logarithmic':
                        result_text += f"  截距: {params[0]:.6f}\n"
                        result_text += f"  斜率: {params[1]:.6f}\n"
                        result_text += f"  方程: RT = {params[0]:.4f} + {params[1]:.4f} * ln(n)\n"
                    else:
                        result_text += f"  截距: {params[0]:.6f}\n"
                        result_text += f"  斜率: {params[1]:.6f}\n"
                        result_text += f"  方程: RT = {params[0]:.4f} + {params[1]:.4f} * n\n"

                self.root.after(0, lambda: self.analysis_message(result_text, "INFO"))
            else:
                self.root.after(0, lambda: self.analysis_message(f"✗ {msg}", "ERROR"))

            self.root.after(0, lambda: self.update_status("就绪"))

        except Exception as e:
            self.root.after(0, lambda: self.analysis_message(f"✗ 拟合PPG标准曲线失败: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("拟合失败"))

    def calculate_ppg_index(self):
        """计算PPG指数"""
        try:
            rt = float(self.calc_rt_var.get().strip())
            condition = self.calc_source_var.get()

            if not condition:
                messagebox.showwarning("警告", "请选择源条件")
                return

            # 计算PPG指数
            ppg_index = self.analyzer.calculate_ppg_index(rt, condition)

            # 显示结果
            self.ppg_index_var.set(f"{ppg_index:.1f}")
            self.analysis_message(f"保留时间 {rt} min → PPG指数: {ppg_index:.1f}", "SUCCESS")

        except ValueError:
            messagebox.showerror("错误", "请输入有效的保留时间数值")
        except Exception as e:
            messagebox.showerror("错误", f"计算PPG指数失败: {str(e)}")

    def calculate_rt_from_index(self):
        """从PPG指数计算保留时间"""
        try:
            index = float(self.rt_from_index_var.get().strip())
            condition = self.calc_target_var.get()
            method = self.calc_method_var.get()

            if not condition:
                messagebox.showwarning("警告", "请选择目标条件")
                return

            # 计算保留时间
            rt = self.analyzer.calculate_rt_from_index(index, condition, method)

            # 显示结果
            self.calc_rt_result_var.set(f"{rt:.2f} min")
            self.analysis_message(f"PPG指数 {index:.1f} → 保留时间: {rt:.2f} min (方法: {method})", "SUCCESS")

        except ValueError:
            messagebox.showerror("错误", "请输入有效的PPG指数数值")
        except Exception as e:
            messagebox.showerror("错误", f"计算保留时间失败: {str(e)}")

    def apply_calibration_methods(self):
        """应用校正方法"""
        source_condition = self.calibration_source_var.get()
        target_condition = self.calibration_target_var.get()

        if not source_condition:
            messagebox.showwarning("警告", "请选择源条件")
            return

        if not target_condition:
            messagebox.showwarning("警告", "请选择目标条件")
            return

        if source_condition == target_condition:
            if not messagebox.askyesno("确认", "源条件和目标条件相同，是否继续？"):
                return

        # 显示分析状态
        self.analysis_message(f"正在应用五种校正方式 ({source_condition} → {target_condition})...", "INFO")
        self.update_status(f"应用校正方法: {source_condition}→{target_condition}")

        # 在线程中应用校正
        self.processing_thread = threading.Thread(
            target=self._apply_calibration_methods_thread,
            args=(source_condition, target_condition)
        )
        self.processing_thread.start()

    def _apply_calibration_methods_thread(self, source_condition, target_condition):
        """应用校正方法的线程函数"""
        try:
            success, msg = self.analyzer.apply_calibration_methods(source_condition, target_condition)

            if success:
                self.root.after(0, lambda: self.analysis_message(f"✓ {msg}", "SUCCESS"))

                # 显示结果摘要
                result_text = "校正方法分析结果:\n"
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
                            result_text += f"  错误: {result['error']}\n"

                        result_text += "\n"

                self.root.after(0, lambda: self.analysis_message(result_text, "INFO"))
            else:
                self.root.after(0, lambda: self.analysis_message(f"✗ {msg}", "ERROR"))

            self.root.after(0, lambda: self.update_status("就绪"))

        except Exception as e:
            self.root.after(0, lambda: self.analysis_message(f"✗ 应用校正方法失败: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("校正失败"))

    def evaluate_model_performance(self):
        """评估模型性能"""
        model_type = self.model_type_var.get()

        # 显示分析状态
        self.analysis_message(f"正在评估{model_type}模型性能...", "INFO")
        self.update_status(f"评估模型性能: {model_type}")

        # 在线程中评估模型
        self.processing_thread = threading.Thread(
            target=self._evaluate_model_performance_thread,
            args=(model_type,)
        )
        self.processing_thread.start()

    def _evaluate_model_performance_thread(self, model_type):
        """评估模型性能的线程函数"""
        try:
            success, msg = self.analyzer.evaluate_model_performance(model_type)

            if success:
                self.root.after(0, lambda: self.analysis_message(f"✓ {msg}", "SUCCESS"))

                # 显示模型性能
                if model_type in self.analyzer.models:
                    model_data = self.analyzer.models[model_type]
                    perf = model_data.get('performance', {})

                    result_text = f"{model_type}模型性能:\n"
                    result_text += f"  样本数: {model_data.get('n_samples', 0)}\n"
                    result_text += f"  MAE: {perf.get('MAE', 0):.4f}\n"
                    result_text += f"  RMSE: {perf.get('RMSE', 0):.4f}\n"
                    result_text += f"  R²: {perf.get('R2', 0):.4f}\n"

                    self.root.after(0, lambda: self.analysis_message(result_text, "INFO"))
            else:
                self.root.after(0, lambda: self.analysis_message(f"✗ {msg}", "ERROR"))

            self.root.after(0, lambda: self.update_status("就绪"))

        except Exception as e:
            self.root.after(0, lambda: self.analysis_message(f"✗ 评估模型性能失败: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("评估失败"))

    def generate_visualization(self):
        """生成可视化图表"""
        viz_type = self.viz_type_var.get()

        if viz_type == "all":
            # 生成所有图表
            viz_types = ["ppg_curves", "calibration_comparison", "index_distribution",
                         "error_analysis", "prediction_vs_actual"]
        else:
            viz_types = [viz_type]

        # 显示生成状态
        self.viz_status_var.set(f"正在生成 {viz_type} 图表...")
        self.update_status(f"生成图表: {viz_type}")

        # 在线程中生成图表
        self.processing_thread = threading.Thread(
            target=self._generate_visualization_thread,
            args=(viz_types,)
        )
        self.processing_thread.start()

    def _generate_visualization_thread(self, viz_types):
        """生成可视化图表的线程函数"""
        try:
            # 先生成所有图表
            success, msg, visualizations = self.analyzer.generate_visualizations()

            if success and visualizations:
                # 显示第一个图表
                if viz_types:
                    first_type = viz_types[0]
                    if first_type in visualizations:
                        self.root.after(0, lambda: self.display_figure(visualizations[first_type], first_type))
                    elif "ppg_curves" in visualizations:
                        self.root.after(0, lambda: self.display_figure(visualizations["ppg_curves"], "ppg_curves"))

                self.root.after(0, lambda: self.viz_status_var.set(f"生成 {len(visualizations)} 个图表"))
                self.root.after(0, lambda: self.log_message(f"✓ {msg}", "SUCCESS"))
            else:
                self.root.after(0, lambda: self.viz_status_var.set("生成图表失败"))
                self.root.after(0, lambda: self.log_message(f"✗ {msg}", "ERROR"))

            self.root.after(0, lambda: self.update_status("就绪"))

        except Exception as e:
            self.root.after(0, lambda: self.viz_status_var.set(f"生成图表失败: {str(e)}"))
            self.root.after(0, lambda: self.log_message(f"✗ 生成图表失败: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("生成失败"))

    def display_figure(self, fig, viz_type):
        """显示图表"""
        # 清除现有图表
        self.clear_visualization()

        # 移除占位标签
        self.viz_placeholder.pack_forget()

        # 创建画布
        canvas = FigureCanvasTkAgg(fig, master=self.viz_placeholder.master)
        canvas.draw()

        # 创建工具栏
        toolbar = NavigationToolbar2Tk(canvas, self.viz_placeholder.master)
        toolbar.update()

        # 显示画布和工具栏
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        toolbar.pack(fill=tk.X)

        # 保存引用
        self.figure_canvas = canvas
        self.figure_toolbar = toolbar
        self.current_figure = fig
        self.current_viz_type = viz_type

        # 更新状态
        self.viz_status_var.set(f"显示图表: {viz_type}")

    def save_visualization(self):
        """保存当前图表"""
        if self.current_figure is None:
            messagebox.showwarning("警告", "没有可保存的图表")
            return

        file_types = [("PNG文件", "*.png"), ("PDF文件", "*.pdf"),
                      ("SVG文件", "*.svg"), ("所有文件", "*.*")]

        default_name = f"PPG_{self.current_viz_type if self.current_viz_type else 'chart'}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

        file_path = filedialog.asksaveasfilename(
            title="保存图表",
            filetypes=file_types,
            defaultextension=".png",
            initialfile=default_name
        )

        if file_path:
            try:
                self.current_figure.savefig(file_path, dpi=300, bbox_inches='tight')
                self.log_message(f"✓ 图表已保存到: {file_path}", "SUCCESS")
                self.viz_status_var.set(f"图表已保存: {Path(file_path).name}")
            except Exception as e:
                self.log_message(f"✗ 保存图表失败: {str(e)}", "ERROR")

    def save_all_visualizations(self):
        """保存所有图表"""
        if not self.analyzer.visualizations:
            messagebox.showwarning("警告", "没有可保存的图表")
            return

        output_dir = filedialog.askdirectory(title="选择保存目录")

        if output_dir:
            try:
                # 生成所有图表
                success, msg, visualizations = self.analyzer.generate_visualizations(output_dir)

                if success:
                    self.log_message(f"✓ 所有图表已保存到: {output_dir}", "SUCCESS")
                    self.viz_status_var.set(f"保存 {len(visualizations)} 个图表")
                else:
                    self.log_message(f"✗ {msg}", "ERROR")

            except Exception as e:
                self.log_message(f"✗ 保存图表失败: {str(e)}", "ERROR")

    def clear_visualization(self):
        """清除图表"""
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

        # 显示占位标签
        self.viz_placeholder.pack(expand=True)

    def generate_report(self):
        """生成报告"""
        # 显示生成状态
        self.update_status("生成报告...")

        # 在线程中生成报告
        self.processing_thread = threading.Thread(
            target=self._generate_report_thread
        )
        self.processing_thread.start()

    def _generate_report_thread(self):
        """生成报告的线程函数"""
        try:
            # 生成报告
            success, msg, report_content = self.analyzer.generate_report()

            if success:
                # 在GUI中显示报告
                if isinstance(report_content, str) and os.path.exists(report_content):
                    # 从文件读取
                    with open(report_content, 'r', encoding='utf-8') as f:
                        report_text = f.read()
                else:
                    # 直接使用报告内容
                    report_text = report_content

                self.root.after(0, lambda: self.report_text.delete(1.0, tk.END))
                self.root.after(0, lambda: self.report_text.insert(tk.END, report_text))

                self.root.after(0, lambda: self.output_status_var.set("报告生成成功"))
                self.root.after(0, lambda: self.log_message(f"✓ {msg}", "SUCCESS"))
            else:
                self.root.after(0, lambda: self.output_status_var.set("报告生成失败"))
                self.root.after(0, lambda: self.log_message(f"✗ {msg}", "ERROR"))

            self.root.after(0, lambda: self.update_status("就绪"))

        except Exception as e:
            self.root.after(0, lambda: self.output_status_var.set("报告生成失败"))
            self.root.after(0, lambda: self.log_message(f"✗ 生成报告失败: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("生成失败"))

    def save_all_results(self):
        """保存所有结果"""
        output_dir = self.output_dir_var.get().strip()

        if not output_dir:
            messagebox.showwarning("警告", "请选择输出目录")
            return

        # 显示保存状态
        self.output_status_var.set("正在保存结果...")
        self.update_status("保存结果")

        # 在线程中保存
        self.processing_thread = threading.Thread(
            target=self._save_all_results_thread,
            args=(output_dir,)
        )
        self.processing_thread.start()

    def _save_all_results_thread(self, output_dir):
        """保存所有结果的线程函数"""
        try:
            success, msg, saved_files = self.analyzer.save_all_results(output_dir)

            if success:
                self.root.after(0, lambda: self.output_status_var.set(f"结果保存完成: {len(saved_files)} 个文件"))
                self.root.after(0, lambda: self.log_message(f"✓ {msg}", "SUCCESS"))

                # 在报告标签页显示保存的文件
                self.root.after(0, lambda: self.report_text.delete(1.0, tk.END))
                self.root.after(0, lambda: self.report_text.insert(tk.END, "保存的文件列表:\n"))
                self.root.after(0, lambda: self.report_text.insert(tk.END, "=" * 60 + "\n\n"))

                for file_path in saved_files:
                    self.root.after(0, lambda fp=file_path: self.report_text.insert(tk.END, f"• {Path(fp).name}\n"))

                self.root.after(0, lambda: self.report_text.insert(tk.END, f"\n所有文件已保存到: {output_dir}"))
            else:
                self.root.after(0, lambda: self.output_status_var.set("保存失败"))
                self.root.after(0, lambda: self.log_message(f"✗ {msg}", "ERROR"))

            self.root.after(0, lambda: self.update_status("就绪"))

        except Exception as e:
            self.root.after(0, lambda: self.output_status_var.set("保存失败"))
            self.root.after(0, lambda: self.log_message(f"✗ 保存结果失败: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("保存失败"))

    def export_analysis_results(self):
        """导出分析结果"""
        if not self.analyzer.results:
            messagebox.showwarning("警告", "没有可导出的分析结果")
            return

        file_path = filedialog.asksaveasfilename(
            title="导出分析结果",
            defaultextension=".csv",
            filetypes=[("CSV文件", "*.csv"), ("Excel文件", "*.xlsx"), ("所有文件", "*.*")]
        )

        if file_path:
            try:
                # 转换结果为DataFrame
                results_list = []
                for method_key, result in self.analyzer.results.items():
                    if method_key.startswith('method'):
                        row = {
                            '方法': result.get('description', method_key),
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

                    self.log_message(f"✓ 分析结果已导出到: {file_path}", "SUCCESS")

            except Exception as e:
                self.log_message(f"✗ 导出分析结果失败: {str(e)}", "ERROR")

    def open_output_dir(self):
        """打开输出目录"""
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
                messagebox.showerror("错误", f"无法打开目录: {str(e)}")
        else:
            messagebox.showwarning("警告", "输出目录不存在或未设置")

    def clear_analysis_text(self):
        """清空分析文本"""
        self.analysis_text.delete(1.0, tk.END)

    def clear_log(self):
        """清空日志"""
        self.log_text.delete(1.0, tk.END)

    def save_log(self):
        """保存日志"""
        file_path = filedialog.asksaveasfilename(
            title="保存日志",
            defaultextension=".txt",
            filetypes=[("文本文件", "*.txt"), ("所有文件", "*.*")]
        )

        if file_path:
            try:
                with open(file_path, 'w', encoding='utf-8') as f:
                    f.write(self.log_text.get(1.0, tk.END))
                self.log_message(f"✓ 日志已保存到: {file_path}", "SUCCESS")
            except Exception as e:
                self.log_message(f"✗ 保存日志失败: {str(e)}", "ERROR")

    def log_message(self, message: str, level: str = "INFO"):
        """添加日志消息"""
        timestamp = datetime.now().strftime("%H:%M:%S")
        log_message = f"[{timestamp}] {message}"

        self.log_text.insert(tk.END, log_message + "\n", level)
        self.log_text.see(tk.END)
        self.root.update()

    def analysis_message(self, message: str, level: str = "INFO"):
        """添加分析消息"""
        self.analysis_text.insert(tk.END, message + "\n", level)
        self.analysis_text.see(tk.END)
        self.root.update()

    def update_status(self, message: str):
        """更新状态栏"""
        self.status_var.set(message)
        self.root.update()

    def on_closing(self):
        """关闭窗口事件"""
        if messagebox.askyesno("确认", "确定要退出程序吗？"):
            # 清理资源
            if self.current_figure:
                plt.close(self.current_figure)
            self.root.destroy()


def check_dependencies():
    """检查依赖库"""
    print("=" * 70)
    print("检查依赖库...")
    print("=" * 70)

    dependencies = {
        'pandas': '数据处理库',
        'numpy': '数值计算库',
        'scipy': '科学计算库',
        'matplotlib': '绘图库',
        'seaborn': '统计可视化库',
        'scikit-learn': '机器学习库',
        'tkinter': 'GUI库'
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
            print(f"✗ {lib}: {desc} - 未安装")
            missing.append(lib)

    if missing:
        print(f"\n缺少以下库:")
        for lib in missing:
            if lib == 'tkinter':
                print("  注意: tkinter通常随Python安装，如未安装请:")
                print("    Ubuntu/Debian: sudo apt-get install python3-tk")
                print("    Windows/macOS: 重新安装Python并选择安装tkinter")
            else:
                print(f"  pip install {lib}")

    print("\n" + "=" * 70)
    return len(missing) == 0


def main():
    """主函数"""

    # 检查GUI支持
    if not GUI_AVAILABLE:
        print("错误: tkinter未安装，无法启动GUI")
        print("请安装tkinter:")
        print("  Ubuntu/Debian: sudo apt-get install python3-tk")
        print("  Windows/macOS: 重新安装Python并选择安装tkinter")
        return

    # 检查依赖
    if not check_dependencies():
        response = input("\n缺少依赖库，是否继续? (y/n): ")
        if response.lower() != 'y':
            return

    # 创建主窗口
    root = tk.Tk()

    # 设置窗口标题和大小
    root.title("PPG保留指数分析系统 - 试验方案实现")

    # 获取屏幕尺寸
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()

    # 设置窗口大小（屏幕的85%）
    window_width = int(screen_width * 0.85)
    window_height = int(screen_height * 0.85)

    # 计算窗口位置（居中）
    x = (screen_width - window_width) // 2
    y = (screen_height - window_height) // 2

    root.geometry(f"{window_width}x{window_height}+{x}+{y}")

    # 创建GUI
    app = PPGExperimentGUI(root)

    # 运行主循环
    root.mainloop()


if __name__ == "__main__":
    print("PPG保留指数分析系统 - 试验方案实现")
    print("=" * 70)
    print("功能:")
    print("  1. 数据加载与管理（PPG、SMRT、验证集、校正化合物）")
    print("  2. PPG保留指数计算与转换")
    print("  3. 五种校正方法分析比较")
    print("  4. 模型性能评估")
    print("  5. 可视化图表生成")
    print("  6. 实验报告导出")
    print("=" * 70)

    main()