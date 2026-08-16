
# !/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PPG保留指数计算与可视化分析程序 - 增强版（含跨条件转换功能）
功能：
1. 加载PPG标准品和化合物的保留时间数据
2. 拟合PPG标准曲线并计算线性关系
3. 计算化合物的PPG保留指数
4. 比较不同色谱条件下的PPG指数
5. 跨条件PPG指数转换与验证
6. 可视化分析结果并生成报告

基于实验方案，实现统一的液相色谱保留指数框架
"""

import os
import sys
import threading
import traceback
import warnings
from datetime import datetime
from typing import Dict, List, Tuple, Optional, Union, Any
from pathlib import Path

# 抑制警告
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
    print("错误: 请先安装必要的库")
    print("安装命令: pip install pandas numpy scipy matplotlib seaborn")
    sys.exit(1)

# 尝试导入GUI库
try:
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox, scrolledtext, Toplevel

    GUI_AVAILABLE = True
except ImportError:
    GUI_AVAILABLE = False
    print("警告: tkinter未安装，GUI不可用")


class PPGIndexCalculator:
    """PPG保留指数计算器核心类"""

    def __init__(self):
        """初始化计算器"""
        self.ppg_data = {}  # 存储不同条件下的PPG数据
        self.compound_data = {}  # 存储化合物数据
        self.standard_curves = {}  # 存储标准曲线参数
        self.ppg_indices = {}  # 存储计算的PPG指数
        self.results_summary = {}  # 存储结果汇总
        self.conversion_results = {}  # 存储跨条件转换结果

    def load_ppg_data(self, file_path: str, condition: str = "default") -> Tuple[bool, str]:
        """
        加载PPG标准品数据

        参数:
            file_path: 数据文件路径（支持Excel和CSV）
            condition: 色谱条件标识

        返回:
            (成功标志, 消息)
        """
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

    def load_compound_data(self, file_path: str, category: str = "validation",
                           condition: str = "default") -> Tuple[bool, str]:
        """
        加载化合物数据

        参数:
            file_path: 数据文件路径
            category: 数据类别（'validation'验证集，'smrt'训练集等）
            condition: 色谱条件标识

        返回:
            (成功标志, 消息)
        """
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
                'CAS': ['CAS', 'CAS号', 'CAS No.', 'CAS号']
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
            key = f"{category}_{condition}"
            self.compound_data[key] = df

            return True, f"成功加载 {len(df)} 个化合物数据（类别: {category}, 条件: {condition}）"

        except Exception as e:
            return False, f"加载化合物数据失败: {str(e)}"

    def fit_standard_curve(self, condition: str = "default",
                           model_type: str = "linear") -> Tuple[bool, str]:
        """
        拟合PPG标准曲线

        参数:
            condition: 色谱条件标识
            model_type: 模型类型 ('linear'线性, 'logarithmic'对数)

        返回:
            (成功标志, 消息)
        """
        try:
            if condition not in self.ppg_data:
                return False, f"未找到条件 {condition} 的PPG数据"

            df = self.ppg_data[condition]

            if len(df) < 3:
                return False, "PPG数据点不足，至少需要3个点建立标准曲线"

            x = df['聚合度'].values
            y = df['保留时间'].values

            if model_type == "logarithmic":
                # 对数模型：RT = a + b * ln(n)
                x_fit = np.log(x)
                model_name = "对数模型 (RT = a + b * ln(n))"
            elif model_type == "linear":
                # 线性模型：RT = a + b * n
                x_fit = x
                model_name = "线性模型 (RT = a + b * n)"
            else:
                return False, f"不支持的模型类型: {model_type}"

            # 线性回归
            slope, intercept, r_value, p_value, std_err = stats.linregress(x_fit, y)

            # 计算预测值和残差
            y_pred = intercept + slope * x_fit
            residuals = y - y_pred

            # 存储标准曲线参数
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

            return True, f"标准曲线拟合成功: {model_name}, R² = {r_value ** 2:.4f}"

        except Exception as e:
            return False, f"拟合标准曲线失败: {str(e)}"

    def calculate_ppg_index(self, condition: str = "default",
                            method: str = "interpolation") -> Tuple[bool, str]:
        """
        计算PPG保留指数

        参数:
            condition: 色谱条件标识
            method: 计算方法 ('interpolation'线性插值, 'regression'回归)

        返回:
            (成功标志, 消息)
        """
        try:
            if condition not in self.ppg_data:
                return False, f"未找到条件 {condition} 的PPG数据"

            df_ppg = self.ppg_data[condition]

            if method == "interpolation":
                # 线性插值法
                ppg_rt = df_ppg['保留时间'].values
                ppg_n = df_ppg['聚合度'].values

                # 对所有化合物数据计算PPG指数
                indices = {}

                for key in self.compound_data:
                    if condition in key:
                        df_comp = self.compound_data[key]
                        results = []

                        for _, row in df_comp.iterrows():
                            rt = row['保留时间']
                            compound_name = row['化合物名称']

                            # 边界处理
                            if rt < ppg_rt[0]:
                                # 外推
                                n_calc = ppg_n[0] - (ppg_rt[0] - rt) / (ppg_rt[1] - ppg_rt[0]) * (ppg_n[1] - ppg_n[0])
                                if n_calc < 0:
                                    n_calc = 0
                                method_used = "外推（低于最小PPG RT）"
                            elif rt > ppg_rt[-1]:
                                # 外推
                                n_calc = ppg_n[-1] + (rt - ppg_rt[-1]) / (ppg_rt[-1] - ppg_rt[-2]) * (
                                        ppg_n[-1] - ppg_n[-2])
                                method_used = "外推（高于最大PPG RT）"
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
                                method_used = "线性插值"

                            # 标准化为保留指数（乘以100）
                            ppg_index = n_calc * 100

                            result = {
                                '化合物名称': compound_name,
                                '保留时间': rt,
                                '计算PPG指数': ppg_index,
                                '计算方法': method_used,
                            }

                            # 添加其他列
                            for col in df_comp.columns:
                                if col not in result:
                                    result[col] = row[col]

                            results.append(result)

                        indices[key] = pd.DataFrame(results)

            elif method == "regression":
                # 回归法（使用标准曲线）
                if condition not in self.standard_curves:
                    success, msg = self.fit_standard_curve(condition, model_type="linear")
                    if not success:
                        return False, f"无法使用回归法: {msg}"

                curve = self.standard_curves[condition]

                # 对每个化合物计算PPG指数
                indices = {}

                for key in self.compound_data:
                    if condition in key:
                        df_comp = self.compound_data[key]
                        results = []

                        for _, row in df_comp.iterrows():
                            rt = row['保留时间']
                            compound_name = row['化合物名称']

                            # 从RT反算n
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
                                '化合物名称': compound_name,
                                '保留时间': rt,
                                '计算PPG指数': ppg_index,
                                '计算方法': "回归法",
                            }

                            # 添加其他列
                            for col in df_comp.columns:
                                if col not in result:
                                    result[col] = row[col]

                            results.append(result)

                        indices[key] = pd.DataFrame(results)
            else:
                return False, f"不支持的计算方法: {method}"

            # 存储结果
            self.ppg_indices[condition] = {
                'method': method,
                'indices': indices
            }

            return True, f"PPG指数计算完成（条件: {condition}, 方法: {method}）"

        except Exception as e:
            return False, f"计算PPG指数失败: {str(e)}"

    def compare_conditions(self, conditions: List[str]) -> pd.DataFrame:
        """
        比较不同色谱条件下的PPG指数

        参数:
            conditions: 要比较的色谱条件列表

        返回:
            比较结果的DataFrame
        """
        comparison_results = []

        for key in self.compound_data:
            # 检查是否包含要比较的条件
            if any(cond in key for cond in conditions):
                category = key.split('_')[0]
                condition = key.split('_')[1] if '_' in key else "default"

                # 获取该条件下的化合物数据
                df_comp = self.compound_data[key]

                for _, row in df_comp.iterrows():
                    compound_name = row['化合物名称']
                    rt = row['保留时间']

                    # 查找该化合物在其他条件下的数据
                    compound_data = {
                        '化合物名称': compound_name,
                        '数据类别': category,
                        f'{condition}_RT': rt
                    }

                    # 添加PPG指数（如果已计算）
                    for cond in conditions:
                        if cond in self.ppg_indices:
                            for data_key, indices_df in self.ppg_indices[cond]['indices'].items():
                                if cond in data_key:
                                    match = indices_df[indices_df['化合物名称'] == compound_name]
                                    if not match.empty:
                                        compound_data[f'{cond}_PPG指数'] = match.iloc[0]['计算PPG指数']
                                        break

                    comparison_results.append(compound_data)

        return pd.DataFrame(comparison_results)

    def calculate_conversion_error(self, from_condition: str, to_condition: str) -> pd.DataFrame:
        """
        计算条件转换误差

        参数:
            from_condition: 源色谱条件
            to_condition: 目标色谱条件

        返回:
            转换误差分析结果
        """
        error_results = []

        # 获取两个条件下的PPG指数
        if from_condition not in self.ppg_indices or to_condition not in self.ppg_indices:
            return pd.DataFrame()

        # 查找两个条件下都有的化合物
        compounds_in_both = set()

        for key in self.compound_data:
            if from_condition in key:
                df = self.compound_data[key]
                compounds_in_both.update(df['化合物名称'].tolist())

        for key in self.compound_data:
            if to_condition in key:
                df = self.compound_data[key]
                compounds_in_both.intersection_update(set(df['化合物名称'].tolist()))

        # 计算误差
        for compound in compounds_in_both:
            # 获取两个条件下的PPG指数
            from_ppg = None
            to_ppg = None

            for key, indices_dict in self.ppg_indices.items():
                if from_condition in key:
                    for data_key, df_indices in indices_dict['indices'].items():
                        match = df_indices[df_indices['化合物名称'] == compound]
                        if not match.empty:
                            from_ppg = match.iloc[0]['计算PPG指数']
                            break

                if to_condition in key:
                    for data_key, df_indices in indices_dict['indices'].items():
                        match = df_indices[df_indices['化合物名称'] == compound]
                        if not match.empty:
                            to_ppg = match.iloc[0]['计算PPG指数']
                            break

            if from_ppg is not None and to_ppg is not None:
                absolute_error = abs(from_ppg - to_ppg)
                relative_error = (absolute_error / from_ppg * 100) if from_ppg != 0 else np.inf

                error_results.append({
                    '化合物名称': compound,
                    f'{from_condition}_PPG指数': from_ppg,
                    f'{to_condition}_PPG指数': to_ppg,
                    '绝对误差': absolute_error,
                    '相对误差(%)': relative_error
                })

        return pd.DataFrame(error_results)

    def convert_ppg_index_to_rt(self, from_condition: str, to_condition: str,
                                compound_names: List[str] = None) -> Tuple[pd.DataFrame, Union[Dict, str]]:
        """
        将源条件下的PPG指数转换为目标条件下的保留时间

        参数:
            from_condition: 源色谱条件（PPG指数来源）
            to_condition: 目标色谱条件（要转换到的条件）
            compound_names: 要转换的化合物名称列表（None表示所有化合物）

        返回:
            (转换结果DataFrame, 统计信息字典或错误消息)
        """
        try:
            # 检查是否有源条件的PPG指数
            if from_condition not in self.ppg_indices:
                return pd.DataFrame(), f"源条件 {from_condition} 没有PPG指数数据"

            # 检查是否有目标条件的标准曲线
            if to_condition not in self.standard_curves:
                # 尝试拟合标准曲线
                success, msg = self.fit_standard_curve(to_condition, "linear")
                if not success:
                    return pd.DataFrame(), f"无法为目标条件 {to_condition} 拟合标准曲线: {msg}"

            # 获取目标条件标准曲线
            curve = self.standard_curves[to_condition]

            # 收集转换结果
            conversion_results = []

            # 遍历源条件的PPG指数
            for key, indices_df in self.ppg_indices[from_condition]['indices'].items():
                if from_condition in key:
                    for _, row in indices_df.iterrows():
                        compound_name = row['化合物名称']

                        # 如果指定了化合物列表，则只处理指定的化合物
                        if compound_names and compound_name not in compound_names:
                            continue

                        # 获取源条件的PPG指数
                        ppg_index = row['计算PPG指数']
                        if pd.isna(ppg_index):
                            continue

                        # 将PPG指数转换为聚合度n (除以100)
                        n_calc = ppg_index / 100

                        # 使用目标条件标准曲线计算保留时间
                        if curve['model_type'] == "logarithmic":
                            # 对数模型: RT = a + b * ln(n)
                            rt_pred = curve['intercept'] + curve['slope'] * np.log(n_calc)
                        else:  # linear
                            # 线性模型: RT = a + b * n
                            rt_pred = curve['intercept'] + curve['slope'] * n_calc

                        # 查找目标条件下的实际保留时间
                        rt_actual = None
                        for comp_key, comp_df in self.compound_data.items():
                            if to_condition in comp_key:
                                match = comp_df[comp_df['化合物名称'] == compound_name]
                                if not match.empty and '保留时间' in match.columns:
                                    rt_actual = match.iloc[0]['保留时间']
                                    break

                        # 计算误差
                        if rt_actual is not None:
                            absolute_error = abs(rt_pred - rt_actual)
                            relative_error = (absolute_error / rt_actual * 100) if rt_actual != 0 else np.nan
                        else:
                            absolute_error = np.nan
                            relative_error = np.nan

                        result = {
                            '化合物名称': compound_name,
                            f'{from_condition}_PPG指数': ppg_index,
                            f'{from_condition}_聚合度': n_calc,
                            f'{to_condition}_预测RT': rt_pred,
                            f'{to_condition}_实际RT': rt_actual,
                            '绝对误差(min)': absolute_error,
                            '相对误差(%)': relative_error,
                            '源条件': from_condition,
                            '目标条件': to_condition
                        }

                        conversion_results.append(result)

            conversion_df = pd.DataFrame(conversion_results)

            # 计算整体统计
            if not conversion_df.empty:
                valid_errors = conversion_df['绝对误差(min)'].dropna()
                if len(valid_errors) > 0:
                    stats = {
                        '平均绝对误差': valid_errors.mean(),
                        '绝对误差标准差': valid_errors.std(),
                        '最大绝对误差': valid_errors.max(),
                        '最小绝对误差': valid_errors.min(),
                        '中位绝对误差': valid_errors.median(),
                        '样本数': len(valid_errors)
                    }

                    # 计算相对误差统计
                    valid_rel_errors = conversion_df['相对误差(%)'].dropna()
                    if len(valid_rel_errors) > 0:
                        stats.update({
                            '平均相对误差(%)': valid_rel_errors.mean(),
                            '相对误差标准差(%)': valid_rel_errors.std(),
                            '最大相对误差(%)': valid_rel_errors.max(),
                            '最小相对误差(%)': valid_rel_errors.min()
                        })

                    return conversion_df, stats
                else:
                    return conversion_df, "没有有效的误差数据"
            else:
                return conversion_df, "没有找到匹配的化合物数据"

        except Exception as e:
            return pd.DataFrame(), f"转换失败: {str(e)}"

    def cross_condition_analysis(self, from_condition: str, to_condition: str,
                                 threshold: float = 0.5) -> Dict[str, Any]:
        """
        跨条件转换的全面分析

        参数:
            from_condition: 源色谱条件
            to_condition: 目标色谱条件
            threshold: 误差阈值（分钟）

        返回:
            分析结果字典
        """
        try:
            # 执行转换
            conversion_df, stats = self.convert_ppg_index_to_rt(from_condition, to_condition)

            if conversion_df.empty:
                return {"error": "没有转换数据"}

            # 基础统计
            analysis_results = {
                '源条件': from_condition,
                '目标条件': to_condition,
                '转换时间': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                '总化合物数': len(conversion_df),
                '有效转换数': len(conversion_df['绝对误差(min)'].dropna()),
                '误差统计': stats if isinstance(stats, dict) else stats,
                '详细数据': conversion_df.to_dict('records')
            }

            # 添加误差分类
            if '绝对误差(min)' in conversion_df.columns:
                errors = conversion_df['绝对误差(min)'].dropna()

                # 误差分布统计
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

                analysis_results['误差分布'] = error_distribution

                # 通过率分析（基于阈值）
                passed = len(errors[errors <= threshold])
                pass_rate = (passed / len(errors) * 100) if len(errors) > 0 else 0

                analysis_results['通过率分析'] = {
                    '阈值(min)': threshold,
                    '通过数': passed,
                    '总数': len(errors),
                    '通过率(%)': pass_rate
                }

            # 存储结果
            key = f"{from_condition}_to_{to_condition}"
            self.conversion_results[key] = analysis_results

            return analysis_results

        except Exception as e:
            return {"error": f"分析失败: {str(e)}"}

    def generate_summary_report(self) -> Dict[str, Any]:
        """
        生成分析报告摘要

        返回:
            报告摘要字典
        """
        summary = {
            '生成时间': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'PPG数据条件数': len(self.ppg_data),
            '化合物数据集数': len(self.compound_data),
            '标准曲线数': len(self.standard_curves),
            'PPG指数计算结果数': len(self.ppg_indices),
            '跨条件转换结果数': len(self.conversion_results),
            '标准曲线性能': {},
            'PPG指数统计': {},
            '跨条件转换统计': {}
        }

        # 标准曲线性能
        for condition, curve in self.standard_curves.items():
            summary['标准曲线性能'][condition] = {
                '模型类型': curve['model_type'],
                'R²': curve['r_squared'],
                '斜率': curve['slope'],
                '截距': curve['intercept'],
                '标准误差': curve['std_err'],
                '数据点数': curve['n_points']
            }

        # PPG指数统计
        for condition, indices_data in self.ppg_indices.items():
            all_indices = []
            for key, df in indices_data['indices'].items():
                if '计算PPG指数' in df.columns:
                    valid_indices = df['计算PPG指数'].dropna()
                    all_indices.extend(valid_indices.tolist())

            if all_indices:
                indices_array = np.array(all_indices)
                summary['PPG指数统计'][condition] = {
                    '计算方法': indices_data['method'],
                    '样本数': len(all_indices),
                    '平均值': np.mean(indices_array),
                    '标准差': np.std(indices_array),
                    '最小值': np.min(indices_array),
                    '最大值': np.max(indices_array),
                    '中位数': np.median(indices_array)
                }

        # 跨条件转换统计
        for key, conversion in self.conversion_results.items():
            summary['跨条件转换统计'][key] = {
                '源条件': conversion.get('源条件', ''),
                '目标条件': conversion.get('目标条件', ''),
                '总化合物数': conversion.get('总化合物数', 0),
                '有效转换数': conversion.get('有效转换数', 0),
                '平均绝对误差': conversion.get('误差统计', {}).get('平均绝对误差', 0) if isinstance(
                    conversion.get('误差统计'), dict) else 0,
                '通过率(%)': conversion.get('通过率分析', {}).get('通过率(%)', 0) if isinstance(
                    conversion.get('通过率分析'), dict) else 0
            }

        self.results_summary = summary
        return summary

    def save_results(self, output_dir: str) -> Tuple[bool, str, List[str]]:
        """
        保存所有结果到文件

        参数:
            output_dir: 输出目录

        返回:
            (成功标志, 消息, 保存的文件列表)
        """
        try:
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)

            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            saved_files = []

            # 1. 保存PPG标准曲线数据
            if self.standard_curves:
                curves_data = []
                for condition, curve in self.standard_curves.items():
                    curves_data.append({
                        '色谱条件': condition,
                        '模型类型': curve['model_type'],
                        '斜率': curve['slope'],
                        '截距': curve['intercept'],
                        'R²': curve['r_squared'],
                        'p值': curve['p_value'],
                        '标准误差': curve['std_err'],
                        '数据点数': curve['n_points']
                    })

                curves_df = pd.DataFrame(curves_data)
                curves_file = output_path / f"PPG标准曲线_{timestamp}.xlsx"
                with pd.ExcelWriter(curves_file, engine='openpyxl') as writer:
                    curves_df.to_excel(writer, sheet_name='标准曲线汇总', index=False)

                    # 每个条件的详细数据
                    for condition, curve in self.standard_curves.items():
                        detail_df = pd.DataFrame({
                            '聚合度': curve['x'],
                            '实测RT': curve['y'],
                            '预测RT': curve['y_pred'],
                            '残差': curve['residuals']
                        })
                        detail_df.to_excel(writer, sheet_name=f'{condition}_详细数据', index=False)

                saved_files.append(str(curves_file))

            # 2. 保存PPG指数计算结果
            if self.ppg_indices:
                for condition, indices_data in self.ppg_indices.items():
                    indices_file = output_path / f"PPG指数_{condition}_{timestamp}.xlsx"

                    with pd.ExcelWriter(indices_file, engine='openpyxl') as writer:
                        for key, df in indices_data['indices'].items():
                            # 简化sheet名称
                            sheet_name = key.replace('_', '-')[:30]
                            if sheet_name in writer.sheets:
                                sheet_name = f"{sheet_name[:25]}_{hash(key) % 1000:03d}"

                            df.to_excel(writer, sheet_name=sheet_name, index=False)

                    saved_files.append(str(indices_file))

            # 3. 保存跨条件转换结果
            if self.conversion_results:
                conversion_file = output_path / f"跨条件转换结果_{timestamp}.xlsx"
                with pd.ExcelWriter(conversion_file, engine='openpyxl') as writer:
                    for key, conversion in self.conversion_results.items():
                        if '详细数据' in conversion:
                            df = pd.DataFrame(conversion['详细数据'])
                            # 简化sheet名称
                            sheet_name = key[:30]
                            if sheet_name in writer.sheets:
                                sheet_name = f"{sheet_name[:25]}_{hash(key) % 1000:03d}"
                            df.to_excel(writer, sheet_name=sheet_name, index=False)

                    # 保存转换统计
                    conversion_stats = []
                    for key, conversion in self.conversion_results.items():
                        stats = {
                            '转换方向': key,
                            '源条件': conversion.get('源条件', ''),
                            '目标条件': conversion.get('目标条件', ''),
                            '总化合物数': conversion.get('总化合物数', 0),
                            '有效转换数': conversion.get('有效转换数', 0)
                        }

                        if isinstance(conversion.get('误差统计'), dict):
                            stats.update({
                                '平均绝对误差(min)': conversion['误差统计'].get('平均绝对误差', 0),
                                '绝对误差标准差': conversion['误差统计'].get('绝对误差标准差', 0),
                                '最大绝对误差(min)': conversion['误差统计'].get('最大绝对误差', 0)
                            })

                        if isinstance(conversion.get('通过率分析'), dict):
                            stats.update({
                                '通过率(%)': conversion['通过率分析'].get('通过率(%)', 0),
                                '通过数': conversion['通过率分析'].get('通过数', 0),
                                '阈值(min)': conversion['通过率分析'].get('阈值(min)', 0.5)
                            })

                        conversion_stats.append(stats)

                    if conversion_stats:
                        stats_df = pd.DataFrame(conversion_stats)
                        stats_df.to_excel(writer, sheet_name='转换统计汇总', index=False)

                saved_files.append(str(conversion_file))

            # 4. 保存分析报告
            if self.results_summary:
                report_file = output_path / f"分析报告_{timestamp}.txt"
                with open(report_file, 'w', encoding='utf-8') as f:
                    f.write("=" * 70 + "\n")
                    f.write("PPG保留指数分析报告\n")
                    f.write("=" * 70 + "\n\n")

                    f.write(f"报告生成时间: {self.results_summary['生成时间']}\n\n")

                    f.write("数据概览:\n")
                    f.write(f"  - PPG数据条件数: {self.results_summary['PPG数据条件数']}\n")
                    f.write(f"  - 化合物数据集数: {self.results_summary['化合物数据集数']}\n")
                    f.write(f"  - 标准曲线数: {self.results_summary['标准曲线数']}\n")
                    f.write(f"  - PPG指数计算结果数: {self.results_summary['PPG指数计算结果数']}\n")
                    f.write(f"  - 跨条件转换结果数: {self.results_summary['跨条件转换结果数']}\n\n")

                    if self.results_summary['标准曲线性能']:
                        f.write("标准曲线性能:\n")
                        for condition, perf in self.results_summary['标准曲线性能'].items():
                            f.write(f"  {condition}:\n")
                            f.write(f"    - 模型类型: {perf['模型类型']}\n")
                            f.write(f"    - R²: {perf['R²']:.4f}\n")
                            f.write(f"    - 斜率: {perf['斜率']:.4f}\n")
                            f.write(f"    - 截距: {perf['截距']:.4f}\n")
                            f.write(f"    - 标准误差: {perf['标准误差']:.4f}\n")
                            f.write(f"    - 数据点数: {perf['数据点数']}\n")
                        f.write("\n")

                    if self.results_summary['PPG指数统计']:
                        f.write("PPG指数统计:\n")
                        for condition, stats in self.results_summary['PPG指数统计'].items():
                            f.write(f"  {condition}:\n")
                            f.write(f"    - 计算方法: {stats['计算方法']}\n")
                            f.write(f"    - 样本数: {stats['样本数']}\n")
                            f.write(f"    - 平均值: {stats['平均值']:.2f}\n")
                            f.write(f"    - 标准差: {stats['标准差']:.2f}\n")
                            f.write(f"    - 范围: {stats['最小值']:.2f} - {stats['最大值']:.2f}\n")
                            f.write(f"    - 中位数: {stats['中位数']:.2f}\n")
                        f.write("\n")

                    if self.results_summary['跨条件转换统计']:
                        f.write("跨条件转换统计:\n")
                        for key, stats in self.results_summary['跨条件转换统计'].items():
                            f.write(f"  {key}:\n")
                            f.write(f"    - 源条件: {stats['源条件']}\n")
                            f.write(f"    - 目标条件: {stats['目标条件']}\n")
                            f.write(f"    - 总化合物数: {stats['总化合物数']}\n")
                            f.write(f"    - 有效转换数: {stats['有效转换数']}\n")
                            f.write(f"    - 平均绝对误差: {stats['平均绝对误差']:.3f} min\n")
                            f.write(f"    - 通过率: {stats['通过率(%)']:.1f}%\n")
                        f.write("\n")

                saved_files.append(str(report_file))

            # 5. 保存化合物数据汇总
            if self.compound_data:
                compounds_file = output_path / f"化合物数据汇总_{timestamp}.xlsx"
                with pd.ExcelWriter(compounds_file, engine='openpyxl') as writer:
                    for key, df in self.compound_data.items():
                        sheet_name = key.replace('_', '-')[:30]
                        if sheet_name in writer.sheets:
                            sheet_name = f"{sheet_name[:25]}_{hash(key) % 1000:03d}"
                        df.to_excel(writer, sheet_name=sheet_name, index=False)

                saved_files.append(str(compounds_file))

            return True, f"结果已保存到 {output_dir}", saved_files

        except Exception as e:
            return False, f"保存结果失败: {str(e)}", []


class PPGVisualizer:
    """PPG数据可视化类"""

    def __init__(self, calculator: PPGIndexCalculator):
        """初始化可视化器"""
        self.calculator = calculator
        self.figures = {}

        # 设置中文字体
        plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
        plt.rcParams['axes.unicode_minus'] = False

    def plot_standard_curves(self, conditions: List[str] = None,
                             save_path: str = None) -> plt.Figure:
        """
        绘制PPG标准曲线

        参数:
            conditions: 要绘制的色谱条件列表（None表示全部）
            save_path: 保存路径（可选）

        返回:
            matplotlib Figure对象
        """
        if conditions is None:
            conditions = list(self.calculator.standard_curves.keys())

        if not conditions:
            print("警告: 没有可绘制的标准曲线数据")
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

            # 原始数据点
            ax.scatter(curve['x'], curve['y'], color='blue', s=50, label='实测数据', zorder=3)

            # 拟合曲线
            if curve['model_type'] == 'logarithmic':
                x_fit = np.log(curve['x'])
                x_range = np.linspace(min(curve['x']), max(curve['x']), 100)
                x_fit_range = np.log(x_range)
            else:
                x_fit = curve['x']
                x_range = np.linspace(min(curve['x']), max(curve['x']), 100)
                x_fit_range = x_range

            y_fit_range = curve['intercept'] + curve['slope'] * x_fit_range
            ax.plot(x_range, y_fit_range, 'r-', label='拟合曲线', linewidth=2)

            # 添加回归信息
            info_text = f"模型: {curve['model_name']}\n"
            info_text += f"R² = {curve['r_squared']:.4f}\n"
            info_text += f"斜率 = {curve['slope']:.4f}\n"
            info_text += f"截距 = {curve['intercept']:.4f}"

            ax.text(0.05, 0.95, info_text, transform=ax.transAxes,
                    verticalalignment='top', fontsize=10,
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

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

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"标准曲线图已保存到: {save_path}")

        self.figures['standard_curves'] = fig
        return fig

    def plot_residuals(self, conditions: List[str] = None,
                       save_path: str = None) -> plt.Figure:
        """
        绘制残差图

        参数:
            conditions: 要绘制的色谱条件列表（None表示全部）
            save_path: 保存路径（可选）

        返回:
            matplotlib Figure对象
        """
        if conditions is None:
            conditions = list(self.calculator.standard_curves.keys())

        if not conditions:
            print("警告: 没有可绘制的残差数据")
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

            # 残差图
            residuals = curve['residuals']
            predicted = curve['y_pred']

            ax.scatter(predicted, residuals, color='blue', s=50, alpha=0.7)
            ax.axhline(y=0, color='red', linestyle='--', linewidth=1)

            # 添加残差统计
            mean_residual = np.mean(residuals)
            std_residual = np.std(residuals)

            info_text = f"残差统计:\n"
            info_text += f"均值 = {mean_residual:.4f}\n"
            info_text += f"标准差 = {std_residual:.4f}\n"
            info_text += f"最大残差 = {max(abs(residuals)):.4f}"

            ax.text(0.05, 0.95, info_text, transform=ax.transAxes,
                    verticalalignment='top', fontsize=9,
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

            ax.set_xlabel('预测保留时间 (min)')
            ax.set_ylabel('残差 (min)')
            ax.set_title(f'残差图 - {condition}')
            ax.grid(True, alpha=0.3)

        # 隐藏多余的子图
        for idx in range(len(conditions), n_rows * n_cols):
            row = idx // n_cols
            col = idx % n_cols
            axes[row, col].set_visible(False)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"残差图已保存到: {save_path}")

        self.figures['residuals'] = fig
        return fig

    def plot_ppg_index_distribution(self, condition: str,
                                    save_path: str = None) -> plt.Figure:
        """
        绘制PPG指数分布图

        参数:
            condition: 色谱条件
            save_path: 保存路径（可选）

        返回:
            matplotlib Figure对象
        """
        if condition not in self.calculator.ppg_indices:
            print(f"警告: 条件 {condition} 没有PPG指数数据")
            return None

        # 收集所有PPG指数
        all_indices = []
        for key, df in self.calculator.ppg_indices[condition]['indices'].items():
            if '计算PPG指数' in df.columns:
                indices = df['计算PPG指数'].dropna().tolist()
                all_indices.extend(indices)

        if not all_indices:
            print(f"警告: 条件 {condition} 没有有效的PPG指数")
            return None

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # 直方图
        ax1.hist(all_indices, bins=30, color='steelblue', edgecolor='black', alpha=0.7)
        ax1.axvline(x=np.mean(all_indices), color='red', linestyle='--',
                    linewidth=2, label=f'均值: {np.mean(all_indices):.1f}')
        ax1.axvline(x=np.median(all_indices), color='green', linestyle='--',
                    linewidth=2, label=f'中位数: {np.median(all_indices):.1f}')

        ax1.set_xlabel('PPG保留指数')
        ax1.set_ylabel('频数')
        ax1.set_title(f'PPG指数分布 - {condition}')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # 箱线图
        ax2.boxplot(all_indices, vert=True, patch_artist=True,
                    boxprops=dict(facecolor='lightblue', color='blue'),
                    medianprops=dict(color='red', linewidth=2))

        # 添加统计信息
        stats_text = f"统计信息:\n"
        stats_text += f"样本数: {len(all_indices)}\n"
        stats_text += f"均值: {np.mean(all_indices):.1f}\n"
        stats_text += f"标准差: {np.std(all_indices):.1f}\n"
        stats_text += f"最小值: {np.min(all_indices):.1f}\n"
        stats_text += f"最大值: {np.max(all_indices):.1f}\n"
        stats_text += f"中位数: {np.median(all_indices):.1f}"

        ax2.text(0.05, 0.95, stats_text, transform=ax2.transAxes,
                 verticalalignment='top', fontsize=10,
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

        ax2.set_ylabel('PPG保留指数')
        ax2.set_title(f'PPG指数箱线图 - {condition}')
        ax2.set_xticks([1])
        ax2.set_xticklabels([condition])
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"PPG指数分布图已保存到: {save_path}")

        self.figures['index_distribution'] = fig
        return fig

    def plot_condition_comparison(self, conditions: List[str],
                                  save_path: str = None) -> plt.Figure:
        """
        绘制不同条件下的PPG指数比较

        参数:
            conditions: 要比较的色谱条件列表
            save_path: 保存路径（可选）

        返回:
            matplotlib Figure对象
        """
        # 获取比较数据
        comparison_df = self.calculator.compare_conditions(conditions)

        if comparison_df.empty:
            print("警告: 没有可用于比较的数据")
            return None

        # 提取每个条件下的PPG指数
        fig, axes = plt.subplots(2, 2, figsize=(14, 12))
        axes = axes.flatten()

        # 1. 散点图比较（两两比较）
        if len(conditions) >= 2:
            ax = axes[0]
            cond1, cond2 = conditions[0], conditions[1]

            # 获取两个条件下都有的化合物
            compounds = set()
            for cond in [cond1, cond2]:
                col_name = f"{cond}_PPG指数"
                if col_name in comparison_df.columns:
                    valid_data = comparison_df[comparison_df[col_name].notna()]
                    compounds.update(valid_data['化合物名称'].tolist())

            # 提取数据
            data_points = []
            for compound in compounds:
                row = comparison_df[comparison_df['化合物名称'] == compound]
                if not row.empty:
                    val1 = row[f"{cond1}_PPG指数"].iloc[0] if f"{cond1}_PPG指数" in row.columns else np.nan
                    val2 = row[f"{cond2}_PPG指数"].iloc[0] if f"{cond2}_PPG指数" in row.columns else np.nan
                    if not (np.isnan(val1) or np.isnan(val2)):
                        data_points.append((val1, val2))

            if data_points:
                val1_vals, val2_vals = zip(*data_points)
                ax.scatter(val1_vals, val2_vals, alpha=0.6, s=50)

                # 添加对角线
                min_val = min(min(val1_vals), min(val2_vals))
                max_val = max(max(val1_vals), max(val2_vals))
                ax.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.5)

                ax.set_xlabel(f'{cond1} PPG指数')
                ax.set_ylabel(f'{cond2} PPG指数')
                ax.set_title(f'{cond1} vs {cond2} PPG指数比较')
                ax.grid(True, alpha=0.3)

        # 2. 多条件箱线图比较
        ax = axes[1]
        box_data = []
        labels = []

        for cond in conditions:
            col_name = f"{cond}_PPG指数"
            if col_name in comparison_df.columns:
                data = comparison_df[col_name].dropna().tolist()
                if data:
                    box_data.append(data)
                    labels.append(cond)

        if box_data:
            bp = ax.boxplot(box_data, labels=labels, patch_artist=True)

            # 设置颜色
            colors = ['lightblue', 'lightgreen', 'lightcoral', 'lightsalmon']
            for patch, color in zip(bp['boxes'], colors[:len(box_data)]):
                patch.set_facecolor(color)

            ax.set_ylabel('PPG保留指数')
            ax.set_title('不同条件下PPG指数分布比较')
            ax.grid(True, alpha=0.3)

        # 3. 误差分布图（如果有误差数据）
        ax = axes[2]
        if len(conditions) >= 2:
            error_df = self.calculator.calculate_conversion_error(conditions[0], conditions[1])
            if not error_df.empty and '相对误差(%)' in error_df.columns:
                errors = error_df['相对误差(%)'].dropna()
                if len(errors) > 0:
                    ax.hist(errors, bins=20, color='coral', edgecolor='black', alpha=0.7)
                    ax.axvline(x=np.mean(errors), color='red', linestyle='--',
                               linewidth=2, label=f'平均误差: {np.mean(errors):.2f}%')

                    ax.set_xlabel('相对误差 (%)')
                    ax.set_ylabel('频数')
                    ax.set_title(f'{conditions[0]} → {conditions[1]} 转换误差分布')
                    ax.legend()
                    ax.grid(True, alpha=0.3)

        # 4. 相关系数热图
        ax = axes[3]
        if len(conditions) >= 2:
            # 构建相关系数矩阵
            corr_data = []
            for cond in conditions:
                col_name = f"{cond}_PPG指数"
                if col_name in comparison_df.columns:
                    corr_data.append(comparison_df[col_name])

            if corr_data and len(corr_data) > 1:
                corr_df = pd.DataFrame(corr_data).T
                corr_df.columns = conditions[:len(corr_data)]
                corr_matrix = corr_df.corr()

                # 绘制热图
                im = ax.imshow(corr_matrix, cmap='RdYlBu', vmin=0, vmax=1)

                # 添加数值
                for i in range(len(corr_matrix)):
                    for j in range(len(corr_matrix)):
                        ax.text(j, i, f'{corr_matrix.iloc[i, j]:.3f}',
                                ha='center', va='center',
                                color='black' if abs(corr_matrix.iloc[i, j]) < 0.7 else 'white')

                ax.set_xticks(range(len(corr_matrix)))
                ax.set_yticks(range(len(corr_matrix)))
                ax.set_xticklabels(corr_matrix.columns)
                ax.set_yticklabels(corr_matrix.columns)
                ax.set_title('PPG指数相关系数矩阵')

                # 添加颜色条
                plt.colorbar(im, ax=ax)

        # 隐藏未使用的子图
        for i in range(len(conditions), 4):
            axes[i].set_visible(False)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"条件比较图已保存到: {save_path}")

        self.figures['condition_comparison'] = fig
        return fig

    def plot_conversion_analysis(self, from_condition: str, to_condition: str,
                                 save_path: str = None) -> plt.Figure:
        """
        绘制跨条件转换分析图

        参数:
            from_condition: 源色谱条件
            to_condition: 目标色谱条件
            save_path: 保存路径（可选）

        返回:
            matplotlib Figure对象
        """
        try:
            # 执行转换分析
            conversion_df, stats = self.calculator.convert_ppg_index_to_rt(from_condition, to_condition)

            if conversion_df.empty or '绝对误差(min)' not in conversion_df.columns:
                print("警告: 没有有效的转换数据")
                return None

            # 过滤有效数据
            valid_df = conversion_df.dropna(
                subset=['绝对误差(min)', f'{to_condition}_预测RT', f'{to_condition}_实际RT'])

            if valid_df.empty:
                print("警告: 没有有效的数据用于绘图")
                return None

            # 创建多子图
            fig, axes = plt.subplots(2, 3, figsize=(18, 12))
            axes = axes.flatten()

            # 1. 预测RT vs 实际RT 散点图
            ax = axes[0]
            predicted = valid_df[f'{to_condition}_预测RT']
            actual = valid_df[f'{to_condition}_实际RT']

            ax.scatter(actual, predicted, alpha=0.6, s=50, color='steelblue')

            # 添加对角线
            min_val = min(predicted.min(), actual.min())
            max_val = max(predicted.max(), actual.max())
            ax.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.5, label='y=x')

            # 添加回归线
            if len(predicted) > 1:
                slope, intercept, r_value, p_value, std_err = stats.linregress(actual, predicted)
                x_range = np.linspace(min_val, max_val, 100)
                y_pred = intercept + slope * x_range
                ax.plot(x_range, y_pred, 'g-', alpha=0.7,
                        label=f'回归线: R²={r_value ** 2:.3f}')

            ax.set_xlabel(f'{to_condition} 实际保留时间 (min)')
            ax.set_ylabel(f'{to_condition} 预测保留时间 (min)')
            ax.set_title(f'{from_condition}→{to_condition}: 预测 vs 实际保留时间')
            ax.legend()
            ax.grid(True, alpha=0.3)

            # 2. 误差分布直方图
            ax = axes[1]
            errors = valid_df['绝对误差(min)']

            ax.hist(errors, bins=20, color='coral', edgecolor='black', alpha=0.7)
            ax.axvline(x=errors.mean(), color='red', linestyle='--',
                       linewidth=2, label=f'均值: {errors.mean():.3f}min')

            # 添加统计信息
            stats_text = f"误差统计:\n"
            stats_text += f"样本数: {len(errors)}\n"
            stats_text += f"均值: {errors.mean():.3f}min\n"
            stats_text += f"标准差: {errors.std():.3f}min\n"
            stats_text += f"最大值: {errors.max():.3f}min\n"
            stats_text += f"中位数: {errors.median():.3f}min"

            ax.text(0.65, 0.95, stats_text, transform=ax.transAxes,
                    verticalalignment='top', fontsize=10,
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

            ax.set_xlabel('绝对误差 (min)')
            ax.set_ylabel('频数')
            ax.set_title(f'{from_condition}→{to_condition}: 误差分布')
            ax.legend()
            ax.grid(True, alpha=0.3)

            # 3. 相对误差分布
            ax = axes[2]
            if '相对误差(%)' in valid_df.columns:
                rel_errors = valid_df['相对误差(%)'].dropna()
                if len(rel_errors) > 0:
                    ax.hist(rel_errors, bins=20, color='lightgreen', edgecolor='black', alpha=0.7)
                    ax.axvline(x=rel_errors.mean(), color='green', linestyle='--',
                               linewidth=2, label=f'均值: {rel_errors.mean():.2f}%')

                    ax.set_xlabel('相对误差 (%)')
                    ax.set_ylabel('频数')
                    ax.set_title(f'{from_condition}→{to_condition}: 相对误差分布')
                    ax.legend()
                    ax.grid(True, alpha=0.3)

            # 4. 误差 vs RT 散点图
            ax = axes[3]
            ax.scatter(actual, errors, alpha=0.6, s=50, color='purple')
            ax.axhline(y=errors.mean(), color='red', linestyle='--',
                       linewidth=1, label=f'平均误差: {errors.mean():.3f}min')

            # 添加趋势线
            if len(actual) > 1:
                z = np.polyfit(actual, errors, 1)
                p = np.poly1d(z)
                ax.plot(actual, p(actual), "b--", alpha=0.5, label='趋势线')

            ax.set_xlabel(f'{to_condition} 实际保留时间 (min)')
            ax.set_ylabel('绝对误差 (min)')
            ax.set_title(f'{from_condition}→{to_condition}: 误差 vs 保留时间')
            ax.legend()
            ax.grid(True, alpha=0.3)

            # 5. 误差排名图（误差最大的化合物）
            ax = axes[4]
            top_n = min(15, len(valid_df))
            top_errors = valid_df.nlargest(top_n, '绝对误差(min)')

            y_pos = np.arange(top_n)
            ax.barh(y_pos, top_errors['绝对误差(min)'], color='tomato', alpha=0.7)

            # 添加化合物名称
            compound_names = []
            for name in top_errors['化合物名称']:
                if len(name) > 20:
                    compound_names.append(name[:17] + '...')
                else:
                    compound_names.append(name)

            ax.set_yticks(y_pos)
            ax.set_yticklabels(compound_names)
            ax.invert_yaxis()  # 最大的在顶部
            ax.set_xlabel('绝对误差 (min)')
            ax.set_title(f'{from_condition}→{to_condition}: 误差最大的{top_n}个化合物')
            ax.grid(True, alpha=0.3, axis='x')

            # 6. Bland-Altman图（一致性分析）
            ax = axes[5]
            mean_values = (predicted + actual) / 2
            differences = predicted - actual

            ax.scatter(mean_values, differences, alpha=0.6, s=50, color='orange')
            ax.axhline(y=differences.mean(), color='red', linestyle='-',
                       linewidth=2, label=f'均值: {differences.mean():.3f}')

            # 添加95%一致性界限
            mean_diff = differences.mean()
            std_diff = differences.std()
            upper_limit = mean_diff + 1.96 * std_diff
            lower_limit = mean_diff - 1.96 * std_diff

            ax.axhline(y=upper_limit, color='red', linestyle='--',
                       linewidth=1, label=f'+1.96SD: {upper_limit:.3f}')
            ax.axhline(y=lower_limit, color='red', linestyle='--',
                       linewidth=1, label=f'-1.96SD: {lower_limit:.3f}')

            ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.5, alpha=0.5)

            ax.set_xlabel('平均保留时间 (min)')
            ax.set_ylabel('预测 - 实际 (min)')
            ax.set_title(f'{from_condition}→{to_condition}: Bland-Altman图')
            ax.legend(fontsize=9)
            ax.grid(True, alpha=0.3)

            plt.tight_layout()

            if save_path:
                plt.savefig(save_path, dpi=300, bbox_inches='tight')
                print(f"转换分析图已保存到: {save_path}")

            self.figures['conversion_analysis'] = fig
            return fig

        except Exception as e:
            print(f"绘制转换分析图失败: {str(e)}")
            return None

    def plot_multiple_conversion_comparison(self, conversions: List[Tuple[str, str]],
                                            save_path: str = None) -> plt.Figure:
        """
        绘制多个转换方案的比较图

        参数:
            conversions: 转换方案列表，每个元素为(from_condition, to_condition)
            save_path: 保存路径（可选）

        返回:
            matplotlib Figure对象
        """
        try:
            # 收集所有转换的误差数据
            all_errors = []
            labels = []

            for from_cond, to_cond in conversions:
                conversion_df, _ = self.calculator.convert_ppg_index_to_rt(from_cond, to_cond)

                if not conversion_df.empty and '绝对误差(min)' in conversion_df.columns:
                    errors = conversion_df['绝对误差(min)'].dropna()
                    if len(errors) > 0:
                        all_errors.append(errors)
                        labels.append(f'{from_cond}→{to_cond}')

            if not all_errors:
                print("警告: 没有有效的转换数据用于比较")
                return None

            # 创建比较图
            fig, axes = plt.subplots(2, 2, figsize=(14, 12))
            axes = axes.flatten()

            # 1. 误差分布箱线图比较
            ax = axes[0]
            bp = ax.boxplot(all_errors, labels=labels, patch_artist=True, showfliers=False)

            # 设置颜色
            colors = ['lightblue', 'lightgreen', 'lightcoral', 'lightsalmon', 'lightyellow']
            for patch, color in zip(bp['boxes'], colors[:len(all_errors)]):
                patch.set_facecolor(color)

            ax.set_ylabel('绝对误差 (min)')
            ax.set_title('不同转换方案的误差分布比较')
            ax.grid(True, alpha=0.3)

            # 添加统计数据
            for i, errors in enumerate(all_errors):
                ax.text(i + 1, errors.max() + 0.05, f'n={len(errors)}',
                        ha='center', va='bottom', fontsize=9)

            # 2. 平均误差柱状图
            ax = axes[1]
            means = [err.mean() for err in all_errors]
            stds = [err.std() for err in all_errors]

            x_pos = np.arange(len(labels))
            bars = ax.bar(x_pos, means, yerr=stds, capsize=5, alpha=0.7, color='steelblue')

            # 添加数值标签
            for bar, mean in zip(bars, means):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width() / 2., height + 0.01,
                        f'{mean:.3f}', ha='center', va='bottom', fontsize=9)

            ax.set_xlabel('转换方案')
            ax.set_ylabel('平均绝对误差 (min)')
            ax.set_title('不同转换方案的平均误差比较')
            ax.set_xticks(x_pos)
            ax.set_xticklabels(labels, rotation=45, ha='right')
            ax.grid(True, alpha=0.3, axis='y')

            # 3. 通过率比较（基于0.5min阈值）
            ax = axes[2]
            pass_rates = []
            for errors in all_errors:
                passed = len(errors[errors <= 0.5])
                pass_rate = (passed / len(errors) * 100) if len(errors) > 0 else 0
                pass_rates.append(pass_rate)

            bars = ax.bar(x_pos, pass_rates, alpha=0.7, color='lightgreen')

            # 添加数值标签
            for bar, rate in zip(bars, pass_rates):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width() / 2., height + 1,
                        f'{rate:.1f}%', ha='center', va='bottom', fontsize=9)

            ax.set_xlabel('转换方案')
            ax.set_ylabel('通过率 (%)')
            ax.set_title('不同转换方案的通过率比较（阈值: 0.5min）')
            ax.set_xticks(x_pos)
            ax.set_xticklabels(labels, rotation=45, ha='right')
            ax.set_ylim([0, 105])
            ax.grid(True, alpha=0.3, axis='y')

            # 4. 误差累积分布图
            ax = axes[3]
            for errors, label in zip(all_errors, labels):
                sorted_errors = np.sort(errors)
                y_vals = np.arange(1, len(sorted_errors) + 1) / len(sorted_errors) * 100

                ax.plot(sorted_errors, y_vals, marker='.', label=label, linewidth=2)

            ax.set_xlabel('绝对误差 (min)')
            ax.set_ylabel('累积百分比 (%)')
            ax.set_title('不同转换方案的误差累积分布')
            ax.legend()
            ax.grid(True, alpha=0.3)

            plt.tight_layout()

            if save_path:
                plt.savefig(save_path, dpi=300, bbox_inches='tight')
                print(f"多转换方案比较图已保存到: {save_path}")

            self.figures['multiple_conversion_comparison'] = fig
            return fig

        except Exception as e:
            print(f"绘制多转换方案比较图失败: {str(e)}")
            return None


class PPGIndexAnalyzerGUI:
    """PPG保留指数分析GUI"""

    def __init__(self, root):
        """初始化GUI"""
        self.root = root
        self.root.title("PPG保留指数计算与可视化分析系统 - 增强版")
        self.root.geometry("1300x950")

        # 设置图标
        try:
            self.root.iconbitmap(default='icon.ico')
        except:
            pass

        # 处理器
        self.calculator = PPGIndexCalculator()
        self.visualizer = None
        self.processing_thread = None
        self.is_processing = False

        # 存储加载的条件信息
        self.loaded_conditions = set()
        self.loaded_compound_datasets = set()

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
        self.setup_conversion_tab()
        self.setup_visualization_tab()
        self.setup_results_tab()
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
        self.condition_var = tk.StringVar(value="condition1")
        condition_entry = ttk.Entry(ppg_frame, textvariable=self.condition_var, width=20)
        condition_entry.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Label(ppg_frame, text="PPG数据文件:").grid(row=1, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.ppg_file_var = tk.StringVar()
        ppg_file_entry = ttk.Entry(ppg_frame, textvariable=self.ppg_file_var, width=60)
        ppg_file_entry.grid(row=1, column=1, sticky=(tk.W, tk.E), padx=(0, 10), pady=5)

        ttk.Button(ppg_frame, text="浏览...", command=self.browse_ppg_file).grid(row=1, column=2, pady=5)

        ttk.Button(ppg_frame, text="加载PPG数据", command=self.load_ppg_data).grid(row=2, column=0, columnspan=3,
                                                                                   pady=10)

        # ==================== 化合物数据加载 ====================
        compound_frame = ttk.LabelFrame(data_frame, text="化合物数据加载", padding=10)
        compound_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(compound_frame, text="数据类别:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.category_var = tk.StringVar(value="validation")
        category_combo = ttk.Combobox(compound_frame, textvariable=self.category_var,
                                      values=["validation", "smrt", "training", "test"],
                                      width=15, state="readonly")
        category_combo.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Label(compound_frame, text="条件名称:").grid(row=0, column=2, sticky=tk.W, padx=(0, 5), pady=5)
        self.compound_condition_var = tk.StringVar(value="condition1")
        compound_condition_entry = ttk.Entry(compound_frame, textvariable=self.compound_condition_var, width=20)
        compound_condition_entry.grid(row=0, column=3, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Label(compound_frame, text="化合物数据文件:").grid(row=1, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.compound_file_var = tk.StringVar()
        compound_file_entry = ttk.Entry(compound_frame, textvariable=self.compound_file_var, width=60)
        compound_file_entry.grid(row=1, column=1, columnspan=3, sticky=(tk.W, tk.E), padx=(0, 10), pady=5)

        ttk.Button(compound_frame, text="浏览...", command=self.browse_compound_file).grid(row=1, column=4, pady=5)

        ttk.Button(compound_frame, text="加载化合物数据", command=self.load_compound_data).grid(row=2, column=0,
                                                                                                columnspan=5, pady=10)

        # ==================== 加载的数据概览 ====================
        overview_frame = ttk.LabelFrame(data_frame, text="已加载数据概览", padding=10)
        overview_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 10))

        # 创建Treeview显示加载的数据
        columns = ("类型", "条件", "数据类别", "数据点数", "状态")
        self.data_tree = ttk.Treeview(overview_frame, columns=columns, show="headings", height=8)

        # 设置列标题
        for col in columns:
            self.data_tree.heading(col, text=col)
            self.data_tree.column(col, width=100)

        # 添加滚动条
        scrollbar = ttk.Scrollbar(overview_frame, orient=tk.VERTICAL, command=self.data_tree.yview)
        self.data_tree.configure(yscrollcommand=scrollbar.set)

        self.data_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # 数据管理按钮
        button_frame = ttk.Frame(data_frame)
        button_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Button(button_frame, text="清空所有数据", command=self.clear_all_data).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="预览PPG数据", command=self.preview_ppg_data).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="预览化合物数据", command=self.preview_compound_data).pack(side=tk.LEFT, padx=5)

        # 底部信息
        self.data_status_var = tk.StringVar(value="等待加载数据...")
        ttk.Label(data_frame, textvariable=self.data_status_var).pack(anchor=tk.W)

    def setup_analysis_tab(self):
        """设置分析标签页"""
        analysis_tab = ttk.Frame(self.notebook)
        self.notebook.add(analysis_tab, text="数据分析")

        # 主框架
        analysis_frame = ttk.LabelFrame(analysis_tab, text="PPG保留指数分析", padding=15)
        analysis_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # ==================== 标准曲线拟合 ====================
        curve_frame = ttk.LabelFrame(analysis_frame, text="PPG标准曲线拟合", padding=10)
        curve_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(curve_frame, text="选择条件:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.curve_condition_var = tk.StringVar()
        self.curve_condition_combo = ttk.Combobox(curve_frame, textvariable=self.curve_condition_var,
                                                  width=25, state="readonly")
        self.curve_condition_combo.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Label(curve_frame, text="模型类型:").grid(row=0, column=2, sticky=tk.W, padx=(0, 5), pady=5)
        self.model_type_var = tk.StringVar(value="logarithmic")
        model_combo = ttk.Combobox(curve_frame, textvariable=self.model_type_var,
                                   values=["logarithmic", "linear"],
                                   width=15, state="readonly")
        model_combo.grid(row=0, column=3, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Button(curve_frame, text="拟合标准曲线", command=self.fit_standard_curve).grid(row=0, column=4,
                                                                                           padx=(20, 0), pady=5)

        # ==================== PPG指数计算 ====================
        calc_frame = ttk.LabelFrame(analysis_frame, text="PPG指数计算", padding=10)
        calc_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(calc_frame, text="选择条件:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.calc_condition_var = tk.StringVar()
        self.calc_condition_combo = ttk.Combobox(calc_frame, textvariable=self.calc_condition_var,
                                                 width=25, state="readonly")
        self.calc_condition_combo.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Label(calc_frame, text="计算方法:").grid(row=0, column=2, sticky=tk.W, padx=(0, 5), pady=5)
        self.calc_method_var = tk.StringVar(value="interpolation")
        method_combo = ttk.Combobox(calc_frame, textvariable=self.calc_method_var,
                                    values=["interpolation", "regression"],
                                    width=15, state="readonly")
        method_combo.grid(row=0, column=3, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Button(calc_frame, text="计算PPG指数", command=self.calculate_ppg_index).grid(row=0, column=4, padx=(20, 0),
                                                                                          pady=5)

        # ==================== 条件比较 ====================
        compare_frame = ttk.LabelFrame(analysis_frame, text="条件比较分析", padding=10)
        compare_frame.pack(fill=tk.X, pady=(0, 15))

        # 条件选择
        ttk.Label(compare_frame, text="选择要比较的条件:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)

        # 创建条件选择框架
        self.condition_checkboxes = {}
        self.condition_checkboxes_frame = ttk.Frame(compare_frame)
        self.condition_checkboxes_frame.grid(row=1, column=0, columnspan=5, sticky=tk.W, padx=5, pady=5)

        # 比较按钮
        ttk.Button(compare_frame, text="比较选定条件", command=self.compare_conditions).grid(row=2, column=0,
                                                                                             sticky=tk.W, padx=5,
                                                                                             pady=10)

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

        ttk.Button(button_frame, text="生成分析报告", command=self.generate_report).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="清空结果显示", command=self.clear_analysis_text).pack(side=tk.LEFT, padx=5)

        # 分析状态
        self.analysis_status_var = tk.StringVar(value="等待分析...")
        ttk.Label(analysis_frame, textvariable=self.analysis_status_var).pack(anchor=tk.W)

    def setup_conversion_tab(self):
        """设置跨条件转换标签页"""
        conversion_tab = ttk.Frame(self.notebook)
        self.notebook.add(conversion_tab, text="跨条件转换")

        # 主框架
        conversion_frame = ttk.LabelFrame(conversion_tab, text="跨条件PPG指数转换与分析", padding=15)
        conversion_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # ==================== 转换设置 ====================
        settings_frame = ttk.LabelFrame(conversion_frame, text="转换设置", padding=10)
        settings_frame.pack(fill=tk.X, pady=(0, 15))

        # 源条件和目标条件选择
        ttk.Label(settings_frame, text="源条件 (PPG指数来源):").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.from_condition_var = tk.StringVar()
        self.from_condition_combo = ttk.Combobox(settings_frame, textvariable=self.from_condition_var,
                                                 width=25, state="readonly")
        self.from_condition_combo.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        ttk.Label(settings_frame, text="目标条件 (转换目标):").grid(row=0, column=2, sticky=tk.W, padx=(0, 5), pady=5)
        self.to_condition_var = tk.StringVar()
        self.to_condition_combo = ttk.Combobox(settings_frame, textvariable=self.to_condition_var,
                                               width=25, state="readonly")
        self.to_condition_combo.grid(row=0, column=3, sticky=tk.W, padx=(0, 10), pady=5)

        # 误差阈值
        ttk.Label(settings_frame, text="误差阈值(min):").grid(row=1, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.threshold_var = tk.DoubleVar(value=0.5)
        threshold_entry = ttk.Entry(settings_frame, textvariable=self.threshold_var, width=10)
        threshold_entry.grid(row=1, column=1, sticky=tk.W, padx=(0, 10), pady=5)

        # 转换按钮
        ttk.Button(settings_frame, text="执行转换分析", command=self.perform_conversion_analysis).grid(
            row=0, column=4, rowspan=2, padx=(20, 0), pady=5, sticky=tk.NS)

        # 多转换方案比较
        multi_compare_frame = ttk.LabelFrame(settings_frame, text="多转换方案比较", padding=5)
        multi_compare_frame.grid(row=2, column=0, columnspan=5, sticky=tk.W, padx=5, pady=10)

        ttk.Label(multi_compare_frame, text="选择转换方案:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=2)

        # 创建多转换方案选择框架
        self.multi_conversion_frame = ttk.Frame(multi_compare_frame)
        self.multi_conversion_frame.grid(row=1, column=0, columnspan=3, sticky=tk.W, padx=5, pady=2)

        ttk.Button(multi_compare_frame, text="比较选定方案", command=self.compare_multiple_conversions).grid(
            row=1, column=3, padx=(10, 0), pady=2)

        # ==================== 转换结果显示 ====================
        results_frame = ttk.LabelFrame(conversion_frame, text="转换结果", padding=10)
        results_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 15))

        # 创建Treeview显示结果
        columns = ("化合物名称", "源PPG指数", "预测RT", "实际RT", "绝对误差", "相对误差(%)")
        self.conversion_tree = ttk.Treeview(results_frame, columns=columns, show="headings", height=10)

        # 设置列标题
        for col in columns:
            self.conversion_tree.heading(col, text=col)
            if col == "化合物名称":
                self.conversion_tree.column(col, width=150)
            else:
                self.conversion_tree.column(col, width=100)

        # 添加滚动条
        scrollbar = ttk.Scrollbar(results_frame, orient=tk.VERTICAL, command=self.conversion_tree.yview)
        self.conversion_tree.configure(yscrollcommand=scrollbar.set)

        self.conversion_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # 统计信息显示
        stats_frame = ttk.LabelFrame(conversion_frame, text="转换统计", padding=10)
        stats_frame.pack(fill=tk.X, pady=(0, 10))

        self.conversion_stats_text = scrolledtext.ScrolledText(stats_frame, width=80, height=6,
                                                               wrap=tk.WORD, font=("Consolas", 9))
        self.conversion_stats_text.pack(fill=tk.BOTH, expand=True)

        # 配置标签颜色
        self.conversion_stats_text.tag_config("INFO", foreground="black")
        self.conversion_stats_text.tag_config("SUCCESS", foreground="green")
        self.conversion_stats_text.tag_config("WARNING", foreground="orange")
        self.conversion_stats_text.tag_config("ERROR", foreground="red")

        # 转换控制按钮
        button_frame = ttk.Frame(conversion_frame)
        button_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Button(button_frame, text="导出转换结果", command=self.export_conversion_results).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="清空结果", command=self.clear_conversion_results).pack(side=tk.LEFT, padx=5)

        # 转换状态
        self.conversion_status_var = tk.StringVar(value="等待转换分析...")
        ttk.Label(conversion_frame, textvariable=self.conversion_status_var).pack(anchor=tk.W)

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
        self.viz_type_var = tk.StringVar(value="standard_curve")
        viz_type_combo = ttk.Combobox(options_frame, textvariable=self.viz_type_var,
                                      values=["standard_curve", "residuals", "index_distribution",
                                              "condition_comparison", "conversion_analysis",
                                              "multiple_conversion_comparison"],
                                      width=25, state="readonly")
        viz_type_combo.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)
        viz_type_combo.bind("<<ComboboxSelected>>", self.on_viz_type_change)

        # 条件选择（动态更新）
        ttk.Label(options_frame, text="选择条件:").grid(row=0, column=2, sticky=tk.W, padx=(0, 5), pady=5)
        self.viz_condition_var = tk.StringVar()
        self.viz_condition_combo = ttk.Combobox(options_frame, textvariable=self.viz_condition_var,
                                                width=25, state="readonly")
        self.viz_condition_combo.grid(row=0, column=3, sticky=tk.W, padx=(0, 10), pady=5)

        # 多条件选择（用于条件比较和转换分析）
        self.viz_conditions_frame = ttk.Frame(options_frame)
        self.viz_conditions_frame.grid(row=1, column=0, columnspan=4, sticky=tk.W, padx=5, pady=5)
        self.viz_conditions_frame.grid_remove()  # 默认隐藏

        # 转换分析参数框架
        self.conversion_viz_frame = ttk.Frame(options_frame)
        self.conversion_viz_frame.grid(row=1, column=0, columnspan=4, sticky=tk.W, padx=5, pady=5)
        self.conversion_viz_frame.grid_remove()  # 默认隐藏

        # 生成图表按钮
        ttk.Button(options_frame, text="生成图表", command=self.generate_visualization).grid(row=0, column=4,
                                                                                             padx=(20, 0), pady=5)

        # ==================== 图表显示区域 ====================
        display_frame = ttk.LabelFrame(viz_frame, text="图表显示", padding=10)
        display_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 10))

        # 创建画布框架
        self.figure_canvas = None
        self.figure_toolbar = None
        self.current_figure = None

        # 创建占位标签
        self.viz_placeholder = ttk.Label(display_frame, text="图表将在此处显示",
                                         font=("Arial", 14), foreground="gray")
        self.viz_placeholder.pack(expand=True)

        # 可视化控制按钮
        button_frame = ttk.Frame(viz_frame)
        button_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Button(button_frame, text="保存图表", command=self.save_figure).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="清除图表", command=self.clear_figure).pack(side=tk.LEFT, padx=5)

        # 可视化状态
        self.viz_status_var = tk.StringVar(value="等待生成图表...")
        ttk.Label(viz_frame, textvariable=self.viz_status_var).pack(anchor=tk.W)

    def setup_results_tab(self):
        """设置结果标签页"""
        results_tab = ttk.Frame(self.notebook)
        self.notebook.add(results_tab, text="结果输出")

        # 主框架
        results_frame = ttk.LabelFrame(results_tab, text="结果输出与管理", padding=15)
        results_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # ==================== 输出设置 ====================
        output_frame = ttk.LabelFrame(results_frame, text="输出设置", padding=10)
        output_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(output_frame, text="输出目录:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.output_dir_var = tk.StringVar(value=os.getcwd())
        output_dir_entry = ttk.Entry(output_frame, textvariable=self.output_dir_var, width=60)
        output_dir_entry.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=(0, 10), pady=5)

        ttk.Button(output_frame, text="浏览...", command=self.browse_output_dir).grid(row=0, column=2, pady=5)

        ttk.Label(output_frame, text="输出文件名前缀:").grid(row=1, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.output_prefix_var = tk.StringVar(value="PPG_Analysis")
        ttk.Entry(output_frame, textvariable=self.output_prefix_var, width=30).grid(row=1, column=1, sticky=tk.W,
                                                                                    padx=(0, 10), pady=5)

        # ==================== 结果预览 ====================
        preview_frame = ttk.LabelFrame(results_frame, text="结果预览", padding=10)
        preview_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 15))

        # 创建文本显示区域
        self.results_text = scrolledtext.ScrolledText(preview_frame, width=80, height=15,
                                                      wrap=tk.WORD, font=("Consolas", 10))
        self.results_text.pack(fill=tk.BOTH, expand=True)

        # ==================== 输出控制 ====================
        control_frame = ttk.Frame(results_frame)
        control_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Button(control_frame, text="保存所有结果", command=self.save_all_results).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_frame, text="预览报告", command=self.preview_report).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_frame, text="打开输出目录", command=self.open_output_dir).pack(side=tk.LEFT, padx=5)

        # 输出状态
        self.output_status_var = tk.StringVar(value="等待输出结果...")
        ttk.Label(results_frame, textvariable=self.output_status_var).pack(anchor=tk.W)

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

    def browse_compound_file(self):
        """浏览化合物数据文件"""
        file_types = [("数据文件", "*.csv *.xlsx *.xls"), ("CSV文件", "*.csv"),
                      ("Excel文件", "*.xlsx *.xls"), ("所有文件", "*.*")]
        file_path = filedialog.askopenfilename(title="选择化合物数据文件", filetypes=file_types)

        if file_path:
            self.compound_file_var.set(file_path)

    def browse_output_dir(self):
        """浏览输出目录"""
        dir_path = filedialog.askdirectory(title="选择输出目录")

        if dir_path:
            self.output_dir_var.set(dir_path)

    def load_ppg_data(self):
        """加载PPG数据"""
        ppg_file = self.ppg_file_var.get().strip()
        condition = self.condition_var.get().strip()

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
            success, msg = self.calculator.load_ppg_data(ppg_file, condition)

            if success:
                self.root.after(0, lambda: self.log_message(f"✓ {msg}", "SUCCESS"))
                self.root.after(0, self.update_data_tree)
                self.root.after(0, lambda: self.update_condition_comboboxes())
                self.loaded_conditions.add(condition)
            else:
                self.root.after(0, lambda: self.log_message(f"✗ {msg}", "ERROR"))

            self.root.after(0, lambda: self.update_status("就绪"))

        except Exception as e:
            self.root.after(0, lambda: self.log_message(f"✗ 加载PPG数据失败: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("加载失败"))

    def load_compound_data(self):
        """加载化合物数据"""
        compound_file = self.compound_file_var.get().strip()
        category = self.category_var.get().strip()
        condition = self.compound_condition_var.get().strip()

        if not compound_file:
            messagebox.showwarning("警告", "请选择化合物数据文件")
            return

        if not category:
            messagebox.showwarning("警告", "请选择数据类别")
            return

        if not condition:
            messagebox.showwarning("警告", "请输入条件名称")
            return

        # 显示加载状态
        self.log_message(f"正在加载化合物数据: {compound_file} (类别: {category}, 条件: {condition})", "INFO")
        self.update_status(f"加载化合物数据: {Path(compound_file).name}")

        # 在线程中加载数据
        self.processing_thread = threading.Thread(
            target=self._load_compound_data_thread,
            args=(compound_file, category, condition)
        )
        self.processing_thread.start()

    def _load_compound_data_thread(self, compound_file, category, condition):
        """加载化合物数据的线程函数"""
        try:
            success, msg = self.calculator.load_compound_data(compound_file, category, condition)

            if success:
                self.root.after(0, lambda: self.log_message(f"✓ {msg}", "SUCCESS"))
                self.root.after(0, self.update_data_tree)
                self.root.after(0, lambda: self.update_condition_comboboxes())
                self.loaded_compound_datasets.add(f"{category}_{condition}")
            else:
                self.root.after(0, lambda: self.log_message(f"✗ {msg}", "ERROR"))

            self.root.after(0, lambda: self.update_status("就绪"))

        except Exception as e:
            self.root.after(0, lambda: self.log_message(f"✗ 加载化合物数据失败: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("加载失败"))

    def update_data_tree(self):
        """更新数据树显示"""
        # 清空现有数据
        for item in self.data_tree.get_children():
            self.data_tree.delete(item)

        # 添加PPG数据
        for condition, df in self.calculator.ppg_data.items():
            self.data_tree.insert("", tk.END, values=("PPG标准品", condition, "", len(df), "已加载"))

        # 添加化合物数据
        for key, df in self.calculator.compound_data.items():
            parts = key.split('_')
            if len(parts) >= 2:
                category, condition = parts[0], parts[1]
                self.data_tree.insert("", tk.END, values=("化合物", condition, category, len(df), "已加载"))

        # 更新状态
        self.data_status_var.set(
            f"已加载 {len(self.calculator.ppg_data)} 个PPG数据集, {len(self.calculator.compound_data)} 个化合物数据集")

    def update_condition_comboboxes(self):
        """更新条件选择下拉框"""
        conditions = list(self.calculator.ppg_data.keys())

        # 更新曲线拟合条件
        self.curve_condition_combo['values'] = conditions
        if conditions and not self.curve_condition_var.get():
            self.curve_condition_var.set(conditions[0])

        # 更新计算条件
        self.calc_condition_combo['values'] = conditions
        if conditions and not self.calc_condition_var.get():
            self.calc_condition_var.set(conditions[0])

        # 更新可视化条件
        self.viz_condition_combo['values'] = conditions
        if conditions and not self.viz_condition_var.get():
            self.viz_condition_var.set(conditions[0])

        # 更新转换条件
        self.from_condition_combo['values'] = conditions
        self.to_condition_combo['values'] = conditions
        if conditions:
            if not self.from_condition_var.get():
                self.from_condition_var.set(conditions[0])
            if not self.to_condition_var.get():
                self.to_condition_var.set(conditions[-1] if len(conditions) > 1 else conditions[0])

        # 更新条件比较复选框
        self.update_condition_checkboxes(conditions)

        # 更新多转换方案复选框
        self.update_multi_conversion_checkboxes(conditions)

    def update_condition_checkboxes(self, conditions):
        """更新条件比较复选框"""
        # 清除现有复选框
        for widget in self.condition_checkboxes_frame.winfo_children():
            widget.destroy()
        self.condition_checkboxes.clear()

        # 创建新的复选框
        for i, condition in enumerate(conditions):
            var = tk.BooleanVar(value=(i < 2))  # 默认选择前两个条件
            cb = ttk.Checkbutton(self.condition_checkboxes_frame, text=condition, variable=var)
            cb.grid(row=i // 4, column=i % 4, sticky=tk.W, padx=5, pady=2)
            self.condition_checkboxes[condition] = var

    def update_multi_conversion_checkboxes(self, conditions):
        """更新多转换方案复选框"""
        # 清除现有复选框
        for widget in self.multi_conversion_frame.winfo_children():
            widget.destroy()

        # 创建转换方案复选框
        self.conversion_scheme_vars = []

        if len(conditions) >= 2:
            for i in range(len(conditions)):
                for j in range(len(conditions)):
                    if i != j:
                        from_cond = conditions[i]
                        to_cond = conditions[j]
                        var = tk.BooleanVar(value=(i == 0 and j == 1))  # 默认选择第一个转换方案
                        cb = ttk.Checkbutton(self.multi_conversion_frame,
                                             text=f"{from_cond}→{to_cond}",
                                             variable=var)
                        row = (i * len(conditions) + j) // 4
                        col = (i * len(conditions) + j) % 4
                        cb.grid(row=row, column=col, sticky=tk.W, padx=2, pady=1)
                        self.conversion_scheme_vars.append((var, from_cond, to_cond))
        else:
            ttk.Label(self.multi_conversion_frame, text="需要至少2个条件才能比较",
                      foreground="gray").grid(row=0, column=0, sticky=tk.W)

    def clear_all_data(self):
        """清空所有数据"""
        if messagebox.askyesno("确认", "确定要清空所有已加载的数据吗？"):
            self.calculator = PPGIndexCalculator()
            self.visualizer = None
            self.loaded_conditions.clear()
            self.loaded_compound_datasets.clear()
            self.update_data_tree()
            self.update_condition_comboboxes()
            self.clear_conversion_results()
            self.log_message("已清空所有数据", "INFO")

    def preview_ppg_data(self):
        """预览PPG数据"""
        conditions = list(self.calculator.ppg_data.keys())
        if not conditions:
            messagebox.showinfo("信息", "没有已加载的PPG数据")
            return

        # 创建预览窗口
        preview_window = Toplevel(self.root)
        preview_window.title("PPG数据预览")
        preview_window.geometry("800x600")

        # 创建Notebook
        notebook = ttk.Notebook(preview_window)
        notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # 为每个条件创建标签页
        for condition in conditions:
            df = self.calculator.ppg_data[condition]

            # 创建框架
            frame = ttk.Frame(notebook)
            notebook.add(frame, text=condition)

            # 创建Treeview
            tree_frame = ttk.Frame(frame)
            tree_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

            # 滚动条
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

            # 定义列
            tree["columns"] = list(df.columns)
            tree["show"] = "headings"

            # 设置列标题
            for col in df.columns:
                tree.heading(col, text=col)
                tree.column(col, width=100)

            # 添加数据
            for _, row in df.iterrows():
                tree.insert("", tk.END, values=list(row))

            # 统计信息
            info_frame = ttk.Frame(frame)
            info_frame.pack(fill=tk.X, padx=5, pady=(0, 5))

            ttk.Label(info_frame, text=f"数据点数: {len(df)}").pack(side=tk.LEFT, padx=10)

            if '聚合度' in df.columns and '保留时间' in df.columns:
                n_range = f"聚合度范围: {df['聚合度'].min()} - {df['聚合度'].max()}"
                rt_range = f"保留时间范围: {df['保留时间'].min():.2f} - {df['保留时间'].max():.2f}"
                ttk.Label(info_frame, text=n_range).pack(side=tk.LEFT, padx=10)
                ttk.Label(info_frame, text=rt_range).pack(side=tk.LEFT, padx=10)

    def preview_compound_data(self):
        """预览化合物数据"""
        if not self.calculator.compound_data:
            messagebox.showinfo("信息", "没有已加载的化合物数据")
            return

        # 创建预览窗口
        preview_window = Toplevel(self.root)
        preview_window.title("化合物数据预览")
        preview_window.geometry("900x700")

        # 创建Notebook
        notebook = ttk.Notebook(preview_window)
        notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # 为每个数据集创建标签页
        for key, df in self.calculator.compound_data.items():
            # 创建框架
            frame = ttk.Frame(notebook)
            notebook.add(frame, text=key)

            # 创建Treeview
            tree_frame = ttk.Frame(frame)
            tree_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

            # 滚动条
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

            # 定义列
            tree["columns"] = list(df.columns)
            tree["show"] = "headings"

            # 设置列标题
            for col in df.columns:
                tree.heading(col, text=col)
                tree.column(col, width=120 if len(col) > 10 else 100)

            # 添加数据
            for _, row in df.iterrows():
                tree.insert("", tk.END, values=list(row))

            # 统计信息
            info_frame = ttk.Frame(frame)
            info_frame.pack(fill=tk.X, padx=5, pady=(0, 5))

            ttk.Label(info_frame, text=f"化合物数量: {len(df)}").pack(side=tk.LEFT, padx=10)

            if '保留时间' in df.columns:
                rt_range = f"保留时间范围: {df['保留时间'].min():.2f} - {df['保留时间'].max():.2f}"
                ttk.Label(info_frame, text=rt_range).pack(side=tk.LEFT, padx=10)

    def fit_standard_curve(self):
        """拟合标准曲线"""
        condition = self.curve_condition_var.get()
        model_type = self.model_type_var.get()

        if not condition:
            messagebox.showwarning("警告", "请选择条件")
            return

        # 显示分析状态
        self.analysis_message(f"正在拟合标准曲线 (条件: {condition}, 模型: {model_type})", "INFO")
        self.update_status(f"拟合标准曲线: {condition}")

        # 在线程中拟合曲线
        self.processing_thread = threading.Thread(
            target=self._fit_standard_curve_thread,
            args=(condition, model_type)
        )
        self.processing_thread.start()

    def _fit_standard_curve_thread(self, condition, model_type):
        """拟合标准曲线的线程函数"""
        try:
            success, msg = self.calculator.fit_standard_curve(condition, model_type)

            if success:
                self.root.after(0, lambda: self.analysis_message(f"✓ {msg}", "SUCCESS"))

                # 显示拟合结果
                curve = self.calculator.standard_curves[condition]
                result_text = f"标准曲线拟合结果 - {condition}:\n"
                result_text += f"  模型类型: {curve['model_name']}\n"
                result_text += f"  R²: {curve['r_squared']:.6f}\n"
                result_text += f"  斜率: {curve['slope']:.6f}\n"
                result_text += f"  截距: {curve['intercept']:.6f}\n"
                result_text += f"  标准误差: {curve['std_err']:.6f}\n"
                result_text += f"  p值: {curve['p_value']:.6f}\n"
                result_text += f"  数据点数: {curve['n_points']}\n"

                self.root.after(0, lambda: self.analysis_message(result_text, "INFO"))

                # 更新可视化器
                self.visualizer = PPGVisualizer(self.calculator)
            else:
                self.root.after(0, lambda: self.analysis_message(f"✗ {msg}", "ERROR"))

            self.root.after(0, lambda: self.update_status("就绪"))

        except Exception as e:
            self.root.after(0, lambda: self.analysis_message(f"✗ 拟合标准曲线失败: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("拟合失败"))

    def calculate_ppg_index(self):
        """计算PPG指数"""
        condition = self.calc_condition_var.get()
        method = self.calc_method_var.get()

        if not condition:
            messagebox.showwarning("警告", "请选择条件")
            return

        # 显示分析状态
        self.analysis_message(f"正在计算PPG指数 (条件: {condition}, 方法: {method})", "INFO")
        self.update_status(f"计算PPG指数: {condition}")

        # 在线程中计算
        self.processing_thread = threading.Thread(
            target=self._calculate_ppg_index_thread,
            args=(condition, method)
        )
        self.processing_thread.start()

    def _calculate_ppg_index_thread(self, condition, method):
        """计算PPG指数的线程函数"""
        try:
            success, msg = self.calculator.calculate_ppg_index(condition, method)

            if success:
                self.root.after(0, lambda: self.analysis_message(f"✓ {msg}", "SUCCESS"))

                # 显示统计结果
                if condition in self.calculator.ppg_indices:
                    indices_data = self.calculator.ppg_indices[condition]
                    all_indices = []

                    for key, df in indices_data['indices'].items():
                        if '计算PPG指数' in df.columns:
                            indices = df['计算PPG指数'].dropna().tolist()
                            all_indices.extend(indices)

                    if all_indices:
                        indices_array = np.array(all_indices)
                        result_text = f"PPG指数统计结果 - {condition}:\n"
                        result_text += f"  计算方法: {method}\n"
                        result_text += f"  总化合物数: {len(all_indices)}\n"
                        result_text += f"  平均值: {np.mean(indices_array):.2f}\n"
                        result_text += f"  标准差: {np.std(indices_array):.2f}\n"
                        result_text += f"  最小值: {np.min(indices_array):.2f}\n"
                        result_text += f"  最大值: {np.max(indices_array):.2f}\n"
                        result_text += f"  中位数: {np.median(indices_array):.2f}\n"

                        self.root.after(0, lambda: self.analysis_message(result_text, "INFO"))
            else:
                self.root.after(0, lambda: self.analysis_message(f"✗ {msg}", "ERROR"))

            self.root.after(0, lambda: self.update_status("就绪"))

        except Exception as e:
            self.root.after(0, lambda: self.analysis_message(f"✗ 计算PPG指数失败: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("计算失败"))

    def compare_conditions(self):
        """比较不同条件"""
        # 获取选定的条件
        selected_conditions = []
        for condition, var in self.condition_checkboxes.items():
            if var.get():
                selected_conditions.append(condition)

        if len(selected_conditions) < 2:
            messagebox.showwarning("警告", "请至少选择两个条件进行比较")
            return

        # 显示分析状态
        self.analysis_message(f"正在比较条件: {', '.join(selected_conditions)}", "INFO")
        self.update_status(f"比较条件")

        # 在线程中比较
        self.processing_thread = threading.Thread(
            target=self._compare_conditions_thread,
            args=(selected_conditions,)
        )
        self.processing_thread.start()

    def _compare_conditions_thread(self, conditions):
        """比较条件的线程函数"""
        try:
            comparison_df = self.calculator.compare_conditions(conditions)

            if not comparison_df.empty:
                self.root.after(0, lambda: self.analysis_message(f"✓ 条件比较完成，共比较 {len(comparison_df)} 个化合物",
                                                                 "SUCCESS"))

                # 显示比较结果
                result_text = f"条件比较结果:\n"
                result_text += f"  比较条件: {', '.join(conditions)}\n"
                result_text += f"  总化合物数: {len(comparison_df)}\n\n"

                # 计算每个条件下的PPG指数统计
                for condition in conditions:
                    col_name = f"{condition}_PPG指数"
                    if col_name in comparison_df.columns:
                        data = comparison_df[col_name].dropna()
                        if len(data) > 0:
                            result_text += f"  {condition}:\n"
                            result_text += f"    有效数据数: {len(data)}\n"
                            result_text += f"    平均值: {np.mean(data):.2f}\n"
                            result_text += f"    标准差: {np.std(data):.2f}\n"

                # 计算条件间的相关系数
                if len(conditions) >= 2:
                    result_text += f"\n  条件间相关系数:\n"

                    # 提取数据
                    corr_data = {}
                    for condition in conditions:
                        col_name = f"{condition}_PPG指数"
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
                self.root.after(0, lambda: self.analysis_message("✗ 没有可用于比较的数据", "WARNING"))

            self.root.after(0, lambda: self.update_status("就绪"))

        except Exception as e:
            self.root.after(0, lambda: self.analysis_message(f"✗ 条件比较失败: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("比较失败"))

    def perform_conversion_analysis(self):
        """执行跨条件转换分析"""
        from_condition = self.from_condition_var.get()
        to_condition = self.to_condition_var.get()
        threshold = self.threshold_var.get()

        if not from_condition or not to_condition:
            messagebox.showwarning("警告", "请选择源条件和目标条件")
            return

        if from_condition == to_condition:
            messagebox.showwarning("警告", "源条件和目标条件不能相同")
            return

        # 显示转换状态
        self.conversion_status_var.set(f"正在执行转换: {from_condition} → {to_condition}")
        self.update_status(f"执行转换分析")

        # 在线程中执行转换
        self.processing_thread = threading.Thread(
            target=self._perform_conversion_analysis_thread,
            args=(from_condition, to_condition, threshold)
        )
        self.processing_thread.start()

    def _perform_conversion_analysis_thread(self, from_condition, to_condition, threshold):
        """执行转换分析的线程函数"""
        try:
            # 执行转换分析
            analysis_results = self.calculator.cross_condition_analysis(from_condition, to_condition, threshold)

            if "error" in analysis_results:
                self.root.after(0, lambda: self.conversion_status_var.set("转换失败"))
                self.root.after(0, lambda: self.log_message(f"✗ 转换分析失败: {analysis_results['error']}", "ERROR"))
                return

            # 显示转换结果
            self.root.after(0, lambda: self.display_conversion_results(analysis_results))

            # 更新日志
            self.root.after(0, lambda: self.log_message(
                f"✓ 转换分析完成: {from_condition} → {to_condition}, " +
                f"有效转换数: {analysis_results['有效转换数']}, " +
                f"平均误差: {analysis_results['误差统计'].get('平均绝对误差', 0):.3f}min",
                "SUCCESS"))

            self.root.after(0, lambda: self.conversion_status_var.set("转换分析完成"))
            self.root.after(0, lambda: self.update_status("就绪"))

        except Exception as e:
            self.root.after(0, lambda: self.conversion_status_var.set("转换失败"))
            self.root.after(0, lambda: self.log_message(f"✗ 转换分析失败: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("转换失败"))

    def display_conversion_results(self, analysis_results):
        """显示转换结果"""
        # 清空现有结果
        for item in self.conversion_tree.get_children():
            self.conversion_tree.delete(item)

        # 清空统计文本
        self.conversion_stats_text.delete(1.0, tk.END)

        # 显示详细数据
        if '详细数据' in analysis_results:
            for item in analysis_results['详细数据']:
                compound_name = item['化合物名称']
                source_ppg = item.get(f"{analysis_results['源条件']}_PPG指数", "")
                pred_rt = item.get(f"{analysis_results['目标条件']}_预测RT", "")
                actual_rt = item.get(f"{analysis_results['目标条件']}_实际RT", "")
                abs_error = item.get('绝对误差(min)', "")
                rel_error = item.get('相对误差(%)', "")

                # 格式化显示
                values = (
                    compound_name,
                    f"{source_ppg:.2f}" if isinstance(source_ppg, (int, float)) else source_ppg,
                    f"{pred_rt:.3f}" if isinstance(pred_rt, (int, float)) else pred_rt,
                    f"{actual_rt:.3f}" if isinstance(actual_rt, (int, float)) else actual_rt,
                    f"{abs_error:.3f}" if isinstance(abs_error, (int, float)) else abs_error,
                    f"{rel_error:.2f}" if isinstance(rel_error, (int, float)) else rel_error
                )

                self.conversion_tree.insert("", tk.END, values=values)

        # 显示统计信息
        stats_text = f"转换分析结果: {analysis_results['源条件']} → {analysis_results['目标条件']}\n"
        stats_text += "=" * 50 + "\n\n"

        stats_text += f"转换时间: {analysis_results['转换时间']}\n"
        stats_text += f"总化合物数: {analysis_results['总化合物数']}\n"
        stats_text += f"有效转换数: {analysis_results['有效转换数']}\n\n"

        if isinstance(analysis_results['误差统计'], dict):
            stats_text += "误差统计:\n"
            for key, value in analysis_results['误差统计'].items():
                if isinstance(value, float):
                    if '误差' in key or '标准差' in key:
                        stats_text += f"  {key}: {value:.3f}\n"
                    else:
                        stats_text += f"  {key}: {value:.2f}\n"
                else:
                    stats_text += f"  {key}: {value}\n"

        if '误差分布' in analysis_results:
            stats_text += "\n误差分布:\n"
            for bin_name, count in analysis_results['误差分布'].items():
                percentage = (count / analysis_results['有效转换数'] * 100) if analysis_results['有效转换数'] > 0 else 0
                stats_text += f"  {bin_name}: {count} ({percentage:.1f}%)\n"

        if '通过率分析' in analysis_results:
            pass_analysis = analysis_results['通过率分析']
            stats_text += f"\n通过率分析 (阈值: {pass_analysis['阈值(min)']}min):\n"
            stats_text += f"  通过数: {pass_analysis['通过数']}\n"
            stats_text += f"  总数: {pass_analysis['总数']}\n"
            stats_text += f"  通过率: {pass_analysis['通过率(%)']:.1f}%\n"

        self.conversion_stats_text.insert(tk.END, stats_text, "INFO")

    def compare_multiple_conversions(self):
        """比较多个转换方案"""
        # 获取选定的转换方案
        selected_conversions = []
        for var, from_cond, to_cond in self.conversion_scheme_vars:
            if var.get():
                selected_conversions.append((from_cond, to_cond))

        if len(selected_conversions) < 2:
            messagebox.showwarning("警告", "请至少选择两个转换方案进行比较")
            return

        # 显示分析状态
        self.log_message(f"正在比较转换方案: {len(selected_conversions)} 个方案", "INFO")
        self.update_status(f"比较转换方案")

        # 在线程中执行比较
        self.processing_thread = threading.Thread(
            target=self._compare_multiple_conversions_thread,
            args=(selected_conversions,)
        )
        self.processing_thread.start()

    def _compare_multiple_conversions_thread(self, conversions):
        """比较多个转换方案的线程函数"""
        try:
            # 确保可视化器已初始化
            if self.visualizer is None:
                self.visualizer = PPGVisualizer(self.calculator)

            # 生成比较图表
            fig = self.visualizer.plot_multiple_conversion_comparison(conversions)

            if fig is not None:
                # 在主线程中显示图表
                self.root.after(0, lambda: self.display_comparison_figure(fig, conversions))

                self.root.after(0, lambda: self.log_message(
                    f"✓ 转换方案比较完成，共比较 {len(conversions)} 个方案", "SUCCESS"))
            else:
                self.root.after(0, lambda: self.log_message("✗ 无法生成比较图表", "ERROR"))

            self.root.after(0, lambda: self.update_status("就绪"))

        except Exception as e:
            self.root.after(0, lambda: self.log_message(f"✗ 比较转换方案失败: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("比较失败"))

    def display_comparison_figure(self, fig, conversions):
        """显示比较图表"""
        # 清除现有图表
        self.clear_figure()

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

        # 更新状态
        self.viz_status_var.set(f"多转换方案比较图生成完成")

    def export_conversion_results(self):
        """导出转换结果"""
        if not self.calculator.conversion_results:
            messagebox.showinfo("信息", "没有可导出的转换结果")
            return

        file_path = filedialog.asksaveasfilename(
            title="导出转换结果",
            defaultextension=".xlsx",
            filetypes=[("Excel文件", "*.xlsx"), ("所有文件", "*.*")]
        )

        if file_path:
            try:
                with pd.ExcelWriter(file_path, engine='openpyxl') as writer:
                    # 保存每个转换结果
                    for key, conversion in self.calculator.conversion_results.items():
                        if '详细数据' in conversion:
                            df = pd.DataFrame(conversion['详细数据'])
                            sheet_name = key[:30]  # Excel sheet名称长度限制
                            df.to_excel(writer, sheet_name=sheet_name, index=False)

                    # 保存转换统计
                    conversion_stats = []
                    for key, conversion in self.calculator.conversion_results.items():
                        stats = {
                            '转换方向': key,
                            '源条件': conversion.get('源条件', ''),
                            '目标条件': conversion.get('目标条件', ''),
                            '总化合物数': conversion.get('总化合物数', 0),
                            '有效转换数': conversion.get('有效转换数', 0)
                        }

                        if isinstance(conversion.get('误差统计'), dict):
                            stats.update({
                                '平均绝对误差(min)': conversion['误差统计'].get('平均绝对误差', 0),
                                '绝对误差标准差': conversion['误差统计'].get('绝对误差标准差', 0),
                                '最大绝对误差(min)': conversion['误差统计'].get('最大绝对误差', 0)
                            })

                        if isinstance(conversion.get('通过率分析'), dict):
                            stats.update({
                                '通过率(%)': conversion['通过率分析'].get('通过率(%)', 0),
                                '通过数': conversion['通过率分析'].get('通过数', 0),
                                '阈值(min)': conversion['通过率分析'].get('阈值(min)', 0.5)
                            })

                        conversion_stats.append(stats)

                    if conversion_stats:
                        stats_df = pd.DataFrame(conversion_stats)
                        stats_df.to_excel(writer, sheet_name='转换统计汇总', index=False)

                self.log_message(f"✓ 转换结果已导出到: {file_path}", "SUCCESS")
                self.conversion_status_var.set(f"结果已导出: {Path(file_path).name}")

            except Exception as e:
                self.log_message(f"✗ 导出转换结果失败: {str(e)}", "ERROR")

    def clear_conversion_results(self):
        """清空转换结果"""
        # 清空Treeview
        for item in self.conversion_tree.get_children():
            self.conversion_tree.delete(item)

        # 清空统计文本
        self.conversion_stats_text.delete(1.0, tk.END)

        # 清空计算器中的转换结果
        self.calculator.conversion_results.clear()

        self.conversion_status_var.set("转换结果已清空")
        self.log_message("已清空转换结果", "INFO")

    def generate_report(self):
        """生成分析报告"""
        try:
            summary = self.calculator.generate_summary_report()

            # 在分析结果中显示报告
            self.analysis_message("=" * 60, "INFO")
            self.analysis_message("PPG保留指数分析报告", "INFO")
            self.analysis_message("=" * 60, "INFO")
            self.analysis_message(f"生成时间: {summary['生成时间']}", "INFO")
            self.analysis_message("", "INFO")

            self.analysis_message("数据概览:", "INFO")
            self.analysis_message(f"  - PPG数据条件数: {summary['PPG数据条件数']}", "INFO")
            self.analysis_message(f"  - 化合物数据集数: {summary['化合物数据集数']}", "INFO")
            self.analysis_message(f"  - 标准曲线数: {summary['标准曲线数']}", "INFO")
            self.analysis_message(f"  - PPG指数计算结果数: {summary['PPG指数计算结果数']}", "INFO")
            self.analysis_message(f"  - 跨条件转换结果数: {summary['跨条件转换结果数']}", "INFO")
            self.analysis_message("", "INFO")

            if summary['标准曲线性能']:
                self.analysis_message("标准曲线性能:", "INFO")
                for condition, perf in summary['标准曲线性能'].items():
                    self.analysis_message(f"  {condition}:", "INFO")
                    self.analysis_message(f"    - 模型类型: {perf['模型类型']}", "INFO")
                    self.analysis_message(f"    - R²: {perf['R²']:.4f}", "INFO")
                    self.analysis_message(f"    - 斜率: {perf['斜率']:.4f}", "INFO")
                    self.analysis_message(f"    - 截距: {perf['截距']:.4f}", "INFO")
                    self.analysis_message(f"    - 标准误差: {perf['标准误差']:.4f}", "INFO")
                self.analysis_message("", "INFO")

            if summary['PPG指数统计']:
                self.analysis_message("PPG指数统计:", "INFO")
                for condition, stats in summary['PPG指数统计'].items():
                    self.analysis_message(f"  {condition}:", "INFO")
                    self.analysis_message(f"    - 计算方法: {stats['计算方法']}", "INFO")
                    self.analysis_message(f"    - 样本数: {stats['样本数']}", "INFO")
                    self.analysis_message(f"    - 平均值: {stats['平均值']:.2f}", "INFO")
                    self.analysis_message(f"    - 标准差: {stats['标准差']:.2f}", "INFO")
                    self.analysis_message(f"    - 范围: {stats['最小值']:.2f} - {stats['最大值']:.2f}", "INFO")
                self.analysis_message("", "INFO")

            if summary['跨条件转换统计']:
                self.analysis_message("跨条件转换统计:", "INFO")
                for key, stats in summary['跨条件转换统计'].items():
                    self.analysis_message(f"  {key}:", "INFO")
                    self.analysis_message(f"    - 源条件: {stats['源条件']}", "INFO")
                    self.analysis_message(f"    - 目标条件: {stats['目标条件']}", "INFO")
                    self.analysis_message(f"    - 总化合物数: {stats['总化合物数']}", "INFO")
                    self.analysis_message(f"    - 有效转换数: {stats['有效转换数']}", "INFO")
                    self.analysis_message(f"    - 平均绝对误差: {stats['平均绝对误差']:.3f} min", "INFO")
                    self.analysis_message(f"    - 通过率: {stats['通过率(%)']:.1f}%", "INFO")
                self.analysis_message("", "INFO")

            self.analysis_message("报告生成完成", "SUCCESS")

            # 在结果标签页显示报告
            self.results_text.delete(1.0, tk.END)
            self.results_text.insert(tk.END, "PPG保留指数分析报告\n")
            self.results_text.insert(tk.END, "=" * 60 + "\n\n")
            self.results_text.insert(tk.END, f"生成时间: {summary['生成时间']}\n\n")

            self.results_text.insert(tk.END, "数据概览:\n")
            self.results_text.insert(tk.END, f"  - PPG数据条件数: {summary['PPG数据条件数']}\n")
            self.results_text.insert(tk.END, f"  - 化合物数据集数: {summary['化合物数据集数']}\n")
            self.results_text.insert(tk.END, f"  - 标准曲线数: {summary['标准曲线数']}\n")
            self.results_text.insert(tk.END, f"  - PPG指数计算结果数: {summary['PPG指数计算结果数']}\n")
            self.results_text.insert(tk.END, f"  - 跨条件转换结果数: {summary['跨条件转换结果数']}\n\n")

            self.analysis_status_var.set("报告生成完成")

        except Exception as e:
            self.analysis_message(f"✗ 生成报告失败: {str(e)}", "ERROR")

    def on_viz_type_change(self, event=None):
        """可视化类型改变事件"""
        viz_type = self.viz_type_var.get()

        # 隐藏所有参数框架
        self.viz_conditions_frame.grid_remove()
        self.conversion_viz_frame.grid_remove()
        self.viz_condition_combo.grid()

        if viz_type in ["condition_comparison", "multiple_conversion_comparison"]:
            # 显示条件选择框架
            self.viz_conditions_frame.grid()
            self.viz_condition_combo.grid_remove()

            # 更新条件复选框
            if viz_type == "condition_comparison":
                self.update_viz_condition_checkboxes()
            else:
                self.update_viz_conversion_checkboxes()

        elif viz_type == "conversion_analysis":
            # 显示转换分析参数框架
            self.conversion_viz_frame.grid()
            self.viz_condition_combo.grid_remove()

            # 更新转换分析参数
            self.update_conversion_viz_params()

    def update_viz_condition_checkboxes(self):
        """更新可视化条件复选框"""
        # 清除现有复选框
        for widget in self.viz_conditions_frame.winfo_children():
            widget.destroy()

        # 获取所有条件
        conditions = list(self.calculator.ppg_data.keys())

        # 创建新的复选框
        ttk.Label(self.viz_conditions_frame, text="选择要比较的条件:").grid(row=0, column=0, columnspan=4, sticky=tk.W,
                                                                            pady=(0, 5))

        self.viz_condition_vars = {}
        for i, condition in enumerate(conditions):
            var = tk.BooleanVar(value=(i < 2))  # 默认选择前两个条件
            cb = ttk.Checkbutton(self.viz_conditions_frame, text=condition, variable=var)
            cb.grid(row=1 + i // 4, column=i % 4, sticky=tk.W, padx=5, pady=2)
            self.viz_condition_vars[condition] = var

    def update_viz_conversion_checkboxes(self):
        """更新可视化转换方案复选框"""
        # 清除现有复选框
        for widget in self.viz_conditions_frame.winfo_children():
            widget.destroy()

        # 获取所有条件
        conditions = list(self.calculator.ppg_data.keys())

        # 创建新的复选框
        ttk.Label(self.viz_conditions_frame, text="选择要比较的转换方案:").grid(row=0, column=0, columnspan=4,
                                                                                sticky=tk.W,
                                                                                pady=(0, 5))

        self.viz_conversion_vars = []

        if len(conditions) >= 2:
            for i in range(len(conditions)):
                for j in range(len(conditions)):
                    if i != j:
                        from_cond = conditions[i]
                        to_cond = conditions[j]
                        var = tk.BooleanVar(value=(i == 0 and j == 1))  # 默认选择第一个转换方案
                        cb = ttk.Checkbutton(self.viz_conditions_frame,
                                             text=f"{from_cond}→{to_cond}",
                                             variable=var)
                        row = (i * len(conditions) + j) // 4 + 1
                        col = (i * len(conditions) + j) % 4
                        cb.grid(row=row, column=col, sticky=tk.W, padx=2, pady=1)
                        self.viz_conversion_vars.append((var, from_cond, to_cond))
        else:
            ttk.Label(self.viz_conditions_frame, text="需要至少2个条件才能比较",
                      foreground="gray").grid(row=1, column=0, sticky=tk.W)

    def update_conversion_viz_params(self):
        """更新转换可视化参数"""
        # 清除现有控件
        for widget in self.conversion_viz_frame.winfo_children():
            widget.destroy()

        # 获取所有条件
        conditions = list(self.calculator.ppg_data.keys())

        # 创建源条件和目标条件选择
        ttk.Label(self.conversion_viz_frame, text="源条件:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=5)
        self.viz_from_condition_var = tk.StringVar()
        viz_from_combo = ttk.Combobox(self.conversion_viz_frame, textvariable=self.viz_from_condition_var,
                                      values=conditions, width=15, state="readonly")
        viz_from_combo.grid(row=0, column=1, sticky=tk.W, padx=(0, 10), pady=5)
        if conditions:
            self.viz_from_condition_var.set(conditions[0])

        ttk.Label(self.conversion_viz_frame, text="目标条件:").grid(row=0, column=2, sticky=tk.W, padx=(0, 5), pady=5)
        self.viz_to_condition_var = tk.StringVar()
        viz_to_combo = ttk.Combobox(self.conversion_viz_frame, textvariable=self.viz_to_condition_var,
                                    values=conditions, width=15, state="readonly")
        viz_to_combo.grid(row=0, column=3, sticky=tk.W, padx=(0, 10), pady=5)
        if len(conditions) > 1:
            self.viz_to_condition_var.set(conditions[1])

    def generate_visualization(self):
        """生成可视化图表"""
        viz_type = self.viz_type_var.get()

        if viz_type == "condition_comparison":
            # 获取选定的条件
            selected_conditions = []
            if hasattr(self, 'viz_condition_vars'):
                for condition, var in self.viz_condition_vars.items():
                    if var.get():
                        selected_conditions.append(condition)

            if len(selected_conditions) < 2:
                messagebox.showwarning("警告", "请至少选择两个条件进行比较")
                return

            self.create_visualization(viz_type, selected_conditions)

        elif viz_type == "multiple_conversion_comparison":
            # 获取选定的转换方案
            selected_conversions = []
            if hasattr(self, 'viz_conversion_vars'):
                for var, from_cond, to_cond in self.viz_conversion_vars:
                    if var.get():
                        selected_conversions.append((from_cond, to_cond))

            if len(selected_conversions) < 2:
                messagebox.showwarning("警告", "请至少选择两个转换方案进行比较")
                return

            self.create_visualization(viz_type, selected_conversions)

        elif viz_type == "conversion_analysis":
            # 获取转换参数
            from_condition = self.viz_from_condition_var.get()
            to_condition = self.viz_to_condition_var.get()

            if not from_condition or not to_condition:
                messagebox.showwarning("警告", "请选择源条件和目标条件")
                return

            if from_condition == to_condition:
                messagebox.showwarning("警告", "源条件和目标条件不能相同")
                return

            self.create_visualization(viz_type, (from_condition, to_condition))

        else:
            condition = self.viz_condition_var.get()
            if not condition:
                messagebox.showwarning("警告", "请选择条件")
                return

            self.create_visualization(viz_type, condition)

    def create_visualization(self, viz_type, conditions):
        """创建可视化图表"""
        try:
            # 确保可视化器已初始化
            if self.visualizer is None and self.calculator.standard_curves:
                self.visualizer = PPGVisualizer(self.calculator)
            elif self.visualizer is None:
                messagebox.showwarning("警告", "请先拟合标准曲线")
                return

            # 显示生成状态
            self.viz_status_var.set(f"正在生成 {viz_type} 图表...")
            self.update_status(f"生成图表: {viz_type}")

            # 生成图表
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
                self.viz_status_var.set("不支持的可视化类型")
                return

            if fig is None:
                self.viz_status_var.set("无法生成图表，请检查数据")
                return

            # 清除现有图表
            self.clear_figure()

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

            self.viz_status_var.set(f"{viz_type} 图表生成完成")
            self.update_status("就绪")

        except Exception as e:
            self.viz_status_var.set(f"生成图表失败: {str(e)}")
            self.log_message(f"✗ 生成图表失败: {str(e)}", "ERROR")

    def save_figure(self):
        """保存图表"""
        if self.current_figure is None:
            messagebox.showwarning("警告", "没有可保存的图表")
            return

        file_types = [("PNG文件", "*.png"), ("PDF文件", "*.pdf"),
                      ("SVG文件", "*.svg"), ("所有文件", "*.*")]

        file_path = filedialog.asksaveasfilename(
            title="保存图表",
            filetypes=file_types,
            defaultextension=".png"
        )

        if file_path:
            try:
                self.current_figure.savefig(file_path, dpi=300, bbox_inches='tight')
                self.log_message(f"✓ 图表已保存到: {file_path}", "SUCCESS")
                self.viz_status_var.set(f"图表已保存: {Path(file_path).name}")
            except Exception as e:
                self.log_message(f"✗ 保存图表失败: {str(e)}", "ERROR")

    def clear_figure(self):
        """清除图表"""
        if self.figure_canvas:
            self.figure_canvas.get_tk_widget().destroy()
            self.figure_toolbar.destroy()
            self.figure_canvas = None
            self.figure_toolbar = None

        if self.current_figure:
            plt.close(self.current_figure)
            self.current_figure = None

        # 显示占位标签
        self.viz_placeholder.pack(expand=True)

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
            success, msg, saved_files = self.calculator.save_results(output_dir)

            if success:
                self.root.after(0, lambda: self.output_status_var.set(f"结果保存完成: {len(saved_files)} 个文件"))
                self.root.after(0, lambda: self.log_message(f"✓ {msg}", "SUCCESS"))

                # 在结果标签页显示保存的文件
                self.results_text.delete(1.0, tk.END)
                self.results_text.insert(tk.END, "保存的文件列表:\n")
                self.results_text.insert(tk.END, "=" * 60 + "\n\n")

                for file_path in saved_files:
                    self.results_text.insert(tk.END, f"• {Path(file_path).name}\n")

                self.results_text.insert(tk.END, f"\n所有文件已保存到: {output_dir}")
            else:
                self.root.after(0, lambda: self.output_status_var.set("保存失败"))
                self.root.after(0, lambda: self.log_message(f"✗ {msg}", "ERROR"))

            self.root.after(0, lambda: self.update_status("就绪"))

        except Exception as e:
            self.root.after(0, lambda: self.output_status_var.set("保存失败"))
            self.root.after(0, lambda: self.log_message(f"✗ 保存结果失败: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.update_status("保存失败"))

    def preview_report(self):
        """预览报告"""
        # 生成报告
        self.generate_report()

        # 切换到结果标签页
        self.notebook.select(4)  # 结果标签页是第5个（索引4）

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

    def clear_log(self):
        """清空日志"""
        self.log_text.delete(1.0, tk.END)

    def clear_analysis_text(self):
        """清空分析文本"""
        self.analysis_text.delete(1.0, tk.END)

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
        'openpyxl': 'Excel文件支持',
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
    root.title("PPG保留指数计算与可视化分析系统 - 增强版")

    # 获取屏幕尺寸
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()

    # 设置窗口大小（屏幕的80%）
    window_width = int(screen_width * 0.8)
    window_height = int(screen_height * 0.8)

    # 计算窗口位置（居中）
    x = (screen_width - window_width) // 2
    y = (screen_height - window_height) // 2

    root.geometry(f"{window_width}x{window_height}+{x}+{y}")

    # 创建GUI
    app = PPGIndexAnalyzerGUI(root)

    # 运行主循环
    root.mainloop()


if __name__ == "__main__":
    print("PPG保留指数计算与可视化分析系统 - 增强版")
    print("版本: 3.0 (含跨条件转换功能)")
    print("功能:")
    print("  1. 加载PPG标准品和化合物保留时间数据")
    print("  2. 拟合PPG标准曲线并计算线性关系")
    print("  3. 计算化合物的PPG保留指数")
    print("  4. 比较不同色谱条件下的PPG指数")
    print("  5. 跨条件PPG指数转换与验证")
    print("  6. 多转换方案比较分析")
    print("  7. 可视化分析结果并生成报告")
    print("=" * 70)

    main()
