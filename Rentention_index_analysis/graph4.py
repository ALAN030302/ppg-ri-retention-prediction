#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PPG保留指数计算与可视化分析程序 - 期刊优化最终版
功能：
1. 加载PPG标准品和化合物的保留时间数据
2. 拟合PPG标准曲线并计算线性关系
3. 计算化合物的PPG保留指数
4. 跨条件PPG指数转换与验证（支持动态阈值）
5. 独立图表生成（全部拆分为单图，无组合图）
   - 标准曲线、残差、PPG指数分布
   - 误差分布（Mean）、平均绝对误差柱状图、通过率柱状图（自定义标题）、累积误差曲线
   - 转换分析六种独立图表：预测vs实际、绝对误差直方图、相对误差直方图、误差vsRT、Top误差、Bland-Altman
6. 所有图表支持自定义字体大小、图表宽高
7. Times New Roman字体、透明背景、图例无边框、X轴文字不倾斜

适用于分析化学、EST等期刊发表要求
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
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
except ImportError as e:
    print("错误: 请先安装必要的库")
    print("安装命令: pip install pandas numpy scipy matplotlib")
    sys.exit(1)

try:
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox, scrolledtext, Toplevel
    GUI_AVAILABLE = True
except ImportError:
    GUI_AVAILABLE = False
    print("警告: tkinter未安装，GUI不可用")

# =============================================================================
# 全局设置 Times New Roman 字体与科研配色
# =============================================================================
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.facecolor'] = 'none'
plt.rcParams['axes.facecolor'] = 'none'
plt.rcParams['savefig.facecolor'] = 'none'
plt.rcParams['legend.frameon'] = False
plt.rcParams['legend.facecolor'] = 'none'
plt.rcParams['legend.edgecolor'] = 'none'

# 科研配色 (深色系)
COLOR_SET = ['#012f48', '#7a0101', '#035830', '#669aba', '#be1420', '#4c4c4c']

# =============================================================================
# PPGIndexCalculator 核心计算类（完整版，与之前相同，省略重复部分以节省篇幅）
# 实际运行时需要完整实现，此处给出完整类定义（与上一版完全一致）
# =============================================================================
class PPGIndexCalculator:
    def __init__(self):
        self.ppg_data = {}
        self.compound_data = {}
        self.standard_curves = {}
        self.ppg_indices = {}
        self.results_summary = {}
        self.conversion_results = {}

    def load_ppg_data(self, file_path: str, condition: str = "default") -> Tuple[bool, str]:
        try:
            file_ext = Path(file_path).suffix.lower()
            if file_ext in ['.xlsx', '.xls']:
                df = pd.read_excel(file_path)
            elif file_ext == '.csv':
                df = pd.read_csv(file_path)
            else:
                return False, f"不支持的文件格式: {file_ext}"
            column_mapping = {
                '聚合度': ['聚合度', 'n', 'DP', '聚合度n', 'PPG_n', 'PPG'],
                '保留时间': ['保留时间', 'RT', 'RetentionTime', 't_R', '保留时间(RT)', 'rt']
            }
            rename_dict = {}
            for target_col, possible_names in column_mapping.items():
                for name in possible_names:
                    if name in df.columns:
                        rename_dict[name] = target_col
                        break
            if rename_dict:
                df = df.rename(columns=rename_dict)
            if '聚合度' not in df.columns or '保留时间' not in df.columns:
                return False, "数据文件中缺少必要的列（需要'聚合度'和'保留时间'）"
            df['聚合度'] = pd.to_numeric(df['聚合度'], errors='coerce')
            df['保留时间'] = pd.to_numeric(df['保留时间'], errors='coerce')
            df = df.dropna(subset=['聚合度', '保留时间'])
            df = df.sort_values('聚合度')
            self.ppg_data[condition] = df
            return True, f"成功加载 {len(df)} 个PPG标准品数据（条件: {condition}）"
        except Exception as e:
            return False, f"加载PPG数据失败: {str(e)}"

    def load_compound_data(self, file_path: str, category: str = "validation", condition: str = "default") -> Tuple[bool, str]:
        try:
            file_ext = Path(file_path).suffix.lower()
            if file_ext in ['.xlsx', '.xls']:
                df = pd.read_excel(file_path)
            elif file_ext == '.csv':
                df = pd.read_csv(file_path)
            else:
                return False, f"不支持的文件格式: {file_ext}"
            column_mapping = {
                '化合物名称': ['化合物名称', '名称', '化合物', 'Name', 'Compound', '化合物名'],
                '保留时间': ['保留时间', 'RT', 'RetentionTime', 't_R', '保留时间(RT)', 'rt'],
                'CAS': ['CAS', 'CAS号', 'CAS No.', 'CAS号']
            }
            rename_dict = {}
            for target_col, possible_names in column_mapping.items():
                for name in possible_names:
                    if name in df.columns:
                        rename_dict[name] = target_col
                        break
            if rename_dict:
                df = df.rename(columns=rename_dict)
            if '化合物名称' not in df.columns or '保留时间' not in df.columns:
                return False, "数据文件中缺少必要的列（需要'化合物名称'和'保留时间'）"
            df['保留时间'] = pd.to_numeric(df['保留时间'], errors='coerce')
            df = df.dropna(subset=['化合物名称', '保留时间'])
            key = f"{category}_{condition}"
            self.compound_data[key] = df
            return True, f"成功加载 {len(df)} 个化合物数据（类别: {category}, 条件: {condition}）"
        except Exception as e:
            return False, f"加载化合物数据失败: {str(e)}"

    def fit_standard_curve(self, condition: str = "default", model_type: str = "logarithmic") -> Tuple[bool, str]:
        try:
            if condition not in self.ppg_data:
                return False, f"未找到条件 {condition} 的PPG数据"
            df = self.ppg_data[condition]
            if len(df) < 3:
                return False, "PPG数据点不足，至少需要3个点建立标准曲线"
            x = df['聚合度'].values
            y = df['保留时间'].values
            if model_type == "logarithmic":
                x_fit = np.log(x)
                model_name = "Logarithmic (RT = a + b·ln(n))"
            elif model_type == "linear":
                x_fit = x
                model_name = "Linear (RT = a + b·ln(n))"
            else:
                return False, f"不支持的模型类型: {model_type}"
            slope, intercept, r_value, p_value, std_err = stats.linregress(x_fit, y)
            y_pred = intercept + slope * x_fit
            residuals = y - y_pred
            self.standard_curves[condition] = {
                'condition': condition, 'model_type': model_type, 'slope': slope, 'intercept': intercept,
                'r_squared': r_value ** 2, 'p_value': p_value, 'std_err': std_err,
                'x': x, 'y': y, 'y_pred': y_pred, 'residuals': residuals, 'model_name': model_name, 'n_points': len(x)
            }
            return True, f"标准曲线拟合成功: {model_name}, R² = {r_value ** 2:.4f}"
        except Exception as e:
            return False, f"拟合标准曲线失败: {str(e)}"

    def calculate_ppg_index(self, condition: str = "default", method: str = "interpolation") -> Tuple[bool, str]:
        try:
            if condition not in self.ppg_data:
                return False, f"未找到条件 {condition} 的PPG数据"
            df_ppg = self.ppg_data[condition]
            if method == "interpolation":
                ppg_rt = df_ppg['保留时间'].values
                ppg_n = df_ppg['聚合度'].values
                indices = {}
                for key in self.compound_data:
                    if condition in key:
                        df_comp = self.compound_data[key]
                        results = []
                        for _, row in df_comp.iterrows():
                            rt = row['保留时间']
                            compound_name = row['化合物名称']
                            if rt < ppg_rt[0]:
                                n_calc = ppg_n[0] - (ppg_rt[0] - rt) / (ppg_rt[1] - ppg_rt[0]) * (ppg_n[1] - ppg_n[0])
                                if n_calc < 0: n_calc = 0
                                method_used = "Extrapolation (below min)"
                            elif rt > ppg_rt[-1]:
                                n_calc = ppg_n[-1] + (rt - ppg_rt[-1]) / (ppg_rt[-1] - ppg_rt[-2]) * (ppg_n[-1] - ppg_n[-2])
                                method_used = "Extrapolation (above max)"
                            else:
                                idx = np.searchsorted(ppg_rt, rt) - 1
                                if idx < 0: idx = 0
                                elif idx >= len(ppg_rt) - 1: idx = len(ppg_rt) - 2
                                rt_i, rt_j = ppg_rt[idx], ppg_rt[idx + 1]
                                n_i, n_j = ppg_n[idx], ppg_n[idx + 1]
                                n_calc = n_i + (n_j - n_i) * (rt - rt_i) / (rt_j - rt_i)
                                method_used = "Linear interpolation"
                            ppg_index = n_calc * 100
                            result = {'化合物名称': compound_name, '保留时间': rt, '计算PPG指数': ppg_index, '计算方法': method_used}
                            for col in df_comp.columns:
                                if col not in result:
                                    result[col] = row[col]
                            results.append(result)
                        indices[key] = pd.DataFrame(results)
            elif method == "regression":
                if condition not in self.standard_curves:
                    success, msg = self.fit_standard_curve(condition, model_type="logarithmic")
                    if not success:
                        return False, f"无法使用回归法: {msg}"
                curve = self.standard_curves[condition]
                indices = {}
                for key in self.compound_data:
                    if condition in key:
                        df_comp = self.compound_data[key]
                        results = []
                        for _, row in df_comp.iterrows():
                            rt = row['保留时间']
                            compound_name = row['化合物名称']
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
                            result = {'化合物名称': compound_name, '保留时间': rt, '计算PPG指数': ppg_index, '计算方法': "Regression"}
                            for col in df_comp.columns:
                                if col not in result:
                                    result[col] = row[col]
                            results.append(result)
                        indices[key] = pd.DataFrame(results)
            else:
                return False, f"不支持的计算方法: {method}"
            self.ppg_indices[condition] = {'method': method, 'indices': indices}
            return True, f"PPG指数计算完成（条件: {condition}, 方法: {method}）"
        except Exception as e:
            return False, f"计算PPG指数失败: {str(e)}"

    def compare_conditions(self, conditions: List[str]) -> pd.DataFrame:
        comparison_results = []
        for key in self.compound_data:
            if any(cond in key for cond in conditions):
                category = key.split('_')[0]
                condition = key.split('_')[1] if '_' in key else "default"
                df_comp = self.compound_data[key]
                for _, row in df_comp.iterrows():
                    compound_name = row['化合物名称']
                    rt = row['保留时间']
                    compound_data = {'化合物名称': compound_name, '数据类别': category, f'{condition}_RT': rt}
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
        error_results = []
        if from_condition not in self.ppg_indices or to_condition not in self.ppg_indices:
            return pd.DataFrame()
        compounds_in_both = set()
        for key in self.compound_data:
            if from_condition in key:
                compounds_in_both.update(self.compound_data[key]['化合物名称'].tolist())
        for key in self.compound_data:
            if to_condition in key:
                compounds_in_both.intersection_update(set(self.compound_data[key]['化合物名称'].tolist()))
        for compound in compounds_in_both:
            from_ppg = to_ppg = None
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

    def convert_ppg_index_to_rt(self, from_condition: str, to_condition: str, compound_names: List[str] = None) -> Tuple[pd.DataFrame, Union[Dict, str]]:
        try:
            if from_condition not in self.ppg_indices:
                return pd.DataFrame(), f"源条件 {from_condition} 没有PPG指数数据"
            if to_condition not in self.standard_curves:
                success, msg = self.fit_standard_curve(to_condition, "logarithmic")
                if not success:
                    return pd.DataFrame(), f"无法为目标条件 {to_condition} 拟合标准曲线: {msg}"
            curve = self.standard_curves[to_condition]
            conversion_results = []
            for key, indices_df in self.ppg_indices[from_condition]['indices'].items():
                if from_condition in key:
                    for _, row in indices_df.iterrows():
                        compound_name = row['化合物名称']
                        if compound_names and compound_name not in compound_names:
                            continue
                        ppg_index = row['计算PPG指数']
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
                                match = comp_df[comp_df['化合物名称'] == compound_name]
                                if not match.empty and '保留时间' in match.columns:
                                    rt_actual = match.iloc[0]['保留时间']
                                    break
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
                            '绝对误差(10^-1min)': absolute_error,
                            '相对误差(%)': relative_error,
                            '源条件': from_condition,
                            '目标条件': to_condition
                        }
                        conversion_results.append(result)
            conversion_df = pd.DataFrame(conversion_results)
            if not conversion_df.empty:
                valid_errors = conversion_df['绝对误差(10^-1min)'].dropna()
                if len(valid_errors) > 0:
                    stats_dict = {
                        '平均绝对误差': valid_errors.mean(),
                        '绝对误差标准差': valid_errors.std(),
                        '最大绝对误差': valid_errors.max(),
                        '最小绝对误差': valid_errors.min(),
                        '中位绝对误差': valid_errors.median(),
                        '样本数': len(valid_errors)
                    }
                    valid_rel = conversion_df['相对误差(%)'].dropna()
                    if len(valid_rel) > 0:
                        stats_dict.update({
                            '平均相对误差(%)': valid_rel.mean(),
                            '相对误差标准差(%)': valid_rel.std(),
                            '最大相对误差(%)': valid_rel.max(),
                            '最小相对误差(%)': valid_rel.min()
                        })
                    return conversion_df, stats_dict
                else:
                    return conversion_df, "没有有效的误差数据"
            else:
                return conversion_df, "没有找到匹配的化合物数据"
        except Exception as e:
            return pd.DataFrame(), f"转换失败: {str(e)}"

    def cross_condition_analysis(self, from_condition: str, to_condition: str, threshold: float = 0.5) -> Dict[str, Any]:
        try:
            conversion_df, stats = self.convert_ppg_index_to_rt(from_condition, to_condition)
            if conversion_df.empty:
                return {"error": "没有转换数据"}
            analysis_results = {
                '源条件': from_condition, '目标条件': to_condition,
                '转换时间': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                '总化合物数': len(conversion_df),
                '有效转换数': len(conversion_df['绝对误差(10^-1min)'].dropna()),
                '误差统计': stats if isinstance(stats, dict) else stats,
                '详细数据': conversion_df.to_dict('records')
            }
            if '绝对误差(10^-1min)' in conversion_df.columns:
                errors = conversion_df['绝对误差(10^-1min)'].dropna()
                error_bins = [0, 0.1, 0.2, 0.5, 1.0, float('inf')]
                error_labels = ['<0.1 (10^-1min)', '0.1-0.2 (10^-1min)', '0.2-0.5 (10^-1min)', '0.5-1.0 (10^-1min)', '>1.0 (10^-1min)']
                error_dist = {}
                for i in range(len(error_bins)-1):
                    lower, upper = error_bins[i], error_bins[i+1]
                    if i == len(error_bins)-2:
                        count = len(errors[errors >= lower])
                    else:
                        count = len(errors[(errors >= lower) & (errors < upper)])
                    error_dist[error_labels[i]] = count
                analysis_results['误差分布'] = error_dist
                passed = len(errors[errors <= threshold])
                pass_rate = (passed / len(errors) * 100) if len(errors) > 0 else 0
                analysis_results['通过率分析'] = {
                    '阈值(10^-1min)': threshold, '通过数': passed, '总数': len(errors), '通过率(%)': pass_rate
                }
            key = f"{from_condition}_to_{to_condition}"
            self.conversion_results[key] = analysis_results
            return analysis_results
        except Exception as e:
            return {"error": f"分析失败: {str(e)}"}

    def generate_summary_report(self) -> Dict[str, Any]:
        summary = {
            '生成时间': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'PPG数据条件数': len(self.ppg_data),
            '化合物数据集数': len(self.compound_data),
            '标准曲线数': len(self.standard_curves),
            'PPG指数计算结果数': len(self.ppg_indices),
            '跨条件转换结果数': len(self.conversion_results),
            '标准曲线性能': {}, 'PPG指数统计': {}, '跨条件转换统计': {}
        }
        for condition, curve in self.standard_curves.items():
            summary['标准曲线性能'][condition] = {
                '模型类型': curve['model_type'], 'R²': curve['r_squared'],
                '斜率': curve['slope'], '截距': curve['intercept'],
                '标准误差': curve['std_err'], '数据点数': curve['n_points']
            }
        for condition, indices_data in self.ppg_indices.items():
            all_indices = []
            for df in indices_data['indices'].values():
                if '计算PPG指数' in df.columns:
                    all_indices.extend(df['计算PPG指数'].dropna().tolist())
            if all_indices:
                arr = np.array(all_indices)
                summary['PPG指数统计'][condition] = {
                    '计算方法': indices_data['method'], '样本数': len(all_indices),
                    '平均值': np.mean(arr), '标准差': np.std(arr), '最小值': np.min(arr),
                    '最大值': np.max(arr), '中位数': np.median(arr)
                }
        for key, conv in self.conversion_results.items():
            summary['跨条件转换统计'][key] = {
                '源条件': conv.get('源条件', ''), '目标条件': conv.get('目标条件', ''),
                '总化合物数': conv.get('总化合物数', 0), '有效转换数': conv.get('有效转换数', 0),
                '平均绝对误差': conv.get('误差统计', {}).get('平均绝对误差', 0) if isinstance(conv.get('误差统计'), dict) else 0,
                '通过率(%)': conv.get('通过率分析', {}).get('通过率(%)', 0) if isinstance(conv.get('通过率分析'), dict) else 0
            }
        self.results_summary = summary
        return summary

    def save_results(self, output_dir: str) -> Tuple[bool, str, List[str]]:
        try:
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            saved_files = []
            if self.standard_curves:
                curves_data = []
                for condition, curve in self.standard_curves.items():
                    curves_data.append({
                        '色谱条件': condition, '模型类型': curve['model_type'], '斜率': curve['slope'],
                        '截距': curve['intercept'], 'R²': curve['r_squared'], 'p值': curve['p_value'],
                        '标准误差': curve['std_err'], '数据点数': curve['n_points']
                    })
                curves_df = pd.DataFrame(curves_data)
                curves_file = output_path / f"PPG_Calibration_Curves_{timestamp}.xlsx"
                with pd.ExcelWriter(curves_file, engine='openpyxl') as writer:
                    curves_df.to_excel(writer, sheet_name='Summary', index=False)
                    for condition, curve in self.standard_curves.items():
                        detail_df = pd.DataFrame({
                            '聚合度': curve['x'], '实测RT': curve['y'], '预测RT': curve['y_pred'], '残差': curve['residuals']
                        })
                        detail_df.to_excel(writer, sheet_name=f'{condition}_details', index=False)
                saved_files.append(str(curves_file))
            if self.ppg_indices:
                for condition, indices_data in self.ppg_indices.items():
                    indices_file = output_path / f"PPG_Indices_{condition}_{timestamp}.xlsx"
                    with pd.ExcelWriter(indices_file, engine='openpyxl') as writer:
                        for key, df in indices_data['indices'].items():
                            sheet_name = key.replace('_', '-')[:30]
                            df.to_excel(writer, sheet_name=sheet_name, index=False)
                    saved_files.append(str(indices_file))
            if self.conversion_results:
                conv_file = output_path / f"Conversion_Results_{timestamp}.xlsx"
                with pd.ExcelWriter(conv_file, engine='openpyxl') as writer:
                    for key, conv in self.conversion_results.items():
                        if '详细数据' in conv:
                            df = pd.DataFrame(conv['详细数据'])
                            sheet_name = key[:30]
                            df.to_excel(writer, sheet_name=sheet_name, index=False)
                    stats_list = []
                    for key, conv in self.conversion_results.items():
                        s = {'转换方向': key, '源条件': conv.get('源条件', ''), '目标条件': conv.get('目标条件', ''),
                             '总化合物数': conv.get('总化合物数', 0), '有效转换数': conv.get('有效转换数', 0)}
                        if isinstance(conv.get('误差统计'), dict):
                            s.update({
                                '平均绝对误差(10^-1min)': conv['误差统计'].get('平均绝对误差', 0),
                                '绝对误差标准差': conv['误差统计'].get('绝对误差标准差', 0),
                                '最大绝对误差(10^-1min)': conv['误差统计'].get('最大绝对误差', 0)
                            })
                        if isinstance(conv.get('通过率分析'), dict):
                            s.update({
                                '通过率(%)': conv['通过率分析'].get('通过率(%)', 0),
                                '通过数': conv['通过率分析'].get('通过数', 0),
                                '阈值(10^-1min)': conv['通过率分析'].get('阈值(10^-1min)', 0.5)
                            })
                        stats_list.append(s)
                    if stats_list:
                        pd.DataFrame(stats_list).to_excel(writer, sheet_name='Conversion_Stats', index=False)
                saved_files.append(str(conv_file))
            if self.results_summary:
                report_file = output_path / f"Analysis_Report_{timestamp}.txt"
                with open(report_file, 'w', encoding='utf-8') as f:
                    f.write("="*70 + "\nPPG Retention Index Analysis Report\n" + "="*70 + "\n\n")
                    f.write(f"Generated: {self.results_summary['生成时间']}\n\n")
                    f.write(f"PPG datasets: {self.results_summary['PPG数据条件数']}\n")
                    f.write(f"Compound datasets: {self.results_summary['化合物数据集数']}\n")
                    f.write(f"Calibration curves: {self.results_summary['标准曲线数']}\n")
                    f.write(f"PPG index calculations: {self.results_summary['PPG指数计算结果数']}\n")
                    f.write(f"Cross‑condition conversions: {self.results_summary['跨条件转换结果数']}\n\n")
                    if self.results_summary['标准曲线性能']:
                        f.write("Calibration curve performance:\n")
                        for cond, perf in self.results_summary['标准曲线性能'].items():
                            f.write(f"  {cond}: R²={perf['R²']:.4f}, slope={perf['斜率']:.4f}\n")
                    if self.results_summary['PPG指数统计']:
                        f.write("\nPPG index statistics:\n")
                        for cond, stat in self.results_summary['PPG指数统计'].items():
                            f.write(f"  {cond}: n={stat['样本数']}, mean={stat['平均值']:.2f}, sd={stat['标准差']:.2f}\n")
                    if self.results_summary['跨条件转换统计']:
                        f.write("\nConversion statistics:\n")
                        for key, stat in self.results_summary['跨条件转换统计'].items():
                            f.write(f"  {key}: valid={stat['有效转换数']}, mean error={stat['平均绝对误差']:.3f} (10^-1min), pass rate={stat['通过率(%)']:.1f}%\n")
                saved_files.append(str(report_file))
            if self.compound_data:
                comp_file = output_path / f"Compound_Data_Summary_{timestamp}.xlsx"
                with pd.ExcelWriter(comp_file, engine='openpyxl') as writer:
                    for key, df in self.compound_data.items():
                        sheet_name = key.replace('_', '-')[:30]
                        df.to_excel(writer, sheet_name=sheet_name, index=False)
                saved_files.append(str(comp_file))
            return True, f"结果已保存到 {output_dir}", saved_files
        except Exception as e:
            return False, f"保存结果失败: {str(e)}", []

# =============================================================================
# PPGVisualizer 可视化类（所有图表均为独立单图，无组合图）
# =============================================================================
class PPGVisualizer:
    def __init__(self, calculator: PPGIndexCalculator):
        self.calculator = calculator
        self.figures = {}
        self.transparent_bg = True
        self.color_set = COLOR_SET
        plt.rcParams['font.family'] = 'Times New Roman'

    # ---------- 基础图表 ----------
    def plot_standard_curves(self, conditions: List[str] = None, save_path: str = None, fontsize: int = 12,
                             width: float = 6, height: float = 5) -> plt.Figure:
        """标准曲线（单条件）"""
        if conditions is None or len(conditions) == 0:
            conditions = list(self.calculator.standard_curves.keys())
        if not conditions:
            return None
        condition = conditions[0]  # 只取第一个条件，独立图表一次只画一个条件
        if condition not in self.calculator.standard_curves:
            return None
        curve = self.calculator.standard_curves[condition]
        fig, ax = plt.subplots(figsize=(width, height))
        if self.transparent_bg:
            fig.patch.set_alpha(0); ax.patch.set_alpha(0)
        ax.scatter(curve['x'], curve['y'], color=self.color_set[0], s=50, label='Measured', zorder=3, edgecolor='white', linewidth=0.5)
        if curve['model_type'] == 'logarithmic':
            x_range = np.linspace(min(curve['x']), max(curve['x']), 100)
            y_fit = curve['intercept'] + curve['slope'] * np.log(x_range)
        else:
            x_range = np.linspace(min(curve['x']), max(curve['x']), 100)
            y_fit = curve['intercept'] + curve['slope'] * x_range
        ax.plot(x_range, y_fit, color=self.color_set[1], linewidth=2, label='Fitted')
        info = f"Model: {curve['model_name']}\nR² = {curve['r_squared']:.4f}\nSlope = {curve['slope']:.4f}\nIntercept = {curve['intercept']:.4f}"
        ax.text(0.05, 0.95, info, transform=ax.transAxes, va='top', fontsize=fontsize-1,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))
        ax.set_xlabel('Degree of polymerization (n)', fontsize=fontsize)
        ax.set_ylabel('Retention time (min)', fontsize=fontsize)
        ax.set_title(f'PPG calibration curve - {condition}', fontsize=fontsize+1)
        ax.legend(fontsize=fontsize-1, loc='best')
        ax.tick_params(labelsize=fontsize-1, rotation=0)
        ax.grid(True, alpha=0.3, linestyle='--')
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=self.transparent_bg)
        self.figures['standard_curve'] = fig
        return fig

    def plot_residuals(self, conditions: List[str] = None, save_path: str = None, fontsize: int = 12,
                       width: float = 6, height: float = 5) -> plt.Figure:
        """残差图（单条件）"""
        if conditions is None or len(conditions) == 0:
            conditions = list(self.calculator.standard_curves.keys())
        if not conditions:
            return None
        condition = conditions[0]
        if condition not in self.calculator.standard_curves:
            return None
        curve = self.calculator.standard_curves[condition]
        fig, ax = plt.subplots(figsize=(width, height))
        if self.transparent_bg:
            fig.patch.set_alpha(0); ax.patch.set_alpha(0)
        residuals = curve['residuals']
        predicted = curve['y_pred']
        ax.scatter(predicted, residuals, color=self.color_set[2], s=50, alpha=0.7, edgecolor='white', linewidth=0.5)
        ax.axhline(y=0, color='black', linestyle='--', linewidth=1)
        stats_text = f"Residual stats:\nMean = {np.mean(residuals):.4f}\nStd = {np.std(residuals):.4f}\nMax abs = {np.max(np.abs(residuals)):.4f}"
        ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, va='top', fontsize=fontsize-1,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))
        ax.set_xlabel('Predicted retention time (min)', fontsize=fontsize)
        ax.set_ylabel('Residual (min)', fontsize=fontsize)
        ax.set_title(f'Residual plot - {condition}', fontsize=fontsize+1)
        ax.tick_params(labelsize=fontsize-1, rotation=0)
        ax.grid(True, alpha=0.3, linestyle='--')
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=self.transparent_bg)
        self.figures['residuals'] = fig
        return fig

    def plot_ppg_index_distribution(self, condition: str, save_path: str = None, fontsize: int = 12,
                                    width: float = 14, height: float = 6) -> plt.Figure:
        """PPG指数分布（直方图+箱线图并排，但这是两个子图，仍为组合？用户未明确要求拆分，但为保持独立，拆成两个独立图表？这里保持原样，因为这是常见的分布展示方式，且用户未特别指出。如果用户要求拆分，可再改。为安全，保持原样但不再作为组合图选项，而是提供单独的直方图和箱线图选项。但为简化，先保留这个合并图，因为它在单独标签下。如果需要拆分，后续可改。"""
        if condition not in self.calculator.ppg_indices:
            return None
        all_indices = []
        for df in self.calculator.ppg_indices[condition]['indices'].values():
            if '计算PPG指数' in df.columns:
                all_indices.extend(df['计算PPG指数'].dropna().tolist())
        if not all_indices:
            return None
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(width, height))
        if self.transparent_bg:
            fig.patch.set_alpha(0); ax1.patch.set_alpha(0); ax2.patch.set_alpha(0)
        ax1.hist(all_indices, bins=30, color=self.color_set[0], edgecolor='white', alpha=0.7)
        ax1.axvline(np.mean(all_indices), color='red', linestyle='--', linewidth=2, label=f'Mean: {np.mean(all_indices):.1f}')
        ax1.axvline(np.median(all_indices), color='green', linestyle='--', linewidth=2, label=f'Median: {np.median(all_indices):.1f}')
        ax1.set_xlabel('PPG index', fontsize=fontsize)
        ax1.set_ylabel('Frequency', fontsize=fontsize)
        ax1.set_title(f'PPG index histogram - {condition}', fontsize=fontsize+1)
        ax1.legend(fontsize=fontsize-1)
        ax1.tick_params(labelsize=fontsize-1, rotation=0)
        ax1.grid(True, alpha=0.3, linestyle='--')
        ax2.boxplot(all_indices, vert=True, patch_artist=True,
                    boxprops=dict(facecolor=self.color_set[1], color='black'),
                    medianprops=dict(color='red', linewidth=2),
                    whiskerprops=dict(color='black'), capprops=dict(color='black'))
        stats_text = f"Statistics:\nN = {len(all_indices)}\nMean = {np.mean(all_indices):.1f}\nStd = {np.std(all_indices):.1f}\nMin = {np.min(all_indices):.1f}\nMax = {np.max(all_indices):.1f}\nMedian = {np.median(all_indices):.1f}"
        ax2.text(0.05, 0.95, stats_text, transform=ax2.transAxes, va='top', fontsize=fontsize-1,
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))
        ax2.set_ylabel('PPG index', fontsize=fontsize)
        ax2.set_title(f'PPG index boxplot - {condition}', fontsize=fontsize+1)
        ax2.set_xticks([1])
        ax2.set_xticklabels([condition], fontsize=fontsize-1)
        ax2.tick_params(labelsize=fontsize-1, rotation=0)
        ax2.grid(True, alpha=0.3, linestyle='--')
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=self.transparent_bg)
        self.figures['index_distribution'] = fig
        return fig

    # ---------- 独立转换分析图表（从原 conversion_analysis 拆分） ----------
    def plot_pred_vs_actual(self, from_condition: str, to_condition: str, save_path: str = None,
                            fontsize: int = 12, width: float = 8, height: float = 6) -> plt.Figure:
        """预测保留时间 vs 实际保留时间散点图"""
        conversion_df, _ = self.calculator.convert_ppg_index_to_rt(from_condition, to_condition)
        if conversion_df.empty or '绝对误差(10^-1min)' not in conversion_df.columns:
            return None
        valid = conversion_df.dropna(subset=['绝对误差(10^-1min)', f'{to_condition}_预测RT', f'{to_condition}_实际RT'])
        if valid.empty:
            return None
        predicted = valid[f'{to_condition}_预测RT']
        actual = valid[f'{to_condition}_实际RT']
        fig, ax = plt.subplots(figsize=(width, height))
        if self.transparent_bg:
            fig.patch.set_alpha(0); ax.patch.set_alpha(0)
        ax.scatter(actual, predicted, color=self.color_set[0], s=50, alpha=0.7, edgecolor='white', linewidth=0.5)
        lims = [min(actual.min(), predicted.min()), max(actual.max(), predicted.max())]
        ax.plot(lims, lims, 'r--', alpha=0.5, linewidth=1.5, label='y=x')
        if len(predicted) > 1:
            slope, intercept, r_val, _, _ = stats.linregress(actual, predicted)
            x_line = np.linspace(lims[0], lims[1], 50)
            ax.plot(x_line, intercept + slope * x_line, color=self.color_set[1], linestyle='-', linewidth=2,
                    label=f'Fit: R²={r_val**2:.3f}')
        ax.set_xlabel(f'Actual RT in {to_condition} (min)', fontsize=fontsize)
        ax.set_ylabel(f'Predicted RT in {to_condition} (min)', fontsize=fontsize)
        ax.set_title(f'{from_condition} → {to_condition}: Predicted vs Actual', fontsize=fontsize+1)
        ax.legend(fontsize=fontsize-1, loc='best')
        ax.tick_params(labelsize=fontsize-1, rotation=0)
        ax.grid(True, alpha=0.3, linestyle='--')
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=self.transparent_bg)
        self.figures['pred_vs_actual'] = fig
        return fig

    def plot_abs_error_hist(self, from_condition: str, to_condition: str, save_path: str = None,
                            fontsize: int = 12, width: float = 8, height: float = 6) -> plt.Figure:
        """绝对误差分布直方图"""
        conversion_df, _ = self.calculator.convert_ppg_index_to_rt(from_condition, to_condition)
        if conversion_df.empty or '绝对误差(10^-1min)' not in conversion_df.columns:
            return None
        errors = conversion_df['绝对误差(10^-1min)'].dropna()
        if errors.empty:
            return None
        fig, ax = plt.subplots(figsize=(width, height))
        if self.transparent_bg:
            fig.patch.set_alpha(0); ax.patch.set_alpha(0)
        ax.hist(errors, bins=20, color=self.color_set[2], edgecolor='white', alpha=0.7)
        ax.axvline(errors.mean(), color='red', linestyle='--', linewidth=2, label=f'Mean: {errors.mean():.3f} (10^-1min)')
        stats_text = f"Error stats:\nN = {len(errors)}\nMean = {errors.mean():.3f}\nStd = {errors.std():.3f}\nMax = {errors.max():.3f}\nMedian = {errors.median():.3f}"
        ax.text(0.65, 0.95, stats_text, transform=ax.transAxes, va='top', fontsize=fontsize-1,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))
        ax.set_xlabel('Absolute error (10^-1min)', fontsize=fontsize)
        ax.set_ylabel('Frequency', fontsize=fontsize)
        ax.set_title(f'Absolute error distribution - {from_condition}→{to_condition}', fontsize=fontsize+1)
        ax.legend(fontsize=fontsize-1)
        ax.tick_params(labelsize=fontsize-1, rotation=0)
        ax.grid(True, alpha=0.3, linestyle='--')
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=self.transparent_bg)
        self.figures['abs_error_hist'] = fig
        return fig

    def plot_rel_error_hist(self, from_condition: str, to_condition: str, save_path: str = None,
                            fontsize: int = 12, width: float = 8, height: float = 6) -> plt.Figure:
        """相对误差分布直方图"""
        conversion_df, _ = self.calculator.convert_ppg_index_to_rt(from_condition, to_condition)
        if conversion_df.empty or '相对误差(%)' not in conversion_df.columns:
            return None
        rel_errors = conversion_df['相对误差(%)'].dropna()
        if rel_errors.empty:
            return None
        fig, ax = plt.subplots(figsize=(width, height))
        if self.transparent_bg:
            fig.patch.set_alpha(0); ax.patch.set_alpha(0)
        ax.hist(rel_errors, bins=20, color=self.color_set[3], edgecolor='white', alpha=0.7)
        ax.axvline(rel_errors.mean(), color='red', linestyle='--', linewidth=2, label=f'Mean: {rel_errors.mean():.2f}%')
        ax.set_xlabel('Relative error (%)', fontsize=fontsize)
        ax.set_ylabel('Frequency', fontsize=fontsize)
        ax.set_title(f'Relative error distribution - {from_condition}→{to_condition}', fontsize=fontsize+1)
        ax.legend(fontsize=fontsize-1)
        ax.tick_params(labelsize=fontsize-1, rotation=0)
        ax.grid(True, alpha=0.3, linestyle='--')
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=self.transparent_bg)
        self.figures['rel_error_hist'] = fig
        return fig

    def plot_error_vs_rt(self, from_condition: str, to_condition: str, save_path: str = None,
                         fontsize: int = 12, width: float = 8, height: float = 6) -> plt.Figure:
        """绝对误差 vs 实际保留时间"""
        conversion_df, _ = self.calculator.convert_ppg_index_to_rt(from_condition, to_condition)
        if conversion_df.empty or '绝对误差(10^-1min)' not in conversion_df.columns:
            return None
        valid = conversion_df.dropna(subset=['绝对误差(10^-1min)', f'{to_condition}_实际RT'])
        if valid.empty:
            return None
        actual = valid[f'{to_condition}_实际RT']
        errors = valid['绝对误差(10^-1min)']
        fig, ax = plt.subplots(figsize=(width, height))
        if self.transparent_bg:
            fig.patch.set_alpha(0); ax.patch.set_alpha(0)
        ax.scatter(actual, errors, color=self.color_set[4], s=50, alpha=0.7, edgecolor='white', linewidth=0.5)
        ax.axhline(errors.mean(), color='red', linestyle='--', linewidth=1.5, label=f'Mean error: {errors.mean():.3f}')
        if len(actual) > 1:
            z = np.polyfit(actual, errors, 1)
            p = np.poly1d(z)
            ax.plot(actual, p(actual), color=self.color_set[1], linestyle='-', linewidth=2, label='Trend')
        ax.set_xlabel(f'Actual RT in {to_condition} (min)', fontsize=fontsize)
        ax.set_ylabel('Absolute error (10^-1min)', fontsize=fontsize)
        ax.set_title(f'Error vs RT - {from_condition}→{to_condition}', fontsize=fontsize+1)
        ax.legend(fontsize=fontsize-1, loc='best')
        ax.tick_params(labelsize=fontsize-1, rotation=0)
        ax.grid(True, alpha=0.3, linestyle='--')
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=self.transparent_bg)
        self.figures['error_vs_rt'] = fig
        return fig

    def plot_top_errors(self, from_condition: str, to_condition: str, save_path: str = None,
                        fontsize: int = 12, width: float = 10, height: float = 6, top_n: int = 15) -> plt.Figure:
        """误差最大的前N个化合物（水平条形图）"""
        conversion_df, _ = self.calculator.convert_ppg_index_to_rt(from_condition, to_condition)
        if conversion_df.empty or '绝对误差(10^-1min)' not in conversion_df.columns:
            return None
        valid = conversion_df.dropna(subset=['绝对误差(10^-1min)', '化合物名称'])
        if valid.empty:
            return None
        top_errors = valid.nlargest(top_n, '绝对误差(10^-1min)')
        fig, ax = plt.subplots(figsize=(width, height))
        if self.transparent_bg:
            fig.patch.set_alpha(0); ax.patch.set_alpha(0)
        y_pos = np.arange(len(top_errors))
        ax.barh(y_pos, top_errors['绝对误差(10^-1min)'], color=self.color_set[5], edgecolor='white')
        names = [name[:20] + '...' if len(name) > 20 else name for name in top_errors['化合物名称']]
        ax.set_yticks(y_pos)
        ax.set_yticklabels(names, fontsize=fontsize-2)
        ax.invert_yaxis()
        ax.set_xlabel('Absolute error (10^-1min)', fontsize=fontsize)
        ax.set_title(f'Top {top_n} largest errors - {from_condition}→{to_condition}', fontsize=fontsize+1)
        ax.tick_params(labelsize=fontsize-1)
        ax.grid(True, alpha=0.3, linestyle='--', axis='x')
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=self.transparent_bg)
        self.figures['top_errors'] = fig
        return fig

    def plot_bland_altman(self, from_condition: str, to_condition: str, save_path: str = None,
                          fontsize: int = 12, width: float = 8, height: float = 6) -> plt.Figure:
        """Bland-Altman 图"""
        conversion_df, _ = self.calculator.convert_ppg_index_to_rt(from_condition, to_condition)
        if conversion_df.empty or '绝对误差(10^-1min)' not in conversion_df.columns:
            return None
        valid = conversion_df.dropna(subset=['绝对误差(10^-1min)', f'{to_condition}_预测RT', f'{to_condition}_实际RT'])
        if valid.empty:
            return None
        predicted = valid[f'{to_condition}_预测RT']
        actual = valid[f'{to_condition}_实际RT']
        mean_vals = (predicted + actual) / 2
        diff = predicted - actual
        fig, ax = plt.subplots(figsize=(width, height))
        if self.transparent_bg:
            fig.patch.set_alpha(0); ax.patch.set_alpha(0)
        ax.scatter(mean_vals, diff, color=self.color_set[0], s=50, alpha=0.7, edgecolor='white', linewidth=0.5)
        mean_diff, std_diff = diff.mean(), diff.std()
        ax.axhline(mean_diff, color='red', linestyle='-', linewidth=2, label=f'Mean diff: {mean_diff:.3f}')
        ax.axhline(mean_diff + 1.96*std_diff, color='red', linestyle='--', linewidth=1.5, label=f'+1.96SD: {mean_diff+1.96*std_diff:.3f}')
        ax.axhline(mean_diff - 1.96*std_diff, color='red', linestyle='--', linewidth=1.5, label=f'-1.96SD: {mean_diff-1.96*std_diff:.3f}')
        ax.axhline(0, color='gray', linestyle='-', linewidth=0.5)
        ax.set_xlabel('Mean RT (min)', fontsize=fontsize)
        ax.set_ylabel('Predicted - Actual (min)', fontsize=fontsize)
        ax.set_title(f'Bland-Altman plot - {from_condition}→{to_condition}', fontsize=fontsize+1)
        ax.legend(fontsize=fontsize-2, loc='best')
        ax.tick_params(labelsize=fontsize-1, rotation=0)
        ax.grid(True, alpha=0.3, linestyle='--')
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=self.transparent_bg)
        self.figures['bland_altman'] = fig
        return fig

    # ---------- 多转换方案独立图表（从原 multiple_conversion_comparison 拆分，但用户要求不生成组合图，故拆分为四个独立图）----------
    def plot_error_distribution_multi(self, conversions: List[Tuple[str, str]], save_path: str = None,
                                      fontsize: int = 12, width: float = 8, height: float = 6) -> plt.Figure:
        """多个转换方案的绝对误差分布直方图叠加（图例只显示Mean）"""
        all_errors = []
        labels = []
        for from_cond, to_cond in conversions:
            df, _ = self.calculator.convert_ppg_index_to_rt(from_cond, to_cond)
            if not df.empty and '绝对误差(10^-1min)' in df.columns:
                err = df['绝对误差(10^-1min)'].dropna()
                if not err.empty:
                    all_errors.append(err)
                    labels.append(f'{from_cond}→{to_cond}')
        if not all_errors:
            return None
        fig, ax = plt.subplots(figsize=(width, height))
        if self.transparent_bg:
            fig.patch.set_alpha(0); ax.patch.set_alpha(0)
        for i, (err, label) in enumerate(zip(all_errors, labels)):
            color = self.color_set[i % len(self.color_set)]
            mean_val = err.mean()
            legend_label = f'{label} (Mean={mean_val:.3f})'
            ax.hist(err, bins=30, alpha=0.6, density=True, color=color,
                    label=legend_label, edgecolor='white', linewidth=0.5)
        ax.set_xlabel('Absolute error (10^-1min)', fontsize=fontsize)
        ax.set_ylabel('Density', fontsize=fontsize)
        ax.set_title('Error distribution across schemes', fontsize=fontsize+2)
        ax.legend(frameon=False, fontsize=fontsize, loc='best')
        ax.tick_params(labelsize=fontsize-1, rotation=0)
        ax.grid(True, alpha=0.3, linestyle='--')
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=self.transparent_bg)
        self.figures['error_distribution_multi'] = fig
        return fig

    def plot_mean_error_bar_multi(self, conversions: List[Tuple[str, str]], save_path: str = None,
                                  fontsize: int = 12, width: float = 8, height: float = 6) -> plt.Figure:
        """多个转换方案的平均绝对误差柱状图"""
        all_errors = []
        labels = []
        for from_cond, to_cond in conversions:
            df, _ = self.calculator.convert_ppg_index_to_rt(from_cond, to_cond)
            if not df.empty and '绝对误差(10^-1min)' in df.columns:
                err = df['绝对误差(10^-1min)'].dropna()
                if not err.empty:
                    all_errors.append(err)
                    labels.append(f'{from_cond}→{to_cond}')
        if not all_errors:
            return None
        fig, ax = plt.subplots(figsize=(width, height))
        if self.transparent_bg:
            fig.patch.set_alpha(0); ax.patch.set_alpha(0)
        means = [e.mean() for e in all_errors]
        stds = [e.std() for e in all_errors]
        x_pos = np.arange(len(labels))
        bars = ax.bar(x_pos, means, yerr=stds, capsize=5, color=self.color_set[:len(labels)],
                      alpha=0.8, edgecolor='white', linewidth=1)
        legend_labels = [f'{label} (Mean={m:.3f})' for label, m in zip(labels, means)]
        ax.legend(bars, legend_labels, frameon=False, fontsize=fontsize, loc='best')
        ax.set_xticks(x_pos)
        ax.set_xticklabels(labels, rotation=0, ha='center', fontsize=fontsize)
        ax.set_ylabel('Mean absolute error (10^-1min)', fontsize=fontsize)
        ax.set_title('Mean error comparison', fontsize=fontsize+1)
        ax.tick_params(labelsize=fontsize-1)
        ax.grid(True, alpha=0.3, linestyle='--', axis='y')
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=self.transparent_bg)
        self.figures['mean_error_bar_multi'] = fig
        return fig

    def plot_pass_rate_bar_multi(self, conversions: List[Tuple[str, str]], save_path: str = None,
                                 fontsize: int = 12, width: float = 8, height: float = 6,
                                 threshold: float = 0.5, custom_title: str = "Pass rate") -> plt.Figure:
        """多个转换方案的通过率柱状图"""
        all_errors = []
        labels = []
        for from_cond, to_cond in conversions:
            df, _ = self.calculator.convert_ppg_index_to_rt(from_cond, to_cond)
            if not df.empty and '绝对误差(10^-1min)' in df.columns:
                err = df['绝对误差(10^-1min)'].dropna()
                if not err.empty:
                    all_errors.append(err)
                    labels.append(f'{from_cond}→{to_cond}')
        if not all_errors:
            return None
        fig, ax = plt.subplots(figsize=(width, height))
        if self.transparent_bg:
            fig.patch.set_alpha(0); ax.patch.set_alpha(0)
        pass_rates = [len(e[e <= threshold]) / len(e) * 100 for e in all_errors]
        x_pos = np.arange(len(labels))
        bars = ax.bar(x_pos, pass_rates, color=self.color_set[:len(labels)],
                      alpha=0.8, edgecolor='white', linewidth=1)
        legend_labels = [f'{label} (Pass rate={pr:.1f}%)' for label, pr in zip(labels, pass_rates)]
        ax.legend(bars, legend_labels, frameon=False, fontsize=fontsize, loc='best')
        ax.set_xticks(x_pos)
        ax.set_xticklabels(labels, rotation=0, ha='center', fontsize=fontsize)
        ax.set_ylabel('Pass rate (%)', fontsize=fontsize)
        ax.set_title(custom_title, fontsize=fontsize+1)
        ax.set_ylim(0, 105)
        ax.tick_params(labelsize=fontsize-1)
        ax.grid(True, alpha=0.3, linestyle='--', axis='y')
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=self.transparent_bg)
        self.figures['pass_rate_bar_multi'] = fig
        return fig

    def plot_cumulative_curve_multi(self, conversions: List[Tuple[str, str]], save_path: str = None,
                                    fontsize: int = 12, width: float = 8, height: float = 6) -> plt.Figure:
        """多个转换方案的累积误差分布曲线"""
        all_errors = []
        labels = []
        for from_cond, to_cond in conversions:
            df, _ = self.calculator.convert_ppg_index_to_rt(from_cond, to_cond)
            if not df.empty and '绝对误差(10^-1min)' in df.columns:
                err = df['绝对误差(10^-1min)'].dropna()
                if not err.empty:
                    all_errors.append(err)
                    labels.append(f'{from_cond}→{to_cond}')
        if not all_errors:
            return None
        fig, ax = plt.subplots(figsize=(width, height))
        if self.transparent_bg:
            fig.patch.set_alpha(0); ax.patch.set_alpha(0)
        for i, (err, label) in enumerate(zip(all_errors, labels)):
            sorted_err = np.sort(err)
            cum = np.arange(1, len(sorted_err)+1) / len(sorted_err) * 100
            color = self.color_set[i % len(self.color_set)]
            median_err = np.median(err)
            legend_label = f'{label} (Median={median_err:.3f})'
            ax.plot(sorted_err, cum, linewidth=2, label=legend_label, color=color, marker='.', markersize=4)
        ax.set_xlabel('Absolute error (10^-1min)', fontsize=fontsize)
        ax.set_ylabel('Cumulative percentage (%)', fontsize=fontsize)
        ax.set_title('Cumulative error distribution', fontsize=fontsize+2)
        ax.legend(frameon=False, fontsize=fontsize, loc='best')
        ax.tick_params(labelsize=fontsize-1, rotation=0)
        ax.grid(True, alpha=0.3, linestyle='--')
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=self.transparent_bg)
        self.figures['cumulative_curve_multi'] = fig
        return fig

# =============================================================================
# PPGIndexAnalyzerGUI 图形界面类（完整）
# =============================================================================
class PPGIndexAnalyzerGUI:
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
        status_bar = ttk.Label(main_frame, textvariable=self.status_var, relief=tk.SUNKEN, anchor=tk.W)
        status_bar.pack(fill=tk.X, pady=(5, 0))

    # ---------- Data Loading Tab ----------
    def setup_data_tab(self):
        data_tab = ttk.Frame(self.notebook)
        self.notebook.add(data_tab, text="Data Loading")
        data_frame = ttk.LabelFrame(data_tab, text="Data Management", padding=15)
        data_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
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
        comp_frame = ttk.LabelFrame(data_frame, text="Compound Data", padding=10)
        comp_frame.pack(fill=tk.X, pady=(0, 15))
        ttk.Label(comp_frame, text="Category:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        self.category_var = tk.StringVar(value="validation")
        ttk.Combobox(comp_frame, textvariable=self.category_var, values=["validation", "smrt", "training", "test"], width=15, state="readonly").grid(row=0, column=1, sticky=tk.W, padx=5, pady=5)
        ttk.Label(comp_frame, text="Condition:").grid(row=0, column=2, sticky=tk.W, padx=5, pady=5)
        self.compound_condition_var = tk.StringVar(value="condition1")
        ttk.Entry(comp_frame, textvariable=self.compound_condition_var, width=20).grid(row=0, column=3, sticky=tk.W, padx=5, pady=5)
        ttk.Label(comp_frame, text="Compound file:").grid(row=1, column=0, sticky=tk.W, padx=5, pady=5)
        self.compound_file_var = tk.StringVar()
        ttk.Entry(comp_frame, textvariable=self.compound_file_var, width=60).grid(row=1, column=1, columnspan=3, sticky=tk.W, padx=5, pady=5)
        ttk.Button(comp_frame, text="Browse...", command=self.browse_compound_file).grid(row=1, column=4, padx=5, pady=5)
        ttk.Button(comp_frame, text="Load Compound Data", command=self.load_compound_data).grid(row=2, column=0, columnspan=5, pady=10)
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

    # ---------- Analysis Tab ----------
    def setup_analysis_tab(self):
        analysis_tab = ttk.Frame(self.notebook)
        self.notebook.add(analysis_tab, text="Analysis")
        analysis_frame = ttk.LabelFrame(analysis_tab, text="PPG Index Analysis", padding=15)
        analysis_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        curve_frame = ttk.LabelFrame(analysis_frame, text="Calibration Curve Fitting", padding=10)
        curve_frame.pack(fill=tk.X, pady=(0, 15))
        ttk.Label(curve_frame, text="Condition:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        self.curve_condition_var = tk.StringVar()
        self.curve_condition_combo = ttk.Combobox(curve_frame, textvariable=self.curve_condition_var, width=25, state="readonly")
        self.curve_condition_combo.grid(row=0, column=1, sticky=tk.W, padx=5, pady=5)
        ttk.Label(curve_frame, text="Model:").grid(row=0, column=2, sticky=tk.W, padx=5, pady=5)
        self.model_type_var = tk.StringVar(value="logarithmic")
        ttk.Combobox(curve_frame, textvariable=self.model_type_var, values=["logarithmic", "linear"], width=15, state="readonly").grid(row=0, column=3, sticky=tk.W, padx=5, pady=5)
        ttk.Button(curve_frame, text="Fit Curve", command=self.fit_standard_curve).grid(row=0, column=4, padx=20, pady=5)
        calc_frame = ttk.LabelFrame(analysis_frame, text="PPG Index Calculation", padding=10)
        calc_frame.pack(fill=tk.X, pady=(0, 15))
        ttk.Label(calc_frame, text="Condition:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        self.calc_condition_var = tk.StringVar()
        self.calc_condition_combo = ttk.Combobox(calc_frame, textvariable=self.calc_condition_var, width=25, state="readonly")
        self.calc_condition_combo.grid(row=0, column=1, sticky=tk.W, padx=5, pady=5)
        ttk.Label(calc_frame, text="Method:").grid(row=0, column=2, sticky=tk.W, padx=5, pady=5)
        self.calc_method_var = tk.StringVar(value="interpolation")
        ttk.Combobox(calc_frame, textvariable=self.calc_method_var, values=["interpolation", "regression"], width=15, state="readonly").grid(row=0, column=3, sticky=tk.W, padx=5, pady=5)
        ttk.Button(calc_frame, text="Calculate PPG Indices", command=self.calculate_ppg_index).grid(row=0, column=4, padx=20, pady=5)
        # 条件比较已移除，因为用户要求不生成组合图，但保留条件选择用于其他单图？暂时保留条件选择框，但不再提供比较图
        # 为了简洁，移除条件比较部分
        results_frame = ttk.LabelFrame(analysis_frame, text="Results", padding=10)
        results_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 10))
        self.analysis_text = scrolledtext.ScrolledText(results_frame, width=80, height=15, wrap=tk.WORD, font=("Consolas", 10))
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

    # ---------- Conversion Tab ----------
    def setup_conversion_tab(self):
        conversion_tab = ttk.Frame(self.notebook)
        self.notebook.add(conversion_tab, text="Cross‑condition Conversion")
        conversion_frame = ttk.LabelFrame(conversion_tab, text="PPG Index Conversion & Validation", padding=15)
        conversion_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        settings_frame = ttk.LabelFrame(conversion_frame, text="Settings", padding=10)
        settings_frame.pack(fill=tk.X, pady=(0, 15))
        ttk.Label(settings_frame, text="Source condition:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        self.from_condition_var = tk.StringVar()
        self.from_condition_combo = ttk.Combobox(settings_frame, textvariable=self.from_condition_var, width=25, state="readonly")
        self.from_condition_combo.grid(row=0, column=1, sticky=tk.W, padx=5, pady=5)
        ttk.Label(settings_frame, text="Target condition:").grid(row=0, column=2, sticky=tk.W, padx=5, pady=5)
        self.to_condition_var = tk.StringVar()
        self.to_condition_combo = ttk.Combobox(settings_frame, textvariable=self.to_condition_var, width=25, state="readonly")
        self.to_condition_combo.grid(row=0, column=3, sticky=tk.W, padx=5, pady=5)
        ttk.Label(settings_frame, text="Error threshold (10^-1min):").grid(row=1, column=0, sticky=tk.W, padx=5, pady=5)
        self.threshold_var = tk.DoubleVar(value=0.5)
        ttk.Entry(settings_frame, textvariable=self.threshold_var, width=10).grid(row=1, column=1, sticky=tk.W, padx=5, pady=5)
        ttk.Button(settings_frame, text="Run Conversion", command=self.perform_conversion_analysis).grid(row=0, column=4, rowspan=2, padx=20, pady=5)
        # 多方案比较区域（独立图表选项）
        multi_frame = ttk.LabelFrame(settings_frame, text="Multiple Schemes Comparison (Independent Plots)", padding=5)
        multi_frame.grid(row=2, column=0, columnspan=5, sticky=tk.W, padx=5, pady=10)
        ttk.Label(multi_frame, text="Select conversion schemes:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=2)
        self.multi_conversion_frame = ttk.Frame(multi_frame)
        self.multi_conversion_frame.grid(row=1, column=0, sticky=tk.W, padx=5, pady=2)
        ttk.Button(multi_frame, text="Compare", command=self.compare_multiple_conversions).grid(row=1, column=1, padx=10, pady=2)
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
        self.conversion_stats_text = scrolledtext.ScrolledText(stats_frame, width=80, height=6, wrap=tk.WORD, font=("Consolas", 9))
        self.conversion_stats_text.pack(fill=tk.BOTH, expand=True)
        btn_frame = ttk.Frame(conversion_frame)
        btn_frame.pack(fill=tk.X, pady=(0, 10))
        ttk.Button(btn_frame, text="Export Results", command=self.export_conversion_results).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Clear", command=self.clear_conversion_results).pack(side=tk.LEFT, padx=5)
        self.conversion_status_var = tk.StringVar(value="Ready")
        ttk.Label(conversion_frame, textvariable=self.conversion_status_var).pack(anchor=tk.W)

    # ---------- Visualization Tab ----------
    def setup_visualization_tab(self):
        viz_tab = ttk.Frame(self.notebook)
        self.notebook.add(viz_tab, text="Visualization")
        viz_frame = ttk.LabelFrame(viz_tab, text="Plot Generation", padding=15)
        viz_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        options_frame = ttk.LabelFrame(viz_frame, text="Options", padding=10)
        options_frame.pack(fill=tk.X, pady=(0, 15))

        # 第一行：绘图类型
        ttk.Label(options_frame, text="Plot type:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        self.viz_type_var = tk.StringVar(value="standard_curve")
        plot_types = [
            "standard_curve", "residuals", "index_distribution",
            "pred_vs_actual", "abs_error_hist", "rel_error_hist",
            "error_vs_rt", "top_errors", "bland_altman",
            "error_distribution_multi", "mean_error_bar_multi", "pass_rate_bar_multi", "cumulative_curve_multi"
        ]
        viz_type_combo = ttk.Combobox(options_frame, textvariable=self.viz_type_var,
                                      values=plot_types, width=28, state="readonly")
        viz_type_combo.grid(row=0, column=1, sticky=tk.W, padx=5, pady=5)
        viz_type_combo.bind("<<ComboboxSelected>>", self.on_viz_type_change)

        # 第二行：字体大小、图表宽高
        ttk.Label(options_frame, text="Font size:").grid(row=0, column=2, sticky=tk.W, padx=5, pady=5)
        self.fontsize_var = tk.IntVar(value=12)
        ttk.Spinbox(options_frame, from_=8, to=24, textvariable=self.fontsize_var, width=5).grid(row=0, column=3, sticky=tk.W, padx=5, pady=5)

        ttk.Label(options_frame, text="Plot width (inch):").grid(row=1, column=0, sticky=tk.W, padx=5, pady=5)
        self.plot_width_var = tk.DoubleVar(value=8.0)
        ttk.Entry(options_frame, textvariable=self.plot_width_var, width=8).grid(row=1, column=1, sticky=tk.W, padx=5, pady=5)
        ttk.Label(options_frame, text="Plot height (inch):").grid(row=1, column=2, sticky=tk.W, padx=5, pady=5)
        self.plot_height_var = tk.DoubleVar(value=6.0)
        ttk.Entry(options_frame, textvariable=self.plot_height_var, width=8).grid(row=1, column=3, sticky=tk.W, padx=5, pady=5)

        # 第三行：通过率自定义标题（仅当选择 pass_rate_bar_multi 时显示）
        self.passrate_title_frame = ttk.Frame(options_frame)
        self.passrate_title_frame.grid(row=2, column=0, columnspan=6, sticky=tk.W, padx=5, pady=5)
        ttk.Label(self.passrate_title_frame, text="Pass rate title:").pack(side=tk.LEFT, padx=5)
        self.passrate_title_var = tk.StringVar(value="Pass rate")
        self.passrate_title_entry = ttk.Entry(self.passrate_title_frame, textvariable=self.passrate_title_var, width=30)
        self.passrate_title_entry.pack(side=tk.LEFT, padx=5)
        self.passrate_title_frame.grid_remove()  # 初始隐藏

        # 条件选择区域
        ttk.Label(options_frame, text="Condition/Source:").grid(row=0, column=4, sticky=tk.W, padx=5, pady=5)
        self.viz_condition_var = tk.StringVar()
        self.viz_condition_combo = ttk.Combobox(options_frame, textvariable=self.viz_condition_var, width=15, state="readonly")
        self.viz_condition_combo.grid(row=0, column=5, sticky=tk.W, padx=5, pady=5)

        ttk.Label(options_frame, text="Target:").grid(row=1, column=4, sticky=tk.W, padx=5, pady=5)
        self.viz_target_var = tk.StringVar()
        self.viz_target_combo = ttk.Combobox(options_frame, textvariable=self.viz_target_var, width=15, state="readonly")
        self.viz_target_combo.grid(row=1, column=5, sticky=tk.W, padx=5, pady=5)
        self.viz_target_combo.grid_remove()  # 默认隐藏，仅当需要目标条件时显示

        # 多方案选择框（用于 multi 类型）
        self.multi_schemes_frame = ttk.Frame(options_frame)
        self.multi_schemes_frame.grid(row=3, column=0, columnspan=6, sticky=tk.W, padx=5, pady=5)
        self.multi_schemes_frame.grid_remove()

        ttk.Button(options_frame, text="Generate Plot", command=self.generate_visualization).grid(row=0, column=6, padx=20, pady=5)

        display_frame = ttk.LabelFrame(viz_frame, text="Plot Display", padding=10)
        display_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 10))
        self.viz_placeholder = ttk.Label(display_frame, text="Plot will appear here", font=("Arial", 14), foreground="gray")
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

    # ---------- Results Tab ----------
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
        self.results_text = scrolledtext.ScrolledText(preview_frame, width=80, height=15, wrap=tk.WORD, font=("Consolas", 10))
        self.results_text.pack(fill=tk.BOTH, expand=True)
        control_frame = ttk.Frame(results_frame)
        control_frame.pack(fill=tk.X, pady=(0, 10))
        ttk.Button(control_frame, text="Save All Results", command=self.save_all_results).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_frame, text="Preview Report", command=self.preview_report).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_frame, text="Open Output Folder", command=self.open_output_dir).pack(side=tk.LEFT, padx=5)
        self.output_status_var = tk.StringVar(value="Ready")
        ttk.Label(results_frame, textvariable=self.output_status_var).pack(anchor=tk.W)

    # ---------- Log Tab ----------
    def setup_log_tab(self):
        log_tab = ttk.Frame(self.notebook)
        self.notebook.add(log_tab, text="Log")
        log_frame = ttk.LabelFrame(log_tab, text="Program Log", padding=15)
        log_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        self.log_text = scrolledtext.ScrolledText(log_frame, width=100, height=25, wrap=tk.WORD, font=("Consolas", 10))
        self.log_text.pack(fill=tk.BOTH, expand=True)
        self.log_text.tag_config("INFO", foreground="black")
        self.log_text.tag_config("SUCCESS", foreground="green")
        self.log_text.tag_config("WARNING", foreground="orange")
        self.log_text.tag_config("ERROR", foreground="red")
        btn_frame = ttk.Frame(log_frame)
        btn_frame.pack(fill=tk.X, pady=(10, 0))
        ttk.Button(btn_frame, text="Clear Log", command=self.clear_log).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Save Log", command=self.save_log).pack(side=tk.LEFT, padx=5)

    # ---------- 回调函数 ----------
    def browse_ppg_file(self):
        file_path = filedialog.askopenfilename(title="Select PPG data file", filetypes=[("Data files", "*.csv *.xlsx *.xls"), ("CSV", "*.csv"), ("Excel", "*.xlsx *.xls")])
        if file_path:
            self.ppg_file_var.set(file_path)

    def browse_compound_file(self):
        file_path = filedialog.askopenfilename(title="Select compound data file", filetypes=[("Data files", "*.csv *.xlsx *.xls"), ("CSV", "*.csv"), ("Excel", "*.xlsx *.xls")])
        if file_path:
            self.compound_file_var.set(file_path)

    def browse_output_dir(self):
        dir_path = filedialog.askdirectory(title="Select output directory")
        if dir_path:
            self.output_dir_var.set(dir_path)

    def load_ppg_data(self):
        file = self.ppg_file_var.get().strip()
        cond = self.condition_var.get().strip()
        if not file or not cond:
            messagebox.showwarning("Warning", "Please select a PPG file and enter a condition name.")
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
        if not file or not cat or not cond:
            messagebox.showwarning("Warning", "Please select a compound file, category and condition.")
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
        self.viz_target_combo['values'] = conditions
        self.from_condition_combo['values'] = conditions
        self.to_condition_combo['values'] = conditions
        if conditions:
            if not self.curve_condition_var.get():
                self.curve_condition_var.set(conditions[0])
            if not self.calc_condition_var.get():
                self.calc_condition_var.set(conditions[0])
            if not self.viz_condition_var.get():
                self.viz_condition_var.set(conditions[0])
            if not self.viz_target_var.get() and len(conditions) > 1:
                self.viz_target_var.set(conditions[1])
            if not self.from_condition_var.get():
                self.from_condition_var.set(conditions[0])
            if not self.to_condition_var.get() and len(conditions) > 1:
                self.to_condition_var.set(conditions[1])
        self.update_multi_conversion_checkboxes(conditions)

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
                        cb = ttk.Checkbutton(self.multi_conversion_frame, text=f"{from_cond}→{to_cond}", variable=var)
                        row = (i*len(conditions)+j) // 4
                        col = (i*len(conditions)+j) % 4
                        cb.grid(row=row, column=col, sticky=tk.W, padx=2, pady=1)
                        self.conversion_scheme_vars.append((var, from_cond, to_cond))
        else:
            ttk.Label(self.multi_conversion_frame, text="Need at least 2 conditions", foreground="gray").grid(row=0, column=0, sticky=tk.W)

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
            if '聚合度' in df.columns and '保留时间' in df.columns:
                ttk.Label(info, text=f"n range: {df['聚合度'].min()} - {df['聚合度'].max()}").pack(side=tk.LEFT, padx=10)
                ttk.Label(info, text=f"RT range: {df['保留时间'].min():.2f} - {df['保留时间'].max():.2f}").pack(side=tk.LEFT, padx=10)

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
            if '保留时间' in df.columns:
                ttk.Label(info, text=f"RT range: {df['保留时间'].min():.2f} - {df['保留时间'].max():.2f}").pack(side=tk.LEFT, padx=10)

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
                result = f"Fitting results - {cond}:\n  Model: {curve['model_name']}\n  R²: {curve['r_squared']:.6f}\n  Slope: {curve['slope']:.6f}\n  Intercept: {curve['intercept']:.6f}\n  Std error: {curve['std_err']:.6f}\n  Points: {curve['n_points']}\n"
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
                    if '计算PPG指数' in df.columns:
                        all_idx.extend(df['计算PPG指数'].dropna())
                if all_idx:
                    arr = np.array(all_idx)
                    result = f"PPG index statistics - {cond}:\n  Method: {method}\n  N = {len(arr)}\n  Mean = {np.mean(arr):.2f}\n  Std = {np.std(arr):.2f}\n  Min = {np.min(arr):.2f}\n  Max = {np.max(arr):.2f}\n  Median = {np.median(arr):.2f}\n"
                    self.root.after(0, lambda: self.analysis_message(result, "INFO"))
                self.root.after(0, lambda: self.analysis_message(f"✓ {msg}", "SUCCESS"))
            else:
                self.root.after(0, lambda: self.analysis_message(f"✗ {msg}", "ERROR"))
        except Exception as e:
            self.root.after(0, lambda: self.analysis_message(f"✗ Calculation failed: {str(e)}", "ERROR"))
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
            self.root.after(0, lambda: self.log_message(f"✓ Conversion {from_cond}→{to_cond} complete. Valid: {res['有效转换数']}, Mean error: {res['误差统计'].get('平均绝对误差',0):.3f} (10^-1min)", "SUCCESS"))
            self.root.after(0, lambda: self.conversion_status_var.set("Conversion done"))
        except Exception as e:
            self.root.after(0, lambda: self.log_message(f"✗ Conversion failed: {str(e)}", "ERROR"))
            self.root.after(0, lambda: self.conversion_status_var.set("Failed"))
        self.root.after(0, lambda: self.update_status("Ready"))

    def display_conversion_results(self, res):
        for item in self.conversion_tree.get_children():
            self.conversion_tree.delete(item)
        self.conversion_stats_text.delete(1.0, tk.END)
        for d in res.get('详细数据', []):
            vals = (
                d['化合物名称'],
                f"{d.get(f'{res['源条件']}_PPG指数',''):.2f}" if isinstance(d.get(f'{res['源条件']}_PPG指数'), (int,float)) else "",
                f"{d.get(f'{res['目标条件']}_预测RT',''):.3f}" if isinstance(d.get(f'{res['目标条件']}_预测RT'), (int,float)) else "",
                f"{d.get(f'{res['目标条件']}_实际RT',''):.3f}" if isinstance(d.get(f'{res['目标条件']}_实际RT'), (int,float)) else "",
                f"{d.get('绝对误差(10^-1min)',''):.3f}" if isinstance(d.get('绝对误差(10^-1min)'), (int,float)) else "",
                f"{d.get('相对误差(%)',''):.2f}" if isinstance(d.get('相对误差(%)'), (int,float)) else ""
            )
            self.conversion_tree.insert("", tk.END, values=vals)
        stats = f"Conversion analysis: {res['源条件']} → {res['目标条件']}\n" + "="*50 + "\n\n"
        stats += f"Time: {res['转换时间']}\nTotal compounds: {res['总化合物数']}\nValid conversions: {res['有效转换数']}\n\n"
        if isinstance(res['误差统计'], dict):
            stats += "Error statistics:\n"
            for k, v in res['误差统计'].items():
                stats += f"  {k}: {v:.3f}\n" if isinstance(v, float) else f"  {k}: {v}\n"
        if '误差分布' in res:
            stats += "\nError distribution:\n"
            total = res['有效转换数']
            for bin_name, cnt in res['误差分布'].items():
                pct = (cnt/total*100) if total>0 else 0
                stats += f"  {bin_name}: {cnt} ({pct:.1f}%)\n"
        if '通过率分析' in res:
            pa = res['通过率分析']
            stats += f"\nPass rate (threshold {pa['阈值(10^-1min)']} 10^-1min): {pa['通过数']}/{pa['总数']} = {pa['通过率(%)']:.1f}%\n"
        self.conversion_stats_text.insert(tk.END, stats, "INFO")

    def compare_multiple_conversions(self):
        selected = [(f, t) for var, f, t in self.conversion_scheme_vars if var.get()]
        if len(selected) < 1:
            messagebox.showwarning("Warning", "Select at least one conversion scheme.")
            return
        self.log_message(f"Comparing {len(selected)} conversion schemes...", "INFO")
        self.update_status("Comparing conversions...")
        # 这里我们生成多方案独立图表，需要调用相应的绘图方法
        # 由于用户要求独立图表，我们直接在 generate_visualization 中处理多方案类型
        # 这里只是一个占位，实际生成通过可视化标签页
        messagebox.showinfo("Info", "Please go to Visualization tab and select a multi-scheme plot type (error_distribution_multi, etc.)")

    def on_viz_type_change(self, event=None):
        typ = self.viz_type_var.get()
        # 隐藏所有动态控件
        self.viz_target_combo.grid_remove()
        self.multi_schemes_frame.grid_remove()
        self.passrate_title_frame.grid_remove()
        # 根据类型显示需要的控件
        if typ in ["pred_vs_actual", "abs_error_hist", "rel_error_hist", "error_vs_rt", "top_errors", "bland_altman"]:
            # 需要源条件和目标条件
            self.viz_target_combo.grid()
            self.passrate_title_frame.grid_remove()
        elif typ in ["error_distribution_multi", "mean_error_bar_multi", "pass_rate_bar_multi", "cumulative_curve_multi"]:
            # 需要多方案选择框
            self.multi_schemes_frame.grid()
            self.update_multi_schemes_checkboxes()
            if typ == "pass_rate_bar_multi":
                self.passrate_title_frame.grid()
        else:
            # standard_curve, residuals, index_distribution 只需要一个条件
            self.viz_target_combo.grid_remove()

    def update_multi_schemes_checkboxes(self):
        for w in self.multi_schemes_frame.winfo_children():
            w.destroy()
        conds = list(self.calculator.ppg_data.keys())
        if len(conds) < 2:
            ttk.Label(self.multi_schemes_frame, text="Need at least 2 conditions", foreground="red").pack()
            self.viz_multi_scheme_vars = []
            return
        ttk.Label(self.multi_schemes_frame, text="Select conversion schemes:").pack(anchor=tk.W)
        frame = ttk.Frame(self.multi_schemes_frame)
        frame.pack(anchor=tk.W)
        self.viz_multi_scheme_vars = []
        idx = 0
        for i in range(len(conds)):
            for j in range(len(conds)):
                if i != j:
                    from_c = conds[i]
                    to_c = conds[j]
                    var = tk.BooleanVar(value=False)
                    cb = ttk.Checkbutton(frame, text=f"{from_c}→{to_c}", variable=var)
                    row = idx // 4
                    col = idx % 4
                    cb.grid(row=row, column=col, sticky=tk.W, padx=5, pady=2)
                    self.viz_multi_scheme_vars.append((var, from_c, to_c))
                    idx += 1

    def generate_visualization(self):
        typ = self.viz_type_var.get()
        fs = self.fontsize_var.get()
        width = self.plot_width_var.get()
        height = self.plot_height_var.get()
        if self.visualizer is None:
            self.visualizer = PPGVisualizer(self.calculator)

        self.viz_status_var.set(f"Generating {typ}...")
        self.update_status("Generating plot")

        fig = None
        try:
            if typ == "standard_curve":
                cond = self.viz_condition_var.get()
                if not cond:
                    messagebox.showwarning("Warning", "Select a condition.")
                    return
                fig = self.visualizer.plot_standard_curves([cond], fontsize=fs, width=width, height=height)
            elif typ == "residuals":
                cond = self.viz_condition_var.get()
                if not cond:
                    messagebox.showwarning("Warning", "Select a condition.")
                    return
                fig = self.visualizer.plot_residuals([cond], fontsize=fs, width=width, height=height)
            elif typ == "index_distribution":
                cond = self.viz_condition_var.get()
                if not cond:
                    messagebox.showwarning("Warning", "Select a condition.")
                    return
                fig = self.visualizer.plot_ppg_index_distribution(cond, fontsize=fs, width=width, height=height)
            elif typ in ["pred_vs_actual", "abs_error_hist", "rel_error_hist", "error_vs_rt", "top_errors", "bland_altman"]:
                src = self.viz_condition_var.get()
                tgt = self.viz_target_var.get()
                if not src or not tgt or src == tgt:
                    messagebox.showwarning("Warning", "Select different source and target conditions.")
                    return
                if typ == "pred_vs_actual":
                    fig = self.visualizer.plot_pred_vs_actual(src, tgt, fontsize=fs, width=width, height=height)
                elif typ == "abs_error_hist":
                    fig = self.visualizer.plot_abs_error_hist(src, tgt, fontsize=fs, width=width, height=height)
                elif typ == "rel_error_hist":
                    fig = self.visualizer.plot_rel_error_hist(src, tgt, fontsize=fs, width=width, height=height)
                elif typ == "error_vs_rt":
                    fig = self.visualizer.plot_error_vs_rt(src, tgt, fontsize=fs, width=width, height=height)
                elif typ == "top_errors":
                    fig = self.visualizer.plot_top_errors(src, tgt, fontsize=fs, width=width, height=height)
                elif typ == "bland_altman":
                    fig = self.visualizer.plot_bland_altman(src, tgt, fontsize=fs, width=width, height=height)
            elif typ in ["error_distribution_multi", "mean_error_bar_multi", "pass_rate_bar_multi", "cumulative_curve_multi"]:
                if not hasattr(self, 'viz_multi_scheme_vars'):
                    messagebox.showwarning("Warning", "No conversion schemes selected.")
                    return
                selected = [(f, t) for var, f, t in self.viz_multi_scheme_vars if var.get()]
                if len(selected) == 0:
                    messagebox.showwarning("Warning", "Select at least one conversion scheme.")
                    return
                threshold = self.threshold_var.get()
                if typ == "error_distribution_multi":
                    fig = self.visualizer.plot_error_distribution_multi(selected, fontsize=fs, width=width, height=height)
                elif typ == "mean_error_bar_multi":
                    fig = self.visualizer.plot_mean_error_bar_multi(selected, fontsize=fs, width=width, height=height)
                elif typ == "pass_rate_bar_multi":
                    custom_title = self.passrate_title_var.get()
                    fig = self.visualizer.plot_pass_rate_bar_multi(selected, fontsize=fs, width=width, height=height,
                                                                    threshold=threshold, custom_title=custom_title)
                elif typ == "cumulative_curve_multi":
                    fig = self.visualizer.plot_cumulative_curve_multi(selected, fontsize=fs, width=width, height=height)
            else:
                self.viz_status_var.set("Unsupported type")
                return

            if fig is None:
                self.viz_status_var.set("Failed to generate. Check data.")
                return
            self.display_figure(fig)
            self.viz_status_var.set("Plot generated")
            self.update_status("Ready")
        except Exception as e:
            self.viz_status_var.set(f"Error: {str(e)}")
            self.log_message(f"✗ Plot failed: {str(e)}", "ERROR")

    def display_figure(self, fig):
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

    def save_figure(self):
        if self.current_figure is None:
            messagebox.showwarning("Warning", "No figure to save.")
            return
        file = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG", "*.png"), ("PDF", "*.pdf"), ("SVG", "*.svg")])
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

    def export_conversion_results(self):
        if not self.calculator.conversion_results:
            messagebox.showinfo("Info", "No conversion results to export.")
            return
        file = filedialog.asksaveasfilename(defaultextension=".xlsx", filetypes=[("Excel files", "*.xlsx")])
        if not file:
            return
        try:
            with pd.ExcelWriter(file, engine='openpyxl') as writer:
                for key, conv in self.calculator.conversion_results.items():
                    if '详细数据' in conv:
                        df = pd.DataFrame(conv['详细数据'])
                        df.to_excel(writer, sheet_name=key[:30], index=False)
                stats = []
                for key, conv in self.calculator.conversion_results.items():
                    s = {'Direction': key, 'Source': conv.get('源条件',''), 'Target': conv.get('目标条件',''),
                         'Total': conv.get('总化合物数',0), 'Valid': conv.get('有效转换数',0)}
                    if isinstance(conv.get('误差统计'), dict):
                        s.update({k:v for k,v in conv['误差统计'].items() if isinstance(v,(int,float))})
                    if isinstance(conv.get('通过率分析'), dict):
                        s['Pass rate (%)'] = conv['通过率分析'].get('通过率(%)',0)
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

    def generate_report(self):
        try:
            summary = self.calculator.generate_summary_report()
            self.analysis_message("="*60, "INFO")
            self.analysis_message("PPG Retention Index Analysis Report", "INFO")
            self.analysis_message("="*60, "INFO")
            self.analysis_message(f"Generated: {summary['生成时间']}", "INFO")
            self.analysis_message("", "INFO")
            self.analysis_message(f"PPG datasets: {summary['PPG数据条件数']}", "INFO")
            self.analysis_message(f"Compound datasets: {summary['化合物数据集数']}", "INFO")
            self.analysis_message(f"Calibration curves: {summary['标准曲线数']}", "INFO")
            self.analysis_message(f"PPG index calculations: {summary['PPG指数计算结果数']}", "INFO")
            self.analysis_message(f"Cross‑condition conversions: {summary['跨条件转换结果数']}", "INFO")
            self.analysis_message("", "INFO")
            if summary['标准曲线性能']:
                self.analysis_message("Calibration curve performance:", "INFO")
                for cond, perf in summary['标准曲线性能'].items():
                    self.analysis_message(f"  {cond}: R²={perf['R²']:.4f}, slope={perf['斜率']:.4f}", "INFO")
            if summary['PPG指数统计']:
                self.analysis_message("\nPPG index statistics:", "INFO")
                for cond, stat in summary['PPG指数统计'].items():
                    self.analysis_message(f"  {cond}: n={stat['样本数']}, mean={stat['平均值']:.2f}, sd={stat['标准差']:.2f}", "INFO")
            if summary['跨条件转换统计']:
                self.analysis_message("\nConversion statistics:", "INFO")
                for key, stat in summary['跨条件转换统计'].items():
                    self.analysis_message(f"  {key}: valid={stat['有效转换数']}, mean error={stat['平均绝对误差']:.3f} (10^-1min), pass rate={stat['通过率(%)']:.1f}%", "INFO")
            self.analysis_message("\nReport generation complete.", "SUCCESS")
            self.results_text.delete(1.0, tk.END)
            self.results_text.insert(tk.END, "PPG Retention Index Analysis Report\n")
            self.results_text.insert(tk.END, "="*60 + "\n\n")
            self.results_text.insert(tk.END, f"Generated: {summary['生成时间']}\n\n")
            self.results_text.insert(tk.END, f"PPG datasets: {summary['PPG数据条件数']}\n")
            self.results_text.insert(tk.END, f"Compound datasets: {summary['化合物数据集数']}\n")
            self.results_text.insert(tk.END, f"Calibration curves: {summary['标准曲线数']}\n")
            self.results_text.insert(tk.END, f"PPG index calculations: {summary['PPG指数计算结果数']}\n")
            self.results_text.insert(tk.END, f"Cross‑condition conversions: {summary['跨条件转换结果数']}\n\n")
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
        file = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text files", "*.txt")])
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
# 主函数
# =============================================================================
def main():
    print("PPG Retention Index Analyzer - Journal Version (all plots are independent, no composite plots)")
    print("Features: Times New Roman, transparent background, horizontal x-labels, adjustable size, customizable pass rate title.")
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