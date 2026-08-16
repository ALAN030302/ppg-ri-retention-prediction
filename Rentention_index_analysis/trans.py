#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PPG保留指数文件转换工具

功能：
1. 加载PPG标准品数据
2. 加载包含PPG保留指数的文件
3. 批量将保留指数转换为保留时间
4. 保存转换结果
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
    """保留指数转保留时间转换器"""

    def __init__(self):
        self.ppg_data = {}  # PPG标准品数据
        self.standard_curves = {}  # 标准曲线参数
        self.index_data = None  # 保留指数数据
        self.converted_data = None  # 转换后的数据

    def load_ppg_data(self, file_path, condition="default"):
        """加载PPG标准品数据"""
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

            return True, f"成功加载 {len(df)} 个PPG标准品数据"

        except Exception as e:
            return False, f"加载PPG数据失败: {str(e)}"

    def load_index_data(self, file_path):
        """加载保留指数数据"""
        try:
            file_ext = Path(file_path).suffix.lower()

            if file_ext in ['.xlsx', '.xls']:
                self.index_data = pd.read_excel(file_path)
            elif file_ext == '.csv':
                self.index_data = pd.read_csv(file_path)
            else:
                return False, f"不支持的文件格式: {file_ext}"

            # 查找可能的PPG指数列
            index_columns = []
            for col in self.index_data.columns:
                col_lower = col.lower()
                if any(keyword in col_lower for keyword in ['ppg', '指数', 'index', '保留指数']):
                    index_columns.append(col)

            if not index_columns:
                # 如果没有明确的PPG指数列，让用户选择
                return True, "未检测到PPG指数列，请在转换时手动选择"

            return True, f"成功加载 {len(self.index_data)} 行数据，检测到PPG指数列: {', '.join(index_columns)}"

        except Exception as e:
            return False, f"加载保留指数数据失败: {str(e)}"

    def fit_standard_curve(self, condition="default", model_type="logarithmic"):
        """拟合标准曲线"""
        try:
            if condition not in self.ppg_data:
                return False, f"未找到条件 {condition} 的PPG数据"

            df = self.ppg_data[condition]
            n_values = df['聚合度'].values
            rt_values = df['保留时间'].values

            if model_type == "logarithmic":
                # 对数模型: RT = a + b * ln(n)
                x = np.log(n_values)
                model_name = "对数模型 (RT = a + b * ln(n))"
            elif model_type == "linear":
                # 线性模型: RT = a + b * n
                x = n_values
                model_name = "线性模型 (RT = a + b * n)"
            else:
                return False, f"不支持的模型类型: {model_type}"

            # 线性回归
            slope, intercept, r_value, p_value, std_err = stats.linregress(x, rt_values)

            # 存储标准曲线参数
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

            return True, f"标准曲线拟合成功: {model_name}, R² = {r_value ** 2:.6f}"

        except Exception as e:
            return False, f"拟合标准曲线失败: {str(e)}"

    def convert_index_to_rt(self, index, condition="default", method="regression"):
        """单个保留指数转保留时间"""
        try:
            if condition not in self.ppg_data:
                return None, f"未找到条件 {condition} 的PPG数据"

            if condition not in self.standard_curves:
                success, msg = self.fit_standard_curve(condition)
                if not success:
                    return None, f"无法转换: {msg}"

            df_ppg = self.ppg_data[condition]
            ppg_rt = df_ppg['保留时间'].values
            ppg_n = df_ppg['聚合度'].values

            # 将PPG指数转换为聚合度
            n_value = index / 100

            if method == "interpolation":
                # 线性插值法
                if n_value < ppg_n[0]:
                    # 外推
                    if len(ppg_n) >= 2:
                        rt_calc = ppg_rt[0] - (ppg_n[0] - n_value) / (ppg_n[1] - ppg_n[0]) * (ppg_rt[1] - ppg_rt[0])
                    else:
                        rt_calc = ppg_rt[0]
                elif n_value > ppg_n[-1]:
                    # 外推
                    if len(ppg_n) >= 2:
                        rt_calc = ppg_rt[-1] + (n_value - ppg_n[-1]) / (ppg_n[-1] - ppg_n[-2]) * (
                                ppg_rt[-1] - ppg_rt[-2])
                    else:
                        rt_calc = ppg_rt[-1]
                else:
                    # 线性插值
                    idx = np.searchsorted(ppg_n, n_value) - 1
                    if idx < 0:
                        idx = 0
                    elif idx >= len(ppg_n) - 1:
                        idx = len(ppg_n) - 2

                    n_i, n_j = ppg_n[idx], ppg_n[idx + 1]
                    rt_i, rt_j = ppg_rt[idx], ppg_rt[idx + 1]

                    rt_calc = rt_i + (rt_j - rt_i) * (n_value - n_i) / (n_j - n_i)

                method_used = "线性插值"

            elif method == "regression":
                # 回归法
                curve = self.standard_curves[condition]

                if curve['model_type'] == "logarithmic":
                    # RT = a + b * ln(n)
                    rt_calc = curve['intercept'] + curve['slope'] * np.log(n_value)
                else:  # linear
                    # RT = a + b * n
                    rt_calc = curve['intercept'] + curve['slope'] * n_value

                method_used = "回归法"
            else:
                return None, f"不支持的计算方法: {method}"

            return rt_calc, method_used

        except Exception as e:
            return None, f"转换失败: {str(e)}"

    def batch_convert_file(self, index_column, condition="default", method="regression"):
        """批量转换文件中的保留指数"""
        try:
            if self.index_data is None:
                return False, "未加载保留指数数据"

            if index_column not in self.index_data.columns:
                return False, f"列 '{index_column}' 不存在于数据中"

            # 创建副本
            self.converted_data = self.index_data.copy()

            # 准备结果列
            rt_column_name = f"保留时间_{condition}"
            method_column_name = f"计算方法_{condition}"

            # 转换每一行
            rt_values = []
            method_values = []

            for idx, value in enumerate(self.index_data[index_column]):
                if pd.isna(value):
                    rt_values.append(np.nan)
                    method_values.append("数据缺失")
                else:
                    try:
                        rt, method_used = self.convert_index_to_rt(float(value), condition, method)
                        rt_values.append(rt)
                        method_values.append(method_used)
                    except Exception as e:
                        rt_values.append(np.nan)
                        method_values.append(f"转换失败: {str(e)}")

            # 添加结果列
            self.converted_data[rt_column_name] = rt_values
            self.converted_data[method_column_name] = method_values

            # 统计信息
            successful_conversions = sum(1 for m in method_values if m in ["线性插值", "回归法"])
            total_conversions = len(self.index_data)

            return True, f"转换完成: {successful_conversions}/{total_conversions} 个数据转换成功"

        except Exception as e:
            return False, f"批量转换失败: {str(e)}"

    def save_converted_data(self, file_path):
        """保存转换后的数据"""
        try:
            if self.converted_data is None:
                return False, "没有可保存的转换数据"

            file_ext = Path(file_path).suffix.lower()

            if file_ext in ['.xlsx', '.xls']:
                self.converted_data.to_excel(file_path, index=False)
            elif file_ext == '.csv':
                self.converted_data.to_csv(file_path, index=False, encoding='utf-8')
            else:
                return False, f"不支持的文件格式: {file_ext}"

            return True, f"数据已保存到: {file_path}"

        except Exception as e:
            return False, f"保存数据失败: {str(e)}"


class IndexFileConverterGUI:
    """保留指数文件转换GUI"""

    def __init__(self, root):
        self.root = root
        self.root.title("PPG保留指数文件转换工具")
        self.root.geometry("1200x800")

        # 初始化转换器
        self.converter = IndexToRTConverter()

        # 创建UI
        self.setup_ui()

    def setup_ui(self):
        """设置UI界面"""
        # 创建主框架
        main_frame = ttk.Frame(self.root, padding=10)
        main_frame.pack(fill=tk.BOTH, expand=True)

        # 创建Notebook（标签页）
        self.notebook = ttk.Notebook(main_frame)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        # 创建标签页
        self.setup_data_tab()
        self.setup_conversion_tab()
        self.setup_result_tab()

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
        data_frame = ttk.LabelFrame(data_tab, text="数据加载", padding=15)
        data_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # ==================== PPG数据加载 ====================
        ppg_frame = ttk.LabelFrame(data_frame, text="1. 加载PPG标准品数据", padding=10)
        ppg_frame.pack(fill=tk.X, pady=(0, 15))

        # 文件选择
        file_frame1 = ttk.Frame(ppg_frame)
        file_frame1.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(file_frame1, text="PPG数据文件:").pack(side=tk.LEFT, padx=(0, 5))

        self.ppg_file_var = tk.StringVar()
        ppg_file_entry = ttk.Entry(file_frame1, textvariable=self.ppg_file_var, width=60)
        ppg_file_entry.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 10))

        ttk.Button(file_frame1, text="浏览...", command=self.browse_ppg_file).pack(side=tk.LEFT)

        # 条件设置
        condition_frame = ttk.Frame(ppg_frame)
        condition_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(condition_frame, text="条件名称:").pack(side=tk.LEFT, padx=(0, 5))

        self.ppg_condition_var = tk.StringVar(value="default")
        ppg_condition_entry = ttk.Entry(condition_frame, textvariable=self.ppg_condition_var, width=20)
        ppg_condition_entry.pack(side=tk.LEFT, padx=(0, 20))

        # 加载按钮
        ttk.Button(ppg_frame, text="加载PPG数据", command=self.load_ppg_data).pack(pady=10)

        # PPG数据信息
        self.ppg_info_var = tk.StringVar(value="等待加载PPG数据...")
        ttk.Label(ppg_frame, textvariable=self.ppg_info_var, wraplength=800).pack(anchor=tk.W)

        # ==================== 保留指数数据加载 ====================
        index_frame = ttk.LabelFrame(data_frame, text="2. 加载保留指数数据", padding=10)
        index_frame.pack(fill=tk.X, pady=(0, 15))

        # 文件选择
        file_frame2 = ttk.Frame(index_frame)
        file_frame2.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(file_frame2, text="保留指数文件:").pack(side=tk.LEFT, padx=(0, 5))

        self.index_file_var = tk.StringVar()
        index_file_entry = ttk.Entry(file_frame2, textvariable=self.index_file_var, width=60)
        index_file_entry.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 10))

        ttk.Button(file_frame2, text="浏览...", command=self.browse_index_file).pack(side=tk.LEFT)

        # 加载按钮
        ttk.Button(index_frame, text="加载保留指数数据", command=self.load_index_data).pack(pady=10)

        # 保留指数数据信息
        self.index_info_var = tk.StringVar(value="等待加载保留指数数据...")
        ttk.Label(index_frame, textvariable=self.index_info_var, wraplength=800).pack(anchor=tk.W)

        # ==================== 数据预览 ====================
        preview_frame = ttk.LabelFrame(data_frame, text="数据预览", padding=10)
        preview_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 10))

        # 创建文本显示区域
        self.preview_text = scrolledtext.ScrolledText(preview_frame, width=100, height=15,
                                                      wrap=tk.WORD, font=("Consolas", 9))
        self.preview_text.pack(fill=tk.BOTH, expand=True)

        # 预览控制
        preview_controls = ttk.Frame(preview_frame)
        preview_controls.pack(fill=tk.X, pady=(5, 0))

        ttk.Button(preview_controls, text="预览PPG数据",
                   command=self.preview_ppg_data).pack(side=tk.LEFT, padx=5)

        ttk.Button(preview_controls, text="预览保留指数数据",
                   command=self.preview_index_data).pack(side=tk.LEFT, padx=5)

        ttk.Button(preview_controls, text="清空预览",
                   command=self.clear_preview).pack(side=tk.LEFT, padx=5)

    def setup_conversion_tab(self):
        """设置转换标签页"""
        conversion_tab = ttk.Frame(self.notebook)
        self.notebook.add(conversion_tab, text="转换设置")

        # 主框架
        conversion_frame = ttk.LabelFrame(conversion_tab, text="转换设置", padding=15)
        conversion_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # ==================== 标准曲线设置 ====================
        curve_frame = ttk.LabelFrame(conversion_frame, text="1. PPG标准曲线设置", padding=10)
        curve_frame.pack(fill=tk.X, pady=(0, 15))

        # 条件选择
        condition_frame = ttk.Frame(curve_frame)
        condition_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(condition_frame, text="选择条件:").pack(side=tk.LEFT, padx=(0, 5))

        self.curve_condition_var = tk.StringVar()
        self.curve_condition_combo = ttk.Combobox(condition_frame,
                                                  textvariable=self.curve_condition_var,
                                                  width=25, state="readonly")
        self.curve_condition_combo.pack(side=tk.LEFT, padx=(0, 20))

        # 模型选择
        ttk.Label(condition_frame, text="模型类型:").pack(side=tk.LEFT, padx=(0, 5))

        self.curve_model_var = tk.StringVar(value="logarithmic")
        curve_model_combo = ttk.Combobox(condition_frame, textvariable=self.curve_model_var,
                                         values=["logarithmic", "linear"],
                                         width=15, state="readonly")
        curve_model_combo.pack(side=tk.LEFT, padx=(0, 20))

        # 拟合按钮
        ttk.Button(curve_frame, text="拟合标准曲线", command=self.fit_standard_curve).pack(pady=10)

        # 曲线结果显示
        self.curve_result_var = tk.StringVar(value="")
        ttk.Label(curve_frame, textvariable=self.curve_result_var,
                  font=("Arial", 10), wraplength=800).pack(anchor=tk.W, pady=(0, 10))

        # ==================== 转换设置 ====================
        settings_frame = ttk.LabelFrame(conversion_frame, text="2. 转换设置", padding=10)
        settings_frame.pack(fill=tk.X, pady=(0, 15))

        # 保留指数列选择
        column_frame = ttk.Frame(settings_frame)
        column_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(column_frame, text="保留指数列:").pack(side=tk.LEFT, padx=(0, 5))

        self.index_column_var = tk.StringVar()
        self.index_column_combo = ttk.Combobox(column_frame,
                                               textvariable=self.index_column_var,
                                               width=30, state="readonly")
        self.index_column_combo.pack(side=tk.LEFT, padx=(0, 20))

        # 刷新列按钮
        ttk.Button(column_frame, text="刷新列列表",
                   command=self.refresh_columns).pack(side=tk.LEFT)

        # 转换条件选择
        conv_condition_frame = ttk.Frame(settings_frame)
        conv_condition_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(conv_condition_frame, text="转换条件:").pack(side=tk.LEFT, padx=(0, 5))

        self.conv_condition_var = tk.StringVar()
        self.conv_condition_combo = ttk.Combobox(conv_condition_frame,
                                                 textvariable=self.conv_condition_var,
                                                 width=25, state="readonly")
        self.conv_condition_combo.pack(side=tk.LEFT, padx=(0, 20))

        # 转换方法选择
        ttk.Label(conv_condition_frame, text="计算方法:").pack(side=tk.LEFT, padx=(0, 5))

        self.conv_method_var = tk.StringVar(value="regression")
        conv_method_combo = ttk.Combobox(conv_condition_frame, textvariable=self.conv_method_var,
                                         values=["interpolation", "regression"],
                                         width=15, state="readonly")
        conv_method_combo.pack(side=tk.LEFT, padx=(0, 20))

        # ==================== 转换执行 ====================
        execute_frame = ttk.LabelFrame(conversion_frame, text="3. 执行转换", padding=10)
        execute_frame.pack(fill=tk.X, pady=(0, 15))

        # 测试转换
        test_frame = ttk.Frame(execute_frame)
        test_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(test_frame, text="测试转换 (输入PPG指数):").pack(side=tk.LEFT, padx=(0, 5))

        self.test_index_var = tk.StringVar(value="550")
        test_index_entry = ttk.Entry(test_frame, textvariable=self.test_index_var, width=15)
        test_index_entry.pack(side=tk.LEFT, padx=(0, 10))

        ttk.Button(test_frame, text="测试转换", command=self.test_conversion).pack(side=tk.LEFT)

        # 测试结果
        self.test_result_var = tk.StringVar(value="")
        ttk.Label(execute_frame, textvariable=self.test_result_var,
                  font=("Arial", 10), wraplength=800).pack(anchor=tk.W, pady=(0, 10))

        # 批量转换按钮
        ttk.Button(execute_frame, text="开始批量转换",
                   command=self.batch_convert, style="Accent.TButton").pack(pady=10)

        # 转换状态
        self.conversion_status_var = tk.StringVar(value="等待转换...")
        ttk.Label(execute_frame, textvariable=self.conversion_status_var,
                  font=("Arial", 10, "bold"), wraplength=800).pack(anchor=tk.W)

        # 创建强调样式
        style = ttk.Style()
        style.configure("Accent.TButton", font=("Arial", 10, "bold"), foreground="blue")

    def setup_result_tab(self):
        """设置结果标签页"""
        result_tab = ttk.Frame(self.notebook)
        self.notebook.add(result_tab, text="结果与导出")

        # 主框架
        result_frame = ttk.LabelFrame(result_tab, text="转换结果与导出", padding=15)
        result_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # ==================== 结果预览 ====================
        preview_frame = ttk.LabelFrame(result_frame, text="转换结果预览", padding=10)
        preview_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 15))

        # 创建文本显示区域
        self.result_text = scrolledtext.ScrolledText(preview_frame, width=100, height=20,
                                                     wrap=tk.WORD, font=("Consolas", 9))
        self.result_text.pack(fill=tk.BOTH, expand=True)

        # ==================== 结果统计 ====================
        stats_frame = ttk.LabelFrame(result_frame, text="转换统计", padding=10)
        stats_frame.pack(fill=tk.X, pady=(0, 15))

        self.stats_var = tk.StringVar(value="暂无转换结果")
        ttk.Label(stats_frame, textvariable=self.stats_var,
                  font=("Arial", 10), wraplength=800).pack(anchor=tk.W)

        # ==================== 导出设置 ====================
        export_frame = ttk.LabelFrame(result_frame, text="导出结果", padding=10)
        export_frame.pack(fill=tk.X, pady=(0, 15))

        # 文件选择
        export_file_frame = ttk.Frame(export_frame)
        export_file_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(export_file_frame, text="输出文件:").pack(side=tk.LEFT, padx=(0, 5))

        self.export_file_var = tk.StringVar(value="")
        export_file_entry = ttk.Entry(export_file_frame, textvariable=self.export_file_var, width=60)
        export_file_entry.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 10))

        ttk.Button(export_file_frame, text="浏览...", command=self.browse_export_file).pack(side=tk.LEFT)

        # 文件格式选择
        format_frame = ttk.Frame(export_frame)
        format_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(format_frame, text="文件格式:").pack(side=tk.LEFT, padx=(0, 5))

        self.export_format_var = tk.StringVar(value="csv")
        ttk.Radiobutton(format_frame, text="CSV (.csv)",
                        variable=self.export_format_var, value="csv").pack(side=tk.LEFT, padx=(0, 20))
        ttk.Radiobutton(format_frame, text="Excel (.xlsx)",
                        variable=self.export_format_var, value="excel").pack(side=tk.LEFT)

        # 导出按钮
        ttk.Button(export_frame, text="导出转换结果",
                   command=self.export_results, style="Accent.TButton").pack(pady=10)

        # 导出状态
        self.export_status_var = tk.StringVar(value="等待导出...")
        ttk.Label(export_frame, textvariable=self.export_status_var,
                  font=("Arial", 10), wraplength=800).pack(anchor=tk.W)

    def browse_ppg_file(self):
        """浏览PPG文件"""
        file_types = [("数据文件", "*.csv *.xlsx *.xls"), ("CSV文件", "*.csv"),
                      ("Excel文件", "*.xlsx *.xls"), ("所有文件", "*.*")]
        file_path = filedialog.askopenfilename(title="选择PPG标准品数据文件", filetypes=file_types)

        if file_path:
            self.ppg_file_var.set(file_path)

    def browse_index_file(self):
        """浏览保留指数文件"""
        file_types = [("数据文件", "*.csv *.xlsx *.xls"), ("CSV文件", "*.csv"),
                      ("Excel文件", "*.xlsx *.xls"), ("所有文件", "*.*")]
        file_path = filedialog.askopenfilename(title="选择保留指数文件", filetypes=file_types)

        if file_path:
            self.index_file_var.set(file_path)

    def browse_export_file(self):
        """浏览导出文件"""
        default_name = "转换结果"
        if self.export_format_var.get() == "csv":
            default_name += ".csv"
            file_types = [("CSV文件", "*.csv"), ("所有文件", "*.*")]
        else:
            default_name += ".xlsx"
            file_types = [("Excel文件", "*.xlsx"), ("所有文件", "*.*")]

        file_path = filedialog.asksaveasfilename(
            title="保存转换结果",
            filetypes=file_types,
            defaultextension=".csv" if self.export_format_var.get() == "csv" else ".xlsx",
            initialfile=default_name
        )

        if file_path:
            self.export_file_var.set(file_path)

    def load_ppg_data(self):
        """加载PPG数据"""
        file_path = self.ppg_file_var.get().strip()
        condition = self.ppg_condition_var.get().strip()

        if not file_path:
            messagebox.showwarning("警告", "请选择PPG数据文件")
            return

        if not condition:
            messagebox.showwarning("警告", "请输入条件名称")
            return

        # 在线程中加载数据
        threading.Thread(target=self._load_ppg_data_thread,
                         args=(file_path, condition)).start()

    def _load_ppg_data_thread(self, file_path, condition):
        """加载PPG数据的线程函数"""
        self.update_status("正在加载PPG数据...")

        success, msg = self.converter.load_ppg_data(file_path, condition)

        if success:
            # 更新信息
            df = self.converter.ppg_data[condition]
            info_text = f"✓ {msg}\n"
            info_text += f"聚合度范围: {df['聚合度'].min()} - {df['聚合度'].max()}, "
            info_text += f"保留时间范围: {df['保留时间'].min():.2f} - {df['保留时间'].max():.2f} min"

            self.ppg_info_var.set(info_text)

            # 更新条件下拉框
            conditions = list(self.converter.ppg_data.keys())
            self.curve_condition_combo['values'] = conditions
            self.conv_condition_combo['values'] = conditions

            if conditions:
                self.curve_condition_combo.set(conditions[0])
                self.conv_condition_combo.set(conditions[0])

            self.log_message(f"✓ PPG数据加载成功: {condition}")
        else:
            self.ppg_info_var.set(f"✗ {msg}")
            self.log_message(f"✗ PPG数据加载失败: {msg}")

        self.update_status("就绪")

    def load_index_data(self):
        """加载保留指数数据"""
        file_path = self.index_file_var.get().strip()

        if not file_path:
            messagebox.showwarning("警告", "请选择保留指数文件")
            return

        # 在线程中加载数据
        threading.Thread(target=self._load_index_data_thread,
                         args=(file_path,)).start()

    def _load_index_data_thread(self, file_path):
        """加载保留指数数据的线程函数"""
        self.update_status("正在加载保留指数数据...")

        success, msg = self.converter.load_index_data(file_path)

        if success:
            # 更新信息
            self.index_info_var.set(f"✓ {msg}")

            # 更新列下拉框
            self.refresh_columns()

            self.log_message(f"✓ 保留指数数据加载成功")
        else:
            self.index_info_var.set(f"✗ {msg}")
            self.log_message(f"✗ 保留指数数据加载失败: {msg}")

        self.update_status("就绪")

    def refresh_columns(self):
        """刷新列下拉框"""
        if self.converter.index_data is not None:
            columns = list(self.converter.index_data.columns)
            self.index_column_combo['values'] = columns

            # 尝试自动选择PPG指数列
            found_column = False
            for col in columns:
                col_lower = col.lower()
                if any(keyword in col_lower for keyword in ['ppg', '指数', 'index', '保留指数']):
                    self.index_column_combo.set(col)
                    found_column = True
                    break

            # 如果没有找到合适的列，使用第一列
            if not found_column and columns:
                self.index_column_combo.set(columns[0])

    def preview_ppg_data(self):
        """预览PPG数据"""
        condition = self.ppg_condition_var.get()

        if condition not in self.converter.ppg_data:
            messagebox.showwarning("警告", "请先加载PPG数据")
            return

        df = self.converter.ppg_data[condition]

        self.preview_text.delete(1.0, tk.END)
        self.preview_text.insert(tk.END, f"PPG标准品数据 - {condition}\n")
        self.preview_text.insert(tk.END, "=" * 60 + "\n\n")
        self.preview_text.insert(tk.END, df.to_string())
        self.preview_text.insert(tk.END, f"\n\n共 {len(df)} 个数据点")

    def preview_index_data(self):
        """预览保留指数数据"""
        if self.converter.index_data is None:
            messagebox.showwarning("警告", "请先加载保留指数数据")
            return

        df = self.converter.index_data

        self.preview_text.delete(1.0, tk.END)
        self.preview_text.insert(tk.END, "保留指数数据\n")
        self.preview_text.insert(tk.END, "=" * 60 + "\n\n")

        # 只显示前20行
        preview_df = df.head(20)
        self.preview_text.insert(tk.END, preview_df.to_string())

        if len(df) > 20:
            self.preview_text.insert(tk.END, f"\n\n... (只显示前20行，共 {len(df)} 行)")
        else:
            self.preview_text.insert(tk.END, f"\n\n共 {len(df)} 行")

        # 显示列信息
        self.preview_text.insert(tk.END, "\n\n列信息:\n")
        for i, col in enumerate(df.columns):
            self.preview_text.insert(tk.END, f"  {i + 1}. {col}\n")

    def clear_preview(self):
        """清空预览"""
        self.preview_text.delete(1.0, tk.END)

    def fit_standard_curve(self):
        """拟合标准曲线"""
        condition = self.curve_condition_var.get()
        model_type = self.curve_model_var.get()

        if not condition:
            messagebox.showwarning("警告", "请选择条件")
            return

        # 在线程中拟合曲线
        threading.Thread(target=self._fit_standard_curve_thread,
                         args=(condition, model_type)).start()

    def _fit_standard_curve_thread(self, condition, model_type):
        """拟合标准曲线的线程函数"""
        self.update_status("正在拟合标准曲线...")

        success, msg = self.converter.fit_standard_curve(condition, model_type)

        if success:
            # 显示拟合结果
            curve = self.converter.standard_curves[condition]
            result_text = f"✓ {msg}\n"
            result_text += f"模型: {curve['model_name']}\n"
            result_text += f"R² = {curve['r_squared']:.6f}\n"
            result_text += f"斜率 = {curve['slope']:.4f}\n"
            result_text += f"截距 = {curve['intercept']:.4f}"

            self.curve_result_var.set(result_text)
            self.log_message(f"✓ 标准曲线拟合成功: {condition}")
        else:
            self.curve_result_var.set(f"✗ {msg}")
            self.log_message(f"✗ 标准曲线拟合失败: {msg}")

        self.update_status("就绪")

    def test_conversion(self):
        """测试转换"""
        try:
            index = float(self.test_index_var.get().strip())
            condition = self.conv_condition_var.get()
            method = self.conv_method_var.get()

            if not condition:
                messagebox.showwarning("警告", "请选择转换条件")
                return

            # 执行转换
            rt, method_used = self.converter.convert_index_to_rt(index, condition, method)

            if rt is not None:
                result_text = f"测试转换成功!\n"
                result_text += f"PPG指数: {index:.2f} → 保留时间: {rt:.4f} min\n"
                result_text += f"计算方法: {method_used}"

                self.test_result_var.set(result_text)
                self.log_message(f"✓ 测试转换成功: {index} → {rt:.4f} min")
            else:
                self.test_result_var.set(f"✗ 转换失败: {method_used}")
                self.log_message(f"✗ 测试转换失败")

        except ValueError:
            messagebox.showerror("错误", "请输入有效的PPG指数数值")
        except Exception as e:
            messagebox.showerror("错误", f"测试转换失败: {str(e)}")

    def batch_convert(self):
        """批量转换"""
        index_column = self.index_column_var.get()
        condition = self.conv_condition_var.get()
        method = self.conv_method_var.get()

        if not index_column:
            messagebox.showwarning("警告", "请选择保留指数列")
            return

        if not condition:
            messagebox.showwarning("警告", "请选择转换条件")
            return

        # 检查是否已加载数据
        if self.converter.index_data is None:
            messagebox.showwarning("警告", "请先加载保留指数数据")
            return

        if condition not in self.converter.ppg_data:
            messagebox.showwarning("警告", f"请先加载条件 '{condition}' 的PPG数据")
            return

        # 在线程中执行批量转换
        threading.Thread(target=self._batch_convert_thread,
                         args=(index_column, condition, method)).start()

    def _batch_convert_thread(self, index_column, condition, method):
        """批量转换的线程函数"""
        self.update_status("正在批量转换...")
        self.conversion_status_var.set("正在转换，请稍候...")

        success, msg = self.converter.batch_convert_file(index_column, condition, method)

        if success:
            # 更新转换状态
            self.conversion_status_var.set(f"✓ {msg}")

            # 显示结果
            self.show_conversion_results()

            # 自动切换到结果标签页
            self.notebook.select(2)

            self.log_message(f"✓ {msg}")
        else:
            self.conversion_status_var.set(f"✗ {msg}")
            self.log_message(f"✗ {msg}")

        self.update_status("就绪")

    def show_conversion_results(self):
        """显示转换结果"""
        if self.converter.converted_data is None:
            return

        df = self.converter.converted_data

        # 清空结果文本
        self.result_text.delete(1.0, tk.END)

        # 显示前50行
        preview_df = df.head(50)
        self.result_text.insert(tk.END, "转换结果预览 (前50行):\n")
        self.result_text.insert(tk.END, "=" * 80 + "\n\n")
        self.result_text.insert(tk.END, preview_df.to_string())

        if len(df) > 50:
            self.result_text.insert(tk.END, f"\n\n... (只显示前50行，共 {len(df)} 行)")
        else:
            self.result_text.insert(tk.END, f"\n\n共 {len(df)} 行")

        # 统计信息
        total_rows = len(df)

        # 找到保留时间列
        rt_columns = [col for col in df.columns if '保留时间' in col]
        if rt_columns:
            rt_col = rt_columns[0]
            successful = df[rt_col].notna().sum()
            failed = total_rows - successful

            stats_text = f"\n\n转换统计:\n"
            stats_text += f"总数据行数: {total_rows}\n"
            stats_text += f"成功转换: {successful} 行\n"
            stats_text += f"转换失败: {failed} 行\n"

            if successful > 0:
                rt_min = df[rt_col].min()
                rt_max = df[rt_col].max()
                rt_mean = df[rt_col].mean()
                stats_text += f"\n保留时间范围: {rt_min:.2f} - {rt_max:.2f} min\n"
                stats_text += f"平均保留时间: {rt_mean:.2f} min"

            self.stats_var.set(stats_text)

    def export_results(self):
        """导出结果"""
        if self.converter.converted_data is None:
            messagebox.showwarning("警告", "没有可导出的转换结果")
            return

        file_path = self.export_file_var.get().strip()

        if not file_path:
            # 使用默认文件名
            default_dir = os.path.dirname(self.index_file_var.get()) if self.index_file_var.get() else ""
            default_name = "转换结果"

            if self.export_format_var.get() == "csv":
                default_path = os.path.join(default_dir, f"{default_name}.csv")
            else:
                default_path = os.path.join(default_dir, f"{default_name}.xlsx")

            self.export_file_var.set(default_path)
            file_path = default_path

        # 确保文件扩展名正确
        if self.export_format_var.get() == "csv" and not file_path.endswith('.csv'):
            file_path += '.csv'
        elif self.export_format_var.get() == "excel" and not file_path.endswith(('.xlsx', '.xls')):
            file_path += '.xlsx'

        # 在线程中保存数据
        threading.Thread(target=self._export_results_thread,
                         args=(file_path,)).start()

    def _export_results_thread(self, file_path):
        """导出结果的线程函数"""
        self.update_status("正在导出结果...")
        self.export_status_var.set("正在导出，请稍候...")

        success, msg = self.converter.save_converted_data(file_path)

        if success:
            self.export_status_var.set(f"✓ {msg}")
            self.log_message(f"✓ 结果导出成功: {file_path}")

            # 询问是否打开文件
            if messagebox.askyesno("导出成功", f"文件已保存到:\n{file_path}\n\n是否要打开文件？"):
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
            self.log_message(f"✗ 结果导出失败: {msg}")

        self.update_status("就绪")

    def log_message(self, message):
        """记录日志消息"""
        print(f"[INFO] {message}")

    def update_status(self, message):
        """更新状态栏"""
        self.status_var.set(message)
        self.root.update()


def main():
    """主函数"""
    root = tk.Tk()
    root.title("PPG保留指数文件转换工具")

    # 设置窗口大小
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()
    window_width = int(screen_width * 0.8)
    window_height = int(screen_height * 0.85)
    x = (screen_width - window_width) // 2
    y = (screen_height - window_height) // 2

    root.geometry(f"{window_width}x{window_height}+{x}+{y}")

    # 创建应用
    app = IndexFileConverterGUI(root)

    # 运行主循环
    root.mainloop()


if __name__ == "__main__":
    print("=" * 70)
    print("PPG保留指数文件转换工具")
    print("=" * 70)
    print("功能:")
    print("  1. 加载PPG标准品数据")
    print("  2. 加载包含PPG保留指数的文件")
    print("  3. 批量将保留指数转换为保留时间")
    print("  4. 保存转换结果")
    print("=" * 70)
    print("使用方法:")
    print("  1. 在'数据加载'标签页加载PPG数据和保留指数文件")
    print("  2. 在'转换设置'标签页设置转换参数并执行转换")
    print("  3. 在'结果与导出'标签页查看和导出结果")
    print("=" * 70)

    main()