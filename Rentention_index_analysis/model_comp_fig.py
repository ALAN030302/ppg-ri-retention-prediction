import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import matplotlib
import matplotlib.font_manager as fm


# 设置中文字体，避免缺少上标2字符的问题
def setup_chinese_font():
    # 尝试不同的中文字体，优先选择包含上标2字符的字体
    chinese_fonts = [
        'Microsoft YaHei',  # 微软雅黑 - 通常包含上标字符
        'Arial Unicode MS',  # 包含完整Unicode字符集
        'DejaVu Sans',  # 开源字体，包含完整字符集
        'SimSun',  # 宋体
        'SimHei',  # 黑体 - 可能缺少上标字符
        'NSimSun',  # 新宋体
        'FangSong',  # 仿宋
        'KaiTi',  # 楷体
        'sans-serif'
    ]

    # 查找系统可用的中文字体
    available_fonts = [f.name for f in fm.fontManager.ttflist]
    selected_font = None

    for font in chinese_fonts:
        if any(font.lower() in f.lower() for f in available_fonts):
            selected_font = font
            break

    if selected_font:
        # 设置matplotlib全局字体
        plt.rcParams['font.sans-serif'] = [selected_font]
        plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题

        # 对于SimHei字体，避免使用上标2字符
        if 'SimHei' in selected_font:
            print(f"使用字体: {selected_font} (注意: 此字体可能缺少上标字符，将使用R^2代替R²)")
        else:
            print(f"使用字体: {selected_font}")
    else:
        print("警告: 未找到中文字体，可能无法正常显示中文")


# 调用字体设置函数
setup_chinese_font()

matplotlib.use('TkAgg')
from tkinter.colorchooser import askcolor
import json
import os
from datetime import datetime
import sys


class RetentionTimeAnalyzerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("四模型保留时间预测分析器")
        self.root.geometry("1600x1000")

        # Initialize variables
        self.df = None
        self.model_names = {
            'rt_smrt_pred': '6万文献模型',
            'rt_M1_pred': '非文献条件模型',
            'rt_M2_pred': '文献条件模型',
            'rt_M3_pred': '文献保留指数模型'
        }

        self.chart_settings = {
            'scatter_alpha': 0.7,
            'scatter_size': 50,
            'line_width': 2,
            'grid_alpha': 0.3,
            'font_size': 10,
            'title_size': 12,
            'colors': {
                'rt_smrt_pred': '#1f77b4',  # 蓝色
                'rt_M1_pred': '#ff7f0e',  # 橙色
                'rt_M2_pred': '#2ca02c',  # 绿色
                'rt_M3_pred': '#d62728',  # 红色
                'rt_actual': '#9467bd',  # 紫色
                'rti_M3_pred': '#8c564b'  # 棕色
            },
            'show_grid': True,
            'show_legend': True,
            'show_r_squared': True,
            'show_regression_line': True,
            'show_error_bars': True,
            'show_trend_line': True,
            'marker_shape': 'o',
            'bar_width': 0.15,
            'dpi': 100
        }

        # Load settings if exists
        self.load_settings()

        # Setup UI
        self.setup_ui()

        # Bind window close event
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)

    def setup_ui(self):
        # Configure style
        style = ttk.Style()
        style.theme_use('clam')

        # Create main paned window for resizable panels
        main_paned = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        main_paned.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Left panel for controls
        left_panel = ttk.Frame(main_paned)
        main_paned.add(left_panel, weight=1)

        # Right panel for charts
        self.right_panel = ttk.Frame(main_paned)
        main_paned.add(self.right_panel, weight=3)

        # Setup left panel controls
        self.setup_left_panel(left_panel)

        # Setup right panel
        self.setup_right_panel()

    def setup_left_panel(self, parent):
        # Notebook for different control sections
        control_notebook = ttk.Notebook(parent)
        control_notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # File control tab
        file_frame = ttk.Frame(control_notebook)
        control_notebook.add(file_frame, text="文件")
        self.setup_file_controls(file_frame)

        # Analysis control tab
        analysis_frame = ttk.Frame(control_notebook)
        control_notebook.add(analysis_frame, text="分析")
        self.setup_analysis_controls(analysis_frame)

        # Model selection tab
        model_frame = ttk.Frame(control_notebook)
        control_notebook.add(model_frame, text="模型选择")
        self.setup_model_controls(model_frame)

        # Chart settings tab
        settings_frame = ttk.Frame(control_notebook)
        control_notebook.add(settings_frame, text="图表设置")
        self.setup_chart_settings(settings_frame)

        # Results tab
        results_frame = ttk.Frame(control_notebook)
        control_notebook.add(results_frame, text="结果")
        self.setup_results_display(results_frame)

    def setup_file_controls(self, parent):
        # File selection
        file_group = ttk.LabelFrame(parent, text="数据文件", padding=10)
        file_group.pack(fill=tk.X, padx=5, pady=5)

        ttk.Button(file_group, text="浏览Excel文件...",
                   command=self.browse_file, width=20).pack(pady=5)

        self.file_path_var = tk.StringVar()
        ttk.Entry(file_group, textvariable=self.file_path_var,
                  state='readonly', width=50).pack(pady=5)

        ttk.Button(file_group, text="加载数据",
                   command=self.load_data, width=20).pack(pady=5)

        # Data preview
        preview_group = ttk.LabelFrame(parent, text="数据预览", padding=10)
        preview_group.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Treeview for data preview
        columns = ("smiles", "rt", "rt_smrt_pred", "rt_M1_pred",
                   "rt_M2_pred", "rt_M3_pred", "rti_M3_pred")
        self.preview_tree = ttk.Treeview(preview_group, columns=columns,
                                         show="headings", height=10)

        # Configure columns
        for col in columns:
            self.preview_tree.heading(col, text=col)
            self.preview_tree.column(col, width=80)

        # Scrollbars
        v_scroll = ttk.Scrollbar(preview_group, orient=tk.VERTICAL,
                                 command=self.preview_tree.yview)
        h_scroll = ttk.Scrollbar(preview_group, orient=tk.HORIZONTAL,
                                 command=self.preview_tree.xview)
        self.preview_tree.configure(yscrollcommand=v_scroll.set,
                                    xscrollcommand=h_scroll.set)

        self.preview_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        v_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        h_scroll.pack(side=tk.BOTTOM, fill=tk.X)

    def setup_analysis_controls(self, parent):
        # Analysis settings
        settings_group = ttk.LabelFrame(parent, text="分析设置", padding=10)
        settings_group.pack(fill=tk.X, padx=5, pady=5)

        # Metrics to calculate
        ttk.Label(settings_group, text="计算指标:").pack(anchor=tk.W, pady=2)

        self.metrics_vars = {
            'MAE': tk.BooleanVar(value=True),
            'RMSE': tk.BooleanVar(value=True),
            'R2': tk.BooleanVar(value=True),
            'Pearson': tk.BooleanVar(value=True),
            'MAPE': tk.BooleanVar(value=False),
            'MedianAE': tk.BooleanVar(value=False),
            'MaxError': tk.BooleanVar(value=True),
            'StdError': tk.BooleanVar(value=True)
        }

        metrics_frame = ttk.Frame(settings_group)
        metrics_frame.pack(fill=tk.X, pady=5)

        for i, (metric, var) in enumerate(self.metrics_vars.items()):
            cb = ttk.Checkbutton(metrics_frame, text=metric, variable=var)
            cb.grid(row=i // 4, column=i % 4, sticky=tk.W, padx=5, pady=2)

        # Analysis button
        ttk.Button(settings_group, text="运行分析",
                   command=self.run_analysis, width=20).pack(pady=10)

        # Chart type selection
        chart_group = ttk.LabelFrame(parent, text="图表生成", padding=10)
        chart_group.pack(fill=tk.X, padx=5, pady=5)

        self.chart_vars = {
            'scatter_all': tk.BooleanVar(value=True),
            'scatter_individual': tk.BooleanVar(value=False),
            'error_distribution': tk.BooleanVar(value=True),
            'absolute_error': tk.BooleanVar(value=True),
            'trend_comparison': tk.BooleanVar(value=True),
            'residual_plot': tk.BooleanVar(value=True),
            'metrics_bar': tk.BooleanVar(value=True),
            'boxplot': tk.BooleanVar(value=True),
            'correlation_matrix': tk.BooleanVar(value=True),
            'rank_comparison': tk.BooleanVar(value=True)
        }

        for i, (chart, var) in enumerate(self.chart_vars.items()):
            cb = ttk.Checkbutton(chart_group, text=chart.replace('_', ' ').title(),
                                 variable=var)
            cb.grid(row=i // 2, column=i % 2, sticky=tk.W, padx=5, pady=2)

    def setup_model_controls(self, parent):
        """Setup model selection controls"""
        model_group = ttk.LabelFrame(parent, text="选择比较的模型", padding=10)
        model_group.pack(fill=tk.X, padx=5, pady=5)

        # Checkboxes for each model
        self.model_vars = {}
        for model_key, model_name in self.model_names.items():
            var = tk.BooleanVar(value=True)
            cb = ttk.Checkbutton(model_group, text=model_name, variable=var)
            cb.pack(anchor=tk.W, padx=5, pady=2)
            self.model_vars[model_key] = var

        # Include RI checkbox
        self.include_ri_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(model_group, text="包含保留指数分析",
                        variable=self.include_ri_var).pack(anchor=tk.W, padx=5, pady=5)

    def setup_chart_settings(self, parent):
        # Scrollable frame for settings
        canvas = tk.Canvas(parent)
        scrollbar = ttk.Scrollbar(parent, orient=tk.VERTICAL, command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)

        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )

        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # Color settings
        color_group = ttk.LabelFrame(scrollable_frame, text="颜色设置", padding=10)
        color_group.pack(fill=tk.X, padx=5, pady=5)

        colors = [
            ('6万文献模型', 'rt_smrt_pred'),
            ('非文献条件模型', 'rt_M1_pred'),
            ('文献条件模型', 'rt_M2_pred'),
            ('文献保留指数模型', 'rt_M3_pred'),
            ('实际值', 'rt_actual'),
            ('保留指数', 'rti_M3_pred')
        ]

        self.color_buttons = {}
        for i, (label, key) in enumerate(colors):
            ttk.Label(color_group, text=label).grid(row=i, column=0, padx=5, pady=2)
            btn = ttk.Button(color_group, text="更改", width=10,
                             command=lambda k=key: self.change_color(k))
            btn.grid(row=i, column=1, padx=5, pady=2)
            self.color_buttons[key] = btn

        # Alpha and size settings
        style_group = ttk.LabelFrame(scrollable_frame, text="样式设置", padding=10)
        style_group.pack(fill=tk.X, padx=5, pady=5)

        settings = [
            ('散点透明度', 'scatter_alpha', 0.1, 1.0),
            ('散点大小', 'scatter_size', 10, 200),
            ('线宽', 'line_width', 1, 5),
            ('字体大小', 'font_size', 8, 20),
            ('网格透明度', 'grid_alpha', 0.0, 1.0),
            ('柱状图宽度', 'bar_width', 0.05, 0.3),
            ('DPI', 'dpi', 50, 300)
        ]

        self.style_widgets = {}
        for i, (label, key, min_val, max_val) in enumerate(settings):
            ttk.Label(style_group, text=label).grid(row=i, column=0, padx=5, pady=2)
            scale = ttk.Scale(style_group, from_=min_val, to=max_val, orient=tk.HORIZONTAL,
                              command=lambda v, k=key: self.update_setting(k, v))
            scale.set(self.chart_settings[key])
            scale.grid(row=i, column=1, padx=5, pady=2)
            self.style_widgets[key] = scale

        # Toggle settings
        toggle_group = ttk.LabelFrame(scrollable_frame, text="显示选项", padding=10)
        toggle_group.pack(fill=tk.X, padx=5, pady=5)

        toggles = [
            ('显示网格', 'show_grid'),
            ('显示图例', 'show_legend'),
            ('显示R^2值', 'show_r_squared'),
            ('显示回归线', 'show_regression_line'),
            ('显示误差条', 'show_error_bars'),
            ('显示趋势线', 'show_trend_line')
        ]

        self.toggle_vars = {}
        for i, (label, key) in enumerate(toggles):
            var = tk.BooleanVar(value=self.chart_settings[key])
            cb = ttk.Checkbutton(toggle_group, text=label, variable=var,
                                 command=lambda k=key, v=var: self.update_toggle(k, v))
            cb.grid(row=i // 2, column=i % 2, sticky=tk.W, padx=5, pady=2)
            self.toggle_vars[key] = var

        # Marker shape
        marker_group = ttk.LabelFrame(scrollable_frame, text="标记样式", padding=10)
        marker_group.pack(fill=tk.X, padx=5, pady=5)

        shapes = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h', 'H', '+', 'x', 'd', '|', '_']
        self.marker_var = tk.StringVar(value=self.chart_settings['marker_shape'])

        marker_frame = ttk.Frame(marker_group)
        marker_frame.pack(fill=tk.X, pady=5)

        for i, shape in enumerate(shapes):
            rb = ttk.Radiobutton(marker_frame, text=shape, value=shape,
                                 variable=self.marker_var,
                                 command=lambda: self.update_setting('marker_shape', self.marker_var.get()))
            rb.grid(row=i // 8, column=i % 8, padx=2, pady=2)

        # Reset and save settings
        button_frame = ttk.Frame(scrollable_frame)
        button_frame.pack(fill=tk.X, padx=5, pady=10)

        ttk.Button(button_frame, text="恢复默认",
                   command=self.reset_settings).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="保存设置",
                   command=self.save_settings).pack(side=tk.LEFT, padx=5)

    def setup_results_display(self, parent):
        # Results text area
        results_group = ttk.LabelFrame(parent, text="分析结果", padding=10)
        results_group.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        self.results_text = scrolledtext.ScrolledText(results_group,
                                                      wrap=tk.WORD,
                                                      font=('Consolas', 10),
                                                      height=20)
        self.results_text.pack(fill=tk.BOTH, expand=True, pady=5)

        # Export buttons
        button_frame = ttk.Frame(results_group)
        button_frame.pack(fill=tk.X, pady=5)

        ttk.Button(button_frame, text="复制结果",
                   command=self.copy_results).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="保存为文本",
                   command=self.save_results_text).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="保存为CSV",
                   command=self.save_results_csv).pack(side=tk.LEFT, padx=5)

    def setup_right_panel(self):
        # Notebook for different chart views
        self.chart_notebook = ttk.Notebook(self.right_panel)
        self.chart_notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Create initial tabs
        self.chart_tabs = {}
        tab_names = ["综合比较", "单个模型", "误差分析", "指标对比", "相关性矩阵", "排名对比"]
        for tab_name in tab_names:
            tab = ttk.Frame(self.chart_notebook)
            self.chart_notebook.add(tab, text=tab_name)
            self.chart_tabs[tab_name] = tab

        # Chart control buttons
        control_frame = ttk.Frame(self.right_panel)
        control_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Button(control_frame, text="刷新图表",
                   command=self.refresh_charts).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_frame, text="保存当前图表",
                   command=self.save_current_chart).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_frame, text="保存所有图表",
                   command=self.save_all_charts).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_frame, text="导出报告",
                   command=self.export_report).pack(side=tk.LEFT, padx=5)

    def browse_file(self):
        filename = filedialog.askopenfilename(
            title="选择Excel文件",
            filetypes=[
                ("Excel files", "*.xlsx *.xls"),
                ("CSV files", "*.csv"),
                ("All files", "*.*")
            ]
        )
        if filename:
            self.file_path_var.set(filename)

    def load_data(self):
        if not self.file_path_var.get():
            messagebox.showwarning("警告", "请先选择文件！")
            return

        try:
            # Read the file
            if self.file_path_var.get().endswith('.csv'):
                self.df = pd.read_csv(self.file_path_var.get())
            else:
                self.df = pd.read_excel(self.file_path_var.get())

            # Check required columns
            required_cols = ['smiles', 'rt', 'rt_smrt_pred', 'rt_M1_pred',
                             'rt_M2_pred', 'rt_M3_pred', 'rti_M3_pred']

            missing_cols = [col for col in required_cols if col not in self.df.columns]
            if missing_cols:
                messagebox.showerror("错误", f"缺少必要的列:\n{', '.join(missing_cols)}")
                return

            # Update preview tree
            self.update_preview_tree()

            # Enable analysis
            messagebox.showinfo("成功",
                                f"数据加载成功！\n"
                                f"数据行数: {len(self.df)}\n"
                                f"化合物数量: {len(self.df)}\n"
                                f"可用模型: {', '.join(self.model_names.values())}")

        except Exception as e:
            messagebox.showerror("错误", f"加载文件失败:\n{str(e)}")

    def update_preview_tree(self):
        # Clear existing items
        for item in self.preview_tree.get_children():
            self.preview_tree.delete(item)

        # Add new items (first 50 rows)
        for idx, row in self.df.head(50).iterrows():
            values = (
                str(row.get('smiles', ''))[:30],
                f"{row.get('rt', 0):.4f}",
                f"{row.get('rt_smrt_pred', 0):.4f}",
                f"{row.get('rt_M1_pred', 0):.4f}",
                f"{row.get('rt_M2_pred', 0):.4f}",
                f"{row.get('rt_M3_pred', 0):.4f}",
                f"{row.get('rti_M3_pred', 0):.4f}"
            )
            self.preview_tree.insert("", tk.END, values=values)

    def run_analysis(self):
        if self.df is None:
            messagebox.showwarning("警告", "请先加载数据！")
            return

        try:
            # Clear previous results
            self.results_text.delete(1.0, tk.END)

            # Get selected models
            selected_models = [key for key, var in self.model_vars.items() if var.get()]

            if not selected_models:
                messagebox.showwarning("警告", "请至少选择一个模型进行比较！")
                return

            # Calculate metrics for each selected model
            self.metrics_results = {}
            for model_key in selected_models:
                self.metrics_results[model_key] = self.calculate_metrics(
                    self.df['rt'], self.df[model_key],
                    self.model_names[model_key]
                )

            # Include RI analysis if selected
            if self.include_ri_var.get() and 'rti_M3_pred' in self.df.columns:
                # Calculate correlation between predicted RI and RT
                ri_corr = self.df['rti_M3_pred'].corr(self.df['rt'])
                self.metrics_results['rti_M3_pred'] = {
                    'Method': '保留指数预测',
                    'Pearson r': ri_corr,
                    'R^2': ri_corr ** 2
                }

            # Display results
            self.display_results()

            # Generate charts
            self.generate_charts()

            messagebox.showinfo("成功", "分析完成！")

        except Exception as e:
            messagebox.showerror("错误", f"分析失败:\n{str(e)}")

    def calculate_metrics(self, y_true, y_pred, method_name):
        metrics = {'Method': method_name}

        y_true = np.array(y_true)
        y_pred = np.array(y_pred)

        # 检查数据是否有效
        if len(y_true) != len(y_pred):
            raise ValueError(f"数据长度不一致: y_true={len(y_true)}, y_pred={len(y_pred)}")

        # 移除NaN值
        valid_mask = ~np.isnan(y_true) & ~np.isnan(y_pred)
        y_true_valid = y_true[valid_mask]
        y_pred_valid = y_pred[valid_mask]

        if len(y_true_valid) == 0:
            metrics['R^2'] = np.nan
            metrics['MAE'] = np.nan
            metrics['RMSE'] = np.nan
            metrics['Pearson r'] = np.nan
            metrics['Pearson p'] = np.nan
            metrics['IsConstant'] = False
            return metrics

        # 计算R² - 使用公式 R² = 1 - (SS_res / SS_tot)
        ss_res = np.sum((y_true_valid - y_pred_valid) ** 2)  # 残差平方和
        ss_tot = np.sum((y_true_valid - np.mean(y_true_valid)) ** 2)  # 总平方和

        if ss_tot == 0:
            # 如果y_true是常数
            if ss_res == 0:
                r2 = 1.0  # 完美预测
            else:
                r2 = -np.inf  # 比常数预测还差
        else:
            r2 = 1 - (ss_res / ss_tot)
            # 确保R²在合理范围内
            if r2 < -1:
                r2 = -1
            elif r2 > 1:
                r2 = 1

        metrics['R^2'] = r2

        # Calculate selected metrics
        if self.metrics_vars['MAE'].get():
            metrics['MAE'] = mean_absolute_error(y_true_valid, y_pred_valid)

        if self.metrics_vars['RMSE'].get():
            metrics['RMSE'] = np.sqrt(mean_squared_error(y_true_valid, y_pred_valid))

        if self.metrics_vars['Pearson'].get():
            if len(y_true_valid) > 1 and np.std(y_true_valid) > 0 and np.std(y_pred_valid) > 0:
                pearson_r, p_value = stats.pearsonr(y_true_valid, y_pred_valid)
                metrics['Pearson r'] = pearson_r
                metrics['Pearson p'] = p_value
            else:
                metrics['Pearson r'] = np.nan
                metrics['Pearson p'] = np.nan

        if self.metrics_vars['MAPE'].get():
            mask = y_true_valid != 0
            if np.any(mask):
                mape = np.mean(np.abs((y_true_valid[mask] - y_pred_valid[mask]) / y_true_valid[mask])) * 100
            else:
                mape = np.nan
            metrics['MAPE (%)'] = mape

        if self.metrics_vars['MedianAE'].get():
            metrics['MedianAE'] = np.median(np.abs(y_true_valid - y_pred_valid))

        if self.metrics_vars['MaxError'].get():
            metrics['MaxError'] = np.max(np.abs(y_true_valid - y_pred_valid))

        if self.metrics_vars['StdError'].get():
            metrics['StdError'] = np.std(y_true_valid - y_pred_valid)

        # Add constant model flag
        metrics['IsConstant'] = np.all(y_pred_valid == y_pred_valid[0]) if len(y_pred_valid) > 0 else False
        if metrics['IsConstant']:
            metrics['ConstantValue'] = y_pred_valid[0] if len(y_pred_valid) > 0 else np.nan

        return metrics

    def display_results(self):
        results_text = "=" * 100 + "\n"
        results_text += "四模型保留时间预测对比分析\n"
        results_text += "=" * 100 + "\n\n"

        results_text += f"数据集大小: {len(self.df)} 个化合物\n"
        results_text += f"有效数据点: {len(self.df.dropna(subset=['rt']))} 个\n"
        results_text += f"分析日期: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n"

        # Check for constant models
        constant_models = []
        for model_key, metrics in self.metrics_results.items():
            if metrics.get('IsConstant', False):
                constant_models.append(self.model_names.get(model_key, model_key))

        if constant_models:
            results_text += "⚠️ 警告: 以下模型预测为常数:\n"
            for model in constant_models:
                results_text += f"   • {model}\n"
            results_text += "\n"

        # Create comparison table
        all_metrics = set()
        for metrics in self.metrics_results.values():
            all_metrics.update(metrics.keys())

        # Remove unwanted columns
        exclude_cols = {'Method', 'IsConstant', 'ConstantValue', 'Pearson p'}
        all_metrics = all_metrics - exclude_cols

        # Sort metrics in logical order
        metric_order = ['R^2', 'Pearson r', 'MAE', 'RMSE', 'MAPE (%)',
                        'MedianAE', 'MaxError', 'StdError']
        ordered_metrics = [m for m in metric_order if m in all_metrics]
        other_metrics = sorted([m for m in all_metrics if m not in metric_order])
        all_metrics_sorted = ordered_metrics + other_metrics

        # Create table header
        header = f"{'指标':<20}"
        for model_key in self.metrics_results.keys():
            model_name = self.model_names.get(model_key, model_key)
            header += f" {model_name[:15]:>15}"
        header += f" {'最优模型':>15}\n"

        results_text += header
        results_text += "-" * (20 + 15 * (len(self.metrics_results) + 1)) + "\n"

        # Add rows
        for metric in all_metrics_sorted:
            row = f"{metric:<20}"
            values = []

            for model_key in self.metrics_results.keys():
                val = self.metrics_results[model_key].get(metric, np.nan)
                if isinstance(val, (int, float, np.number)) and not np.isnan(val):
                    if metric in ['R^2', 'Pearson r']:
                        val_fmt = f"{val:.4f}"
                    elif metric in ['MAPE (%)']:
                        val_fmt = f"{val:.2f}%"
                    else:
                        val_fmt = f"{val:.4f}"
                else:
                    val_fmt = 'N/A'

                row += f" {val_fmt:>15}"
                values.append(val if val_fmt != 'N/A' else None)

            # Determine best model
            valid_values = [v for v in values if v is not None and not np.isnan(v)]
            if valid_values:
                if metric in ['MAE', 'RMSE', 'MAPE (%)', 'MedianAE', 'MaxError', 'StdError']:
                    # 对于误差指标，越小越好
                    best_idx = np.argmin(valid_values)
                else:  # R^2, Pearson r
                    # 对于相关系数指标，越大越好
                    best_idx = np.argmax(valid_values)

                best_model = list(self.metrics_results.keys())[best_idx]
                best_name = self.model_names.get(best_model, best_model)
                row += f" {best_name[:15]:>15}"
            else:
                row += f" {'N/A':>15}"

            results_text += row + "\n"

        # Add R² interpretation
        results_text += "\n" + "=" * 100 + "\n"
        results_text += "R^2值详细解释:\n"
        results_text += "=" * 100 + "\n"

        for model_key, metrics in self.metrics_results.items():
            if 'R^2' in metrics:
                r2 = metrics['R^2']
                model_name = self.model_names.get(model_key, model_key)
                results_text += f"\n{model_name}: R^2 = {r2:.6f}\n"

                if np.isnan(r2):
                    results_text += "  • R^2 = NaN: 数据存在问题，无法计算R^2\n"
                elif r2 == -np.inf:
                    results_text += "  • R^2 = -∞: 模型预测比使用常数更差\n"
                elif r2 < 0:
                    results_text += f"  • R^2 = {r2:.3f}: 模型预测比使用简单平均值更差\n"
                    results_text += f"    残差平方和(SS_res) = {np.sum((self.df['rt'] - self.df[model_key]) ** 2):.2f}\n"
                    results_text += f"    总平方和(SS_tot) = {np.sum((self.df['rt'] - np.mean(self.df['rt'])) ** 2):.2f}\n"
                elif r2 == 0:
                    results_text += "  • R^2 = 0: 模型预测与使用简单平均值相同\n"
                elif r2 < 0.3:
                    results_text += f"  • R^2 = {r2:.3f}: 预测能力较弱 (解释{round(r2 * 100, 1)}%的方差)\n"
                elif r2 < 0.7:
                    results_text += f"  • R^2 = {r2:.3f}: 预测能力中等 (解释{round(r2 * 100, 1)}%的方差)\n"
                else:
                    results_text += f"  • R^2 = {r2:.3f}: 预测能力强 (解释{round(r2 * 100, 1)}%的方差)\n"

        # Summary conclusions
        results_text += "\n" + "=" * 100 + "\n"
        results_text += "综合分析结论:\n"
        results_text += "=" * 100 + "\n"

        # Count best performances
        performance_counts = {model_key: 0 for model_key in self.metrics_results.keys()}

        for metric in all_metrics_sorted:
            values = []
            valid_models = []

            for model_key in self.metrics_results.keys():
                val = self.metrics_results[model_key].get(metric, np.nan)
                if isinstance(val, (int, float, np.number)) and not np.isnan(val):
                    values.append(val)
                    valid_models.append(model_key)

            if values:
                if metric in ['MAE', 'RMSE', 'MAPE (%)', 'MedianAE', 'MaxError', 'StdError']:
                    best_idx = np.argmin(values)
                else:
                    best_idx = np.argmax(values)

                best_model = valid_models[best_idx]
                performance_counts[best_model] += 1

        # Sort models by performance
        sorted_models = sorted(performance_counts.items(), key=lambda x: x[1], reverse=True)

        results_text += "模型性能排名:\n"
        for i, (model_key, count) in enumerate(sorted_models, 1):
            model_name = self.model_names.get(model_key, model_key)
            results_text += f"{i}. {model_name}: 在{len(all_metrics_sorted)}个指标中，{count}个指标最优\n"

        self.results_text.insert(1.0, results_text)

    def generate_charts(self):
        # Clear all chart frames
        for tab in self.chart_tabs.values():
            for widget in tab.winfo_children():
                widget.destroy()

        # Get selected models
        selected_models = [key for key in self.metrics_results.keys()
                           if key in self.model_names]

        if not selected_models:
            return

        # Generate charts based on selections
        if any(var.get() for var in self.chart_vars.values()):
            self.create_comprehensive_charts(selected_models)
            self.create_individual_charts(selected_models)
            self.create_error_analysis_charts(selected_models)
            self.create_metrics_comparison_charts(selected_models)
            self.create_correlation_matrix(selected_models)
            self.create_rank_comparison_charts(selected_models)

    def create_comprehensive_charts(self, selected_models):
        """Create comprehensive comparison charts"""
        tab = self.chart_tabs["综合比较"]

        # Determine which charts to create
        charts_to_create = []
        for chart_name, var in self.chart_vars.items():
            if var.get() and chart_name in ['scatter_all', 'trend_comparison', 'boxplot']:
                charts_to_create.append(chart_name)

        if not charts_to_create:
            return

        # Create grid
        n_charts = len(charts_to_create)
        n_cols = min(2, n_charts)
        n_rows = (n_charts + n_cols - 1) // n_cols

        fig = Figure(figsize=(7 * n_cols, 5 * n_rows), dpi=self.chart_settings['dpi'])

        for idx, chart_type in enumerate(charts_to_create):
            ax = fig.add_subplot(n_rows, n_cols, idx + 1)

            if chart_type == 'scatter_all':
                self.create_all_models_scatter(ax, selected_models)
            elif chart_type == 'trend_comparison':
                self.create_trend_comparison_plot(ax, selected_models)
            elif chart_type == 'boxplot':
                self.create_error_boxplot(ax, selected_models)

        fig.tight_layout(pad=3.0)

        # Embed in tkinter
        canvas = FigureCanvasTkAgg(fig, master=tab)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Add toolbar
        toolbar = NavigationToolbar2Tk(canvas, tab)
        toolbar.update()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        self.current_comprehensive_fig = fig
        self.current_comprehensive_canvas = canvas

    def create_all_models_scatter(self, ax, selected_models):
        """Create scatter plot for all models"""
        alpha = self.chart_settings['scatter_alpha']
        size = self.chart_settings['scatter_size']
        marker = self.chart_settings['marker_shape']

        # Plot each model
        for model_key in selected_models:
            color = self.chart_settings['colors'][model_key]
            label = self.model_names[model_key]
            r2 = self.metrics_results[model_key].get('R^2', 0)

            if np.isnan(r2):
                r2_label = f'{label} (R^2=NaN)'
            elif r2 < 0:
                r2_label = f'{label} (R^2={r2:.3f}) - 比平均值差'
            else:
                r2_label = f'{label} (R^2={r2:.3f})'

            # 移除NaN值
            valid_mask = ~self.df['rt'].isna() & ~self.df[model_key].isna()
            x_data = self.df.loc[valid_mask, 'rt']
            y_data = self.df.loc[valid_mask, model_key]

            ax.scatter(x_data, y_data,
                       alpha=alpha, s=size, marker=marker,
                       color=color, label=r2_label)

        # Add regression lines
        if self.chart_settings['show_regression_line']:
            for model_key in selected_models:
                color = self.chart_settings['colors'][model_key]
                y_pred = self.df[model_key]
                valid_mask = ~self.df['rt'].isna() & ~y_pred.isna()
                x_data = self.df.loc[valid_mask, 'rt']
                y_data = y_pred.loc[valid_mask]

                if len(x_data) > 1 and np.std(y_data) > 0:
                    z = np.polyfit(x_data, y_data, 1)
                    p = np.poly1d(z)
                    x_range = np.linspace(x_data.min(), x_data.max(), 100)
                    ax.plot(x_range, p(x_range), '--', color=color,
                            linewidth=self.chart_settings['line_width'] * 0.7, alpha=0.7)

        # Add diagonal line
        min_val = min(self.df['rt'].min(), min([self.df[m].min() for m in selected_models]))
        max_val = max(self.df['rt'].max(), max([self.df[m].max() for m in selected_models]))
        ax.plot([min_val, max_val], [min_val, max_val], 'k:', alpha=0.3, label='理想线')

        # Configure plot with explicit font properties
        ax.set_xlabel('实际保留时间', fontsize=self.chart_settings['font_size'])
        ax.set_ylabel('预测保留时间', fontsize=self.chart_settings['font_size'])
        ax.set_title('四模型预测对比散点图', fontsize=self.chart_settings['font_size'] + 2)

        if self.chart_settings['show_grid']:
            ax.grid(True, alpha=self.chart_settings['grid_alpha'])

        if self.chart_settings['show_legend']:
            ax.legend(fontsize=self.chart_settings['font_size'] - 2,
                      loc='upper left', bbox_to_anchor=(1, 1))

    def create_trend_comparison_plot(self, ax, selected_models):
        """Create trend comparison plot"""
        # Remove NaN values
        valid_mask = ~self.df['rt'].isna()
        for model in selected_models:
            valid_mask = valid_mask & ~self.df[model].isna()

        # Sort by actual retention time
        df_valid = self.df[valid_mask]
        sorted_idx = df_valid['rt'].sort_values().index
        x = range(len(df_valid))

        # Plot actual values
        ax.plot(x, df_valid.loc[sorted_idx, 'rt'].values,
                color=self.chart_settings['colors']['rt_actual'],
                linewidth=self.chart_settings['line_width'] * 1.5,
                label='实际值', zorder=5)

        # Plot predictions
        for model_key in selected_models:
            color = self.chart_settings['colors'][model_key]
            label = self.model_names[model_key]

            ax.plot(x, df_valid.loc[sorted_idx, model_key].values,
                    '--', color=color, alpha=0.7,
                    linewidth=self.chart_settings['line_width'],
                    label=label)

        # Configure plot
        ax.set_xlabel('化合物 (按实际RT排序)', fontsize=self.chart_settings['font_size'])
        ax.set_ylabel('保留时间', fontsize=self.chart_settings['font_size'])
        ax.set_title('保留时间趋势对比', fontsize=self.chart_settings['font_size'] + 2)

        if self.chart_settings['show_grid']:
            ax.grid(True, alpha=self.chart_settings['grid_alpha'])

        if self.chart_settings['show_legend']:
            ax.legend(fontsize=self.chart_settings['font_size'] - 2)

    def create_error_boxplot(self, ax, selected_models):
        """Create error distribution boxplot"""
        # Calculate absolute errors for each model
        error_data = []
        labels = []

        for model_key in selected_models:
            # Remove NaN values
            valid_mask = ~self.df['rt'].isna() & ~self.df[model_key].isna()
            errors = np.abs(self.df.loc[valid_mask, model_key] - self.df.loc[valid_mask, 'rt'])
            if len(errors) > 0:
                error_data.append(errors)
                labels.append(self.model_names[model_key])

        if not error_data:
            ax.text(0.5, 0.5, '无有效数据', transform=ax.transAxes,
                    ha='center', va='center', fontsize=self.chart_settings['font_size'])
            return

        # Create boxplot - 修复Matplotlib 3.9+的兼容性问题
        try:
            # Matplotlib 3.9+ 使用 tick_labels 参数
            import matplotlib
            if hasattr(matplotlib, '__version__'):
                version = matplotlib.__version__
                if int(version.split('.')[0]) >= 3 and int(version.split('.')[1]) >= 9:
                    # 使用新的参数名
                    bp = ax.boxplot(error_data, tick_labels=labels, patch_artist=True)
                else:
                    # 使用旧的参数名
                    bp = ax.boxplot(error_data, labels=labels, patch_artist=True)
            else:
                # 默认使用旧参数名
                bp = ax.boxplot(error_data, labels=labels, patch_artist=True)
        except:
            # 如果出错，尝试两种方式
            try:
                bp = ax.boxplot(error_data, tick_labels=labels, patch_artist=True)
            except:
                bp = ax.boxplot(error_data, labels=labels, patch_artist=True)

        # Color boxes
        for i, box in enumerate(bp['boxes']):
            model_key = selected_models[i]
            box.set_facecolor(self.chart_settings['colors'][model_key])
            box.set_alpha(0.7)

        # Configure plot
        ax.set_xlabel('预测模型', fontsize=self.chart_settings['font_size'])
        ax.set_ylabel('绝对误差', fontsize=self.chart_settings['font_size'])
        ax.set_title('绝对误差分布箱线图', fontsize=self.chart_settings['font_size'] + 2)

        if self.chart_settings['show_grid']:
            ax.grid(True, alpha=self.chart_settings['grid_alpha'], axis='y')

        # Add mean values as text
        for i, model_key in enumerate(selected_models):
            valid_mask = ~self.df['rt'].isna() & ~self.df[model_key].isna()
            errors = np.abs(self.df.loc[valid_mask, model_key] - self.df.loc[valid_mask, 'rt'])
            if len(errors) > 0:
                mean_error = errors.mean()
                ax.text(i + 1, ax.get_ylim()[1] * 0.95, f'均值: {mean_error:.3f}',
                        ha='center', va='top', fontsize=self.chart_settings['font_size'] - 2)

    def create_individual_charts(self, selected_models):
        """Create individual charts for each model"""
        if not self.chart_vars['scatter_individual'].get():
            return

        tab = self.chart_tabs["单个模型"]

        n_models = len(selected_models)
        n_cols = min(2, n_models)
        n_rows = (n_models + n_cols - 1) // n_cols

        fig = Figure(figsize=(7 * n_cols, 5 * n_rows), dpi=self.chart_settings['dpi'])

        for idx, model_key in enumerate(selected_models):
            ax = fig.add_subplot(n_rows, n_cols, idx + 1)
            self.create_model_scatter(ax, model_key)

        fig.tight_layout(pad=3.0)

        # Embed in tkinter
        canvas = FigureCanvasTkAgg(fig, master=tab)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Add toolbar
        toolbar = NavigationToolbar2Tk(canvas, tab)
        toolbar.update()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        self.current_individual_fig = fig
        self.current_individual_canvas = canvas

    def create_model_scatter(self, ax, model_key):
        """Create scatter plot for individual model"""
        alpha = self.chart_settings['scatter_alpha']
        size = self.chart_settings['scatter_size']
        marker = self.chart_settings['marker_shape']
        color = self.chart_settings['colors'][model_key]
        model_name = self.model_names[model_key]

        # Remove NaN values
        valid_mask = ~self.df['rt'].isna() & ~self.df[model_key].isna()
        x_data = self.df.loc[valid_mask, 'rt']
        y_data = self.df.loc[valid_mask, model_key]

        if len(x_data) == 0:
            ax.text(0.5, 0.5, '无有效数据', transform=ax.transAxes,
                    ha='center', va='center', fontsize=self.chart_settings['font_size'])
            ax.set_title(f'{model_name}预测结果', fontsize=self.chart_settings['font_size'] + 2)
            return

        # Scatter plot
        ax.scatter(x_data, y_data,
                   alpha=alpha, s=size, marker=marker, color=color)

        # Regression line
        if self.chart_settings['show_regression_line'] and len(x_data) > 1:
            if np.std(y_data) > 0:
                z = np.polyfit(x_data, y_data, 1)
                p = np.poly1d(z)
                x_range = np.linspace(x_data.min(), x_data.max(), 100)
                ax.plot(x_range, p(x_range), 'r--',
                        linewidth=self.chart_settings['line_width'] * 1.5,
                        label=f'y = {z[0]:.3f}x + {z[1]:.3f}')

        # Add diagonal line
        min_val = min(x_data.min(), y_data.min())
        max_val = max(x_data.max(), y_data.max())
        ax.plot([min_val, max_val], [min_val, max_val], 'k:', alpha=0.3, label='理想线')

        # Add R^2 text
        if self.chart_settings['show_r_squared']:
            r2 = self.metrics_results[model_key].get('R^2', 0)
            if np.isnan(r2):
                r2_text = f'R^2 = NaN (数据问题)'
                box_color = 'gray'
            elif r2 == -np.inf:
                r2_text = f'R^2 = -∞ (比常数预测差)'
                box_color = 'red'
            elif r2 < 0:
                r2_text = f'R^2 = {r2:.3f} (比平均值差)'
                box_color = 'red'
            elif r2 < 0.3:
                r2_text = f'R^2 = {r2:.3f} (弱)'
                box_color = 'orange'
            elif r2 < 0.7:
                r2_text = f'R^2 = {r2:.3f} (中等)'
                box_color = 'yellow'
            else:
                r2_text = f'R^2 = {r2:.3f} (强)'
                box_color = 'green'

            ax.text(0.05, 0.95, r2_text,
                    transform=ax.transAxes,
                    fontsize=self.chart_settings['font_size'],
                    verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor=box_color, alpha=0.5))

        # Configure plot
        ax.set_xlabel('实际保留时间', fontsize=self.chart_settings['font_size'])
        ax.set_ylabel('预测保留时间', fontsize=self.chart_settings['font_size'])
        ax.set_title(f'{model_name}预测结果', fontsize=self.chart_settings['font_size'] + 2)

        if self.chart_settings['show_grid']:
            ax.grid(True, alpha=self.chart_settings['grid_alpha'])

        if self.chart_settings['show_legend'] and self.chart_settings['show_regression_line']:
            ax.legend(fontsize=self.chart_settings['font_size'] - 2)

    def create_error_analysis_charts(self, selected_models):
        """Create error analysis charts"""
        if not (self.chart_vars['error_distribution'].get() or
                self.chart_vars['absolute_error'].get() or
                self.chart_vars['residual_plot'].get()):
            return

        tab = self.chart_tabs["误差分析"]

        charts_to_create = []
        for chart_name, var in self.chart_vars.items():
            if var.get() and chart_name in ['error_distribution', 'absolute_error', 'residual_plot']:
                charts_to_create.append(chart_name)

        if not charts_to_create:
            return

        n_charts = len(charts_to_create)
        n_cols = min(2, n_charts)
        n_rows = (n_charts + n_cols - 1) // n_cols

        fig = Figure(figsize=(7 * n_cols, 5 * n_rows), dpi=self.chart_settings['dpi'])

        for idx, chart_type in enumerate(charts_to_create):
            ax = fig.add_subplot(n_rows, n_cols, idx + 1)

            if chart_type == 'error_distribution':
                self.create_error_distribution_plot(ax, selected_models)
            elif chart_type == 'absolute_error':
                self.create_absolute_error_violin(ax, selected_models)
            elif chart_type == 'residual_plot':
                self.create_residual_plot_multi(ax, selected_models)

        fig.tight_layout(pad=3.0)

        # Embed in tkinter
        canvas = FigureCanvasTkAgg(fig, master=tab)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Add toolbar
        toolbar = NavigationToolbar2Tk(canvas, tab)
        toolbar.update()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        self.current_error_fig = fig
        self.current_error_canvas = canvas

    def create_error_distribution_plot(self, ax, selected_models):
        """Create error distribution histogram"""
        bins = 30

        for model_key in selected_models:
            color = self.chart_settings['colors'][model_key]
            label = self.model_names[model_key]

            # Remove NaN values
            valid_mask = ~self.df['rt'].isna() & ~self.df[model_key].isna()
            errors = self.df.loc[valid_mask, model_key] - self.df.loc[valid_mask, 'rt']

            if len(errors) > 0:
                ax.hist(errors, bins=bins, alpha=0.5,
                        color=color, label=label, edgecolor='black', density=True)

        # Add vertical line at zero
        ax.axvline(x=0, color='red', linestyle='--',
                   linewidth=self.chart_settings['line_width'],
                   label='零误差')

        # Configure plot
        ax.set_xlabel('预测误差 (预测值 - 实际值)', fontsize=self.chart_settings['font_size'])
        ax.set_ylabel('密度', fontsize=self.chart_settings['font_size'])
        ax.set_title('误差分布直方图 (归一化)', fontsize=self.chart_settings['font_size'] + 2)

        if self.chart_settings['show_grid']:
            ax.grid(True, alpha=self.chart_settings['grid_alpha'])

        if self.chart_settings['show_legend']:
            ax.legend(fontsize=self.chart_settings['font_size'] - 2)

    def create_absolute_error_violin(self, ax, selected_models):
        """Create absolute error violin plot"""
        # Calculate absolute errors
        error_data = []
        labels = []

        for model_key in selected_models:
            # Remove NaN values
            valid_mask = ~self.df['rt'].isna() & ~self.df[model_key].isna()
            abs_errors = np.abs(self.df.loc[valid_mask, model_key] - self.df.loc[valid_mask, 'rt'])
            if len(abs_errors) > 0:
                error_data.append(abs_errors)
                labels.append(self.model_names[model_key])

        if not error_data:
            ax.text(0.5, 0.5, '无有效数据', transform=ax.transAxes,
                    ha='center', va='center', fontsize=self.chart_settings['font_size'])
            return

        # Create violin plot
        vp = ax.violinplot(error_data, showmeans=True, showmedians=True)

        # Color the violins
        for i, pc in enumerate(vp['bodies']):
            model_key = selected_models[i]
            pc.set_facecolor(self.chart_settings['colors'][model_key])
            pc.set_alpha(0.7)

        # Configure plot
        ax.set_xlabel('预测模型', fontsize=self.chart_settings['font_size'])
        ax.set_ylabel('绝对误差', fontsize=self.chart_settings['font_size'])
        ax.set_title('绝对误差分布小提琴图', fontsize=self.chart_settings['font_size'] + 2)
        ax.set_xticks(range(1, len(labels) + 1))
        ax.set_xticklabels(labels, rotation=15)

        if self.chart_settings['show_grid']:
            ax.grid(True, alpha=self.chart_settings['grid_alpha'], axis='y')

    def create_residual_plot_multi(self, ax, selected_models):
        """Create residual plot for multiple models"""
        for model_key in selected_models:
            color = self.chart_settings['colors'][model_key]
            label = self.model_names[model_key]

            # Remove NaN values
            valid_mask = ~self.df['rt'].isna() & ~self.df[model_key].isna()
            residuals = self.df.loc[valid_mask, model_key] - self.df.loc[valid_mask, 'rt']
            x_data = self.df.loc[valid_mask, 'rt']

            if len(residuals) > 0:
                ax.scatter(x_data, residuals,
                           alpha=self.chart_settings['scatter_alpha'] * 0.7,
                           s=self.chart_settings['scatter_size'] * 0.7,
                           color=color, label=label,
                           marker=self.chart_settings['marker_shape'])

        # Add horizontal line at zero
        ax.axhline(y=0, color='red', linestyle='--',
                   linewidth=self.chart_settings['line_width'])

        # Configure plot
        ax.set_xlabel('实际保留时间', fontsize=self.chart_settings['font_size'])
        ax.set_ylabel('残差 (预测值 - 实际值)', fontsize=self.chart_settings['font_size'])
        ax.set_title('残差图', fontsize=self.chart_settings['font_size'] + 2)

        if self.chart_settings['show_grid']:
            ax.grid(True, alpha=self.chart_settings['grid_alpha'])

        if self.chart_settings['show_legend']:
            ax.legend(fontsize=self.chart_settings['font_size'] - 2,
                      loc='upper left', bbox_to_anchor=(1, 1))

    def create_metrics_comparison_charts(self, selected_models):
        """Create metrics comparison charts"""
        if not self.chart_vars['metrics_bar'].get():
            return

        tab = self.chart_tabs["指标对比"]

        # Get metrics to display
        metrics_to_show = []
        for metric, var in self.metrics_vars.items():
            if var.get():
                if metric == 'R2':
                    metrics_to_show.append('R^2')
                elif metric == 'Pearson':
                    metrics_to_show.append('Pearson r')
                elif metric == 'MAPE':
                    metrics_to_show.append('MAPE (%)')
                elif metric == 'MedianAE':
                    metrics_to_show.append('MedianAE')
                elif metric == 'MaxError':
                    metrics_to_show.append('MaxError')
                elif metric == 'StdError':
                    metrics_to_show.append('StdError')
                else:
                    metrics_to_show.append(metric)

        n_metrics = len(metrics_to_show)
        n_cols = min(2, n_metrics)
        n_rows = (n_metrics + n_cols - 1) // n_cols

        fig = Figure(figsize=(7 * n_cols, 5 * n_rows), dpi=self.chart_settings['dpi'])

        for idx, metric in enumerate(metrics_to_show):
            ax = fig.add_subplot(n_rows, n_cols, idx + 1)
            self.create_metric_bar_chart(ax, selected_models, metric)

        fig.tight_layout(pad=3.0)

        # Embed in tkinter
        canvas = FigureCanvasTkAgg(fig, master=tab)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Add toolbar
        toolbar = NavigationToolbar2Tk(canvas, tab)
        toolbar.update()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        self.current_metrics_fig = fig
        self.current_metrics_canvas = canvas

    def create_metric_bar_chart(self, ax, selected_models, metric):
        """Create bar chart for specific metric"""
        # Get values for each model
        values = []
        labels = []

        # 将显示标签中的R^2转换为内部键名R^2
        internal_metric = 'R^2' if metric == 'R^2' else metric

        for model_key in selected_models:
            val = self.metrics_results[model_key].get(internal_metric, 0)
            values.append(val if not np.isnan(val) else 0)
            labels.append(self.model_names[model_key])

        # Create bar chart
        x = np.arange(len(labels))
        width = self.chart_settings['bar_width']

        bars = ax.bar(x, values, width)

        # Color bars
        for i, bar in enumerate(bars):
            model_key = selected_models[i]
            bar.set_facecolor(self.chart_settings['colors'][model_key])
            bar.set_alpha(0.7)

            # Add value label
            height = bar.get_height()
            if not np.isnan(height):
                if metric == 'MAPE (%)':
                    label = f'{height:.2f}%'
                elif metric == 'R^2':
                    label = f'{height:.4f}'
                else:
                    label = f'{height:.4f}'

                ax.text(bar.get_x() + bar.get_width() / 2., height + 0.01,
                        label, ha='center', va='bottom',
                        fontsize=self.chart_settings['font_size'] - 2)

        # Configure plot
        ax.set_xlabel('预测模型', fontsize=self.chart_settings['font_size'])
        ax.set_ylabel(metric, fontsize=self.chart_settings['font_size'])
        ax.set_title(f'{metric}对比', fontsize=self.chart_settings['font_size'] + 2)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=15)

        if self.chart_settings['show_grid']:
            ax.grid(True, alpha=self.chart_settings['grid_alpha'], axis='y')

        # Add reference lines
        if metric == 'R^2':
            ax.axhline(y=1.0, color='green', linestyle=':', alpha=0.5, label='完美预测')
            ax.axhline(y=0.0, color='red', linestyle=':', alpha=0.5, label='与平均值相同')
        elif metric in ['MAE', 'RMSE', 'MAPE (%)', 'MedianAE', 'MaxError', 'StdError']:
            ax.axhline(y=0.0, color='black', linestyle='-', linewidth=0.5)

    def create_correlation_matrix(self, selected_models):
        """Create correlation matrix heatmap"""
        if not self.chart_vars['correlation_matrix'].get():
            return

        tab = self.chart_tabs["相关性矩阵"]

        # Prepare data for correlation matrix
        columns = ['rt'] + selected_models
        if self.include_ri_var.get() and 'rti_M3_pred' in self.df.columns:
            columns.append('rti_M3_pred')

        # Remove rows with NaN values
        df_valid = self.df[columns].dropna()

        if len(df_valid) == 0:
            fig = Figure(figsize=(8, 6), dpi=self.chart_settings['dpi'])
            ax = fig.add_subplot(111)
            ax.text(0.5, 0.5, '无有效数据计算相关性', transform=ax.transAxes,
                    ha='center', va='center', fontsize=self.chart_settings['font_size'])
            ax.set_title('相关性矩阵', fontsize=self.chart_settings['font_size'] + 2)

            canvas = FigureCanvasTkAgg(fig, master=tab)
            canvas.draw()
            canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
            return

        corr_data = df_valid.corr()

        # Rename columns for display
        display_names = {'rt': '实际值'}
        display_names.update(self.model_names)
        if 'rti_M3_pred' in columns:
            display_names['rti_M3_pred'] = '保留指数'

        corr_data.index = [display_names.get(col, col) for col in corr_data.index]
        corr_data.columns = [display_names.get(col, col) for col in corr_data.columns]

        fig = Figure(figsize=(10, 8), dpi=self.chart_settings['dpi'])
        ax = fig.add_subplot(111)

        # Create heatmap
        im = ax.imshow(corr_data.values, cmap='coolwarm', vmin=-1, vmax=1)

        # Add text annotations
        for i in range(len(corr_data)):
            for j in range(len(corr_data)):
                ax.text(j, i, f'{corr_data.iloc[i, j]:.3f}',
                        ha='center', va='center',
                        color='white' if abs(corr_data.iloc[i, j]) > 0.5 else 'black',
                        fontsize=self.chart_settings['font_size'])

        # Configure plot
        ax.set_xticks(range(len(corr_data)))
        ax.set_yticks(range(len(corr_data)))
        ax.set_xticklabels(corr_data.columns, rotation=45, ha='right')
        ax.set_yticklabels(corr_data.index)
        ax.set_title('相关性矩阵', fontsize=self.chart_settings['font_size'] + 2)

        # Add colorbar
        fig.colorbar(im, ax=ax)

        fig.tight_layout()

        # Embed in tkinter
        canvas = FigureCanvasTkAgg(fig, master=tab)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Add toolbar
        toolbar = NavigationToolbar2Tk(canvas, tab)
        toolbar.update()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        self.current_correlation_fig = fig
        self.current_correlation_canvas = canvas

    def create_rank_comparison_charts(self, selected_models):
        """Create rank comparison charts"""
        if not self.chart_vars['rank_comparison'].get():
            return

        tab = self.chart_tabs["排名对比"]

        # Remove rows with NaN values
        df_valid = self.df[['rt'] + selected_models].dropna()

        if len(df_valid) == 0:
            fig = Figure(figsize=(12, 5), dpi=self.chart_settings['dpi'])
            ax = fig.add_subplot(111)
            ax.text(0.5, 0.5, '无有效数据计算排名', transform=ax.transAxes,
                    ha='center', va='center', fontsize=self.chart_settings['font_size'])

            canvas = FigureCanvasTkAgg(fig, master=tab)
            canvas.draw()
            canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
            return

        # Calculate compound-wise performance ranking
        compound_rankings = []

        for idx, row in df_valid.iterrows():
            actual = row['rt']
            predictions = {model_key: row[model_key] for model_key in selected_models}

            # Calculate absolute errors
            errors = {model: abs(predictions[model] - actual) for model in predictions}

            # Rank models (lower error is better)
            ranked_models = sorted(errors.items(), key=lambda x: x[1])
            ranks = {model: rank + 1 for rank, (model, _) in enumerate(ranked_models)}

            compound_rankings.append(ranks)

        # Calculate average rank for each model
        avg_ranks = {}
        for model_key in selected_models:
            ranks = [ranking[model_key] for ranking in compound_rankings]
            avg_ranks[model_key] = np.mean(ranks)

        # Create bar chart of average ranks
        fig = Figure(figsize=(12, 5), dpi=self.chart_settings['dpi'])

        # Subplot 1: Average rank bar chart
        ax1 = fig.add_subplot(121)

        models_sorted = sorted(avg_ranks.items(), key=lambda x: x[1])
        x = range(len(models_sorted))
        values = [rank for _, rank in models_sorted]
        labels = [self.model_names[model] for model, _ in models_sorted]

        bars = ax1.bar(x, values, width=0.6)

        # Color bars
        for i, bar in enumerate(bars):
            model_key, _ = models_sorted[i]
            bar.set_facecolor(self.chart_settings['colors'][model_key])
            bar.set_alpha(0.7)

            # Add value label
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width() / 2., height + 0.1,
                     f'{height:.2f}', ha='center', va='bottom',
                     fontsize=self.chart_settings['font_size'] - 2)

        ax1.set_xlabel('预测模型', fontsize=self.chart_settings['font_size'])
        ax1.set_ylabel('平均排名 (1=最优)', fontsize=self.chart_settings['font_size'])
        ax1.set_title('模型平均性能排名', fontsize=self.chart_settings['font_size'] + 2)
        ax1.set_xticks(x)
        ax1.set_xticklabels(labels, rotation=15)

        if self.chart_settings['show_grid']:
            ax1.grid(True, alpha=self.chart_settings['grid_alpha'], axis='y')

        # Subplot 2: Rank distribution heatmap
        ax2 = fig.add_subplot(122)

        # Prepare rank matrix
        rank_matrix = np.zeros((len(selected_models), len(selected_models)))

        for compound_ranking in compound_rankings:
            for model_key, rank in compound_ranking.items():
                idx = list(selected_models).index(model_key)
                rank_matrix[idx, rank - 1] += 1

        # Normalize by row
        rank_matrix_norm = rank_matrix / rank_matrix.sum(axis=1, keepdims=True)

        im = ax2.imshow(rank_matrix_norm, cmap='YlOrRd', aspect='auto')

        # Add text annotations
        for i in range(len(selected_models)):
            for j in range(len(selected_models)):
                value = rank_matrix_norm[i, j]
                if value > 0:
                    ax2.text(j, i, f'{value:.2f}', ha='center', va='center',
                             color='white' if value > 0.5 else 'black',
                             fontsize=self.chart_settings['font_size'] - 2)

        ax2.set_xlabel('排名位置', fontsize=self.chart_settings['font_size'])
        ax2.set_ylabel('预测模型', fontsize=self.chart_settings['font_size'])
        ax2.set_title('排名分布热图', fontsize=self.chart_settings['font_size'] + 2)
        ax2.set_xticks(range(len(selected_models)))
        ax2.set_yticks(range(len(selected_models)))
        ax2.set_xticklabels([f'第{i + 1}名' for i in range(len(selected_models))])
        ax2.set_yticklabels([self.model_names[m] for m in selected_models])

        fig.colorbar(im, ax=ax2, label='频率')

        fig.tight_layout()

        # Embed in tkinter
        canvas = FigureCanvasTkAgg(fig, master=tab)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Add toolbar
        toolbar = NavigationToolbar2Tk(canvas, tab)
        toolbar.update()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        self.current_rank_fig = fig
        self.current_rank_canvas = canvas

    def refresh_charts(self):
        """Refresh all charts with current settings"""
        if self.df is not None and hasattr(self, 'metrics_results'):
            self.generate_charts()
            messagebox.showinfo("成功", "图表已使用当前设置刷新！")
        else:
            messagebox.showwarning("警告", "请先运行分析！")

    def save_current_chart(self):
        """Save the currently active chart"""
        current_tab = self.chart_notebook.tab(self.chart_notebook.select(), "text")

        fig_attributes = {
            "综合比较": ('current_comprehensive_fig', 'current_comprehensive_canvas'),
            "单个模型": ('current_individual_fig', 'current_individual_canvas'),
            "误差分析": ('current_error_fig', 'current_error_canvas'),
            "指标对比": ('current_metrics_fig', 'current_metrics_canvas'),
            "相关性矩阵": ('current_correlation_fig', 'current_correlation_canvas'),
            "排名对比": ('current_rank_fig', 'current_rank_canvas')
        }

        if current_tab in fig_attributes:
            fig_attr, canvas_attr = fig_attributes[current_tab]
            if hasattr(self, fig_attr):
                fig = getattr(self, fig_attr)
            else:
                messagebox.showwarning("警告", "当前标签页没有可保存的图表！")
                return
        else:
            messagebox.showwarning("警告", "当前标签页没有可保存的图表！")
            return

        filename = filedialog.asksaveasfilename(
            title="保存图表",
            defaultextension=".png",
            filetypes=[
                ("PNG files", "*.png"),
                ("PDF files", "*.pdf"),
                ("SVG files", "*.svg"),
                ("JPEG files", "*.jpg *.jpeg"),
                ("All files", "*.*")
            ]
        )

        if filename:
            try:
                fig.savefig(filename, dpi=self.chart_settings['dpi'], bbox_inches='tight')
                messagebox.showinfo("成功", f"图表已保存到:\n{filename}")
            except Exception as e:
                messagebox.showerror("错误", f"保存图表失败:\n{str(e)}")

    def save_all_charts(self):
        """Save all charts to a directory"""
        directory = filedialog.askdirectory(title="选择保存图表的目录")

        if not directory:
            return

        try:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            saved_files = []

            # Save all available charts
            chart_types = [
                ("综合比较", 'current_comprehensive_fig'),
                ("单个模型", 'current_individual_fig'),
                ("误差分析", 'current_error_fig'),
                ("指标对比", 'current_metrics_fig'),
                ("相关性矩阵", 'current_correlation_fig'),
                ("排名对比", 'current_rank_fig')
            ]

            for chart_name, fig_attr in chart_types:
                if hasattr(self, fig_attr):
                    fig = getattr(self, fig_attr)
                    filename = os.path.join(directory, f"{chart_name}_{timestamp}.png")
                    fig.savefig(filename, dpi=self.chart_settings['dpi'], bbox_inches='tight')
                    saved_files.append(f"{chart_name}.png")

            if saved_files:
                messagebox.showinfo("成功", f"所有图表已保存到:\n{directory}\n\n"
                                            f"保存的文件:\n" + "\n".join(saved_files))
            else:
                messagebox.showwarning("警告", "没有可保存的图表！")

        except Exception as e:
            messagebox.showerror("错误", f"保存图表失败:\n{str(e)}")

    def export_report(self):
        """Export a complete analysis report"""
        if self.df is None or not hasattr(self, 'metrics_results'):
            messagebox.showwarning("警告", "请先运行分析！")
            return

        filename = filedialog.asksaveasfilename(
            title="保存报告",
            defaultextension=".pdf",
            filetypes=[
                ("PDF files", "*.pdf"),
                ("HTML files", "*.html"),
                ("Text files", "*.txt"),
                ("All files", "*.*")
            ]
        )

        if not filename:
            return

        try:
            # Get results text
            results_text = self.results_text.get(1.0, tk.END)

            # Save based on file type
            if filename.endswith('.txt'):
                with open(filename, 'w', encoding='utf-8') as f:
                    f.write(results_text)
            elif filename.endswith('.html'):
                self.save_html_report(filename, results_text)
            elif filename.endswith('.pdf'):
                self.save_pdf_report(filename, results_text)
            else:
                # Default to text
                with open(filename, 'w', encoding='utf-8') as f:
                    f.write(results_text)

            messagebox.showinfo("成功", f"报告已保存到:\n{filename}")

        except Exception as e:
            messagebox.showerror("错误", f"保存报告失败:\n{str(e)}")

    def save_html_report(self, filename, results_text):
        """Save report as HTML"""
        # Prepare metrics data for HTML table
        metrics_html = ""
        selected_models = list(self.metrics_results.keys())

        # Get all metrics
        all_metrics = set()
        for metrics in self.metrics_results.values():
            all_metrics.update(metrics.keys())

        exclude_cols = {'Method', 'IsConstant', 'ConstantValue', 'Pearson p'}
        all_metrics = [m for m in all_metrics if m not in exclude_cols]

        for metric in sorted(all_metrics):
            metrics_html += "<tr>"
            metrics_html += f"<td>{metric}</td>"

            for model_key in selected_models:
                val = self.metrics_results[model_key].get(metric, 'N/A')
                if isinstance(val, (int, float, np.number)) and not np.isnan(val):
                    if metric in ['R^2', 'Pearson r']:
                        val_fmt = f"{val:.4f}"
                    elif metric in ['MAPE (%)']:
                        val_fmt = f"{val:.2f}%"
                    else:
                        val_fmt = f"{val:.4f}"
                else:
                    val_fmt = 'N/A'

                metrics_html += f"<td>{val_fmt}</td>"

            metrics_html += "</tr>"

        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>四模型保留时间预测分析报告</title>
            <style>
                body {{ font-family: 'Microsoft YaHei', Arial, sans-serif; margin: 40px; line-height: 1.6; }}
                h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
                h2 {{ color: #34495e; border-left: 5px solid #3498db; padding-left: 10px; }}
                h3 {{ color: #7f8c8d; }}
                .container {{ max-width: 1200px; margin: 0 auto; }}
                .metrics-table {{ width: 100%; border-collapse: collapse; margin: 20px 0; }}
                .metrics-table th, .metrics-table td {{ border: 1px solid #ddd; padding: 12px; text-align: center; }}
                .metrics-table th {{ background-color: #3498db; color: white; }}
                .metrics-table tr:nth-child(even) {{ background-color: #f9f9f9; }}
                .metrics-table tr:hover {{ background-color: #f5f5f5; }}
                .summary {{ background-color: #ecf0f1; padding: 20px; border-radius: 5px; margin: 20px 0; }}
                .warning {{ background-color: #ffeaa7; border: 1px solid #fdcb6e; padding: 15px; border-radius: 5px; margin: 15px 0; }}
                .success {{ background-color: #d1ecf1; border: 1px solid #bee5eb; padding: 15px; border-radius: 5px; margin: 15px 0; }}
                .model-card {{ display: inline-block; width: 200px; margin: 10px; padding: 15px; 
                              border-radius: 5px; text-align: center; color: white; }}
                .footer {{ margin-top: 40px; padding-top: 20px; border-top: 1px solid #ddd; 
                          color: #7f8c8d; font-size: 0.9em; }}
            </style>
        </head>
        <body>
            <div class="container">
                <h1>四模型保留时间预测分析报告</h1>
                <p><strong>生成时间:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
                <p><strong>数据集大小:</strong> {len(self.df)} 个化合物</p>
                <p><strong>分析模型数:</strong> {len(self.metrics_results)} 个</p>

                {"<div class='warning'><strong>⚠️ 警告:</strong> 存在常数预测模型！</div>"
        if any(m.get('IsConstant', False) for m in self.metrics_results.values()) else ""}

                <div class="summary">
                    <h2>模型概览</h2>
                    {"".join([f'<div class="model-card" style="background-color: {self.chart_settings["colors"][model_key]};">'
                              f'<h3>{self.model_names[model_key]}</h3>'
                              f'<p>R^2: {self.metrics_results[model_key].get("R^2", "N/A"):.4f}</p>'
                              f'<p>MAE: {self.metrics_results[model_key].get("MAE", "N/A"):.4f}</p>'
                              f'</div>'
                              for model_key in selected_models])}
                </div>

                <h2>详细指标对比</h2>
                <table class="metrics-table">
                    <thead>
                        <tr>
                            <th>指标</th>
                            {"".join([f'<th>{self.model_names[model_key]}</th>' for model_key in selected_models])}
                        </tr>
                    </thead>
                    <tbody>
                        {metrics_html}
                    </tbody>
                </table>

                <h2>分析结果摘要</h2>
                <pre style="background-color: #f8f9fa; padding: 20px; border-radius: 5px; overflow-x: auto;">{results_text}</pre>

                <div class="footer">
                    <p>报告生成工具: 四模型保留时间预测分析器 v1.0</p>
                    <p>© 2024 保留时间分析系统</p>
                </div>
            </div>
        </body>
        </html>
        """

        with open(filename, 'w', encoding='utf-8') as f:
            f.write(html_content)

    def save_pdf_report(self, filename, results_text):
        """Save report as PDF (simplified version)"""
        # For PDF, we'll use a simpler text format
        txt_filename = filename.replace('.pdf', '.txt')

        with open(txt_filename, 'w', encoding='utf-8') as f:
            f.write("=" * 100 + "\n")
            f.write("四模型保留时间预测分析报告\n")
            f.write("=" * 100 + "\n\n")
            f.write(f"生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"数据集大小: {len(self.df)} 个化合物\n")
            f.write(f"分析模型数: {len(self.metrics_results)} 个\n\n")

            # Check for constant models
            constant_models = []
            for model_key, metrics in self.metrics_results.items():
                if metrics.get('IsConstant', False):
                    constant_models.append(self.model_names.get(model_key, model_key))

            if constant_models:
                f.write("⚠️ 警告: 以下模型预测为常数:\n")
                for model in constant_models:
                    f.write(f"   • {model}\n")
                f.write("\n")

            f.write(results_text)

        messagebox.showinfo("提示", "PDF生成需要额外的库(reportlab)。\n"
                                    "已保存文本版本到:\n" + txt_filename)

    def copy_results(self):
        """Copy results to clipboard"""
        results = self.results_text.get(1.0, tk.END)
        self.root.clipboard_clear()
        self.root.clipboard_append(results)
        messagebox.showinfo("成功", "结果已复制到剪贴板！")

    def save_results_text(self):
        """Save results as text file"""
        filename = filedialog.asksaveasfilename(
            title="保存结果",
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )

        if filename:
            try:
                with open(filename, 'w', encoding='utf-8') as f:
                    f.write(self.results_text.get(1.0, tk.END))
                messagebox.showinfo("成功", f"结果已保存到:\n{filename}")
            except Exception as e:
                messagebox.showerror("错误", f"保存结果失败:\n{str(e)}")

    def save_results_csv(self):
        """Save metrics as CSV"""
        if not hasattr(self, 'metrics_results'):
            messagebox.showwarning("警告", "没有可保存的结果！")
            return

        filename = filedialog.asksaveasfilename(
            title="保存指标为CSV",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )

        if filename:
            try:
                # Create DataFrame with all metrics
                metrics_list = []
                for model_key, metrics in self.metrics_results.items():
                    metrics_dict = metrics.copy()
                    metrics_dict['模型名称'] = self.model_names.get(model_key, model_key)
                    metrics_list.append(metrics_dict)

                df_results = pd.DataFrame(metrics_list)
                df_results.to_csv(filename, index=False, encoding='utf-8-sig')
                messagebox.showinfo("成功", f"指标已保存到:\n{filename}")
            except Exception as e:
                messagebox.showerror("错误", f"保存指标失败:\n{str(e)}")

    def change_color(self, color_key):
        """Change color setting"""
        if color_key in self.chart_settings['colors']:
            current_color = self.chart_settings['colors'][color_key]
        else:
            current_color = '#ffffff'

        color = askcolor(color=current_color)[1]
        if color:
            if color_key in self.chart_settings['colors']:
                self.chart_settings['colors'][color_key] = color
            else:
                self.chart_settings[color_key] = color

    def update_setting(self, key, value):
        """Update chart setting"""
        try:
            if key in ['scatter_alpha', 'grid_alpha', 'bar_width']:
                self.chart_settings[key] = float(value)
            elif key in ['scatter_size', 'font_size', 'dpi']:
                self.chart_settings[key] = int(float(value))
            elif key == 'line_width':
                self.chart_settings[key] = float(value)
            elif key == 'marker_shape':
                self.chart_settings[key] = value
        except ValueError:
            pass

    def update_toggle(self, key, var):
        """Update toggle setting"""
        self.chart_settings[key] = var.get()

    def reset_settings(self):
        """Reset all settings to defaults"""
        defaults = {
            'scatter_alpha': 0.7,
            'scatter_size': 50,
            'line_width': 2,
            'grid_alpha': 0.3,
            'font_size': 10,
            'title_size': 12,
            'colors': {
                'rt_smrt_pred': '#1f77b4',
                'rt_M1_pred': '#ff7f0e',
                'rt_M2_pred': '#2ca02c',
                'rt_M3_pred': '#d62728',
                'rt_actual': '#9467bd',
                'rti_M3_pred': '#8c564b'
            },
            'show_grid': True,
            'show_legend': True,
            'show_r_squared': True,
            'show_regression_line': True,
            'show_error_bars': True,
            'show_trend_line': True,
            'marker_shape': 'o',
            'bar_width': 0.15,
            'dpi': 100
        }

        self.chart_settings.update(defaults)

        # Update UI widgets if they exist
        for key, value in defaults.items():
            if key == 'colors':
                continue
            if key in self.style_widgets:
                self.style_widgets[key].set(value)
            if key in self.toggle_vars:
                self.toggle_vars[key].set(value)
            if key == 'marker_shape':
                self.marker_var.set(value)

        messagebox.showinfo("设置已重置", "所有设置已恢复为默认值！")

    def save_settings(self):
        """Save chart settings to file"""
        filename = filedialog.asksaveasfilename(
            title="保存设置",
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
        )

        if filename:
            try:
                with open(filename, 'w') as f:
                    json.dump(self.chart_settings, f, indent=4, ensure_ascii=False)
                messagebox.showinfo("成功", f"设置已保存到:\n{filename}")
            except Exception as e:
                messagebox.showerror("错误", f"保存设置失败:\n{str(e)}")

    def load_settings(self):
        """Load chart settings from file"""
        settings_file = "chart_settings.json"
        if os.path.exists(settings_file):
            try:
                with open(settings_file, 'r') as f:
                    loaded_settings = json.load(f)
                self.chart_settings.update(loaded_settings)
            except:
                pass

    def on_closing(self):
        """Handle window closing"""
        # Save settings before closing
        try:
            with open("chart_settings.json", 'w') as f:
                json.dump(self.chart_settings, f, indent=4, ensure_ascii=False)
        except:
            pass

        self.root.destroy()


def main():
    root = tk.Tk()
    app = RetentionTimeAnalyzerApp(root)

    # Center window
    root.update_idletasks()
    width = root.winfo_width()
    height = root.winfo_height()
    x = (root.winfo_screenwidth() // 2) - (width // 2)
    y = (root.winfo_screenheight() // 2) - (height // 2)
    root.geometry(f'{width}x{height}+{x}+{y}')

    root.mainloop()


if __name__ == "__main__":
    main()