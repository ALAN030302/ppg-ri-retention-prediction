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


# Set fonts, 2
def setup_chinese_font():
    # , 2
    chinese_fonts = [
        'Microsoft YaHei', # -
        'Arial Unicode MS', # Unicode
        'DejaVu Sans', # ,
        'SimSun', #
        'SimHei', # -
        'NSimSun', #
        'FangSong', #
        'KaiTi', #
        'sans-serif'
    ]

    # translated note
    available_fonts = [f.name for f in fm.fontManager.ttflist]
    selected_font = None

    for font in chinese_fonts:
        if any(font.lower() in f.lower() for f in available_fonts):
            selected_font = font
            break

    if selected_font:
        # matplotlib
        plt.rcParams['font.sans-serif'] = [selected_font]
        plt.rcParams['axes.unicode_minus'] = False #

        # SimHei, Use2
        if 'SimHei' in selected_font:
            print(f"Use: {selected_font} (: , UseR^2R²)")
        else:
            print(f"Use: {selected_font}")
    else:
        print("Warning: , ")


# translated note
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
        self.root.title("modelretention_timeanalysis")
        self.root.geometry("1600x1000")

        # Initialize variables
        self.df = None
        self.model_names = {
            'rt_smrt_pred': '6model',
            'rt_M1_pred': 'conditionmodel',
            'rt_M2_pred': 'conditionmodel',
            'rt_M3_pred': 'retention_indexmodel'
        }

        self.chart_settings = {
            'scatter_alpha': 0.7,
            'scatter_size': 50,
            'line_width': 2,
            'grid_alpha': 0.3,
            'font_size': 10,
            'title_size': 12,
            'colors': {
                'rt_smrt_pred': '#1f77b4', #
                'rt_M1_pred': '#ff7f0e', #
                'rt_M2_pred': '#2ca02c', #
                'rt_M3_pred': '#d62728', #
                'rt_actual': '#9467bd', #
                'rti_M3_pred': '#8c564b' #
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
        control_notebook.add(file_frame, text="file")
        self.setup_file_controls(file_frame)

        # Analysis control tab
        analysis_frame = ttk.Frame(control_notebook)
        control_notebook.add(analysis_frame, text="analysis")
        self.setup_analysis_controls(analysis_frame)

        # Model selection tab
        model_frame = ttk.Frame(control_notebook)
        control_notebook.add(model_frame, text="model")
        self.setup_model_controls(model_frame)

        # Chart settings tab
        settings_frame = ttk.Frame(control_notebook)
        control_notebook.add(settings_frame, text="chart")
        self.setup_chart_settings(settings_frame)

        # Results tab
        results_frame = ttk.Frame(control_notebook)
        control_notebook.add(results_frame, text="results")
        self.setup_results_display(results_frame)

    def setup_file_controls(self, parent):
        # File selection
        file_group = ttk.LabelFrame(parent, text="datafile", padding=10)
        file_group.pack(fill=tk.X, padx=5, pady=5)

        ttk.Button(file_group, text="Excelfile...",
                   command=self.browse_file, width=20).pack(pady=5)

        self.file_path_var = tk.StringVar()
        ttk.Entry(file_group, textvariable=self.file_path_var,
                  state='readonly', width=50).pack(pady=5)

        ttk.Button(file_group, text="loaddata",
                   command=self.load_data, width=20).pack(pady=5)

        # Data preview
        preview_group = ttk.LabelFrame(parent, text="data", padding=10)
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
        settings_group = ttk.LabelFrame(parent, text="analysis", padding=10)
        settings_group.pack(fill=tk.X, padx=5, pady=5)

        # Metrics to calculate
        ttk.Label(settings_group, text="calculate:").pack(anchor=tk.W, pady=2)

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
        ttk.Button(settings_group, text="analysis",
                   command=self.run_analysis, width=20).pack(pady=10)

        # Chart type selection
        chart_group = ttk.LabelFrame(parent, text="chart", padding=10)
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
        model_group = ttk.LabelFrame(parent, text="comparemodel", padding=10)
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
        ttk.Checkbutton(model_group, text="retention_indexanalysis",
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
        color_group = ttk.LabelFrame(scrollable_frame, text="", padding=10)
        color_group.pack(fill=tk.X, padx=5, pady=5)

        colors = [
            ('6model', 'rt_smrt_pred'),
            ('conditionmodel', 'rt_M1_pred'),
            ('conditionmodel', 'rt_M2_pred'),
            ('retention_indexmodel', 'rt_M3_pred'),
            ('', 'rt_actual'),
            ('retention_index', 'rti_M3_pred')
        ]

        self.color_buttons = {}
        for i, (label, key) in enumerate(colors):
            ttk.Label(color_group, text=label).grid(row=i, column=0, padx=5, pady=2)
            btn = ttk.Button(color_group, text="", width=10,
                             command=lambda k=key: self.change_color(k))
            btn.grid(row=i, column=1, padx=5, pady=2)
            self.color_buttons[key] = btn

        # Alpha and size settings
        style_group = ttk.LabelFrame(scrollable_frame, text="", padding=10)
        style_group.pack(fill=tk.X, padx=5, pady=5)

        settings = [
            ('', 'scatter_alpha', 0.1, 1.0),
            ('', 'scatter_size', 10, 200),
            ('', 'line_width', 1, 5),
            ('', 'font_size', 8, 20),
            ('', 'grid_alpha', 0.0, 1.0),
            ('', 'bar_width', 0.05, 0.3),
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
        toggle_group = ttk.LabelFrame(scrollable_frame, text="", padding=10)
        toggle_group.pack(fill=tk.X, padx=5, pady=5)

        toggles = [
            ('', 'show_grid'),
            ('', 'show_legend'),
            ('R^2', 'show_r_squared'),
            ('', 'show_regression_line'),
            ('', 'show_error_bars'),
            ('', 'show_trend_line')
        ]

        self.toggle_vars = {}
        for i, (label, key) in enumerate(toggles):
            var = tk.BooleanVar(value=self.chart_settings[key])
            cb = ttk.Checkbutton(toggle_group, text=label, variable=var,
                                 command=lambda k=key, v=var: self.update_toggle(k, v))
            cb.grid(row=i // 2, column=i % 2, sticky=tk.W, padx=5, pady=2)
            self.toggle_vars[key] = var

        # Marker shape
        marker_group = ttk.LabelFrame(scrollable_frame, text="", padding=10)
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

        ttk.Button(button_frame, text="",
                   command=self.reset_settings).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="save",
                   command=self.save_settings).pack(side=tk.LEFT, padx=5)

    def setup_results_display(self, parent):
        # Results text area
        results_group = ttk.LabelFrame(parent, text="analysis results", padding=10)
        results_group.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        self.results_text = scrolledtext.ScrolledText(results_group,
                                                      wrap=tk.WORD,
                                                      font=('Consolas', 10),
                                                      height=20)
        self.results_text.pack(fill=tk.BOTH, expand=True, pady=5)

        # Export buttons
        button_frame = ttk.Frame(results_group)
        button_frame.pack(fill=tk.X, pady=5)

        ttk.Button(button_frame, text="results",
                   command=self.copy_results).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="save",
                   command=self.save_results_text).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Save as CSV",
                   command=self.save_results_csv).pack(side=tk.LEFT, padx=5)

    def setup_right_panel(self):
        # Notebook for different chart views
        self.chart_notebook = ttk.Notebook(self.right_panel)
        self.chart_notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Create initial tabs
        self.chart_tabs = {}
        tab_names = ["compare", "model", "analysis", "", "", "rank"]
        for tab_name in tab_names:
            tab = ttk.Frame(self.chart_notebook)
            self.chart_notebook.add(tab, text=tab_name)
            self.chart_tabs[tab_name] = tab

        # Chart control buttons
        control_frame = ttk.Frame(self.right_panel)
        control_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Button(control_frame, text="chart",
                   command=self.refresh_charts).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_frame, text="savechart",
                   command=self.save_current_chart).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_frame, text="savechart",
                   command=self.save_all_charts).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_frame, text="",
                   command=self.export_report).pack(side=tk.LEFT, padx=5)

    def browse_file(self):
        filename = filedialog.askopenfilename(
            title="Excelfile",
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
            messagebox.showwarning("Warning", "file！")
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
                messagebox.showerror("", f"column:\n{', '.join(missing_cols)}")
                return

            # Update preview tree
            self.update_preview_tree()

            # Enable analysis
            messagebox.showinfo("",
                                f"dataload！\n"
                                f"data: {len(self.df)}\n"
                                f"compound: {len(self.df)}\n"
                                f"model: {', '.join(self.model_names.values())}")

        except Exception as e:
            messagebox.showerror("", f"loadfilefailed:\n{str(e)}")

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
            messagebox.showwarning("Warning", "loaddata！")
            return

        try:
            # Clear previous results
            self.results_text.delete(1.0, tk.END)

            # Get selected models
            selected_models = [key for key, var in self.model_vars.items() if var.get()]

            if not selected_models:
                messagebox.showwarning("Warning", "modelcompare！")
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
                    'Method': 'retention_index',
                    'Pearson r': ri_corr,
                    'R^2': ri_corr ** 2
                }

            # Display results
            self.display_results()

            # Generate charts
            self.generate_charts()

            messagebox.showinfo("", "analysis！")

        except Exception as e:
            messagebox.showerror("", f"analysisfailed:\n{str(e)}")

    def calculate_metrics(self, y_true, y_pred, method_name):
        metrics = {'Method': method_name}

        y_true = np.array(y_true)
        y_pred = np.array(y_pred)

        # data
        if len(y_true) != len(y_pred):
            raise ValueError(f"data: y_true={len(y_true)}, y_pred={len(y_pred)}")

        # NaN
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

        # calculateR² - Use R² = 1 - (SS_res / SS_tot)
        ss_res = np.sum((y_true_valid - y_pred_valid) ** 2) #
        ss_tot = np.sum((y_true_valid - np.mean(y_true_valid)) ** 2) #

        if ss_tot == 0:
            # y_true
            if ss_res == 0:
                r2 = 1.0 #
            else:
                r2 = -np.inf #
        else:
            r2 = 1 - (ss_res / ss_tot)
            # R²
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
        results_text += "modelretention_timeanalysis\n"
        results_text += "=" * 100 + "\n\n"

        results_text += f"data: {len(self.df)} compound\n"
        results_text += f"data: {len(self.df.dropna(subset=['rt']))} \n"
        results_text += f"analysisDate: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n"

        # Check for constant models
        constant_models = []
        for model_key, metrics in self.metrics_results.items():
            if metrics.get('IsConstant', False):
                constant_models.append(self.model_names.get(model_key, model_key))

        if constant_models:
            results_text += "⚠️ Warning: model:\n"
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
        header = f"{'':<20}"
        for model_key in self.metrics_results.keys():
            model_name = self.model_names.get(model_key, model_key)
            header += f" {model_name[:15]:>15}"
        header += f" {'model':>15}\n"

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
                    # ,
                    best_idx = np.argmin(valid_values)
                else:  # R^2, Pearson r
                    # ,
                    best_idx = np.argmax(valid_values)

                best_model = list(self.metrics_results.keys())[best_idx]
                best_name = self.model_names.get(best_model, best_model)
                row += f" {best_name[:15]:>15}"
            else:
                row += f" {'N/A':>15}"

            results_text += row + "\n"

        # Add R² interpretation
        results_text += "\n" + "=" * 100 + "\n"
        results_text += "R^2:\n"
        results_text += "=" * 100 + "\n"

        for model_key, metrics in self.metrics_results.items():
            if 'R^2' in metrics:
                r2 = metrics['R^2']
                model_name = self.model_names.get(model_key, model_key)
                results_text += f"\n{model_name}: R^2 = {r2:.6f}\n"

                if np.isnan(r2):
                    results_text += " • R^2 = NaN: data, calculateR^2\n"
                elif r2 == -np.inf:
                    results_text += " • R^2 = -∞: modelUse\n"
                elif r2 < 0:
                    results_text += f" • R^2 = {r2:.3f}: modelUse\n"
                    results_text += f" (SS_res) = {np.sum((self.df['rt'] - self.df[model_key]) ** 2):.2f}\n"
                    results_text += f" (SS_tot) = {np.sum((self.df['rt'] - np.mean(self.df['rt'])) ** 2):.2f}\n"
                elif r2 == 0:
                    results_text += " • R^2 = 0: modelUse\n"
                elif r2 < 0.3:
                    results_text += f" • R^2 = {r2:.3f}: ({round(r2 * 100, 1)}%)\n"
                elif r2 < 0.7:
                    results_text += f" • R^2 = {r2:.3f}: ({round(r2 * 100, 1)}%)\n"
                else:
                    results_text += f" • R^2 = {r2:.3f}: ({round(r2 * 100, 1)}%)\n"

        # Summary conclusions
        results_text += "\n" + "=" * 100 + "\n"
        results_text += "analysis:\n"
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

        results_text += "modelrank:\n"
        for i, (model_key, count) in enumerate(sorted_models, 1):
            model_name = self.model_names.get(model_key, model_key)
            results_text += f"{i}. {model_name}: {len(all_metrics_sorted)}, {count}\n"

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
        tab = self.chart_tabs["compare"]

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
                r2_label = f'{label} (R^2={r2:.3f}) - '
            else:
                r2_label = f'{label} (R^2={r2:.3f})'

            # NaN
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
        ax.plot([min_val, max_val], [min_val, max_val], 'k:', alpha=0.3, label='')

        # Configure plot with explicit font properties
        ax.set_xlabel('retention_time', fontsize=self.chart_settings['font_size'])
        ax.set_ylabel('retention_time', fontsize=self.chart_settings['font_size'])
        ax.set_title('model', fontsize=self.chart_settings['font_size'] + 2)

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
                label='', zorder=5)

        # Plot predictions
        for model_key in selected_models:
            color = self.chart_settings['colors'][model_key]
            label = self.model_names[model_key]

            ax.plot(x, df_valid.loc[sorted_idx, model_key].values,
                    '--', color=color, alpha=0.7,
                    linewidth=self.chart_settings['line_width'],
                    label=label)

        # Configure plot
        ax.set_xlabel('compound (RT)', fontsize=self.chart_settings['font_size'])
        ax.set_ylabel('retention_time', fontsize=self.chart_settings['font_size'])
        ax.set_title('retention_time', fontsize=self.chart_settings['font_size'] + 2)

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
            ax.text(0.5, 0.5, 'data', transform=ax.transAxes,
                    ha='center', va='center', fontsize=self.chart_settings['font_size'])
            return

        # Create boxplot - Matplotlib 3.9+
        try:
            # Matplotlib 3.9+ Use tick_labels Parameters
            import matplotlib
            if hasattr(matplotlib, '__version__'):
                version = matplotlib.__version__
                if int(version.split('.')[0]) >= 3 and int(version.split('.')[1]) >= 9:
                    # UseParameters
                    bp = ax.boxplot(error_data, tick_labels=labels, patch_artist=True)
                else:
                    # UseParameters
                    bp = ax.boxplot(error_data, labels=labels, patch_artist=True)
            else:
                # UseParameters
                bp = ax.boxplot(error_data, labels=labels, patch_artist=True)
        except:
            # ,
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
        ax.set_xlabel('model', fontsize=self.chart_settings['font_size'])
        ax.set_ylabel('', fontsize=self.chart_settings['font_size'])
        ax.set_title('', fontsize=self.chart_settings['font_size'] + 2)

        if self.chart_settings['show_grid']:
            ax.grid(True, alpha=self.chart_settings['grid_alpha'], axis='y')

        # Add mean values as text
        for i, model_key in enumerate(selected_models):
            valid_mask = ~self.df['rt'].isna() & ~self.df[model_key].isna()
            errors = np.abs(self.df.loc[valid_mask, model_key] - self.df.loc[valid_mask, 'rt'])
            if len(errors) > 0:
                mean_error = errors.mean()
                ax.text(i + 1, ax.get_ylim()[1] * 0.95, f': {mean_error:.3f}',
                        ha='center', va='top', fontsize=self.chart_settings['font_size'] - 2)

    def create_individual_charts(self, selected_models):
        """Create individual charts for each model"""
        if not self.chart_vars['scatter_individual'].get():
            return

        tab = self.chart_tabs["model"]

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
            ax.text(0.5, 0.5, 'data', transform=ax.transAxes,
                    ha='center', va='center', fontsize=self.chart_settings['font_size'])
            ax.set_title(f'{model_name}results', fontsize=self.chart_settings['font_size'] + 2)
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
        ax.plot([min_val, max_val], [min_val, max_val], 'k:', alpha=0.3, label='')

        # Add R^2 text
        if self.chart_settings['show_r_squared']:
            r2 = self.metrics_results[model_key].get('R^2', 0)
            if np.isnan(r2):
                r2_text = f'R^2 = NaN (data)'
                box_color = 'gray'
            elif r2 == -np.inf:
                r2_text = f'R^2 = -∞ ()'
                box_color = 'red'
            elif r2 < 0:
                r2_text = f'R^2 = {r2:.3f} ()'
                box_color = 'red'
            elif r2 < 0.3:
                r2_text = f'R^2 = {r2:.3f} ()'
                box_color = 'orange'
            elif r2 < 0.7:
                r2_text = f'R^2 = {r2:.3f} ()'
                box_color = 'yellow'
            else:
                r2_text = f'R^2 = {r2:.3f} ()'
                box_color = 'green'

            ax.text(0.05, 0.95, r2_text,
                    transform=ax.transAxes,
                    fontsize=self.chart_settings['font_size'],
                    verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor=box_color, alpha=0.5))

        # Configure plot
        ax.set_xlabel('retention_time', fontsize=self.chart_settings['font_size'])
        ax.set_ylabel('retention_time', fontsize=self.chart_settings['font_size'])
        ax.set_title(f'{model_name}results', fontsize=self.chart_settings['font_size'] + 2)

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

        tab = self.chart_tabs["analysis"]

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
                   label='')

        # Configure plot
        ax.set_xlabel(' ( - )', fontsize=self.chart_settings['font_size'])
        ax.set_ylabel('', fontsize=self.chart_settings['font_size'])
        ax.set_title(' ()', fontsize=self.chart_settings['font_size'] + 2)

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
            ax.text(0.5, 0.5, 'data', transform=ax.transAxes,
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
        ax.set_xlabel('model', fontsize=self.chart_settings['font_size'])
        ax.set_ylabel('', fontsize=self.chart_settings['font_size'])
        ax.set_title('', fontsize=self.chart_settings['font_size'] + 2)
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
        ax.set_xlabel('retention_time', fontsize=self.chart_settings['font_size'])
        ax.set_ylabel(' ( - )', fontsize=self.chart_settings['font_size'])
        ax.set_title('', fontsize=self.chart_settings['font_size'] + 2)

        if self.chart_settings['show_grid']:
            ax.grid(True, alpha=self.chart_settings['grid_alpha'])

        if self.chart_settings['show_legend']:
            ax.legend(fontsize=self.chart_settings['font_size'] - 2,
                      loc='upper left', bbox_to_anchor=(1, 1))

    def create_metrics_comparison_charts(self, selected_models):
        """Create metrics comparison charts"""
        if not self.chart_vars['metrics_bar'].get():
            return

        tab = self.chart_tabs[""]

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

        # R^2R^2
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
        ax.set_xlabel('model', fontsize=self.chart_settings['font_size'])
        ax.set_ylabel(metric, fontsize=self.chart_settings['font_size'])
        ax.set_title(f'{metric}', fontsize=self.chart_settings['font_size'] + 2)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=15)

        if self.chart_settings['show_grid']:
            ax.grid(True, alpha=self.chart_settings['grid_alpha'], axis='y')

        # Add reference lines
        if metric == 'R^2':
            ax.axhline(y=1.0, color='green', linestyle=':', alpha=0.5, label='')
            ax.axhline(y=0.0, color='red', linestyle=':', alpha=0.5, label='')
        elif metric in ['MAE', 'RMSE', 'MAPE (%)', 'MedianAE', 'MaxError', 'StdError']:
            ax.axhline(y=0.0, color='black', linestyle='-', linewidth=0.5)

    def create_correlation_matrix(self, selected_models):
        """Create correlation matrix heatmap"""
        if not self.chart_vars['correlation_matrix'].get():
            return

        tab = self.chart_tabs[""]

        # Prepare data for correlation matrix
        columns = ['rt'] + selected_models
        if self.include_ri_var.get() and 'rti_M3_pred' in self.df.columns:
            columns.append('rti_M3_pred')

        # Remove rows with NaN values
        df_valid = self.df[columns].dropna()

        if len(df_valid) == 0:
            fig = Figure(figsize=(8, 6), dpi=self.chart_settings['dpi'])
            ax = fig.add_subplot(111)
            ax.text(0.5, 0.5, 'datacalculate', transform=ax.transAxes,
                    ha='center', va='center', fontsize=self.chart_settings['font_size'])
            ax.set_title('', fontsize=self.chart_settings['font_size'] + 2)

            canvas = FigureCanvasTkAgg(fig, master=tab)
            canvas.draw()
            canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
            return

        corr_data = df_valid.corr()

        # Rename columns for display
        display_names = {'rt': ''}
        display_names.update(self.model_names)
        if 'rti_M3_pred' in columns:
            display_names['rti_M3_pred'] = 'retention_index'

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
        ax.set_title('', fontsize=self.chart_settings['font_size'] + 2)

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

        tab = self.chart_tabs["rank"]

        # Remove rows with NaN values
        df_valid = self.df[['rt'] + selected_models].dropna()

        if len(df_valid) == 0:
            fig = Figure(figsize=(12, 5), dpi=self.chart_settings['dpi'])
            ax = fig.add_subplot(111)
            ax.text(0.5, 0.5, 'datacalculaterank', transform=ax.transAxes,
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

        ax1.set_xlabel('model', fontsize=self.chart_settings['font_size'])
        ax1.set_ylabel('rank (1=)', fontsize=self.chart_settings['font_size'])
        ax1.set_title('modelrank', fontsize=self.chart_settings['font_size'] + 2)
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

        ax2.set_xlabel('rank', fontsize=self.chart_settings['font_size'])
        ax2.set_ylabel('model', fontsize=self.chart_settings['font_size'])
        ax2.set_title('rank', fontsize=self.chart_settings['font_size'] + 2)
        ax2.set_xticks(range(len(selected_models)))
        ax2.set_yticks(range(len(selected_models)))
        ax2.set_xticklabels([f'{i + 1}' for i in range(len(selected_models))])
        ax2.set_yticklabels([self.model_names[m] for m in selected_models])

        fig.colorbar(im, ax=ax2, label='')

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
            messagebox.showinfo("", "chartUse！")
        else:
            messagebox.showwarning("Warning", "analysis！")

    def save_current_chart(self):
        """Save the currently active chart"""
        current_tab = self.chart_notebook.tab(self.chart_notebook.select(), "text")

        fig_attributes = {
            "compare": ('current_comprehensive_fig', 'current_comprehensive_canvas'),
            "model": ('current_individual_fig', 'current_individual_canvas'),
            "analysis": ('current_error_fig', 'current_error_canvas'),
            "": ('current_metrics_fig', 'current_metrics_canvas'),
            "": ('current_correlation_fig', 'current_correlation_canvas'),
            "rank": ('current_rank_fig', 'current_rank_canvas')
        }

        if current_tab in fig_attributes:
            fig_attr, canvas_attr = fig_attributes[current_tab]
            if hasattr(self, fig_attr):
                fig = getattr(self, fig_attr)
            else:
                messagebox.showwarning("Warning", "savechart！")
                return
        else:
            messagebox.showwarning("Warning", "savechart！")
            return

        filename = filedialog.asksaveasfilename(
            title="savechart",
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
                messagebox.showinfo("", f"chartsave:\n{filename}")
            except Exception as e:
                messagebox.showerror("", f"savechartfailed:\n{str(e)}")

    def save_all_charts(self):
        """Save all charts to a directory"""
        directory = filedialog.askdirectory(title="savechartdirectory")

        if not directory:
            return

        try:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            saved_files = []

            # Save all available charts
            chart_types = [
                ("compare", 'current_comprehensive_fig'),
                ("model", 'current_individual_fig'),
                ("analysis", 'current_error_fig'),
                ("", 'current_metrics_fig'),
                ("", 'current_correlation_fig'),
                ("rank", 'current_rank_fig')
            ]

            for chart_name, fig_attr in chart_types:
                if hasattr(self, fig_attr):
                    fig = getattr(self, fig_attr)
                    filename = os.path.join(directory, f"{chart_name}_{timestamp}.png")
                    fig.savefig(filename, dpi=self.chart_settings['dpi'], bbox_inches='tight')
                    saved_files.append(f"{chart_name}.png")

            if saved_files:
                messagebox.showinfo("", f"chartsave:\n{directory}\n\n"
                                            f"savefile:\n" + "\n".join(saved_files))
            else:
                messagebox.showwarning("Warning", "savechart！")

        except Exception as e:
            messagebox.showerror("", f"savechartfailed:\n{str(e)}")

    def export_report(self):
        """Export a complete analysis report"""
        if self.df is None or not hasattr(self, 'metrics_results'):
            messagebox.showwarning("Warning", "analysis！")
            return

        filename = filedialog.asksaveasfilename(
            title="save",
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

            messagebox.showinfo("", f"save:\n{filename}")

        except Exception as e:
            messagebox.showerror("", f"savefailed:\n{str(e)}")

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
            <title>modelretention_timeanalysis</title>
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
                <h1>modelretention_timeanalysis</h1>
                <p><strong>:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
                <p><strong>data:</strong> {len(self.df)} compound</p>
                <p><strong>analysismodel:</strong> {len(self.metrics_results)} </p>

                {"<div class='warning'><strong>⚠️ Warning:</strong> model！</div>"
        if any(m.get('IsConstant', False) for m in self.metrics_results.values()) else ""}

                <div class="summary">
                    <h2>model</h2>
                    {"".join([f'<div class="model-card" style="background-color: {self.chart_settings["colors"][model_key]};">'
                              f'<h3>{self.model_names[model_key]}</h3>'
                              f'<p>R^2: {self.metrics_results[model_key].get("R^2", "N/A"):.4f}</p>'
                              f'<p>MAE: {self.metrics_results[model_key].get("MAE", "N/A"):.4f}</p>'
                              f'</div>'
                              for model_key in selected_models])}
                </div>

                <h2></h2>
                <table class="metrics-table">
                    <thead>
                        <tr>
                            <th></th>
                            {"".join([f'<th>{self.model_names[model_key]}</th>' for model_key in selected_models])}
                        </tr>
                    </thead>
                    <tbody>
                        {metrics_html}
                    </tbody>
                </table>

                <h2>analysis results</h2>
                <pre style="background-color: #f8f9fa; padding: 20px; border-radius: 5px; overflow-x: auto;">{results_text}</pre>

                <div class="footer">
                    <p>: modelretention_timeanalysis v1.0</p>
                    <p>© 2024 retention_timeanalysis</p>
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
            f.write("modelretention_timeanalysis\n")
            f.write("=" * 100 + "\n\n")
            f.write(f": {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"data: {len(self.df)} compound\n")
            f.write(f"analysismodel: {len(self.metrics_results)} \n\n")

            # Check for constant models
            constant_models = []
            for model_key, metrics in self.metrics_results.items():
                if metrics.get('IsConstant', False):
                    constant_models.append(self.model_names.get(model_key, model_key))

            if constant_models:
                f.write("⚠️ Warning: model:\n")
                for model in constant_models:
                    f.write(f"   • {model}\n")
                f.write("\n")

            f.write(results_text)

        messagebox.showinfo("", "PDF(reportlab). \n"
                                    "saveVersion:\n" + txt_filename)

    def copy_results(self):
        """Copy results to clipboard"""
        results = self.results_text.get(1.0, tk.END)
        self.root.clipboard_clear()
        self.root.clipboard_append(results)
        messagebox.showinfo("", "results！")

    def save_results_text(self):
        """Save results as text file"""
        filename = filedialog.asksaveasfilename(
            title="saveresults",
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )

        if filename:
            try:
                with open(filename, 'w', encoding='utf-8') as f:
                    f.write(self.results_text.get(1.0, tk.END))
                messagebox.showinfo("", f"resultssave:\n{filename}")
            except Exception as e:
                messagebox.showerror("", f"saveresultsfailed:\n{str(e)}")

    def save_results_csv(self):
        """Save metrics as CSV"""
        if not hasattr(self, 'metrics_results'):
            messagebox.showwarning("Warning", "saveresults！")
            return

        filename = filedialog.asksaveasfilename(
            title="saveCSV",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )

        if filename:
            try:
                # Create DataFrame with all metrics
                metrics_list = []
                for model_key, metrics in self.metrics_results.items():
                    metrics_dict = metrics.copy()
                    metrics_dict['modelname'] = self.model_names.get(model_key, model_key)
                    metrics_list.append(metrics_dict)

                df_results = pd.DataFrame(metrics_list)
                df_results.to_csv(filename, index=False, encoding='utf-8-sig')
                messagebox.showinfo("", f"save:\n{filename}")
            except Exception as e:
                messagebox.showerror("", f"savefailed:\n{str(e)}")

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

        messagebox.showinfo("", "！")

    def save_settings(self):
        """Save chart settings to file"""
        filename = filedialog.asksaveasfilename(
            title="save",
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
        )

        if filename:
            try:
                with open(filename, 'w') as f:
                    json.dump(self.chart_settings, f, indent=4, ensure_ascii=False)
                messagebox.showinfo("", f"save:\n{filename}")
            except Exception as e:
                messagebox.showerror("", f"savefailed:\n{str(e)}")

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