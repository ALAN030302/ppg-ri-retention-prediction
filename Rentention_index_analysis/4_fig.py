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

matplotlib.use('TkAgg')
from tkinter.colorchooser import askcolor
import json
import os
from datetime import datetime
import sys


class RetentionTimeAnalyzerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Retention Time Prediction Analyzer")
        self.root.geometry("1400x900")

        # Initialize variables
        self.df = None
        self.chart_settings = {
            'scatter_alpha': 0.7,
            'scatter_size': 50,
            'line_width': 2,
            'grid_alpha': 0.3,
            'font_size': 10,
            'title_size': 12,
            'color_direct': '#1f77b4',  # blue
            'color_ri': '#ff7f0e',  # orange
            'color_actual': '#2ca02c',  # green
            'color_error': '#d62728',  # red
            'show_grid': True,
            'show_legend': True,
            'show_r_squared': True,
            'show_regression_line': True,
            'show_error_bars': True,
            'show_trend_line': True,
            'marker_shape': 'o'
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
        control_notebook.add(file_frame, text="File")
        self.setup_file_controls(file_frame)

        # Analysis control tab
        analysis_frame = ttk.Frame(control_notebook)
        control_notebook.add(analysis_frame, text="Analysis")
        self.setup_analysis_controls(analysis_frame)

        # Chart settings tab
        settings_frame = ttk.Frame(control_notebook)
        control_notebook.add(settings_frame, text="Chart Settings")
        self.setup_chart_settings(settings_frame)

        # Results tab
        results_frame = ttk.Frame(control_notebook)
        control_notebook.add(results_frame, text="Results")
        self.setup_results_display(results_frame)

    def setup_file_controls(self, parent):
        # File selection
        file_group = ttk.LabelFrame(parent, text="Data File", padding=10)
        file_group.pack(fill=tk.X, padx=5, pady=5)

        ttk.Button(file_group, text="Browse Excel File...",
                   command=self.browse_file, width=20).pack(pady=5)

        self.file_path_var = tk.StringVar()
        ttk.Entry(file_group, textvariable=self.file_path_var,
                  state='readonly', width=50).pack(pady=5)

        ttk.Button(file_group, text="Load Data",
                   command=self.load_data, width=20).pack(pady=5)

        # Data preview
        preview_group = ttk.LabelFrame(parent, text="Data Preview", padding=10)
        preview_group.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Treeview for data preview
        columns = ("SMILES", "rt_pred", "RTI", "rt_RI_pred", "rt")
        self.preview_tree = ttk.Treeview(preview_group, columns=columns,
                                         show="headings", height=10)

        # Configure columns
        for col in columns:
            self.preview_tree.heading(col, text=col)
            self.preview_tree.column(col, width=100)

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
        settings_group = ttk.LabelFrame(parent, text="Analysis Settings", padding=10)
        settings_group.pack(fill=tk.X, padx=5, pady=5)

        # Metrics to calculate
        ttk.Label(settings_group, text="Metrics to Calculate:").pack(anchor=tk.W, pady=2)

        self.metrics_vars = {
            'MAE': tk.BooleanVar(value=True),
            'RMSE': tk.BooleanVar(value=True),
            'R2': tk.BooleanVar(value=True),
            'Pearson': tk.BooleanVar(value=True),
            'MAPE': tk.BooleanVar(value=False),
            'MedianAE': tk.BooleanVar(value=False)
        }

        metrics_frame = ttk.Frame(settings_group)
        metrics_frame.pack(fill=tk.X, pady=5)

        for i, (metric, var) in enumerate(self.metrics_vars.items()):
            cb = ttk.Checkbutton(metrics_frame, text=metric, variable=var)
            cb.grid(row=i // 3, column=i % 3, sticky=tk.W, padx=5, pady=2)

        # Analysis button
        ttk.Button(settings_group, text="Run Analysis",
                   command=self.run_analysis, width=20).pack(pady=10)

        # Chart type selection
        chart_group = ttk.LabelFrame(parent, text="Chart Generation", padding=10)
        chart_group.pack(fill=tk.X, padx=5, pady=5)

        self.chart_vars = {
            'scatter': tk.BooleanVar(value=True),
            'error_dist': tk.BooleanVar(value=True),
            'abs_error': tk.BooleanVar(value=True),
            'trend': tk.BooleanVar(value=True),
            'residual': tk.BooleanVar(value=True),
            'metrics': tk.BooleanVar(value=True),
            'individual': tk.BooleanVar(value=True),
            'correlation': tk.BooleanVar(value=False)
        }

        for i, (chart, var) in enumerate(self.chart_vars.items()):
            cb = ttk.Checkbutton(chart_group, text=chart.replace('_', ' ').title(),
                                 variable=var)
            cb.grid(row=i // 2, column=i % 2, sticky=tk.W, padx=5, pady=2)

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
        color_group = ttk.LabelFrame(scrollable_frame, text="Colors", padding=10)
        color_group.pack(fill=tk.X, padx=5, pady=5)

        colors = [
            ('Direct Model', 'color_direct'),
            ('RI-based Model', 'color_ri'),
            ('Actual Values', 'color_actual'),
            ('Error', 'color_error')
        ]

        self.color_buttons = {}
        for i, (label, key) in enumerate(colors):
            ttk.Label(color_group, text=label).grid(row=i, column=0, padx=5, pady=2)
            btn = ttk.Button(color_group, text="Change", width=10,
                             command=lambda k=key: self.change_color(k))
            btn.grid(row=i, column=1, padx=5, pady=2)
            self.color_buttons[key] = btn

        # Alpha and size settings
        style_group = ttk.LabelFrame(scrollable_frame, text="Style", padding=10)
        style_group.pack(fill=tk.X, padx=5, pady=5)

        settings = [
            ('Scatter Alpha', 'scatter_alpha', 0.1, 1.0),
            ('Scatter Size', 'scatter_size', 10, 200),
            ('Line Width', 'line_width', 1, 5),
            ('Font Size', 'font_size', 8, 20),
            ('Grid Alpha', 'grid_alpha', 0.0, 1.0)
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
        toggle_group = ttk.LabelFrame(scrollable_frame, text="Display Options", padding=10)
        toggle_group.pack(fill=tk.X, padx=5, pady=5)

        toggles = [
            ('Show Grid', 'show_grid'),
            ('Show Legend', 'show_legend'),
            ('Show R²', 'show_r_squared'),
            ('Show Regression Line', 'show_regression_line'),
            ('Show Error Bars', 'show_error_bars'),
            ('Show Trend Line', 'show_trend_line')
        ]

        self.toggle_vars = {}
        for i, (label, key) in enumerate(toggles):
            var = tk.BooleanVar(value=self.chart_settings[key])
            cb = ttk.Checkbutton(toggle_group, text=label, variable=var,
                                 command=lambda k=key, v=var: self.update_toggle(k, v))
            cb.grid(row=i // 2, column=i % 2, sticky=tk.W, padx=5, pady=2)
            self.toggle_vars[key] = var

        # Marker shape
        marker_group = ttk.LabelFrame(scrollable_frame, text="Marker Style", padding=10)
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

        ttk.Button(button_frame, text="Reset to Defaults",
                   command=self.reset_settings).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Save Settings",
                   command=self.save_settings).pack(side=tk.LEFT, padx=5)

    def setup_results_display(self, parent):
        # Results text area
        results_group = ttk.LabelFrame(parent, text="Analysis Results", padding=10)
        results_group.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        self.results_text = scrolledtext.ScrolledText(results_group,
                                                      wrap=tk.WORD,
                                                      font=('Consolas', 10),
                                                      height=20)
        self.results_text.pack(fill=tk.BOTH, expand=True, pady=5)

        # Export buttons
        button_frame = ttk.Frame(results_group)
        button_frame.pack(fill=tk.X, pady=5)

        ttk.Button(button_frame, text="Copy Results",
                   command=self.copy_results).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Save Results as Text",
                   command=self.save_results_text).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Save Results as CSV",
                   command=self.save_results_csv).pack(side=tk.LEFT, padx=5)

    def setup_right_panel(self):
        # Notebook for different chart views
        self.chart_notebook = ttk.Notebook(self.right_panel)
        self.chart_notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Create initial tabs
        self.chart_tabs = {}
        for tab_name in ["Comparison Charts", "Individual Charts", "Correlation Matrix"]:
            tab = ttk.Frame(self.chart_notebook)
            self.chart_notebook.add(tab, text=tab_name)
            self.chart_tabs[tab_name] = tab

        # Chart control buttons
        control_frame = ttk.Frame(self.right_panel)
        control_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Button(control_frame, text="Refresh Charts",
                   command=self.refresh_charts).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_frame, text="Save Current Chart",
                   command=self.save_current_chart).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_frame, text="Save All Charts",
                   command=self.save_all_charts).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_frame, text="Export Report",
                   command=self.export_report).pack(side=tk.LEFT, padx=5)

    def browse_file(self):
        filename = filedialog.askopenfilename(
            title="Select Excel file",
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
            messagebox.showwarning("Warning", "Please select a file first!")
            return

        try:
            # Read the file
            if self.file_path_var.get().endswith('.csv'):
                self.df = pd.read_csv(self.file_path_var.get())
            else:
                self.df = pd.read_excel(self.file_path_var.get())

            # Rename columns for consistency
            if 'retention_index' in self.df.columns:
                self.df = self.df.rename(columns={'retention_index': 'rt_RI_pred'})

            # Update preview tree
            self.update_preview_tree()

            # Enable analysis
            messagebox.showinfo("Success",
                                f"Data loaded successfully!\n"
                                f"Rows: {len(self.df)}\n"
                                f"Columns: {', '.join(self.df.columns)}")

        except Exception as e:
            messagebox.showerror("Error", f"Failed to load file:\n{str(e)}")

    def update_preview_tree(self):
        # Clear existing items
        for item in self.preview_tree.get_children():
            self.preview_tree.delete(item)

        # Add new items (first 50 rows)
        for idx, row in self.df.head(50).iterrows():
            values = (
                str(row.get('smiles', ''))[:30],
                f"{row.get('rt_pred', 0):.4f}",
                f"{row.get('RTI', 0):.4f}",
                f"{row.get('rt_RI_pred', 0):.4f}",
                f"{row.get('rt', 0):.4f}"
            )
            self.preview_tree.insert("", tk.END, values=values)

    def run_analysis(self):
        if self.df is None:
            messagebox.showwarning("Warning", "Please load data first!")
            return

        try:
            # Clear previous results
            self.results_text.delete(1.0, tk.END)

            # Calculate metrics - Fix the R² calculation here
            self.metrics_rt_pred = self.calculate_metrics(
                self.df['rt'], self.df['rt_pred'],
                'Direct Model Prediction'
            )
            self.metrics_rt_RI = self.calculate_metrics(
                self.df['rt'], self.df['rt_RI_pred'],
                'RI-based Prediction'
            )

            # Display results
            self.display_results()

            # Generate charts
            self.generate_charts()

            messagebox.showinfo("Success", "Analysis completed successfully!")

        except Exception as e:
            messagebox.showerror("Error", f"Analysis failed:\n{str(e)}")

    def calculate_metrics(self, y_true, y_pred, method_name):
        metrics = {'Method': method_name}

        y_true = np.array(y_true)
        y_pred = np.array(y_pred)

        # Check whether y_pred is constant
        if np.all(y_pred == y_pred[0]):
            # If it is constant, the R² calculation can fail; use a custom calculation
            ss_res = np.sum((y_true - y_pred) ** 2)
            ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
            if ss_tot == 0:
                r2 = 1.0 # If y_true is also constant, this is a perfect prediction
            else:
                r2 = 1 - (ss_res / ss_tot)
                # R² can be negative, indicating performance worse than mean prediction
        else:
            # Use sklearn r2_score
            r2 = r2_score(y_true, y_pred)

        metrics['R²'] = r2

        if self.metrics_vars['MAE'].get():
            metrics['MAE'] = mean_absolute_error(y_true, y_pred)

        if self.metrics_vars['RMSE'].get():
            metrics['RMSE'] = np.sqrt(mean_squared_error(y_true, y_pred))

        if self.metrics_vars['Pearson'].get():
            # Check whether the data are suitable for calculating correlation coefficients
            if len(y_true) > 1 and np.std(y_true) > 0 and np.std(y_pred) > 0:
                pearson_r, p_value = stats.pearsonr(y_true, y_pred)
                metrics['Pearson r'] = pearson_r
                metrics['Pearson p'] = p_value
            else:
                metrics['Pearson r'] = np.nan
                metrics['Pearson p'] = np.nan

        if self.metrics_vars['MAPE'].get():
            # Mean Absolute Percentage Error
            # Avoid division by zero
            mask = y_true != 0
            if np.any(mask):
                mape = np.mean(np.abs((y_true[mask] - y_pred[mask]) / y_true[mask])) * 100
            else:
                mape = np.nan
            metrics['MAPE (%)'] = mape

        if self.metrics_vars['MedianAE'].get():
            # Median Absolute Error
            metrics['MedianAE'] = np.median(np.abs(y_true - y_pred))

        # Additional useful metrics
        metrics['Max Error'] = np.max(np.abs(y_true - y_pred))
        metrics['Std Error'] = np.std(y_true - y_pred)

        # Add a flag indicating whether the model is constant
        metrics['Is Constant'] = np.all(y_pred == y_pred[0])
        if metrics['Is Constant']:
            metrics['Constant Value'] = y_pred[0]

        return metrics

    def display_results(self):
        # Create comparison table
        results_text = "=" * 80 + "\n"
        results_text += "RETENTION TIME PREDICTION ANALYSIS\n"
        results_text += "=" * 80 + "\n\n"

        results_text += f"Dataset size: {len(self.df)} compounds\n"
        results_text += f"Analysis date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n"

        # Check whether the direct model is constant
        if self.metrics_rt_pred.get('Is Constant', False):
            results_text += "⚠️ WARNING: Direct Model predictions are CONSTANT!\n"
            results_text += f"   All predictions = {self.metrics_rt_pred.get('Constant Value', 'N/A'):.4f}\n\n"

        # Create comparison table
        methods = ['Direct Model Prediction', 'RI-based Prediction']
        metrics_list = [self.metrics_rt_pred, self.metrics_rt_RI]

        # Get all metric names
        all_metrics = set()
        for metrics in metrics_list:
            all_metrics.update(metrics.keys())
        # Remove columns that do not need to be displayed
        exclude_cols = {'Method', 'Is Constant', 'Constant Value', 'Pearson p'}
        all_metrics = all_metrics - exclude_cols

        # Create table header
        results_text += f"{'Metric':<20} {'Direct Model':>15} {'RI-based':>15} {'Better Method':>15}\n"
        results_text += "-" * 65 + "\n"

        # Add rows
        for metric in sorted(all_metrics):
            val1 = self.metrics_rt_pred.get(metric, 'N/A')
            val2 = self.metrics_rt_RI.get(metric, 'N/A')

            if val1 != 'N/A' and val2 != 'N/A' and not np.isnan(val1) and not np.isnan(val2):
                # Format values
                if isinstance(val1, (int, float, np.number)):
                    val1_fmt = f"{val1:.4f}"
                    val2_fmt = f"{val2:.4f}"

                    # Determine better method
                    if metric in ['MAE', 'RMSE', 'MAPE (%)', 'MedianAE', 'Max Error', 'Std Error']:
                        better = 'Direct' if val1 < val2 else 'RI-based'
                    else:  # R², Pearson r
                        better = 'Direct' if val1 > val2 else 'RI-based'
                else:
                    val1_fmt = str(val1)
                    val2_fmt = str(val2)
                    better = 'N/A'
            else:
                val1_fmt = 'N/A'
                val2_fmt = 'N/A'
                better = 'N/A'

            results_text += f"{metric:<20} {val1_fmt:>15} {val2_fmt:>15} {better:>15}\n"

        # Add p-value if calculated
        if 'Pearson p' in self.metrics_rt_pred and not np.isnan(self.metrics_rt_pred['Pearson p']):
            results_text += "\nPearson correlation p-values:\n"
            results_text += f"  Direct Model: {self.metrics_rt_pred['Pearson p']:.4e}\n"
            results_text += f"  RI-based: {self.metrics_rt_RI['Pearson p']:.4e}\n"

        # R² value interpretation
        results_text += "\n" + "=" * 80 + "\n"
        results_text += "R² VALUE INTERPRETATION:\n"
        results_text += "=" * 80 + "\n"

        r2_direct = self.metrics_rt_pred.get('R²', 0)
        r2_ri = self.metrics_rt_RI.get('R²', 0)

        results_text += f"Direct Model R² = {r2_direct:.6f}\n"
        if r2_direct < 0:
            results_text += "  • R² < 0 indicates model is WORSE than using simple mean\n"
            results_text += "  • Model has NO predictive power\n"
        elif r2_direct == 0:
            results_text += "  • R² = 0 indicates model performs same as using simple mean\n"
        elif 0 < r2_direct < 0.3:
            results_text += f"  • R² = {r2_direct:.3f} indicates VERY WEAK predictive power\n"
            results_text += f"  • Model explains only {r2_direct * 100:.1f}% of variance\n"
        elif 0.3 <= r2_direct < 0.7:
            results_text += f"  • R² = {r2_direct:.3f} indicates MODERATE predictive power\n"
            results_text += f"  • Model explains {r2_direct * 100:.1f}% of variance\n"
        else:
            results_text += f"  • R² = {r2_direct:.3f} indicates STRONG predictive power\n"
            results_text += f"  • Model explains {r2_direct * 100:.1f}% of variance\n"

        results_text += f"\nRI-based Model R² = {r2_ri:.6f}\n"
        if r2_ri < 0:
            results_text += "  • R² < 0 indicates model is WORSE than using simple mean\n"
            results_text += "  • Model has NO predictive power\n"
        elif r2_ri == 0:
            results_text += "  • R² = 0 indicates model performs same as using simple mean\n"
        elif 0 < r2_ri < 0.3:
            results_text += f"  • R² = {r2_ri:.3f} indicates WEAK predictive power\n"
            results_text += f"  • Model explains only {r2_ri * 100:.1f}% of variance\n"
        elif 0.3 <= r2_ri < 0.7:
            results_text += f"  • R² = {r2_ri:.3f} indicates MODERATE predictive power\n"
            results_text += f"  • Model explains {r2_ri * 100:.1f}% of variance\n"
        else:
            results_text += f"  • R² = {r2_ri:.3f} indicates STRONG predictive power\n"
            results_text += f"  • Model explains {r2_ri * 100:.1f}% of variance\n"

        # Summary conclusion
        results_text += "\n" + "=" * 80 + "\n"
        results_text += "SUMMARY CONCLUSIONS:\n"
        results_text += "=" * 80 + "\n"

        # Compare key metrics
        better_count = {'Direct': 0, 'RI-based': 0, 'Tie': 0}

        for metric in ['MAE', 'RMSE', 'R²', 'Pearson r']:
            if metric in self.metrics_rt_pred and metric in self.metrics_rt_RI:
                val1 = self.metrics_rt_pred[metric]
                val2 = self.metrics_rt_RI[metric]

                # Check whether the value is NaN
                if np.isnan(val1) or np.isnan(val2):
                    continue

                if metric in ['MAE', 'RMSE']:
                    if val1 < val2:
                        better_count['Direct'] += 1
                    elif val1 > val2:
                        better_count['RI-based'] += 1
                    else:
                        better_count['Tie'] += 1
                else:  # R², Pearson r
                    if val1 > val2:
                        better_count['Direct'] += 1
                    elif val1 < val2:
                        better_count['RI-based'] += 1
                    else:
                        better_count['Tie'] += 1

        if better_count['Direct'] > better_count['RI-based']:
            results_text += "Overall, the Direct Model Prediction performs better.\n"
        elif better_count['RI-based'] > better_count['Direct']:
            results_text += "Overall, the RI-based Prediction performs better.\n"
        else:
            results_text += "Both methods perform similarly.\n"

        # Add a warning about constant models
        if self.metrics_rt_pred.get('Is Constant', False):
            results_text += "\n⚠️ CRITICAL ISSUE: Direct Model output is constant!\n"
            results_text += "   • This model has not learned any patterns\n"
            results_text += "   • Possible causes:\n"
            results_text += "     - Model not properly trained\n"
            results_text += "     - Data leakage or preprocessing error\n"
            results_text += "     - Features not predictive of target\n"
            results_text += "   • RECOMMENDATION: Retrain or debug the Direct Model\n"

        self.results_text.insert(1.0, results_text)

    def generate_charts(self):
        # Clear all chart frames
        for tab in self.chart_tabs.values():
            for widget in tab.winfo_children():
                widget.destroy()

        # Generate charts based on selections
        if any(var.get() for var in self.chart_vars.values()):
            self.create_comparison_charts()

        if self.chart_vars['individual'].get():
            self.create_individual_charts()

        if self.chart_vars['correlation'].get():
            self.create_correlation_matrix()

    def create_comparison_charts(self):
        tab = self.chart_tabs["Comparison Charts"]

        # Determine which charts to create
        charts_to_create = []
        for chart_name, var in self.chart_vars.items():
            if var.get() and chart_name in ['scatter', 'error_dist', 'abs_error', 'trend', 'residual', 'metrics']:
                charts_to_create.append(chart_name)

        if not charts_to_create:
            return

        # Create appropriate grid
        n_charts = len(charts_to_create)
        n_cols = min(3, n_charts)
        n_rows = (n_charts + n_cols - 1) // n_cols

        fig = Figure(figsize=(5 * n_cols, 4 * n_rows), dpi=100)

        for idx, chart_type in enumerate(charts_to_create):
            ax = fig.add_subplot(n_rows, n_cols, idx + 1)

            if chart_type == 'scatter':
                self.create_scatter_plot(ax)
            elif chart_type == 'error_dist':
                self.create_error_distribution(ax)
            elif chart_type == 'abs_error':
                self.create_absolute_error(ax)
            elif chart_type == 'trend':
                self.create_trend_plot(ax)
            elif chart_type == 'residual':
                self.create_residual_plot(ax)
            elif chart_type == 'metrics':
                self.create_metrics_comparison(ax)

        fig.tight_layout(pad=3.0)

        # Embed in tkinter
        canvas = FigureCanvasTkAgg(fig, master=tab)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Add toolbar
        toolbar = NavigationToolbar2Tk(canvas, tab)
        toolbar.update()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        self.current_comparison_fig = fig
        self.current_comparison_canvas = canvas

    def create_scatter_plot(self, ax):
        """Create scatter plot comparing both predictions"""
        # Get settings
        alpha = self.chart_settings['scatter_alpha']
        size = self.chart_settings['scatter_size']
        marker = self.chart_settings['marker_shape']

        # Get the R² value
        r2_direct = self.metrics_rt_pred.get('R²', 0)
        r2_ri = self.metrics_rt_RI.get('R²', 0)

        # Format the R² display
        if r2_direct < 0:
            r2_direct_label = f'Direct Model (R²={r2_direct:.3f}) - WORSE than mean'
        else:
            r2_direct_label = f'Direct Model (R²={r2_direct:.3f})'

        if r2_ri < 0:
            r2_ri_label = f'RI-based (R²={r2_ri:.3f}) - WORSE than mean'
        else:
            r2_ri_label = f'RI-based (R²={r2_ri:.3f})'

        # Plot direct model predictions
        ax.scatter(self.df['rt'], self.df['rt_pred'],
                   alpha=alpha, s=size, marker=marker,
                   color=self.chart_settings['color_direct'],
                   label=r2_direct_label)

        # Plot RI-based predictions
        ax.scatter(self.df['rt'], self.df['rt_RI_pred'],
                   alpha=alpha, s=size, marker=marker,
                   color=self.chart_settings['color_ri'],
                   label=r2_ri_label)

        # Add regression lines if enabled
        if self.chart_settings['show_regression_line']:
            # Direct model regression
            if np.std(self.df['rt_pred']) > 0: # check whether there is variation
                z1 = np.polyfit(self.df['rt'], self.df['rt_pred'], 1)
                p1 = np.poly1d(z1)
                x_range = np.linspace(self.df['rt'].min(), self.df['rt'].max(), 100)
                ax.plot(x_range, p1(x_range), '--',
                        color=self.chart_settings['color_direct'],
                        linewidth=self.chart_settings['line_width'],
                        label=f'Direct fit: y={z1[0]:.3f}x+{z1[1]:.3f}')
            else:
                # If it is constant, draw a horizontal line
                constant_val = self.df['rt_pred'].iloc[0]
                ax.axhline(y=constant_val, color=self.chart_settings['color_direct'],
                           linestyle='--', linewidth=self.chart_settings['line_width'],
                           label=f'Direct: y={constant_val:.3f}')

            # RI-based regression
            if np.std(self.df['rt_RI_pred']) > 0:
                z2 = np.polyfit(self.df['rt'], self.df['rt_RI_pred'], 1)
                p2 = np.poly1d(z2)
                x_range = np.linspace(self.df['rt'].min(), self.df['rt'].max(), 100)
                ax.plot(x_range, p2(x_range), '--',
                        color=self.chart_settings['color_ri'],
                        linewidth=self.chart_settings['line_width'],
                        label=f'RI-based fit: y={z2[0]:.3f}x+{z2[1]:.3f}')

        # Configure plot
        ax.set_xlabel('Actual Retention Time')
        ax.set_ylabel('Predicted Retention Time')
        ax.set_title('Prediction Scatter Plot')

        if self.chart_settings['show_grid']:
            ax.grid(True, alpha=self.chart_settings['grid_alpha'])

        if self.chart_settings['show_legend']:
            ax.legend(fontsize=self.chart_settings['font_size'] - 2)

        # Add annotation for the ideal line
        ax.text(0.05, 0.95, 'Ideal: Points on diagonal line',
                transform=ax.transAxes, fontsize=9,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.3))

    def create_error_distribution(self, ax):
        """Create error distribution histogram"""
        # Calculate errors
        errors_direct = self.df['rt_pred'] - self.df['rt']
        errors_ri = self.df['rt_RI_pred'] - self.df['rt']

        # Create histogram
        bins = 30
        ax.hist(errors_direct, bins=bins, alpha=0.5,
                color=self.chart_settings['color_direct'],
                label='Direct Model', edgecolor='black', density=True)

        ax.hist(errors_ri, bins=bins, alpha=0.5,
                color=self.chart_settings['color_ri'],
                label='RI-based', edgecolor='black', density=True)

        # Add vertical line at zero
        ax.axvline(x=0, color='red', linestyle='--',
                   linewidth=self.chart_settings['line_width'],
                   label='Zero Error')

        # Add mean lines
        mean_direct = errors_direct.mean()
        mean_ri = errors_ri.mean()

        ax.axvline(x=mean_direct, color=self.chart_settings['color_direct'],
                   linestyle='-', linewidth=2,
                   label=f'Direct mean: {mean_direct:.3f}')

        ax.axvline(x=mean_ri, color=self.chart_settings['color_ri'],
                   linestyle='-', linewidth=2,
                   label=f'RI-based mean: {mean_ri:.3f}')

        # Configure plot
        ax.set_xlabel('Prediction Error (Predicted - Actual)')
        ax.set_ylabel('Density')
        ax.set_title('Error Distribution (Normalized)')

        if self.chart_settings['show_grid']:
            ax.grid(True, alpha=self.chart_settings['grid_alpha'])

        if self.chart_settings['show_legend']:
            ax.legend(fontsize=self.chart_settings['font_size'] - 2)

    def create_absolute_error(self, ax):
        """Create absolute error comparison"""
        # Calculate absolute errors
        abs_errors_direct = np.abs(self.df['rt_pred'] - self.df['rt'])
        abs_errors_ri = np.abs(self.df['rt_RI_pred'] - self.df['rt'])

        # Create violin plot for better visualization
        data = [abs_errors_direct, abs_errors_ri]
        positions = [1, 2]
        colors = [self.chart_settings['color_direct'], self.chart_settings['color_ri']]

        # Create violin plot
        vp = ax.violinplot(data, positions=positions, showmeans=True, showmedians=True)

        # Color the violins
        for i, pc in enumerate(vp['bodies']):
            pc.set_facecolor(colors[i])
            pc.set_alpha(0.7)

        # Color mean and median lines
        vp['cmeans'].set_colors(['black', 'black'])
        vp['cmedians'].set_colors(['white', 'white'])

        # Configure plot
        ax.set_xlabel('Prediction Method')
        ax.set_ylabel('Absolute Error')
        ax.set_title('Absolute Error Distribution')
        ax.set_xticks(positions)
        ax.set_xticklabels(['Direct Model', 'RI-based'])

        if self.chart_settings['show_grid']:
            ax.grid(True, alpha=self.chart_settings['grid_alpha'], axis='y')

        # Add statistics as text
        stats_text = (f"Direct Model:\n"
                      f"Mean: {abs_errors_direct.mean():.3f}\n"
                      f"Median: {np.median(abs_errors_direct):.3f}\n"
                      f"Std: {abs_errors_direct.std():.3f}")

        ax.text(0.7, 0.85, stats_text, transform=ax.transAxes,
                fontsize=9, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))

    def create_trend_plot(self, ax):
        """Create trend comparison plot"""
        # Sort by actual retention time
        sorted_idx = self.df['rt'].sort_values().index
        x = range(len(self.df))

        # Plot actual values
        ax.plot(x, self.df.loc[sorted_idx, 'rt'].values,
                color=self.chart_settings['color_actual'],
                linewidth=self.chart_settings['line_width'] * 1.5,
                label='Actual', zorder=3)

        # Plot predictions
        if self.chart_settings['show_trend_line']:
            ax.plot(x, self.df.loc[sorted_idx, 'rt_pred'].values,
                    '--', color=self.chart_settings['color_direct'],
                    linewidth=self.chart_settings['line_width'],
                    label='Direct Model', alpha=0.8, zorder=2)

            ax.plot(x, self.df.loc[sorted_idx, 'rt_RI_pred'].values,
                    '-.', color=self.chart_settings['color_ri'],
                    linewidth=self.chart_settings['line_width'],
                    label='RI-based', alpha=0.8, zorder=1)

        # Configure plot
        ax.set_xlabel('Compounds (sorted by actual RT)')
        ax.set_ylabel('Retention Time')
        ax.set_title('Retention Time Trend Comparison')

        if self.chart_settings['show_grid']:
            ax.grid(True, alpha=self.chart_settings['grid_alpha'])

        if self.chart_settings['show_legend']:
            ax.legend(fontsize=self.chart_settings['font_size'] - 2)

    def create_residual_plot(self, ax):
        """Create residual plot"""
        # Calculate residuals
        residuals_direct = self.df['rt_pred'] - self.df['rt']
        residuals_ri = self.df['rt_RI_pred'] - self.df['rt']

        # Plot residuals
        ax.scatter(self.df['rt'], residuals_direct,
                   alpha=self.chart_settings['scatter_alpha'],
                   s=self.chart_settings['scatter_size'],
                   color=self.chart_settings['color_direct'],
                   label='Direct Model', marker=self.chart_settings['marker_shape'])

        ax.scatter(self.df['rt'], residuals_ri,
                   alpha=self.chart_settings['scatter_alpha'],
                   s=self.chart_settings['scatter_size'],
                   color=self.chart_settings['color_ri'],
                   label='RI-based', marker=self.chart_settings['marker_shape'])

        # Add horizontal line at zero
        ax.axhline(y=0, color='red', linestyle='--',
                   linewidth=self.chart_settings['line_width'])

        # Add loess fit lines if enabled
        if self.chart_settings['show_regression_line']:
            # Simple moving average for trend
            window = max(3, len(self.df) // 20)

            sorted_idx = self.df['rt'].sort_values().index
            rt_sorted = self.df.loc[sorted_idx, 'rt']
            res_direct_sorted = residuals_direct.loc[sorted_idx]
            res_ri_sorted = residuals_ri.loc[sorted_idx]

            # Calculate moving averages
            ma_direct = res_direct_sorted.rolling(window=window, center=True).mean()
            ma_ri = res_ri_sorted.rolling(window=window, center=True).mean()

            ax.plot(rt_sorted, ma_direct, '--',
                    color=self.chart_settings['color_direct'],
                    linewidth=self.chart_settings['line_width'] * 1.2,
                    label='Direct trend')

            ax.plot(rt_sorted, ma_ri, '-.',
                    color=self.chart_settings['color_ri'],
                    linewidth=self.chart_settings['line_width'] * 1.2,
                    label='RI-based trend')

        # Configure plot
        ax.set_xlabel('Actual Retention Time')
        ax.set_ylabel('Residual (Predicted - Actual)')
        ax.set_title('Residual Plot')

        if self.chart_settings['show_grid']:
            ax.grid(True, alpha=self.chart_settings['grid_alpha'])

        if self.chart_settings['show_legend']:
            ax.legend(fontsize=self.chart_settings['font_size'] - 2)

    def create_metrics_comparison(self, ax):
        """Create metrics comparison bar chart"""
        # Select metrics to display
        metrics_to_show = []
        for metric, var in self.metrics_vars.items():
            if var.get():
                if metric == 'R2':
                    metrics_to_show.append('R²')
                elif metric == 'Pearson':
                    metrics_to_show.append('Pearson r')
                else:
                    metrics_to_show.append(metric)

        # Get values
        direct_vals = []
        ri_vals = []

        for metric in metrics_to_show:
            if metric == 'R²':
                direct_vals.append(self.metrics_rt_pred.get('R²', 0))
                ri_vals.append(self.metrics_rt_RI.get('R²', 0))
            elif metric == 'Pearson r':
                direct_vals.append(self.metrics_rt_pred.get('Pearson r', 0))
                ri_vals.append(self.metrics_rt_RI.get('Pearson r', 0))
            else:
                direct_vals.append(self.metrics_rt_pred.get(metric, 0))
                ri_vals.append(self.metrics_rt_RI.get(metric, 0))

        # Create bar chart
        x = np.arange(len(metrics_to_show))
        width = 0.35

        bars1 = ax.bar(x - width / 2, direct_vals, width,
                       color=self.chart_settings['color_direct'],
                       alpha=0.7, label='Direct Model')

        bars2 = ax.bar(x + width / 2, ri_vals, width,
                       color=self.chart_settings['color_ri'],
                       alpha=0.7, label='RI-based')

        # Add value labels on bars
        for bars in [bars1, bars2]:
            for bar in bars:
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width() / 2., height + 0.01,
                        f'{height:.3f}', ha='center', va='bottom',
                        fontsize=self.chart_settings['font_size'] - 2)

        # Configure plot
        ax.set_xlabel('Evaluation Metric')
        ax.set_ylabel('Value')
        ax.set_title('Performance Metrics Comparison')
        ax.set_xticks(x)
        ax.set_xticklabels(metrics_to_show, rotation=15)

        # Add grid
        if self.chart_settings['show_grid']:
            ax.grid(True, alpha=self.chart_settings['grid_alpha'], axis='y')

        if self.chart_settings['show_legend']:
            ax.legend(fontsize=self.chart_settings['font_size'] - 2)

        # Add a reference line
        ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)

        # For R², add a target line
        if 'R²' in metrics_to_show:
            idx = metrics_to_show.index('R²')
            ax.axhline(y=1.0, color='green', linestyle=':', alpha=0.5,
                       xmin=(idx - 0.4) / len(metrics_to_show),
                       xmax=(idx + 0.4) / len(metrics_to_show))

    def create_individual_charts(self):
        """Create individual charts for each method"""
        tab = self.chart_tabs["Individual Charts"]

        fig = Figure(figsize=(14, 6), dpi=100)

        # Direct model subplot
        ax1 = fig.add_subplot(121)
        self.create_method_scatter(ax1, self.df['rt'], self.df['rt_pred'],
                                   self.metrics_rt_pred, 'Direct Model Prediction',
                                   self.chart_settings['color_direct'])

        # RI-based model subplot
        ax2 = fig.add_subplot(122)
        self.create_method_scatter(ax2, self.df['rt'], self.df['rt_RI_pred'],
                                   self.metrics_rt_RI, 'RI-based Prediction',
                                   self.chart_settings['color_ri'])

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

    def create_method_scatter(self, ax, x_true, y_pred, metrics, title, color):
        """Create scatter plot for individual method"""
        alpha = self.chart_settings['scatter_alpha']
        size = self.chart_settings['scatter_size']
        marker = self.chart_settings['marker_shape']

        # Scatter plot
        ax.scatter(x_true, y_pred, alpha=alpha, s=size, marker=marker, color=color)

        # Regression line
        if self.chart_settings['show_regression_line']:
            if np.std(y_pred) > 0: # check whether there is variation
                z = np.polyfit(x_true, y_pred, 1)
                p = np.poly1d(z)
                x_range = np.linspace(min(x_true), max(x_true), 100)
                ax.plot(x_range, p(x_range), 'r--',
                        linewidth=self.chart_settings['line_width'] * 1.5,
                        label=f'y = {z[0]:.3f}x + {z[1]:.3f}')
            else:
                # If it is constant, draw a horizontal line
                constant_val = y_pred.iloc[0] if hasattr(y_pred, 'iloc') else y_pred[0]
                ax.axhline(y=constant_val, color='red',
                           linestyle='--', linewidth=self.chart_settings['line_width'] * 1.5,
                           label=f'y = {constant_val:.3f}')

        # Add R² text
        if self.chart_settings['show_r_squared']:
            r2 = metrics.get('R²', 0)
            if r2 < 0:
                r2_text = f'R² = {r2:.3f} (Worse than mean!)'
                box_color = 'red'
            elif r2 < 0.3:
                r2_text = f'R² = {r2:.3f} (Weak)'
                box_color = 'orange'
            elif r2 < 0.7:
                r2_text = f'R² = {r2:.3f} (Moderate)'
                box_color = 'yellow'
            else:
                r2_text = f'R² = {r2:.3f} (Strong)'
                box_color = 'green'

            ax.text(0.05, 0.95, r2_text,
                    transform=ax.transAxes,
                    fontsize=self.chart_settings['font_size'],
                    verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor=box_color, alpha=0.5))

        # Configure plot
        ax.set_xlabel('Actual Retention Time')
        ax.set_ylabel('Predicted Retention Time')
        ax.set_title(title)

        if self.chart_settings['show_grid']:
            ax.grid(True, alpha=self.chart_settings['grid_alpha'])

        if self.chart_settings['show_legend'] and self.chart_settings['show_regression_line']:
            ax.legend(fontsize=self.chart_settings['font_size'] - 2)

        # Add a diagonal reference line
        min_val = min(min(x_true), min(y_pred))
        max_val = max(max(x_true), max(y_pred))
        ax.plot([min_val, max_val], [min_val, max_val], 'k:', alpha=0.3, label='Ideal')

    def create_correlation_matrix(self):
        """Create correlation matrix heatmap"""
        tab = self.chart_tabs["Correlation Matrix"]

        # Calculate correlations
        corr_data = self.df[['rt', 'rt_pred', 'rt_RI_pred']].corr()

        fig = Figure(figsize=(8, 6), dpi=100)
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
        ax.set_xticklabels(corr_data.columns, rotation=45)
        ax.set_yticklabels(corr_data.columns)
        ax.set_title('Correlation Matrix')

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

    def refresh_charts(self):
        """Refresh all charts with current settings"""
        if self.df is not None and hasattr(self, 'metrics_rt_pred'):
            self.generate_charts()
            messagebox.showinfo("Success", "Charts refreshed with current settings!")
        else:
            messagebox.showwarning("Warning", "Please run analysis first!")

    def save_current_chart(self):
        """Save the currently active chart"""
        current_tab = self.chart_notebook.tab(self.chart_notebook.select(), "text")

        if current_tab == "Comparison Charts" and hasattr(self, 'current_comparison_fig'):
            fig = self.current_comparison_fig
        elif current_tab == "Individual Charts" and hasattr(self, 'current_individual_fig'):
            fig = self.current_individual_fig
        elif current_tab == "Correlation Matrix" and hasattr(self, 'current_correlation_fig'):
            fig = self.current_correlation_fig
        else:
            messagebox.showwarning("Warning", "No chart to save in current tab!")
            return

        filename = filedialog.asksaveasfilename(
            title="Save Chart",
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
                fig.savefig(filename, dpi=300, bbox_inches='tight')
                messagebox.showinfo("Success", f"Chart saved to:\n{filename}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to save chart:\n{str(e)}")

    def save_all_charts(self):
        """Save all charts to a directory"""
        directory = filedialog.askdirectory(title="Select Directory to Save Charts")

        if not directory:
            return

        try:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

            # Save comparison charts
            if hasattr(self, 'current_comparison_fig'):
                filename = os.path.join(directory, f"comparison_charts_{timestamp}.png")
                self.current_comparison_fig.savefig(filename, dpi=300, bbox_inches='tight')

            # Save individual charts
            if hasattr(self, 'current_individual_fig'):
                filename = os.path.join(directory, f"individual_charts_{timestamp}.png")
                self.current_individual_fig.savefig(filename, dpi=300, bbox_inches='tight')

            # Save correlation matrix
            if hasattr(self, 'current_correlation_fig'):
                filename = os.path.join(directory, f"correlation_matrix_{timestamp}.png")
                self.current_correlation_fig.savefig(filename, dpi=300, bbox_inches='tight')

            messagebox.showinfo("Success", f"All charts saved to:\n{directory}")

        except Exception as e:
            messagebox.showerror("Error", f"Failed to save charts:\n{str(e)}")

    def export_report(self):
        """Export a complete analysis report"""
        if self.df is None or not hasattr(self, 'metrics_rt_pred'):
            messagebox.showwarning("Warning", "Please run analysis first!")
            return

        filename = filedialog.asksaveasfilename(
            title="Save Report",
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

            messagebox.showinfo("Success", f"Report saved to:\n{filename}")

        except Exception as e:
            messagebox.showerror("Error", f"Failed to save report:\n{str(e)}")

    def save_html_report(self, filename, results_text):
        """Save report as HTML"""
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>Retention Time Analysis Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 40px; }}
                h1 {{ color: #333; border-bottom: 2px solid #333; }}
                pre {{ background-color: #f5f5f5; padding: 15px; border-radius: 5px; }}
                .metrics {{ display: flex; justify-content: space-around; margin: 20px 0; }}
                .metric-box {{ background-color: #e9ecef; padding: 15px; border-radius: 5px; text-align: center; }}
                .warning {{ background-color: #fff3cd; border: 1px solid #ffeaa7; padding: 10px; margin: 10px 0; }}
            </style>
        </head>
        <body>
            <h1>Retention Time Prediction Analysis Report</h1>
            <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            <p>Dataset size: {len(self.df)} compounds</p>

            {"<div class='warning'><strong>WARNING:</strong> Direct Model predictions are constant!</div>"
        if self.metrics_rt_pred.get('Is Constant', False) else ""}

            <h2>Results Summary</h2>
            <pre>{results_text}</pre>

            <h2>Key Metrics</h2>
            <div class="metrics">
                <div class="metric-box">
                    <h3>Direct Model</h3>
                    <p>R²: {self.metrics_rt_pred.get('R²', 'N/A'):.3f}</p>
                    <p>MAE: {self.metrics_rt_pred.get('MAE', 'N/A'):.3f}</p>
                </div>
                <div class="metric-box">
                    <h3>RI-based Model</h3>
                    <p>R²: {self.metrics_rt_RI.get('R²', 'N/A'):.3f}</p>
                    <p>MAE: {self.metrics_rt_RI.get('MAE', 'N/A'):.3f}</p>
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
        with open(filename.replace('.pdf', '.txt'), 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n")
            f.write("RETENTION TIME ANALYSIS REPORT\n")
            f.write("=" * 80 + "\n\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Dataset size: {len(self.df)} compounds\n\n")

            if self.metrics_rt_pred.get('Is Constant', False):
                f.write("⚠️ WARNING: Direct Model predictions are CONSTANT!\n")
                f.write(f"   All predictions = {self.metrics_rt_pred.get('Constant Value', 'N/A'):.4f}\n\n")

            f.write(results_text)

        # Note: For actual PDF generation, you would need additional libraries like reportlab
        messagebox.showinfo("Info", "PDF generation requires additional libraries.\n"
                                    "Text version saved instead.")

    def copy_results(self):
        """Copy results to clipboard"""
        results = self.results_text.get(1.0, tk.END)
        self.root.clipboard_clear()
        self.root.clipboard_append(results)
        messagebox.showinfo("Success", "Results copied to clipboard!")

    def save_results_text(self):
        """Save results as text file"""
        filename = filedialog.asksaveasfilename(
            title="Save Results",
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )

        if filename:
            try:
                with open(filename, 'w', encoding='utf-8') as f:
                    f.write(self.results_text.get(1.0, tk.END))
                messagebox.showinfo("Success", f"Results saved to:\n{filename}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to save results:\n{str(e)}")

    def save_results_csv(self):
        """Save metrics as CSV"""
        if not hasattr(self, 'metrics_rt_pred'):
            messagebox.showwarning("Warning", "No results to save!")
            return

        filename = filedialog.asksaveasfilename(
            title="Save Metrics as CSV",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )

        if filename:
            try:
                # Create DataFrame with both methods
                df_results = pd.DataFrame([self.metrics_rt_pred, self.metrics_rt_RI])
                df_results.to_csv(filename, index=False, encoding='utf-8-sig')
                messagebox.showinfo("Success", f"Metrics saved to:\n{filename}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to save metrics:\n{str(e)}")

    def change_color(self, color_key):
        """Change color setting"""
        color = askcolor(color=self.chart_settings[color_key])[1]
        if color:
            self.chart_settings[color_key] = color

    def update_setting(self, key, value):
        """Update chart setting"""
        try:
            if key in ['scatter_alpha', 'grid_alpha']:
                self.chart_settings[key] = float(value)
            elif key in ['scatter_size', 'font_size']:
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
            'color_direct': '#1f77b4',
            'color_ri': '#ff7f0e',
            'color_actual': '#2ca02c',
            'color_error': '#d62728',
            'show_grid': True,
            'show_legend': True,
            'show_r_squared': True,
            'show_regression_line': True,
            'show_error_bars': True,
            'show_trend_line': True,
            'marker_shape': 'o'
        }

        self.chart_settings.update(defaults)

        # Update UI widgets if they exist
        for key, value in defaults.items():
            if key in self.style_widgets:
                self.style_widgets[key].set(value)
            if key in self.toggle_vars:
                self.toggle_vars[key].set(value)
            if key == 'marker_shape':
                self.marker_var.set(value)

        messagebox.showinfo("Settings Reset", "All settings have been reset to defaults!")

    def save_settings(self):
        """Save chart settings to file"""
        filename = filedialog.asksaveasfilename(
            title="Save Settings",
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
        )

        if filename:
            try:
                with open(filename, 'w') as f:
                    json.dump(self.chart_settings, f, indent=4)
                messagebox.showinfo("Success", f"Settings saved to:\n{filename}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to save settings:\n{str(e)}")

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
                json.dump(self.chart_settings, f, indent=4)
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