import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import numpy as np
from matplotlib import gridspec
import os

# 解决中文显示问题
import matplotlib

matplotlib.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']  # 设置中文字体
matplotlib.rcParams['axes.unicode_minus'] = False  # 正确显示负号


class ChempropMetricsVisualizer:
    def __init__(self, root):
        self.root = root
        self.root.title("Chemprop训练指标可视化工具")
        self.root.geometry("1400x900")

        # 设置样式
        self.setup_styles()

        # 数据存储
        self.data1 = None
        self.data2 = None
        self.data1_name = "未加载"
        self.data2_name = "未加载"

        # 创建GUI
        self.create_widgets()

    def setup_styles(self):
        """设置GUI样式"""
        self.style = ttk.Style()
        self.style.theme_use('clam')

        # 自定义颜色
        self.bg_color = '#f0f0f0'
        self.frame_bg = '#ffffff'
        self.accent_color = '#2c6fb7'
        self.secondary_color = '#6c757d'

        self.root.configure(bg=self.bg_color)

    def create_widgets(self):
        """创建所有GUI组件"""
        # 主框架布局
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # 左侧控制面板
        control_panel = ttk.LabelFrame(main_frame, text="控制面板", padding=15)
        control_panel.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10))

        # 右侧显示面板
        display_panel = ttk.Frame(main_frame)
        display_panel.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        # ========== 控制面板内容 ==========
        # 文件1加载部分
        file1_frame = ttk.LabelFrame(control_panel, text="文件1", padding=10)
        file1_frame.pack(fill=tk.X, pady=(0, 10))

        self.file1_label = ttk.Label(file1_frame, text="未加载文件", wraplength=250)
        self.file1_label.pack(pady=(0, 5))

        ttk.Button(file1_frame, text="加载文件1",
                   command=lambda: self.load_file(1)).pack(fill=tk.X)

        # 文件2加载部分
        file2_frame = ttk.LabelFrame(control_panel, text="文件2", padding=10)
        file2_frame.pack(fill=tk.X, pady=(0, 10))

        self.file2_label = ttk.Label(file2_frame, text="未加载文件", wraplength=250)
        self.file2_label.pack(pady=(0, 5))

        ttk.Button(file2_frame, text="加载文件2",
                   command=lambda: self.load_file(2)).pack(fill=tk.X)

        # 标签设置
        label_frame = ttk.LabelFrame(control_panel, text="标签设置", padding=10)
        label_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(label_frame, text="实验1标签:").grid(row=0, column=0, sticky=tk.W, pady=5)
        self.label1_entry = ttk.Entry(label_frame)
        self.label1_entry.insert(0, "实验1")
        self.label1_entry.grid(row=0, column=1, padx=5, pady=5)

        ttk.Label(label_frame, text="实验2标签:").grid(row=1, column=0, sticky=tk.W, pady=5)
        self.label2_entry = ttk.Entry(label_frame)
        self.label2_entry.insert(0, "实验2")
        self.label2_entry.grid(row=1, column=1, padx=5, pady=5)

        # 可视化选项
        viz_frame = ttk.LabelFrame(control_panel, text="可视化选项", padding=10)
        viz_frame.pack(fill=tk.X, pady=(0, 10))

        self.smooth_var = tk.IntVar(value=1)
        ttk.Checkbutton(viz_frame, text="显示平滑曲线",
                        variable=self.smooth_var).pack(anchor=tk.W, pady=2)

        self.show_best_var = tk.IntVar(value=1)
        ttk.Checkbutton(viz_frame, text="标记最佳epoch",
                        variable=self.show_best_var).pack(anchor=tk.W, pady=2)

        self.show_gap_var = tk.IntVar(value=1)
        ttk.Checkbutton(viz_frame, text="显示泛化差距",
                        variable=self.show_gap_var).pack(anchor=tk.W, pady=2)

        ttk.Label(viz_frame, text="平滑窗口大小:").pack(anchor=tk.W, pady=(10, 2))
        self.window_size = ttk.Spinbox(viz_frame, from_=3, to=15, width=10)
        self.window_size.set(5)
        self.window_size.pack(anchor=tk.W, pady=(0, 10))

        # 按钮区域
        button_frame = ttk.Frame(control_panel)
        button_frame.pack(fill=tk.X, pady=(10, 0))

        ttk.Button(button_frame, text="生成可视化",
                   command=self.generate_visualization).pack(fill=tk.X, pady=5)

        ttk.Button(button_frame, text="保存图像",
                   command=self.save_image).pack(fill=tk.X, pady=5)

        ttk.Button(button_frame, text="显示统计摘要",
                   command=self.show_summary).pack(fill=tk.X, pady=5)

        ttk.Button(button_frame, text="清除所有数据",
                   command=self.clear_all).pack(fill=tk.X, pady=5)

        # 字体设置按钮
        ttk.Button(button_frame, text="设置字体",
                   command=self.set_font).pack(fill=tk.X, pady=5)

        # ========== 显示面板内容 ==========
        # 图表显示区域
        self.chart_frame = ttk.LabelFrame(display_panel, text="训练指标可视化", padding=10)
        self.chart_frame.pack(fill=tk.BOTH, expand=True)

        # 初始化图表
        self.fig, self.ax = plt.subplots(figsize=(10, 6))
        # 设置图表标题和标签的字体
        self.fig.suptitle("训练指标可视化", fontsize=14, fontproperties=self.get_font())

        self.canvas = FigureCanvasTkAgg(self.fig, self.chart_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # 添加工具栏
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.chart_frame)
        self.toolbar.update()

        # 状态栏
        self.status_bar = ttk.Label(self.root, text="就绪", relief=tk.SUNKEN, anchor=tk.W)
        self.status_bar.pack(side=tk.BOTTOM, fill=tk.X)

    def get_font(self):
        """获取中文字体"""
        import matplotlib.font_manager as fm

        # 尝试查找系统中可用的中文字体
        chinese_fonts = [
            'SimHei',  # 黑体
            'Microsoft YaHei',  # 微软雅黑
            'FangSong',  # 仿宋
            'KaiTi',  # 楷体
            'STSong',  # 华文宋体
            'STKaiti',  # 华文楷体
            'Arial Unicode MS',
            'DejaVu Sans'
        ]

        for font_name in chinese_fonts:
            try:
                # 检查字体是否可用
                font_path = fm.findfont(font_name, fallback_to_default=False)
                if font_path:
                    return fm.FontProperties(fname=font_path)
            except:
                continue

        # 如果找不到中文字体，使用默认字体
        return fm.FontProperties()

    def load_file(self, file_num):
        """加载CSV文件"""
        file_path = filedialog.askopenfilename(
            title=f"选择文件{file_num}",
            filetypes=[("CSV files", "*.csv"), ("所有文件", "*.*")]
        )

        if file_path:
            try:
                data = self.load_and_process_metrics(file_path)

                if file_num == 1:
                    self.data1 = data
                    self.data1_name = os.path.basename(file_path)
                    self.file1_label.config(text=self.data1_name)
                else:
                    self.data2 = data
                    self.data2_name = os.path.basename(file_path)
                    self.file2_label.config(text=self.data2_name)

                self.update_status(f"文件{file_num}加载成功: {os.path.basename(file_path)}")

            except Exception as e:
                messagebox.showerror("加载错误", f"加载文件时出错:\n{str(e)}")

    def load_and_process_metrics(self, file_path):
        """加载并处理metrics.csv文件"""
        df = pd.read_csv(file_path)

        # 分离有验证损失的行和有训练损失的行
        val_rows = df[df['val/mse'].notna()].copy()
        train_rows = df[df['train_loss_epoch'].notna()].copy()

        # 重置索引以便合并
        val_rows = val_rows.reset_index(drop=True)
        train_rows = train_rows.reset_index(drop=True)

        # 合并数据
        processed_data = []
        for i in range(min(len(val_rows), len(train_rows))):
            epoch_data = {
                'epoch': val_rows.loc[i, 'epoch'],
                'step': val_rows.loc[i, 'step'],
                'val_mse': val_rows.loc[i, 'val/mse'],
                'val_loss': val_rows.loc[i, 'val_loss'],
                'train_loss_epoch': train_rows.loc[i, 'train_loss_epoch'],
            }

            # 检查是否有train_loss_step数据
            if 'train_loss_step' in train_rows.columns and not pd.isna(train_rows.loc[i, 'train_loss_step']):
                epoch_data['train_loss_step'] = train_rows.loc[i, 'train_loss_step']

            processed_data.append(epoch_data)

        return pd.DataFrame(processed_data)

    def generate_visualization(self):
        """生成可视化图表"""
        if self.data1 is None:
            messagebox.showwarning("警告", "请先加载至少一个文件")
            return

        try:
            # 清除之前的图表
            self.fig.clf()

            # 获取标签
            label1 = self.label1_entry.get() or "实验1"
            label2 = self.label2_entry.get() or "实验2"

            # 根据加载的数据数量决定图表布局
            if self.data2 is None:
                # 只显示一个实验
                self.create_single_plot(self.data1, label1)
            else:
                # 比较两个实验
                self.create_comparison_plot(self.data1, self.data2, label1, label2)

            # 设置全局字体
            font_prop = self.get_font()
            plt.rcParams['font.sans-serif'] = [font_prop.get_name()]

            # 更新图表
            self.canvas.draw()
            self.update_status("可视化生成完成")

        except Exception as e:
            messagebox.showerror("生成错误", f"生成可视化时出错:\n{str(e)}")

    def create_single_plot(self, data, label):
        """创建单个实验的图表"""
        gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1])

        # 1. 训练损失
        ax1 = self.fig.add_subplot(gs[0, 0])
        ax1.plot(data['epoch'], data['train_loss_epoch'], 'b-', linewidth=2, label='训练损失', alpha=0.8)
        ax1.set_xlabel('Epoch')
        ax1.set_ylabel('训练损失')
        ax1.set_title(f'{label} - 训练损失')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # 2. 验证损失
        ax2 = self.fig.add_subplot(gs[0, 1])
        ax2.plot(data['epoch'], data['val_loss'], 'r-', linewidth=2, label='验证损失', alpha=0.8)
        ax2.set_xlabel('Epoch')
        ax2.set_ylabel('验证损失')
        ax2.set_title(f'{label} - 验证损失')
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        # 3. 训练与验证对比
        ax3 = self.fig.add_subplot(gs[1, :])
        ax3.plot(data['epoch'], data['train_loss_epoch'], 'b-', linewidth=2, label='训练损失', alpha=0.8)
        ax3.plot(data['epoch'], data['val_loss'], 'r-', linewidth=2, label='验证损失', alpha=0.8)
        ax3.set_xlabel('Epoch')
        ax3.set_ylabel('损失')
        ax3.set_title(f'{label} - 训练 vs 验证损失')
        ax3.legend()
        ax3.grid(True, alpha=0.3)

        # 标记最佳epoch
        if self.show_best_var.get():
            best_epoch = data['val_loss'].idxmin()
            best_val = data.loc[best_epoch, 'val_loss']
            ax2.scatter([data.loc[best_epoch, 'epoch']], [best_val],
                        color='red', s=100, zorder=5, edgecolors='black')
            ax2.text(data.loc[best_epoch, 'epoch'], best_val,
                     f' 最佳\n epoch={int(data.loc[best_epoch, "epoch"])}',
                     verticalalignment='bottom')

            ax3.scatter([data.loc[best_epoch, 'epoch']], [best_val],
                        color='red', s=100, zorder=5, edgecolors='black')

        # 平滑曲线
        if self.smooth_var.get():
            window_size = int(self.window_size.get())
            if window_size > 0 and len(data) > window_size:
                train_smooth = data['train_loss_epoch'].rolling(window=window_size, center=True).mean()
                val_smooth = data['val_loss'].rolling(window=window_size, center=True).mean()

                ax3.plot(data['epoch'], train_smooth, 'b--', linewidth=1.5, label='训练(平滑)', alpha=0.6)
                ax3.plot(data['epoch'], val_smooth, 'r--', linewidth=1.5, label='验证(平滑)', alpha=0.6)
                ax3.legend()

        self.fig.suptitle(f'{label} - Chemprop训练指标', fontsize=14)
        self.fig.tight_layout()

    def create_comparison_plot(self, data1, data2, label1, label2):
        """创建两个实验的比较图表"""
        gs = gridspec.GridSpec(3, 2, height_ratios=[1, 1, 1])

        # 1. 训练损失对比
        ax1 = self.fig.add_subplot(gs[0, 0])
        ax1.plot(data1['epoch'], data1['train_loss_epoch'], 'b-', linewidth=2, label=label1, alpha=0.8)
        ax1.plot(data2['epoch'], data2['train_loss_epoch'], 'r-', linewidth=2, label=label2, alpha=0.8)
        ax1.set_xlabel('Epoch')
        ax1.set_ylabel('训练损失')
        ax1.set_title('训练损失对比')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # 2. 验证损失对比
        ax2 = self.fig.add_subplot(gs[0, 1])
        ax2.plot(data1['epoch'], data1['val_loss'], 'b-', linewidth=2, label=label1, alpha=0.8)
        ax2.plot(data2['epoch'], data2['val_loss'], 'r-', linewidth=2, label=label2, alpha=0.8)
        ax2.set_xlabel('Epoch')
        ax2.set_ylabel('验证损失')
        ax2.set_title('验证损失对比')
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        # 3. 训练与验证对比（每个模型单独）
        ax3 = self.fig.add_subplot(gs[1, 0])
        ax3.plot(data1['epoch'], data1['train_loss_epoch'], 'b-', linewidth=2, label=f'{label1}训练', alpha=0.8)
        ax3.plot(data1['epoch'], data1['val_loss'], 'b--', linewidth=2, label=f'{label1}验证', alpha=0.8)
        ax3.set_xlabel('Epoch')
        ax3.set_ylabel('损失')
        ax3.set_title(f'{label1} - 训练 vs 验证')
        ax3.legend()
        ax3.grid(True, alpha=0.3)

        ax4 = self.fig.add_subplot(gs[1, 1])
        ax4.plot(data2['epoch'], data2['train_loss_epoch'], 'r-', linewidth=2, label=f'{label2}训练', alpha=0.8)
        ax4.plot(data2['epoch'], data2['val_loss'], 'r--', linewidth=2, label=f'{label2}验证', alpha=0.8)
        ax4.set_xlabel('Epoch')
        ax4.set_ylabel('损失')
        ax4.set_title(f'{label2} - 训练 vs 验证')
        ax4.legend()
        ax4.grid(True, alpha=0.3)

        # 4. 泛化差距分析
        if self.show_gap_var.get():
            ax5 = self.fig.add_subplot(gs[2, 0])
            gap1 = data1['val_loss'] - data1['train_loss_epoch']
            gap2 = data2['val_loss'] - data2['train_loss_epoch']
            ax5.plot(data1['epoch'], gap1, 'b-', linewidth=2, label=label1, alpha=0.8)
            ax5.plot(data2['epoch'], gap2, 'r-', linewidth=2, label=label2, alpha=0.8)
            ax5.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
            ax5.set_xlabel('Epoch')
            ax5.set_ylabel('验证 - 训练损失')
            ax5.set_title('泛化差距 (过拟合指标)')
            ax5.legend()
            ax5.grid(True, alpha=0.3)

        # 5. 平滑曲线和最佳点
        ax6 = self.fig.add_subplot(gs[2, 1])

        # 标记最佳点
        best_epoch1 = data1['val_loss'].idxmin()
        best_epoch2 = data2['val_loss'].idxmin()

        ax6.plot(data1['epoch'], data1['val_loss'], 'b-', linewidth=2, label=label1, alpha=0.8)
        ax6.plot(data2['epoch'], data2['val_loss'], 'r-', linewidth=2, label=label2, alpha=0.8)

        if self.show_best_var.get():
            ax6.scatter([data1.loc[best_epoch1, 'epoch']], [data1.loc[best_epoch1, 'val_loss']],
                        color='blue', s=100, zorder=5, edgecolors='black',
                        label=f'最佳{label1}: Epoch {int(data1.loc[best_epoch1, "epoch"])}')
            ax6.scatter([data2.loc[best_epoch2, 'epoch']], [data2.loc[best_epoch2, 'val_loss']],
                        color='red', s=100, zorder=5, edgecolors='black',
                        label=f'最佳{label2}: Epoch {int(data2.loc[best_epoch2, "epoch"])}')

        if self.smooth_var.get():
            window_size = int(self.window_size.get())
            if window_size > 0 and len(data1) > window_size:
                data1_val_smooth = data1['val_loss'].rolling(window=window_size, center=True).mean()
                data2_val_smooth = data2['val_loss'].rolling(window=window_size, center=True).mean()

                ax6.plot(data1['epoch'], data1['val_loss'], 'b-', linewidth=1, alpha=0.3, label=f'{label1}原始')
                ax6.plot(data1['epoch'], data1_val_smooth, 'b-', linewidth=2, label=f'{label1}平滑', alpha=0.8)
                ax6.plot(data2['epoch'], data2['val_loss'], 'r-', linewidth=1, alpha=0.3, label=f'{label2}原始')
                ax6.plot(data2['epoch'], data2_val_smooth, 'r-', linewidth=2, label=f'{label2}平滑', alpha=0.8)

        ax6.set_xlabel('Epoch')
        ax6.set_ylabel('验证损失')
        ax6.set_title('验证损失对比 (最佳点标记)')
        ax6.legend()
        ax6.grid(True, alpha=0.3)

        self.fig.suptitle('Chemprop训练指标对比分析', fontsize=14)
        self.fig.tight_layout()

    def save_image(self):
        """保存当前图像到文件"""
        if not hasattr(self, 'fig') or len(self.fig.axes) == 0:
            messagebox.showwarning("警告", "没有可保存的图像，请先生成可视化")
            return

        file_path = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[("PNG files", "*.png"), ("PDF files", "*.pdf"), ("所有文件", "*.*")],
            title="保存图像"
        )

        if file_path:
            try:
                # 设置字体以防止中文乱码
                font_prop = self.get_font()
                plt.rcParams['font.sans-serif'] = [font_prop.get_name()]

                self.fig.savefig(file_path, dpi=300, bbox_inches='tight')
                self.update_status(f"图像已保存: {os.path.basename(file_path)}")
                messagebox.showinfo("保存成功", f"图像已成功保存到:\n{file_path}")
            except Exception as e:
                messagebox.showerror("保存错误", f"保存图像时出错:\n{str(e)}")

    def show_summary(self):
        """显示统计摘要"""
        if self.data1 is None:
            messagebox.showwarning("警告", "请先加载至少一个文件")
            return

        try:
            # 创建摘要窗口
            summary_window = tk.Toplevel(self.root)
            summary_window.title("统计摘要")
            summary_window.geometry("800x500")

            # 创建文本框显示摘要
            text_widget = scrolledtext.ScrolledText(summary_window, wrap=tk.WORD)
            text_widget.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

            # 生成摘要文本
            summary_text = "=" * 60 + "\n"
            summary_text += "Chemprop训练指标统计摘要\n"
            summary_text += "=" * 60 + "\n\n"

            if self.data1 is not None:
                summary_text += self.create_summary_text(self.data1, self.label1_entry.get() or "实验1")
                summary_text += "\n" + "-" * 60 + "\n\n"

            if self.data2 is not None:
                summary_text += self.create_summary_text(self.data2, self.label2_entry.get() or "实验2")
                summary_text += "\n" + "-" * 60 + "\n\n"

            if self.data1 is not None and self.data2 is not None:
                summary_text += self.create_comparison_text(
                    self.data1, self.data2,
                    self.label1_entry.get() or "实验1",
                    self.label2_entry.get() or "实验2"
                )

            # 插入文本
            text_widget.insert(tk.END, summary_text)
            text_widget.configure(state='disabled')

            # 添加保存按钮
            save_button = ttk.Button(summary_window, text="保存摘要",
                                     command=lambda: self.save_summary(summary_text))
            save_button.pack(pady=(0, 10))

        except Exception as e:
            messagebox.showerror("摘要错误", f"生成摘要时出错:\n{str(e)}")

    def create_summary_text(self, data, label):
        """创建单个实验的摘要文本"""
        best_idx = data['val_loss'].idxmin()
        final_idx = len(data) - 1

        text = f"【{label}】统计摘要:\n"
        text += f"总epoch数: {len(data)}\n"
        text += f"最佳epoch: {int(data.loc[best_idx, 'epoch'])} (验证损失: {data.loc[best_idx, 'val_loss']:.6f})\n"
        text += f"最终epoch: {int(data.loc[final_idx, 'epoch'])} (验证损失: {data.loc[final_idx, 'val_loss']:.6f})\n"
        text += f"训练损失范围: {data['train_loss_epoch'].min():.6f} - {data['train_loss_epoch'].max():.6f}\n"
        text += f"验证损失范围: {data['val_loss'].min():.6f} - {data['val_loss'].max():.6f}\n"
        text += f"平均训练损失: {data['train_loss_epoch'].mean():.6f}\n"
        text += f"平均验证损失: {data['val_loss'].mean():.6f}\n"

        # 计算收敛速度（损失下降50%所需的epoch数）
        initial_loss = data['val_loss'].iloc[0]
        half_loss = initial_loss * 0.5
        faster_loss = data['val_loss'][data['val_loss'] <= half_loss]
        if len(faster_loss) > 0:
            convergence_epoch = faster_loss.index[0]
            text += f"损失下降50%所需epoch: {convergence_epoch}\n"

        return text

    def create_comparison_text(self, data1, data2, label1, label2):
        """创建两个实验的比较文本"""
        best_idx1 = data1['val_loss'].idxmin()
        best_idx2 = data2['val_loss'].idxmin()

        text = f"【{label1} vs {label2}】比较分析:\n"
        text += f"最佳验证损失:\n"
        text += f"  {label1}: {data1.loc[best_idx1, 'val_loss']:.6f} (epoch {int(data1.loc[best_idx1, 'epoch'])})\n"
        text += f"  {label2}: {data2.loc[best_idx2, 'val_loss']:.6f} (epoch {int(data2.loc[best_idx2, 'epoch'])})\n"

        # 判断哪个模型更好
        if data1.loc[best_idx1, 'val_loss'] < data2.loc[best_idx2, 'val_loss']:
            text += f"结论: {label1} 表现更好\n"
        elif data1.loc[best_idx1, 'val_loss'] > data2.loc[best_idx2, 'val_loss']:
            text += f"结论: {label2} 表现更好\n"
        else:
            text += "结论: 两个模型表现相当\n"

        # 稳定性分析（验证损失的波动性）
        std1 = data1['val_loss'].std()
        std2 = data2['val_loss'].std()
        text += f"\n验证损失标准差 (稳定性指标):\n"
        text += f"  {label1}: {std1:.6f}\n"
        text += f"  {label2}: {std2:.6f}\n"

        if std1 < std2:
            text += f"稳定性: {label1} 更稳定\n"
        elif std1 > std2:
            text += f"稳定性: {label2} 更稳定\n"

        return text

    def save_summary(self, summary_text):
        """保存摘要到文件"""
        file_path = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("所有文件", "*.*")],
            title="保存摘要"
        )

        if file_path:
            try:
                with open(file_path, 'w', encoding='utf-8') as f:
                    f.write(summary_text)
                messagebox.showinfo("保存成功", f"摘要已保存到:\n{file_path}")
                self.update_status(f"摘要已保存: {os.path.basename(file_path)}")
            except Exception as e:
                messagebox.showerror("保存错误", f"保存摘要时出错:\n{str(e)}")

    def clear_all(self):
        """清除所有数据"""
        if messagebox.askyesno("确认", "确定要清除所有数据吗？"):
            self.data1 = None
            self.data2 = None
            self.data1_name = "未加载"
            self.data2_name = "未加载"
            self.file1_label.config(text=self.data1_name)
            self.file2_label.config(text=self.data2_name)

            # 清除图表
            self.fig.clf()
            self.canvas.draw()

            self.update_status("所有数据已清除")

    def set_font(self):
        """设置字体对话框"""
        font_window = tk.Toplevel(self.root)
        font_window.title("设置字体")
        font_window.geometry("400x300")

        ttk.Label(font_window, text="选择字体:").pack(pady=10)

        # 字体选择列表
        import matplotlib.font_manager as fm
        font_list = fm.findSystemFonts(fontpaths=None, fontext='ttf')
        font_names = []

        for font_path in font_list:
            try:
                font_prop = fm.FontProperties(fname=font_path)
                font_name = font_prop.get_name()
                if font_name and font_name not in font_names:
                    font_names.append(font_name)
            except:
                pass

        # 创建列表框
        font_listbox = tk.Listbox(font_window, height=10)
        font_listbox.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        for font_name in sorted(font_names)[:50]:  # 只显示前50个字体
            font_listbox.insert(tk.END, font_name)

        def apply_font():
            selection = font_listbox.curselection()
            if selection:
                selected_font = font_listbox.get(selection[0])
                matplotlib.rcParams['font.sans-serif'] = [selected_font]
                messagebox.showinfo("成功", f"字体已设置为: {selected_font}")
                font_window.destroy()
            else:
                messagebox.showwarning("警告", "请选择一个字体")

        ttk.Button(font_window, text="应用字体", command=apply_font).pack(pady=10)

    def update_status(self, message):
        """更新状态栏"""
        self.status_bar.config(text=f"状态: {message}")
        self.root.update_idletasks()


def main():
    root = tk.Tk()
    app = ChempropMetricsVisualizer(root)
    root.mainloop()


if __name__ == "__main__":
    main()