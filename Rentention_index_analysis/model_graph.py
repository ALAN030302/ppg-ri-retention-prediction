import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import numpy as np
from matplotlib import gridspec
import os

# translated note
import matplotlib

matplotlib.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans'] # Set fonts
matplotlib.rcParams['axes.unicode_minus'] = False #


class ChempropMetricsVisualizer:
    def __init__(self, root):
        self.root = root
        self.root.title("Chemproptraining")
        self.root.geometry("1400x900")

        # translated note
        self.setup_styles()

        # data
        self.data1 = None
        self.data2 = None
        self.data1_name = "load"
        self.data2_name = "load"

        # GUI
        self.create_widgets()

    def setup_styles(self):
        """GUI"""
        self.style = ttk.Style()
        self.style.theme_use('clam')

        # translated note
        self.bg_color = '#f0f0f0'
        self.frame_bg = '#ffffff'
        self.accent_color = '#2c6fb7'
        self.secondary_color = '#6c757d'

        self.root.configure(bg=self.bg_color)

    def create_widgets(self):
        """GUI"""
        # translated note
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # translated note
        control_panel = ttk.LabelFrame(main_frame, text="", padding=15)
        control_panel.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10))

        # translated note
        display_panel = ttk.Frame(main_frame)
        display_panel.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        # ========== ==========
        # file1load
        file1_frame = ttk.LabelFrame(control_panel, text="file1", padding=10)
        file1_frame.pack(fill=tk.X, pady=(0, 10))

        self.file1_label = ttk.Label(file1_frame, text="loadfile", wraplength=250)
        self.file1_label.pack(pady=(0, 5))

        ttk.Button(file1_frame, text="loadfile1",
                   command=lambda: self.load_file(1)).pack(fill=tk.X)

        # file2load
        file2_frame = ttk.LabelFrame(control_panel, text="file2", padding=10)
        file2_frame.pack(fill=tk.X, pady=(0, 10))

        self.file2_label = ttk.Label(file2_frame, text="loadfile", wraplength=250)
        self.file2_label.pack(pady=(0, 5))

        ttk.Button(file2_frame, text="loadfile2",
                   command=lambda: self.load_file(2)).pack(fill=tk.X)

        # translated note
        label_frame = ttk.LabelFrame(control_panel, text="", padding=10)
        label_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(label_frame, text="1:").grid(row=0, column=0, sticky=tk.W, pady=5)
        self.label1_entry = ttk.Entry(label_frame)
        self.label1_entry.insert(0, "1")
        self.label1_entry.grid(row=0, column=1, padx=5, pady=5)

        ttk.Label(label_frame, text="2:").grid(row=1, column=0, sticky=tk.W, pady=5)
        self.label2_entry = ttk.Entry(label_frame)
        self.label2_entry.insert(0, "2")
        self.label2_entry.grid(row=1, column=1, padx=5, pady=5)

        # translated note
        viz_frame = ttk.LabelFrame(control_panel, text="", padding=10)
        viz_frame.pack(fill=tk.X, pady=(0, 10))

        self.smooth_var = tk.IntVar(value=1)
        ttk.Checkbutton(viz_frame, text="",
                        variable=self.smooth_var).pack(anchor=tk.W, pady=2)

        self.show_best_var = tk.IntVar(value=1)
        ttk.Checkbutton(viz_frame, text="epoch",
                        variable=self.show_best_var).pack(anchor=tk.W, pady=2)

        self.show_gap_var = tk.IntVar(value=1)
        ttk.Checkbutton(viz_frame, text="",
                        variable=self.show_gap_var).pack(anchor=tk.W, pady=2)

        ttk.Label(viz_frame, text=":").pack(anchor=tk.W, pady=(10, 2))
        self.window_size = ttk.Spinbox(viz_frame, from_=3, to=15, width=10)
        self.window_size.set(5)
        self.window_size.pack(anchor=tk.W, pady=(0, 10))

        # translated note
        button_frame = ttk.Frame(control_panel)
        button_frame.pack(fill=tk.X, pady=(10, 0))

        ttk.Button(button_frame, text="",
                   command=self.generate_visualization).pack(fill=tk.X, pady=5)

        ttk.Button(button_frame, text="save",
                   command=self.save_image).pack(fill=tk.X, pady=5)

        ttk.Button(button_frame, text="",
                   command=self.show_summary).pack(fill=tk.X, pady=5)

        ttk.Button(button_frame, text="data",
                   command=self.clear_all).pack(fill=tk.X, pady=5)

        # translated note
        ttk.Button(button_frame, text="",
                   command=self.set_font).pack(fill=tk.X, pady=5)

        # ========== ==========
        # chart
        self.chart_frame = ttk.LabelFrame(display_panel, text="training", padding=10)
        self.chart_frame.pack(fill=tk.BOTH, expand=True)

        # chart
        self.fig, self.ax = plt.subplots(figsize=(10, 6))
        # chart
        self.fig.suptitle("training", fontsize=14, fontproperties=self.get_font())

        self.canvas = FigureCanvasTkAgg(self.fig, self.chart_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # translated note
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.chart_frame)
        self.toolbar.update()

        # translated note
        self.status_bar = ttk.Label(self.root, text="", relief=tk.SUNKEN, anchor=tk.W)
        self.status_bar.pack(side=tk.BOTTOM, fill=tk.X)

    def get_font(self):
        """"""
        import matplotlib.font_manager as fm

        # translated note
        chinese_fonts = [
            'SimHei', #
            'Microsoft YaHei', #
            'FangSong', #
            'KaiTi', #
            'STSong', #
            'STKaiti', #
            'Arial Unicode MS',
            'DejaVu Sans'
        ]

        for font_name in chinese_fonts:
            try:
                # translated note
                font_path = fm.findfont(font_name, fallback_to_default=False)
                if font_path:
                    return fm.FontProperties(fname=font_path)
            except:
                continue

        # , Use
        return fm.FontProperties()

    def load_file(self, file_num):
        """loadCSVfile"""
        file_path = filedialog.askopenfilename(
            title=f"file{file_num}",
            filetypes=[("CSV files", "*.csv"), ("file", "*.*")]
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

                self.update_status(f"file{file_num}load: {os.path.basename(file_path)}")

            except Exception as e:
                messagebox.showerror("load", f"loadfile:\n{str(e)}")

    def load_and_process_metrics(self, file_path):
        """loadmetrics.csvfile"""
        df = pd.read_csv(file_path)

        # training
        val_rows = df[df['val/mse'].notna()].copy()
        train_rows = df[df['train_loss_epoch'].notna()].copy()

        # translated note
        val_rows = val_rows.reset_index(drop=True)
        train_rows = train_rows.reset_index(drop=True)

        # data
        processed_data = []
        for i in range(min(len(val_rows), len(train_rows))):
            epoch_data = {
                'epoch': val_rows.loc[i, 'epoch'],
                'step': val_rows.loc[i, 'step'],
                'val_mse': val_rows.loc[i, 'val/mse'],
                'val_loss': val_rows.loc[i, 'val_loss'],
                'train_loss_epoch': train_rows.loc[i, 'train_loss_epoch'],
            }

            # train_loss_stepdata
            if 'train_loss_step' in train_rows.columns and not pd.isna(train_rows.loc[i, 'train_loss_step']):
                epoch_data['train_loss_step'] = train_rows.loc[i, 'train_loss_step']

            processed_data.append(epoch_data)

        return pd.DataFrame(processed_data)

    def generate_visualization(self):
        """visualization charts"""
        if self.data1 is None:
            messagebox.showwarning("Warning", "loadfile")
            return

        try:
            # chart
            self.fig.clf()

            # translated note
            label1 = self.label1_entry.get() or "1"
            label2 = self.label2_entry.get() or "2"

            # loaddatachart
            if self.data2 is None:
                # translated note
                self.create_single_plot(self.data1, label1)
            else:
                # compare
                self.create_comparison_plot(self.data1, self.data2, label1, label2)

            # translated note
            font_prop = self.get_font()
            plt.rcParams['font.sans-serif'] = [font_prop.get_name()]

            # chart
            self.canvas.draw()
            self.update_status("")

        except Exception as e:
            messagebox.showerror("", f":\n{str(e)}")

    def create_single_plot(self, data, label):
        """chart"""
        gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1])

        # 1. training
        ax1 = self.fig.add_subplot(gs[0, 0])
        ax1.plot(data['epoch'], data['train_loss_epoch'], 'b-', linewidth=2, label='training', alpha=0.8)
        ax1.set_xlabel('Epoch')
        ax1.set_ylabel('training')
        ax1.set_title(f'{label} - training')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # 2.
        ax2 = self.fig.add_subplot(gs[0, 1])
        ax2.plot(data['epoch'], data['val_loss'], 'r-', linewidth=2, label='', alpha=0.8)
        ax2.set_xlabel('Epoch')
        ax2.set_ylabel('')
        ax2.set_title(f'{label} - ')
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        # 3. training
        ax3 = self.fig.add_subplot(gs[1, :])
        ax3.plot(data['epoch'], data['train_loss_epoch'], 'b-', linewidth=2, label='training', alpha=0.8)
        ax3.plot(data['epoch'], data['val_loss'], 'r-', linewidth=2, label='', alpha=0.8)
        ax3.set_xlabel('Epoch')
        ax3.set_ylabel('')
        ax3.set_title(f'{label} - training vs ')
        ax3.legend()
        ax3.grid(True, alpha=0.3)

        # epoch
        if self.show_best_var.get():
            best_epoch = data['val_loss'].idxmin()
            best_val = data.loc[best_epoch, 'val_loss']
            ax2.scatter([data.loc[best_epoch, 'epoch']], [best_val],
                        color='red', s=100, zorder=5, edgecolors='black')
            ax2.text(data.loc[best_epoch, 'epoch'], best_val,
                     f' \n epoch={int(data.loc[best_epoch, "epoch"])}',
                     verticalalignment='bottom')

            ax3.scatter([data.loc[best_epoch, 'epoch']], [best_val],
                        color='red', s=100, zorder=5, edgecolors='black')

        # translated note
        if self.smooth_var.get():
            window_size = int(self.window_size.get())
            if window_size > 0 and len(data) > window_size:
                train_smooth = data['train_loss_epoch'].rolling(window=window_size, center=True).mean()
                val_smooth = data['val_loss'].rolling(window=window_size, center=True).mean()

                ax3.plot(data['epoch'], train_smooth, 'b--', linewidth=1.5, label='training()', alpha=0.6)
                ax3.plot(data['epoch'], val_smooth, 'r--', linewidth=1.5, label='()', alpha=0.6)
                ax3.legend()

        self.fig.suptitle(f'{label} - Chemproptraining', fontsize=14)
        self.fig.tight_layout()

    def create_comparison_plot(self, data1, data2, label1, label2):
        """comparechart"""
        gs = gridspec.GridSpec(3, 2, height_ratios=[1, 1, 1])

        # 1. training
        ax1 = self.fig.add_subplot(gs[0, 0])
        ax1.plot(data1['epoch'], data1['train_loss_epoch'], 'b-', linewidth=2, label=label1, alpha=0.8)
        ax1.plot(data2['epoch'], data2['train_loss_epoch'], 'r-', linewidth=2, label=label2, alpha=0.8)
        ax1.set_xlabel('Epoch')
        ax1.set_ylabel('training')
        ax1.set_title('training')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # 2.
        ax2 = self.fig.add_subplot(gs[0, 1])
        ax2.plot(data1['epoch'], data1['val_loss'], 'b-', linewidth=2, label=label1, alpha=0.8)
        ax2.plot(data2['epoch'], data2['val_loss'], 'r-', linewidth=2, label=label2, alpha=0.8)
        ax2.set_xlabel('Epoch')
        ax2.set_ylabel('')
        ax2.set_title('')
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        # 3. training (model)
        ax3 = self.fig.add_subplot(gs[1, 0])
        ax3.plot(data1['epoch'], data1['train_loss_epoch'], 'b-', linewidth=2, label=f'{label1}training', alpha=0.8)
        ax3.plot(data1['epoch'], data1['val_loss'], 'b--', linewidth=2, label=f'{label1}', alpha=0.8)
        ax3.set_xlabel('Epoch')
        ax3.set_ylabel('')
        ax3.set_title(f'{label1} - training vs ')
        ax3.legend()
        ax3.grid(True, alpha=0.3)

        ax4 = self.fig.add_subplot(gs[1, 1])
        ax4.plot(data2['epoch'], data2['train_loss_epoch'], 'r-', linewidth=2, label=f'{label2}training', alpha=0.8)
        ax4.plot(data2['epoch'], data2['val_loss'], 'r--', linewidth=2, label=f'{label2}', alpha=0.8)
        ax4.set_xlabel('Epoch')
        ax4.set_ylabel('')
        ax4.set_title(f'{label2} - training vs ')
        ax4.legend()
        ax4.grid(True, alpha=0.3)

        # 4. analysis
        if self.show_gap_var.get():
            ax5 = self.fig.add_subplot(gs[2, 0])
            gap1 = data1['val_loss'] - data1['train_loss_epoch']
            gap2 = data2['val_loss'] - data2['train_loss_epoch']
            ax5.plot(data1['epoch'], gap1, 'b-', linewidth=2, label=label1, alpha=0.8)
            ax5.plot(data2['epoch'], gap2, 'r-', linewidth=2, label=label2, alpha=0.8)
            ax5.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
            ax5.set_xlabel('Epoch')
            ax5.set_ylabel(' - training')
            ax5.set_title(' ()')
            ax5.legend()
            ax5.grid(True, alpha=0.3)

        # 5.
        ax6 = self.fig.add_subplot(gs[2, 1])

        # translated note
        best_epoch1 = data1['val_loss'].idxmin()
        best_epoch2 = data2['val_loss'].idxmin()

        ax6.plot(data1['epoch'], data1['val_loss'], 'b-', linewidth=2, label=label1, alpha=0.8)
        ax6.plot(data2['epoch'], data2['val_loss'], 'r-', linewidth=2, label=label2, alpha=0.8)

        if self.show_best_var.get():
            ax6.scatter([data1.loc[best_epoch1, 'epoch']], [data1.loc[best_epoch1, 'val_loss']],
                        color='blue', s=100, zorder=5, edgecolors='black',
                        label=f'{label1}: Epoch {int(data1.loc[best_epoch1, "epoch"])}')
            ax6.scatter([data2.loc[best_epoch2, 'epoch']], [data2.loc[best_epoch2, 'val_loss']],
                        color='red', s=100, zorder=5, edgecolors='black',
                        label=f'{label2}: Epoch {int(data2.loc[best_epoch2, "epoch"])}')

        if self.smooth_var.get():
            window_size = int(self.window_size.get())
            if window_size > 0 and len(data1) > window_size:
                data1_val_smooth = data1['val_loss'].rolling(window=window_size, center=True).mean()
                data2_val_smooth = data2['val_loss'].rolling(window=window_size, center=True).mean()

                ax6.plot(data1['epoch'], data1['val_loss'], 'b-', linewidth=1, alpha=0.3, label=f'{label1}')
                ax6.plot(data1['epoch'], data1_val_smooth, 'b-', linewidth=2, label=f'{label1}', alpha=0.8)
                ax6.plot(data2['epoch'], data2['val_loss'], 'r-', linewidth=1, alpha=0.3, label=f'{label2}')
                ax6.plot(data2['epoch'], data2_val_smooth, 'r-', linewidth=2, label=f'{label2}', alpha=0.8)

        ax6.set_xlabel('Epoch')
        ax6.set_ylabel('')
        ax6.set_title(' ()')
        ax6.legend()
        ax6.grid(True, alpha=0.3)

        self.fig.suptitle('Chemproptraininganalysis', fontsize=14)
        self.fig.tight_layout()

    def save_image(self):
        """savefile"""
        if not hasattr(self, 'fig') or len(self.fig.axes) == 0:
            messagebox.showwarning("Warning", "save, ")
            return

        file_path = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[("PNG files", "*.png"), ("PDF files", "*.pdf"), ("file", "*.*")],
            title="save"
        )

        if file_path:
            try:
                # translated note
                font_prop = self.get_font()
                plt.rcParams['font.sans-serif'] = [font_prop.get_name()]

                self.fig.savefig(file_path, dpi=300, bbox_inches='tight')
                self.update_status(f"save: {os.path.basename(file_path)}")
                messagebox.showinfo("save", f"save:\n{file_path}")
            except Exception as e:
                messagebox.showerror("save", f"save:\n{str(e)}")

    def show_summary(self):
        """"""
        if self.data1 is None:
            messagebox.showwarning("Warning", "loadfile")
            return

        try:
            # translated note
            summary_window = tk.Toplevel(self.root)
            summary_window.title("")
            summary_window.geometry("800x500")

            # translated note
            text_widget = scrolledtext.ScrolledText(summary_window, wrap=tk.WORD)
            text_widget.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

            # translated note
            summary_text = "=" * 60 + "\n"
            summary_text += "Chemproptraining\n"
            summary_text += "=" * 60 + "\n\n"

            if self.data1 is not None:
                summary_text += self.create_summary_text(self.data1, self.label1_entry.get() or "1")
                summary_text += "\n" + "-" * 60 + "\n\n"

            if self.data2 is not None:
                summary_text += self.create_summary_text(self.data2, self.label2_entry.get() or "2")
                summary_text += "\n" + "-" * 60 + "\n\n"

            if self.data1 is not None and self.data2 is not None:
                summary_text += self.create_comparison_text(
                    self.data1, self.data2,
                    self.label1_entry.get() or "1",
                    self.label2_entry.get() or "2"
                )

            # translated note
            text_widget.insert(tk.END, summary_text)
            text_widget.configure(state='disabled')

            # save
            save_button = ttk.Button(summary_window, text="save",
                                     command=lambda: self.save_summary(summary_text))
            save_button.pack(pady=(0, 10))

        except Exception as e:
            messagebox.showerror("", f":\n{str(e)}")

    def create_summary_text(self, data, label):
        """"""
        best_idx = data['val_loss'].idxmin()
        final_idx = len(data) - 1

        text = f"【{label}】:\n"
        text += f"epoch: {len(data)}\n"
        text += f"epoch: {int(data.loc[best_idx, 'epoch'])} (: {data.loc[best_idx, 'val_loss']:.6f})\n"
        text += f"epoch: {int(data.loc[final_idx, 'epoch'])} (: {data.loc[final_idx, 'val_loss']:.6f})\n"
        text += f"training: {data['train_loss_epoch'].min():.6f} - {data['train_loss_epoch'].max():.6f}\n"
        text += f": {data['val_loss'].min():.6f} - {data['val_loss'].max():.6f}\n"
        text += f"training: {data['train_loss_epoch'].mean():.6f}\n"
        text += f": {data['val_loss'].mean():.6f}\n"

        # calculate (50%epoch)
        initial_loss = data['val_loss'].iloc[0]
        half_loss = initial_loss * 0.5
        faster_loss = data['val_loss'][data['val_loss'] <= half_loss]
        if len(faster_loss) > 0:
            convergence_epoch = faster_loss.index[0]
            text += f"50%epoch: {convergence_epoch}\n"

        return text

    def create_comparison_text(self, data1, data2, label1, label2):
        """compare"""
        best_idx1 = data1['val_loss'].idxmin()
        best_idx2 = data2['val_loss'].idxmin()

        text = f"【{label1} vs {label2}】compareanalysis:\n"
        text += f":\n"
        text += f"  {label1}: {data1.loc[best_idx1, 'val_loss']:.6f} (epoch {int(data1.loc[best_idx1, 'epoch'])})\n"
        text += f"  {label2}: {data2.loc[best_idx2, 'val_loss']:.6f} (epoch {int(data2.loc[best_idx2, 'epoch'])})\n"

        # model
        if data1.loc[best_idx1, 'val_loss'] < data2.loc[best_idx2, 'val_loss']:
            text += f": {label1} \n"
        elif data1.loc[best_idx1, 'val_loss'] > data2.loc[best_idx2, 'val_loss']:
            text += f": {label2} \n"
        else:
            text += ": model\n"

        # analysis ()
        std1 = data1['val_loss'].std()
        std2 = data2['val_loss'].std()
        text += f"\n ():\n"
        text += f"  {label1}: {std1:.6f}\n"
        text += f"  {label2}: {std2:.6f}\n"

        if std1 < std2:
            text += f": {label1} \n"
        elif std1 > std2:
            text += f": {label2} \n"

        return text

    def save_summary(self, summary_text):
        """savefile"""
        file_path = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("file", "*.*")],
            title="save"
        )

        if file_path:
            try:
                with open(file_path, 'w', encoding='utf-8') as f:
                    f.write(summary_text)
                messagebox.showinfo("save", f"save:\n{file_path}")
                self.update_status(f"save: {os.path.basename(file_path)}")
            except Exception as e:
                messagebox.showerror("save", f"save:\n{str(e)}")

    def clear_all(self):
        """data"""
        if messagebox.askyesno("", "data？"):
            self.data1 = None
            self.data2 = None
            self.data1_name = "load"
            self.data2_name = "load"
            self.file1_label.config(text=self.data1_name)
            self.file2_label.config(text=self.data2_name)

            # chart
            self.fig.clf()
            self.canvas.draw()

            self.update_status("data")

    def set_font(self):
        """"""
        font_window = tk.Toplevel(self.root)
        font_window.title("")
        font_window.geometry("400x300")

        ttk.Label(font_window, text=":").pack(pady=10)

        # column
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

        # column
        font_listbox = tk.Listbox(font_window, height=10)
        font_listbox.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        for font_name in sorted(font_names)[:50]: # 50
            font_listbox.insert(tk.END, font_name)

        def apply_font():
            selection = font_listbox.curselection()
            if selection:
                selected_font = font_listbox.get(selection[0])
                matplotlib.rcParams['font.sans-serif'] = [selected_font]
                messagebox.showinfo("", f": {selected_font}")
                font_window.destroy()
            else:
                messagebox.showwarning("Warning", "")

        ttk.Button(font_window, text="", command=apply_font).pack(pady=10)

    def update_status(self, message):
        """"""
        self.status_bar.config(text=f": {message}")
        self.root.update_idletasks()


def main():
    root = tk.Tk()
    app = ChempropMetricsVisualizer(root)
    root.mainloop()


if __name__ == "__main__":
    main()