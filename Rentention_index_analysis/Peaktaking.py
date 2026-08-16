#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
mzML文件峰提取与化合物匹配程序
功能：
1. 从mzML文件提取峰数据
2. 从Excel表格加载化合物m/z信息
3. 根据化合物m/z匹配峰数据
4. 将匹配结果保存到新的Excel表格
作者: 姚钱磊
版本: 6.0
日期: 2025
"""
import os
import sys
import threading
import traceback
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Any
import warnings

# 抑制警告
warnings.filterwarnings('ignore')

try:
    import pandas as pd
    import numpy as np
except ImportError:
    print("错误: 请先安装 pandas 和 numpy")
    print("安装命令: pip install pandas numpy")
    sys.exit(1)

# 尝试导入GUI库
try:
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox, scrolledtext
    from tkinter.font import Font

    GUI_AVAILABLE = True
except ImportError:
    GUI_AVAILABLE = False
    print("警告: tkinter未安装，GUI不可用")


class MzMLProcessor:
    """mzML文件处理器 - 核心处理逻辑"""

    def __init__(self, progress_callback=None, log_callback=None):
        """初始化处理器"""
        self.results_summary = []
        self.temp_files = []
        self.progress_callback = progress_callback
        self.log_callback = log_callback
        self.cancel_requested = False

    def log(self, message: str, level: str = "INFO"):
        """记录日志"""
        timestamp = datetime.now().strftime("%H:%M:%S")
        log_message = f"[{timestamp}] {message}"

        if self.log_callback:
            self.log_callback(log_message, level)
        else:
            print(log_message)

    def update_progress(self, value: int, message: str = ""):
        """更新进度"""
        if self.progress_callback:
            self.progress_callback(value, message)

    def check_mzml_file(self, file_path: Path) -> Tuple[bool, str]:
        """检查mzML文件是否有效"""
        try:
            if not file_path.exists():
                return False, f"文件不存在: {file_path}"

            file_size = file_path.stat().st_size
            if file_size == 0:
                return False, "文件为空"

            # 检查文件格式
            with open(file_path, 'rb') as f:
                header = f.read(500)
                header_str = header.decode('utf-8', errors='ignore')

                if '<mzML' in header_str or '<indexedmzML' in header_str:
                    return True, f"文件有效 ({file_size / 1024 / 1024:.2f} MB)"
                else:
                    return False, "不是有效的mzML文件"

        except Exception as e:
            return False, f"检查文件时出错: {str(e)}"

    def extract_peaks_with_pyopenms(self, mzml_path: str, intensity_threshold: float = 1000) -> Tuple[
        Optional[List[Dict]], str]:
        """使用pyOpenMS提取峰"""
        try:
            from pyopenms import MSExperiment, MzMLFile

            self.log(f"使用pyOpenMS加载文件: {Path(mzml_path).name}")
            exp = MSExperiment()
            MzMLFile().load(mzml_path, exp)

            if exp.size() == 0:
                return None, "加载了0个质谱图"

            peaks = []
            total_spectra = exp.size()

            for i, spectrum in enumerate(exp):
                if self.cancel_requested:
                    return None, "用户取消操作"

                if spectrum.size() > 0:
                    mz_array, intensity_array = spectrum.get_peaks()
                    rt = spectrum.getRT()

                    for j in range(len(mz_array)):
                        if intensity_array[j] > intensity_threshold:
                            peaks.append({
                                '谱图索引': i,
                                '保留时间(RT)': rt,
                                '质荷比(m/z)': float(mz_array[j]),
                                '强度': float(intensity_array[j]),
                                'MS级别': spectrum.getMSLevel(),
                                '特征ID': f"F{i:05d}_{j:04d}"
                            })

                # 更新进度
                if i % 10 == 0 or i == total_spectra - 1:
                    progress = int((i + 1) / total_spectra * 50)  # 前50%用于加载和处理
                    self.update_progress(progress, f"处理谱图 {i + 1}/{total_spectra}")

            return peaks, f"pyOpenMS提取成功: {len(peaks)}个峰, {exp.size()}个谱图"

        except ImportError:
            return None, "pyOpenMS未安装"
        except Exception as e:
            return None, f"pyOpenMS处理失败: {str(e)}"

    def extract_peaks_with_pymzml(self, mzml_path: str, intensity_threshold: float = 1000) -> Tuple[
        Optional[List[Dict]], str]:
        """使用pymzML提取峰"""
        try:
            import pymzml

            self.log(f"使用pymzML加载文件: {Path(mzml_path).name}")
            peaks = []

            run = pymzml.run.Reader(mzml_path)
            spectrum_count = 0
            total_spectra = None

            try:
                # 尝试获取总谱图数
                total_spectra = len(run)
            except:
                pass

            for spectrum in run:
                if self.cancel_requested:
                    return None, "用户取消操作"

                if spectrum.ms_level == 1:
                    mz_array = spectrum.mz
                    intensity_array = spectrum.i

                    if len(mz_array) > 0:
                        for idx, (mz, intensity) in enumerate(zip(mz_array, intensity_array)):
                            if intensity > intensity_threshold:
                                peaks.append({
                                    '谱图索引': spectrum_count,
                                    '保留时间(RT)': spectrum.scan_time[0] if hasattr(spectrum,
                                                                                     'scan_time') and spectrum.scan_time else spectrum_count * 0.5,
                                    '质荷比(m/z)': float(mz),
                                    '强度': float(intensity),
                                    'MS级别': 1,
                                    '特征ID': f"F{spectrum_count:05d}_{idx:04d}"
                                })

                        spectrum_count += 1

                # 更新进度
                if spectrum_count % 10 == 0:
                    progress_msg = f"处理谱图 {spectrum_count}"
                    if total_spectra:
                        progress = int(spectrum_count / total_spectra * 50)
                        progress_msg += f"/{total_spectra}"
                    else:
                        progress = min(50, spectrum_count // 2)

                    self.update_progress(progress, progress_msg)

            if peaks:
                return peaks, f"pymzML提取成功: {len(peaks)}个峰"
            else:
                return None, "pymzML未找到符合条件的峰"

        except ImportError:
            return None, "pymzML未安装"
        except Exception as e:
            return None, f"pymzML处理失败: {str(e)}"

    def save_results(self, peaks: List[Dict], output_path: Path, base_name: str, method_used: str) -> Dict[str, str]:
        """保存结果到CSV和Excel格式"""
        try:
            self.log(f"开始保存结果...")
            self.update_progress(60, "保存CSV文件...")

            # 创建DataFrame
            df = pd.DataFrame(peaks)

            # 保存为CSV
            csv_file = output_path / f"{base_name}_峰数据.csv"
            df.to_csv(csv_file, index=False, encoding='utf-8-sig')
            self.log(f"CSV已保存: {csv_file.name}")

            result_files = {'csv': str(csv_file)}

            # 保存为Excel (XLSX)
            self.update_progress(70, "保存Excel文件...")
            try:
                import openpyxl
                excel_file = output_path / f"{base_name}_峰数据.xlsx"
                df.to_excel(excel_file, index=False, sheet_name='峰数据')
                result_files['excel'] = str(excel_file)
                self.log(f"Excel(XLSX)已保存: {excel_file.name}")
            except ImportError:
                self.log("openpyxl未安装，跳过Excel(XLSX)保存", "WARNING")
            except Exception as e:
                self.log(f"保存Excel(XLSX)失败: {str(e)}", "WARNING")

            # 保存为Excel (XLS) - 只在小数据量时尝试
            self.update_progress(80, "保存统计信息...")
            if len(df) <= 10000:  # XLS格式最多支持65536行
                try:
                    xls_file = output_path / f"{base_name}_峰数据.xls"
                    df.to_excel(xls_file, index=False, sheet_name='峰数据')
                    result_files['xls'] = str(xls_file)
                    self.log(f"Excel(XLS)已保存: {xls_file.name}")
                except Exception as e:
                    self.log(f"保存Excel(XLS)失败: {str(e)}", "WARNING")

            # 保存统计信息
            self.update_progress(90, "生成统计报告...")
            stats_file = output_path / f"{base_name}_统计信息.txt"
            self.save_statistics(df, stats_file, method_used)
            result_files['stats'] = str(stats_file)

            # 保存Excel统计信息
            try:
                import openpyxl
                stats_excel = output_path / f"{base_name}_统计信息.xlsx"
                self.save_statistics_excel(df, stats_excel, method_used)
                result_files['stats_excel'] = str(stats_excel)
            except:
                pass

            self.update_progress(100, "保存完成!")
            return result_files

        except Exception as e:
            self.log(f"保存结果失败: {str(e)}", "ERROR")
            return {}

    def save_statistics(self, df: pd.DataFrame, stats_file: Path, method_used: str):
        """保存统计信息到文本文件"""
        try:
            with open(stats_file, 'w', encoding='utf-8') as f:
                f.write("=" * 70 + "\n")
                f.write("质谱峰提取统计报告\n")
                f.write("=" * 70 + "\n\n")

                f.write(f"文件名称: {stats_file.stem}\n")
                f.write(f"提取时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"提取方法: {method_used}\n")
                f.write(f"总峰数量: {len(df):,}\n\n")

                if len(df) > 0:
                    f.write("基本统计:\n")
                    f.write("-" * 50 + "\n")

                    # m/z统计
                    f.write(f"m/z统计:\n")
                    f.write(f"  最小值: {df['质荷比(m/z)'].min():.6f}\n")
                    f.write(f"  最大值: {df['质荷比(m/z)'].max():.6f}\n")
                    f.write(f"  平均值: {df['质荷比(m/z)'].mean():.6f}\n")
                    f.write(f"  中位数: {df['质荷比(m/z)'].median():.6f}\n")
                    f.write(f"  标准差: {df['质荷比(m/z)'].std():.6f}\n\n")

                    # RT统计
                    f.write(f"RT统计:\n")
                    f.write(f"  最小值: {df['保留时间(RT)'].min():.2f}秒\n")
                    f.write(f"  最大值: {df['保留时间(RT)'].max():.2f}秒\n")
                    f.write(f"  平均值: {df['保留时间(RT)'].mean():.2f}秒\n")
                    f.write(f"  中位数: {df['保留时间(RT)'].median():.2f}秒\n\n")

                    # 强度统计
                    f.write(f"强度统计:\n")
                    f.write(f"  最小值: {df['强度'].min():.2e}\n")
                    f.write(f"  最大值: {df['强度'].max():.2e}\n")
                    f.write(f"  平均值: {df['强度'].mean():.2e}\n")
                    f.write(f"  中位数: {df['强度'].median():.2e}\n")
                    f.write(f"  总和: {df['强度'].sum():.2e}\n\n")

                    # 最强峰
                    f.write(f"最强的前10个峰:\n")
                    f.write("-" * 50 + "\n")
                    top_peaks = df.nlargest(10, '强度')
                    for i, (_, row) in enumerate(top_peaks.iterrows()):
                        f.write(
                            f"{i + 1:2d}. m/z={row['质荷比(m/z)']:9.6f}, RT={row['保留时间(RT)']:7.2f}秒, 强度={row['强度']:12.2e}\n")

            self.log(f"统计信息已保存: {stats_file.name}")

        except Exception as e:
            self.log(f"保存统计信息失败: {str(e)}", "ERROR")

    def save_statistics_excel(self, df: pd.DataFrame, excel_file: Path, method_used: str):
        """保存统计信息到Excel"""
        try:
            import openpyxl

            with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
                # 基本统计表
                basic_stats = pd.DataFrame({
                    '统计项': ['总峰数', 'm/z最小值', 'm/z最大值', 'm/z平均值', 'm/z中位数',
                               'RT最小值(秒)', 'RT最大值(秒)', 'RT平均值(秒)', 'RT中位数(秒)',
                               '强度最小值', '强度最大值', '强度平均值', '强度中位数'],
                    '数值': [len(df),
                             df['质荷比(m/z)'].min() if len(df) > 0 else 0,
                             df['质荷比(m/z)'].max() if len(df) > 0 else 0,
                             df['质荷比(m/z)'].mean() if len(df) > 0 else 0,
                             df['质荷比(m/z)'].median() if len(df) > 0 else 0,
                             df['保留时间(RT)'].min() if len(df) > 0 else 0,
                             df['保留时间(RT)'].max() if len(df) > 0 else 0,
                             df['保留时间(RT)'].mean() if len(df) > 0 else 0,
                             df['保留时间(RT)'].median() if len(df) > 0 else 0,
                             df['强度'].min() if len(df) > 0 else 0,
                             df['强度'].max() if len(df) > 0 else 0,
                             df['强度'].mean() if len(df) > 0 else 0,
                             df['强度'].median() if len(df) > 0 else 0]
                })
                basic_stats.to_excel(writer, sheet_name='基本统计', index=False)

                # 强度百分位
                if len(df) > 0:
                    percentiles = [1, 5, 10, 25, 50, 75, 90, 95, 99]
                    intensity_percentiles = np.percentile(df['强度'], percentiles)

                    percentile_df = pd.DataFrame({
                        '百分位(%)': percentiles,
                        '强度值': intensity_percentiles
                    })
                    percentile_df.to_excel(writer, sheet_name='强度百分位', index=False)

                # 最强峰
                if len(df) > 0:
                    top_50 = df.nlargest(50, '强度').copy()
                    top_50.insert(0, '排名', range(1, len(top_50) + 1))
                    top_50.to_excel(writer, sheet_name='最强峰Top50', index=False)

            self.log(f"Excel统计已保存: {excel_file.name}")

        except Exception as e:
            self.log(f"保存Excel统计失败: {str(e)}", "WARNING")

    def process_file(self, input_file: str, output_dir: str, method: str = 'auto',
                     intensity_threshold: float = 1000) -> Dict[str, Any]:
        """处理单个文件"""
        try:
            input_path = Path(input_file)

            # 检查文件
            self.log(f"检查文件: {input_path.name}")
            self.update_progress(5, "检查文件中...")

            is_valid, message = self.check_mzml_file(input_path)
            if not is_valid:
                self.log(f"文件检查失败: {message}", "ERROR")
                return {
                    'input_file': input_file,
                    'status': 'error',
                    'message': message,
                    'peak_count': 0
                }

            self.log(f"文件检查通过: {message}")

            # 创建输出目录
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)

            # 提取峰数据
            self.log(f"开始提取峰数据...")
            peaks = None
            method_used = ""
            extract_message = ""

            # 根据方法选择提取器
            if method == 'auto':
                # 先尝试pyOpenMS
                self.update_progress(10, "尝试pyOpenMS...")
                peaks, msg = self.extract_peaks_with_pyopenms(input_file, intensity_threshold)
                if peaks is not None:
                    method_used = "pyOpenMS"
                    extract_message = msg
                else:
                    # 再尝试pymzML
                    self.log(f"pyOpenMS失败: {msg}", "WARNING")
                    self.update_progress(10, "尝试pymzML...")
                    peaks, msg = self.extract_peaks_with_pymzml(input_file, intensity_threshold)
                    if peaks is not None:
                        method_used = "pymzML"
                        extract_message = msg
                    else:
                        extract_message = f"所有方法都失败: {msg}"

            elif method == 'pyopenms':
                self.update_progress(10, "使用pyOpenMS...")
                peaks, msg = self.extract_peaks_with_pyopenms(input_file, intensity_threshold)
                method_used = "pyOpenMS"
                extract_message = msg

            elif method == 'pymzml':
                self.update_progress(10, "使用pymzML...")
                peaks, msg = self.extract_peaks_with_pymzml(input_file, intensity_threshold)
                method_used = "pymzML"
                extract_message = msg

            # 检查是否取消
            if self.cancel_requested:
                self.log("处理被用户取消", "WARNING")
                return {
                    'input_file': input_file,
                    'status': 'cancelled',
                    'message': '用户取消操作',
                    'peak_count': 0
                }

            # 保存结果
            if peaks is not None and len(peaks) > 0:
                base_name = input_path.stem
                result_files = self.save_results(peaks, output_path, base_name, method_used)

                return {
                    'input_file': input_file,
                    'output_dir': str(output_path),
                    'peak_count': len(peaks),
                    'method': method_used,
                    'status': 'success',
                    'message': extract_message,
                    'output_files': result_files
                }
            else:
                self.log(f"未提取到峰数据: {extract_message}", "WARNING")
                return {
                    'input_file': input_file,
                    'output_dir': str(output_path),
                    'peak_count': 0,
                    'method': method_used,
                    'status': 'failed',
                    'message': extract_message
                }

        except Exception as e:
            error_msg = f"处理文件时出错: {str(e)}"
            self.log(error_msg, "ERROR")
            traceback.print_exc()
            return {
                'input_file': input_file,
                'status': 'error',
                'message': error_msg,
                'peak_count': 0
            }

    def process_batch(self, input_dir: str, output_dir: str, file_pattern: str = "*.mzML",
                      method: str = 'auto', intensity_threshold: float = 1000) -> List[Dict[str, Any]]:
        """批量处理文件"""
        try:
            input_path = Path(input_dir)

            # 检查输入目录
            if not input_path.exists():
                error_msg = f"输入目录不存在: {input_dir}"
                self.log(error_msg, "ERROR")
                raise FileNotFoundError(error_msg)

            # 查找文件
            self.log(f"在目录中查找文件: {input_dir}")
            self.update_progress(5, "查找文件中...")

            mzml_files = list(input_path.glob(file_pattern))
            if not mzml_files:
                error_msg = f"在 {input_dir} 中未找到 {file_pattern} 文件"
                self.log(error_msg, "ERROR")
                raise FileNotFoundError(error_msg)

            self.log(f"找到 {len(mzml_files)} 个文件")

            # 创建输出目录
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)

            self.results_summary = []

            # 处理每个文件
            for i, mzml_file in enumerate(mzml_files):
                if self.cancel_requested:
                    self.log("批量处理被用户取消", "WARNING")
                    break

                file_num = i + 1
                total_files = len(mzml_files)

                self.log(f"处理文件 {file_num}/{total_files}: {mzml_file.name}")

                # 为每个文件创建单独的输出子目录
                file_output_dir = output_path / mzml_file.stem
                file_output_dir.mkdir(exist_ok=True)

                # 更新进度
                progress = int(file_num / total_files * 100)
                self.update_progress(progress, f"处理文件 {file_num}/{total_files}")

                # 处理文件
                result = self.process_file(str(mzml_file), str(file_output_dir), method, intensity_threshold)

                if result:
                    self.results_summary.append(result)

            # 生成批量报告
            if self.results_summary:
                self.generate_batch_report(output_path)

            return self.results_summary

        except Exception as e:
            error_msg = f"批量处理时出错: {str(e)}"
            self.log(error_msg, "ERROR")
            traceback.print_exc()
            return []

    def generate_batch_report(self, output_path: Path):
        """生成批量处理报告"""
        try:
            self.log("生成批量处理报告...")
            self.update_progress(95, "生成报告中...")

            successful = [r for r in self.results_summary if r.get('status') == 'success']
            failed = [r for r in self.results_summary if r.get('status') == 'failed']
            errors = [r for r in self.results_summary if r.get('status') == 'error']

            # 保存汇总报告
            summary_file = output_path / "批量处理汇总.csv"
            summary_data = []

            for result in self.results_summary:
                summary_data.append({
                    '文件名': Path(result['input_file']).name,
                    '状态': result.get('status', 'unknown'),
                    '峰数量': result.get('peak_count', 0),
                    '提取方法': result.get('method', ''),
                    '消息': result.get('message', '')[:100]
                })

            df_summary = pd.DataFrame(summary_data)
            df_summary.to_csv(summary_file, index=False, encoding='utf-8-sig')
            self.log(f"汇总报告已保存: {summary_file.name}")

            # 保存Excel汇总
            try:
                import openpyxl
                excel_summary = output_path / "批量处理汇总.xlsx"
                df_summary.to_excel(excel_summary, index=False, sheet_name='批量汇总')
                self.log(f"Excel汇总已保存: {excel_summary.name}")
            except:
                pass

            # 生成详细报告
            report_file = output_path / "批量处理报告.txt"
            with open(report_file, 'w', encoding='utf-8') as f:
                f.write("批量处理详细报告\n")
                f.write("=" * 70 + "\n\n")
                f.write(f"生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"输入目录: {output_path.parent}\n")
                f.write(f"输出目录: {output_path}\n\n")

                f.write(f"处理统计:\n")
                f.write(f"  总文件数: {len(self.results_summary)}")
                f.write(f"  成功处理: {len(successful)}")
                f.write(f"  处理失败: {len(failed)}")
                f.write(f"  处理错误: {len(errors)}\n\n")

                if successful:
                    total_peaks = sum(r.get('peak_count', 0) for r in successful)
                    f.write(f"成功处理文件列表 (共提取 {total_peaks:,} 个峰):\n")
                    f.write("-" * 50 + "\n")
                    for result in successful:
                        f.write(f"  {Path(result['input_file']).name}: {result.get('peak_count', 0)}个峰\n")
                    f.write("\n")

                if failed:
                    f.write(f"处理失败文件列表:\n")
                    f.write("-" * 50 + "\n")
                    for result in failed:
                        f.write(f"  {Path(result['input_file']).name}: {result.get('message', '未知原因')}\n")
                    f.write("\n")

            self.log(f"详细报告已保存: {report_file.name}")
            self.update_progress(100, "批量处理完成!")

        except Exception as e:
            self.log(f"生成批量报告失败: {str(e)}", "ERROR")

    def cancel(self):
        """取消处理"""
        self.cancel_requested = True
        self.log("正在取消处理...", "WARNING")


class MzMLGUI:
    """mzML文件处理GUI"""

    def __init__(self, root):
        """初始化GUI"""
        self.root = root
        self.root.title("mzML文件峰提取程序")
        self.root.geometry("900x700")

        # 设置图标
        try:
            self.root.iconbitmap(default='icon.ico')
        except:
            pass

        # 处理器
        self.processor = None
        self.processing_thread = None
        self.is_processing = False

        # 创建UI
        self.setup_ui()

        # 绑定关闭事件
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)

    def setup_ui(self):
        """设置UI界面"""
        # 主框架
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # 标题
        title_label = ttk.Label(main_frame, text="mzML文件峰提取程序",
                                font=("Arial", 16, "bold"))
        title_label.grid(row=0, column=0, columnspan=3, pady=10)

        # 版本信息
        version_label = ttk.Label(main_frame, text="版本 5.1 - 支持CSV和Excel输出",
                                  font=("Arial", 10))
        version_label.grid(row=1, column=0, columnspan=3, pady=(0, 10))

        # 分隔线
        ttk.Separator(main_frame, orient=tk.HORIZONTAL).grid(
            row=2, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=10
        )

        # 处理模式
        mode_frame = ttk.LabelFrame(main_frame, text="处理模式", padding="10")
        mode_frame.grid(row=3, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5)

        self.mode_var = tk.StringVar(value="single")

        ttk.Radiobutton(mode_frame, text="单个文件", variable=self.mode_var,
                        value="single", command=self.on_mode_change).grid(
            row=0, column=0, sticky=tk.W, padx=10)

        ttk.Radiobutton(mode_frame, text="批量处理", variable=self.mode_var,
                        value="batch", command=self.on_mode_change).grid(
            row=0, column=1, sticky=tk.W, padx=10)

        # 文件选择区域
        file_frame = ttk.LabelFrame(main_frame, text="文件选择", padding="10")
        file_frame.grid(row=4, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=10)

        # 单个文件选择
        self.single_file_frame = ttk.Frame(file_frame)
        self.single_file_frame.grid(row=0, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5)

        ttk.Label(self.single_file_frame, text="输入文件:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5))

        self.input_file_var = tk.StringVar()
        input_file_entry = ttk.Entry(self.single_file_frame, textvariable=self.input_file_var, width=60)
        input_file_entry.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=(0, 5))

        ttk.Button(self.single_file_frame, text="浏览...", command=self.browse_input_file).grid(
            row=0, column=2, padx=(0, 5))

        ttk.Button(self.single_file_frame, text="清空", command=self.clear_input_file).grid(
            row=0, column=3)

        # 批量处理选择
        self.batch_file_frame = ttk.Frame(file_frame)
        self.batch_file_frame.grid(row=1, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5)

        ttk.Label(self.batch_file_frame, text="输入目录:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5))

        self.input_dir_var = tk.StringVar()
        input_dir_entry = ttk.Entry(self.batch_file_frame, textvariable=self.input_dir_var, width=60)
        input_dir_entry.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=(0, 5))

        ttk.Button(self.batch_file_frame, text="浏览...", command=self.browse_input_dir).grid(
            row=0, column=2, padx=(0, 5))

        ttk.Label(self.batch_file_frame, text="文件模式:").grid(row=1, column=0, sticky=tk.W, padx=(0, 5), pady=(5, 0))

        self.file_pattern_var = tk.StringVar(value="*.mzML")
        file_pattern_combo = ttk.Combobox(self.batch_file_frame, textvariable=self.file_pattern_var,
                                          values=["*.mzML", "*.mzML.gz", "*.mzXML"], width=20)
        file_pattern_combo.grid(row=1, column=1, sticky=tk.W, padx=(0, 5), pady=(5, 0))

        # 初始隐藏批量框架
        self.batch_file_frame.grid_remove()

        # 输出目录
        ttk.Label(file_frame, text="输出目录:").grid(row=2, column=0, sticky=tk.W, padx=(0, 5), pady=(10, 0))

        self.output_dir_var = tk.StringVar()
        output_dir_entry = ttk.Entry(file_frame, textvariable=self.output_dir_var, width=60)
        output_dir_entry.grid(row=2, column=1, sticky=(tk.W, tk.E), padx=(0, 5), pady=(10, 0))

        ttk.Button(file_frame, text="浏览...", command=self.browse_output_dir).grid(
            row=2, column=2, pady=(10, 0))

        # 拖拽支持提示（移除拖拽相关代码，只保留浏览按钮）
        ttk.Label(file_frame, text="提示: 请使用浏览按钮选择文件或目录",
                  font=("Arial", 9, "italic")).grid(
            row=3, column=0, columnspan=3, sticky=tk.W, pady=(10, 0))

        # 处理参数
        param_frame = ttk.LabelFrame(main_frame, text="处理参数", padding="10")
        param_frame.grid(row=5, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=10)

        ttk.Label(param_frame, text="提取方法:").grid(row=0, column=0, sticky=tk.W, padx=(0, 10))

        self.method_var = tk.StringVar(value="auto")
        method_combo = ttk.Combobox(param_frame, textvariable=self.method_var,
                                    values=["auto", "pyopenms", "pymzml"],
                                    state="readonly", width=15)
        method_combo.grid(row=0, column=1, sticky=tk.W)

        ttk.Label(param_frame, text="强度阈值:").grid(row=0, column=2, sticky=tk.W, padx=(20, 10))

        self.threshold_var = tk.StringVar(value="1000")
        threshold_entry = ttk.Entry(param_frame, textvariable=self.threshold_var, width=15)
        threshold_entry.grid(row=0, column=3, sticky=tk.W)

        # 进度条
        progress_frame = ttk.LabelFrame(main_frame, text="处理进度", padding="10")
        progress_frame.grid(row=6, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=10)

        self.progress_var = tk.StringVar(value="准备就绪...")
        ttk.Label(progress_frame, textvariable=self.progress_var).grid(
            row=0, column=0, sticky=tk.W, pady=(0, 5))

        self.progress_bar = ttk.Progressbar(progress_frame, mode='determinate', length=800)
        self.progress_bar.grid(row=1, column=0, sticky=(tk.W, tk.E), pady=(0, 5))

        self.progress_percent = tk.StringVar(value="0%")
        ttk.Label(progress_frame, textvariable=self.progress_percent).grid(
            row=1, column=1, sticky=tk.E, padx=(5, 0), pady=(0, 5))

        # 日志区域
        log_frame = ttk.LabelFrame(main_frame, text="处理日志", padding="10")
        log_frame.grid(row=7, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), pady=10)

        # 日志文本框
        self.log_text = scrolledtext.ScrolledText(log_frame, width=100, height=15,
                                                  wrap=tk.WORD, font=("Consolas", 10))
        self.log_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # 配置标签颜色
        self.log_text.tag_config("INFO", foreground="black")
        self.log_text.tag_config("SUCCESS", foreground="green")
        self.log_text.tag_config("WARNING", foreground="orange")
        self.log_text.tag_config("ERROR", foreground="red")

        # 按钮区域
        button_frame = ttk.Frame(main_frame)
        button_frame.grid(row=8, column=0, columnspan=3, pady=10)

        self.process_button = ttk.Button(button_frame, text="开始处理",
                                         command=self.start_processing, width=15)
        self.process_button.pack(side=tk.LEFT, padx=5)

        self.cancel_button = ttk.Button(button_frame, text="取消处理",
                                        command=self.cancel_processing, width=15,
                                        state=tk.DISABLED)
        self.cancel_button.pack(side=tk.LEFT, padx=5)

        ttk.Button(button_frame, text="清空日志", command=self.clear_log).pack(side=tk.LEFT, padx=5)

        ttk.Button(button_frame, text="打开输出目录", command=self.open_output_dir).pack(side=tk.LEFT, padx=5)

        ttk.Button(button_frame, text="退出", command=self.on_closing).pack(side=tk.LEFT, padx=5)

        # 状态栏
        self.status_var = tk.StringVar(value="就绪")
        status_bar = ttk.Label(main_frame, textvariable=self.status_var, relief=tk.SUNKEN)
        status_bar.grid(row=9, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(10, 0))

        # 配置网格权重
        main_frame.columnconfigure(1, weight=1)
        log_frame.columnconfigure(0, weight=1)
        log_frame.rowconfigure(0, weight=1)

        # 初始化
        self.on_mode_change()

    def on_mode_change(self):
        """处理模式改变事件"""
        mode = self.mode_var.get()

        if mode == "single":
            self.single_file_frame.grid()
            self.batch_file_frame.grid_remove()
        else:
            self.single_file_frame.grid_remove()
            self.batch_file_frame.grid()

    def browse_input_file(self):
        """浏览输入文件"""
        file_types = [("mzML文件", "*.mzML *.mzML.gz"), ("所有文件", "*.*")]
        file_path = filedialog.askopenfilename(title="选择mzML文件", filetypes=file_types)

        if file_path:
            self.input_file_var.set(file_path)

            # 自动设置输出目录
            if not self.output_dir_var.get():
                input_path = Path(file_path)
                output_dir = input_path.parent / "output"
                self.output_dir_var.set(str(output_dir))

    def browse_input_dir(self):
        """浏览输入目录"""
        dir_path = filedialog.askdirectory(title="选择输入目录")

        if dir_path:
            self.input_dir_var.set(dir_path)

            # 检查目录中是否有mzML文件
            file_pattern = self.file_pattern_var.get()
            input_path = Path(dir_path)
            mzml_files = list(input_path.glob(file_pattern))

            if not mzml_files:
                messagebox.showwarning("警告", f"在目录中未找到 {file_pattern} 文件")
            else:
                messagebox.showinfo("信息", f"找到 {len(mzml_files)} 个文件")

            # 自动设置输出目录
            if not self.output_dir_var.get():
                output_dir = input_path.parent / "batch_output"
                self.output_dir_var.set(str(output_dir))

    def browse_output_dir(self):
        """浏览输出目录"""
        dir_path = filedialog.askdirectory(title="选择输出目录")

        if dir_path:
            self.output_dir_var.set(dir_path)

    def clear_input_file(self):
        """清空输入文件"""
        self.input_file_var.set("")

    def log_message(self, message: str, level: str = "INFO"):
        """添加日志消息"""
        self.log_text.insert(tk.END, message + "\n", level)
        self.log_text.see(tk.END)
        self.root.update()

    def update_progress(self, value: int, message: str = ""):
        """更新进度条"""
        self.progress_bar['value'] = value
        if message:
            self.progress_var.set(message)
        self.progress_percent.set(f"{value}%")
        self.root.update()

    def update_status(self, message: str):
        """更新状态栏"""
        self.status_var.set(message)
        self.root.update()

    def clear_log(self):
        """清空日志"""
        self.log_text.delete(1.0, tk.END)

    def validate_inputs(self) -> Tuple[bool, str]:
        """验证输入参数"""
        mode = self.mode_var.get()

        if mode == "single":
            input_path = self.input_file_var.get().strip()
            if not input_path:
                return False, "请选择输入文件"

            if not os.path.exists(input_path):
                return False, f"输入文件不存在: {input_path}"

        else:  # batch模式
            input_dir = self.input_dir_var.get().strip()
            if not input_dir:
                return False, "请选择输入目录"

            if not os.path.exists(input_dir):
                return False, f"输入目录不存在: {input_dir}"

            # 检查目录中是否有mzML文件
            file_pattern = self.file_pattern_var.get()
            input_path = Path(input_dir)
            mzml_files = list(input_path.glob(file_pattern))
            if not mzml_files:
                return False, f"在目录中未找到 {file_pattern} 文件"

        # 检查输出目录
        output_dir = self.output_dir_var.get().strip()
        if not output_dir:
            return False, "请选择输出目录"

        # 检查强度阈值
        try:
            threshold = float(self.threshold_var.get())
            if threshold <= 0:
                return False, "强度阈值必须大于0"
        except ValueError:
            return False, "强度阈值必须是数字"

        return True, "验证通过"

    def start_processing(self):
        """开始处理"""
        # 验证输入
        is_valid, message = self.validate_inputs()
        if not is_valid:
            messagebox.showerror("输入错误", message)
            return

        # 禁用开始按钮，启用取消按钮
        self.process_button.config(state=tk.DISABLED)
        self.cancel_button.config(state=tk.NORMAL)

        # 清空进度和日志
        self.progress_bar['value'] = 0
        self.progress_var.set("开始处理...")
        self.progress_percent.set("0%")
        self.clear_log()

        # 获取参数
        mode = self.mode_var.get()
        output_dir = self.output_dir_var.get()
        method = self.method_var.get()

        try:
            threshold = float(self.threshold_var.get())
        except ValueError:
            threshold = 1000

        # 创建处理器
        self.processor = MzMLProcessor(
            progress_callback=self.update_progress,
            log_callback=self.log_message
        )

        # 在新线程中处理
        if mode == "single":
            input_file = self.input_file_var.get()
            self.processing_thread = threading.Thread(
                target=self.process_single_file,
                args=(input_file, output_dir, method, threshold)
            )
        else:
            input_dir = self.input_dir_var.get()
            file_pattern = self.file_pattern_var.get()
            self.processing_thread = threading.Thread(
                target=self.process_batch,
                args=(input_dir, output_dir, file_pattern, method, threshold)
            )

        self.is_processing = True
        self.processing_thread.start()

        # 启动进度监视器
        self.root.after(100, self.check_processing_status)

    def process_single_file(self, input_file: str, output_dir: str, method: str, threshold: float):
        """处理单个文件（在线程中）"""
        try:
            result = self.processor.process_file(input_file, output_dir, method, threshold)

            # 在主线程中显示结果
            self.root.after(0, self.show_processing_result, result)

        except Exception as e:
            error_msg = f"处理过程中出错: {str(e)}"
            self.root.after(0, self.log_message, error_msg, "ERROR")
            self.root.after(0, self.processing_finished)

    def process_batch(self, input_dir: str, output_dir: str, file_pattern: str, method: str, threshold: float):
        """批量处理（在线程中）"""
        try:
            results = self.processor.process_batch(input_dir, output_dir, file_pattern, method, threshold)

            # 在主线程中显示结果
            self.root.after(0, self.show_batch_result, results)

        except Exception as e:
            error_msg = f"批量处理过程中出错: {str(e)}"
            self.root.after(0, self.log_message, error_msg, "ERROR")
            self.root.after(0, self.processing_finished)

    def show_processing_result(self, result: Dict[str, Any]):
        """显示单个文件处理结果"""
        if result['status'] == 'success':
            self.log_message(f"✓ 处理成功!", "SUCCESS")
            self.log_message(f"  提取到 {result['peak_count']:,} 个峰")
            self.log_message(f"  输出目录: {result['output_dir']}")

            if 'output_files' in result:
                self.log_message(f"  生成的文件:")
                for file_type, file_path in result['output_files'].items():
                    if file_path:
                        file_name = Path(file_path).name
                        self.log_message(f"    • {file_type}: {file_name}")

            # 显示成功对话框
            self.root.after(0, lambda: messagebox.showinfo(
                "处理成功",
                f"成功提取 {result['peak_count']:,} 个峰\n输出目录: {result['output_dir']}"
            ))

        elif result['status'] == 'cancelled':
            self.log_message("处理被用户取消", "WARNING")

        else:
            self.log_message(f"✗ 处理失败: {result.get('message', '未知错误')}", "ERROR")

        self.processing_finished()

    def show_batch_result(self, results: List[Dict[str, Any]]):
        """显示批量处理结果"""
        if not results:
            self.log_message("没有处理结果", "WARNING")
            self.processing_finished()
            return

        successful = [r for r in results if r.get('status') == 'success']
        failed = [r for r in results if r.get('status') == 'failed']
        errors = [r for r in results if r.get('status') == 'error']

        self.log_message(f"批量处理完成!", "SUCCESS")
        self.log_message(f"  总文件数: {len(results)}")
        self.log_message(f"  成功处理: {len(successful)}")
        self.log_message(f"  处理失败: {len(failed)}")
        self.log_message(f"  处理错误: {len(errors)}")

        if successful:
            total_peaks = sum(r.get('peak_count', 0) for r in successful)
            self.log_message(f"  总提取峰数: {total_peaks:,}")

        # 显示结果对话框
        result_text = f"批量处理完成!\n\n"
        result_text += f"总文件数: {len(results)}\n"
        result_text += f"成功处理: {len(successful)}\n"
        result_text += f"处理失败: {len(failed)}\n"

        if successful:
            total_peaks = sum(r.get('peak_count', 0) for r in successful)
            result_text += f"总提取峰数: {total_peaks:,}\n"

        result_text += f"\n详细结果请查看输出目录中的报告文件。"

        self.root.after(0, lambda: messagebox.showinfo("批量处理完成", result_text))

        self.processing_finished()

    def check_processing_status(self):
        """检查处理状态"""
        if self.is_processing and self.processing_thread.is_alive():
            # 继续检查
            self.root.after(100, self.check_processing_status)
        elif self.is_processing:
            # 线程已结束，但结果可能还未显示
            self.root.after(100, self.check_processing_status)

    def processing_finished(self):
        """处理完成"""
        self.is_processing = False
        self.processor = None

        # 启用开始按钮，禁用取消按钮
        self.process_button.config(state=tk.NORMAL)
        self.cancel_button.config(state=tk.DISABLED)

        # 更新状态
        self.update_status("处理完成")

    def cancel_processing(self):
        """取消处理"""
        if self.processor and self.is_processing:
            self.processor.cancel()
            self.log_message("正在取消处理...", "WARNING")
            self.update_status("正在取消...")

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

    def on_closing(self):
        """关闭窗口事件"""
        if self.is_processing:
            if messagebox.askyesno("确认", "正在处理中，确定要退出吗？"):
                if self.processor:
                    self.processor.cancel()
                self.root.destroy()
        else:
            self.root.destroy()


def check_dependencies():
    """检查依赖库"""
    print("=" * 70)
    print("检查依赖库...")
    print("=" * 70)

    dependencies = {
        'pandas': '数据处理库',
        'numpy': '数值计算库',
        'openpyxl': 'Excel XLSX格式支持',
        'pyopenms': '质谱处理库 (方法1)',
        'pymzml': '质谱处理库 (方法2)'
    }

    missing = []

    for lib, desc in dependencies.items():
        try:
            __import__(lib)
            print(f"✓ {lib}: {desc}")
        except ImportError:
            print(f"✗ {lib}: {desc} - 未安装")
            missing.append(lib)

    if missing:
        print(f"\n缺少以下库:")
        for lib in missing:
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
        print("  Windows/macOS: 通常已预装")
        return

    # 检查依赖
    if not check_dependencies():
        response = input("\n缺少依赖库，是否继续? (y/n): ")
        if response.lower() != 'y':
            return

    # 创建主窗口
    root = tk.Tk()

    # 设置窗口居中
    window_width = 900
    window_height = 700
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()
    x = (screen_width - window_width) // 2
    y = (screen_height - window_height) // 2
    root.geometry(f"{window_width}x{window_height}+{x}+{y}")

    # 创建GUI
    app = MzMLGUI(root)

    # 运行主循环
    root.mainloop()


if __name__ == "__main__":
    print("mzML文件峰提取程序 - GUI版本 (修复版)")
    print("版本: 5.1")
    print("功能: 支持进度显示、批量处理、CSV和Excel输出")
    print("注意: 已移除拖拽功能以避免兼容性问题")
    print("=" * 70)

    main()
