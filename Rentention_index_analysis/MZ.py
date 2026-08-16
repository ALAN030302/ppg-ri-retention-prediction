#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
mzML文件峰提取与化合物匹配及保留指数计算程序
功能：
1. 从SMILES/InChI计算化合物的m/z值
2. 从mzML文件提取峰数据
3. 从Excel表格加载化合物信息
4. 根据化合物m/z匹配峰数据
5. 计算化合物的PPG保留指数
6. 将匹配结果和保留指数保存到新的Excel表格
作者: 姚钱磊
版本: 9.0
日期: 2025
"""

import os
import sys
import threading
import traceback
import tempfile
import shutil
import json
import re
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Any, Union, Callable
import warnings
import math

# 抑制警告
warnings.filterwarnings('ignore')

try:
    import pandas as pd
    import numpy as np
    from scipy import stats
    from scipy.interpolate import interp1d
except ImportError:
    print("错误: 请先安装 pandas, numpy 和 scipy")
    print("安装命令: pip install pandas numpy scipy")
    sys.exit(1)

# 尝试导入化学信息学库
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    from rdkit.Chem import AllChem
    from rdkit.Chem import rdMolDescriptors

    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("警告: RDKit未安装，无法从SMILES/InChI计算m/z")
    print("安装命令: conda install -c conda-forge rdkit")

# 尝试导入GUI库
try:
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox, scrolledtext

    GUI_AVAILABLE = True
except ImportError:
    GUI_AVAILABLE = False
    print("警告: tkinter未安装，GUI不可用")


class MzCalculator:
    """m/z计算器 - 从SMILES/InChI计算化合物的m/z值"""

    def __init__(self):
        """初始化m/z计算器"""
        self.ion_modes = {
            'M+H': ('+', 1.007276),
            'M-H': ('-', -1.007276),
            'M+Na': ('+', 22.989218),
            'M+K': ('+', 38.963158),
            'M+NH4': ('+', 18.033823),
            'M+CH3COO': ('-', 59.013304),
            'M+H-H2O': ('+', -17.003289),
            'M+2H': ('+', 2.014552),
            'M+2Na': ('+', 45.978436),
            'M+2K': ('+', 77.926316),
            'M+H+Na': ('+', 23.996494),
            'M+H+K': ('+', 39.970434),
            'M+Na+K': ('+', 61.952376),
            'M+CH3OH+H': ('+', 33.033489),
            'M+ACN+H': ('+', 42.033823),
            'M+FA-H': ('-', 44.998203),  # 甲酸根
            'M+Cl': ('-', 34.969402)  # 氯离子
        }

        # 常见元素的质量
        self.element_masses = {
            'H': 1.007825,
            'C': 12.0,
            'N': 14.003074,
            'O': 15.994915,
            'P': 30.973762,
            'S': 31.972071,
            'F': 18.998403,
            'Cl': 34.968853,
            'Br': 78.918336,
            'I': 126.904473,
            'Na': 22.989770,
            'K': 38.963707,
            'Ca': 39.962591,
            'Mg': 23.985042,
            'Fe': 55.934939,
            'Zn': 63.929145,
            'Cu': 62.929598,
            'Mn': 54.938047
        }

    def calculate_monoisotopic_mass_from_formula(self, formula: str) -> float:
        """从分子式计算单同位素质量"""
        if not RDKIT_AVAILABLE:
            raise ImportError("RDKit未安装，无法计算分子式质量")

        try:
            # 使用RDKit从分子式计算质量
            mol = Chem.MolFromSmiles('')  # 创建空分子
            mol = Chem.AddHs(mol)  # 添加氢

            # 解析分子式
            pattern = r'([A-Z][a-z]*)(\d*)'
            elements = re.findall(pattern, formula)

            # 计算总质量
            total_mass = 0.0
            for element, count in elements:
                element = element.capitalize()
                if count == '':
                    count = 1
                else:
                    count = int(count)

                # 获取原子序数
                atomic_num = Chem.GetPeriodicTable().GetAtomicNumber(element)
                # 获取同位素质量
                isotopes = Chem.GetPeriodicTable().GetMostCommonIsotopeMass(atomic_num)
                if element in self.element_masses:
                    mass = self.element_masses[element] * count
                else:
                    # 使用RDKit的近似值
                    mass = isotopes * count

                total_mass += mass

            return total_mass

        except Exception as e:
            # 如果RDKit失败，尝试简单计算
            try:
                return self._simple_formula_mass(formula)
            except:
                raise ValueError(f"无法解析分子式 {formula}: {str(e)}")

    def _simple_formula_mass(self, formula: str) -> float:
        """简单分子式计算（不使用RDKit）"""
        pattern = r'([A-Z][a-z]*)(\d*)'
        elements = re.findall(pattern, formula)

        total_mass = 0.0
        for element, count in elements:
            element = element.capitalize()
            if count == '':
                count = 1
            else:
                count = int(count)

            if element in self.element_masses:
                total_mass += self.element_masses[element] * count
            else:
                # 使用近似值
                approx_masses = {
                    'Li': 6.941, 'Be': 9.012, 'B': 10.811, 'Al': 26.982,
                    'Si': 28.086, 'Ar': 39.948, 'Sc': 44.956, 'Ti': 47.867,
                    'V': 50.942, 'Cr': 51.996, 'Co': 58.933, 'Ni': 58.693,
                    'Ga': 69.723, 'Ge': 72.630, 'As': 74.922, 'Se': 78.971,
                    'Kr': 83.798, 'Rb': 85.468, 'Sr': 87.620, 'Y': 88.906,
                    'Zr': 91.224, 'Nb': 92.906, 'Mo': 95.950, 'Tc': 98.000,
                    'Ru': 101.070, 'Rh': 102.906, 'Pd': 106.420, 'Ag': 107.868,
                    'Cd': 112.414, 'In': 114.818, 'Sn': 118.710, 'Sb': 121.760,
                    'Te': 127.600, 'Xe': 131.293, 'Cs': 132.905, 'Ba': 137.327,
                    'La': 138.905, 'Ce': 140.116, 'Pr': 140.908, 'Nd': 144.242,
                    'Pm': 145.000, 'Sm': 150.360, 'Eu': 151.964, 'Gd': 157.250,
                    'Tb': 158.925, 'Dy': 162.500, 'Ho': 164.930, 'Er': 167.259,
                    'Tm': 168.934, 'Yb': 173.054, 'Lu': 174.967, 'Hf': 178.490,
                    'Ta': 180.948, 'W': 183.840, 'Re': 186.207, 'Os': 190.230,
                    'Ir': 192.217, 'Pt': 195.084, 'Au': 196.967, 'Hg': 200.592,
                    'Tl': 204.380, 'Pb': 207.200, 'Bi': 208.980, 'Th': 232.038,
                    'Pa': 231.036, 'U': 238.029
                }

                if element in approx_masses:
                    total_mass += approx_masses[element] * count
                else:
                    raise ValueError(f"未知元素: {element}")

        return total_mass

    def calculate_mz_from_smiles(self, smiles: str, ion_mode: str = 'M+H', charge: int = 1) -> float:
        """从SMILES计算m/z值

        参数:
            smiles: SMILES字符串
            ion_mode: 离子化模式
            charge: 电荷数

        返回:
            m/z值
        """
        if not RDKIT_AVAILABLE:
            raise ImportError("RDKit未安装，无法从SMILES计算m/z")

        try:
            # 从SMILES创建分子
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f"无法解析SMILES: {smiles}")

            # 计算单同位素质量
            mass = rdMolDescriptors.CalcExactMolWt(mol)

            # 应用离子化
            mz = self._apply_ionization(mass, ion_mode, charge)

            return mz

        except Exception as e:
            raise ValueError(f"从SMILES计算m/z失败: {str(e)}")

    def calculate_mz_from_inchi(self, inchi: str, ion_mode: str = 'M+H', charge: int = 1) -> float:
        """从InChI计算m/z值"""
        if not RDKIT_AVAILABLE:
            raise ImportError("RDKit未安装，无法从InChI计算m/z")

        try:
            # 从InChI创建分子
            mol = Chem.MolFromInchi(inchi)
            if mol is None:
                raise ValueError(f"无法解析InChI: {inchi}")

            # 计算单同位素质量
            mass = rdMolDescriptors.CalcExactMolWt(mol)

            # 应用离子化
            mz = self._apply_ionization(mass, ion_mode, charge)

            return mz

        except Exception as e:
            raise ValueError(f"从InChI计算m/z失败: {str(e)}")

    def calculate_mz_from_formula(self, formula: str, ion_mode: str = 'M+H', charge: int = 1) -> float:
        """从分子式计算m/z值"""
        try:
            # 计算单同位素质量
            mass = self.calculate_monoisotopic_mass_from_formula(formula)

            # 应用离子化
            mz = self._apply_ionization(mass, ion_mode, charge)

            return mz

        except Exception as e:
            raise ValueError(f"从分子式计算m/z失败: {str(e)}")

    def _apply_ionization(self, mass: float, ion_mode: str, charge: int) -> float:
        """应用离子化模式"""
        if ion_mode not in self.ion_modes:
            raise ValueError(f"不支持的离子模式: {ion_mode}")

        polarity, offset = self.ion_modes[ion_mode]

        # 计算m/z
        if charge == 0:
            raise ValueError("电荷数不能为0")

        mz = (mass + offset) / abs(charge)

        return mz

    def detect_ion_mode_from_mz(self, observed_mz: float, calculated_mass: float,
                                tolerance_ppm: float = 10) -> List[Tuple[str, float]]:
        """从观察到的m/z和计算的质量推断可能的离子化模式

        返回:
            可能的离子化模式列表，包含模式名称和匹配误差(ppm)
        """
        possible_modes = []

        for mode, (polarity, offset) in self.ion_modes.items():
            for charge in [1, 2, 3]:  # 考虑1, 2, 3个电荷
                expected_mz = (calculated_mass + offset) / charge
                error_ppm = abs((observed_mz - expected_mz) / expected_mz * 1e6)

                if error_ppm <= tolerance_ppm:
                    possible_modes.append((f"{mode} ({charge}+)" if charge > 1 else mode, error_ppm))

        # 按误差排序
        possible_modes.sort(key=lambda x: x[1])

        return possible_modes

    def batch_calculate_mz(self, compounds_df: pd.DataFrame,
                           smiles_col: str = None,
                           inchi_col: str = None,
                           formula_col: str = None,
                           ion_mode: str = 'M+H') -> pd.DataFrame:
        """批量计算m/z值"""
        result_df = compounds_df.copy()

        if 'mz' not in result_df.columns:
            result_df['mz'] = np.nan

        if 'mz_source' not in result_df.columns:
            result_df['mz_source'] = ''

        for idx, row in result_df.iterrows():
            calculated = False

            # 尝试从SMILES计算
            if smiles_col and smiles_col in row and pd.notna(row[smiles_col]):
                try:
                    mz = self.calculate_mz_from_smiles(str(row[smiles_col]), ion_mode)
                    result_df.at[idx, 'mz'] = mz
                    result_df.at[idx, 'mz_source'] = f'SMILES ({ion_mode})'
                    calculated = True
                except Exception as e:
                    pass

            # 尝试从InChI计算
            if not calculated and inchi_col and inchi_col in row and pd.notna(row[inchi_col]):
                try:
                    mz = self.calculate_mz_from_inchi(str(row[inchi_col]), ion_mode)
                    result_df.at[idx, 'mz'] = mz
                    result_df.at[idx, 'mz_source'] = f'InChI ({ion_mode})'
                    calculated = True
                except Exception as e:
                    pass

            # 尝试从分子式计算
            if not calculated and formula_col and formula_col in row and pd.notna(row[formula_col]):
                try:
                    mz = self.calculate_mz_from_formula(str(row[formula_col]), ion_mode)
                    result_df.at[idx, 'mz'] = mz
                    result_df.at[idx, 'mz_source'] = f'Formula ({ion_mode})'
                    calculated = True
                except Exception as e:
                    pass

        return result_df


class PPGIndexCalculator:
    """PPG保留指数计算器"""

    def __init__(self):
        """初始化保留指数计算器"""
        self.ppg_standards = None  # PPG标准品数据
        self.regression_params = None  # 线性回归参数
        self.calibration_curve = None  # 校准曲线

    def load_ppg_standards(self, excel_file: str) -> Tuple[bool, str, pd.DataFrame]:
        """从Excel文件加载PPG标准品数据"""
        try:
            file_path = Path(excel_file)

            if not file_path.exists():
                return False, "PPG标准品文件不存在", None

            # 尝试读取Excel文件
            try:
                # 尝试所有可能的sheet
                xls = pd.ExcelFile(file_path)
                sheet_names = xls.sheet_names

                # 尝试第一个sheet
                df = pd.read_excel(file_path, sheet_name=sheet_names[0])

                # 查找必要的列
                rt_col = self._find_column(df, ['RT', 'rt', '保留时间', 'RetentionTime', 'retention_time'])
                n_col = self._find_column(df, ['n', 'N', '聚合度', '聚合度n', 'degree', 'Degree'])
                name_col = self._find_column(df, ['name', 'Name', '化合物', 'Compound', 'PPG'])

                if rt_col is None or n_col is None:
                    return False, "未找到保留时间或聚合度列，请确保Excel文件包含RT和n信息", None

                # 重命名列
                rename_dict = {}
                if rt_col:
                    rename_dict[rt_col] = 'RT'
                if n_col:
                    rename_dict[n_col] = 'n'
                if name_col:
                    rename_dict[name_col] = 'name'

                if rename_dict:
                    df = df.rename(columns=rename_dict)

                # 确保数据类型正确
                if 'RT' in df.columns:
                    df['RT'] = pd.to_numeric(df['RT'], errors='coerce')
                if 'n' in df.columns:
                    df['n'] = pd.to_numeric(df['n'], errors='coerce')

                # 排序
                df = df.sort_values('n').reset_index(drop=True)

                self.ppg_standards = df

                return True, f"成功加载 {len(df)} 个PPG标准品", df

            except Exception as e:
                return False, f"读取PPG标准品文件失败: {str(e)}", None

        except Exception as e:
            return False, f"加载PPG标准品信息失败: {str(e)}", None

    def _find_column(self, df: pd.DataFrame, possible_names: List[str]) -> Optional[str]:
        """在DataFrame中查找可能的列名"""
        for name in possible_names:
            if name in df.columns:
                return name
        return None

    def build_calibration_curve(self, method: str = 'linear') -> Tuple[bool, str]:
        """建立PPG校准曲线

        参数:
            method: 校准方法，可选 'linear' (线性回归) 或 'interpolation' (插值)
        """
        if self.ppg_standards is None or len(self.ppg_standards) < 2:
            return False, "PPG标准品数据不足，至少需要2个标准品"

        try:
            df = self.ppg_standards.copy()
            df = df.dropna(subset=['RT', 'n'])

            if len(df) < 2:
                return False, "有效的PPG标准品数据不足"

            if method == 'linear':
                # 线性回归方法
                x = df['n'].values
                y = df['RT'].values

                # 线性回归
                slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

                # 计算预测值
                y_pred = intercept + slope * x
                residuals = y - y_pred

                # 计算R²和调整R²
                ss_res = np.sum(residuals ** 2)
                ss_tot = np.sum((y - np.mean(y)) ** 2)
                r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0
                n = len(x)
                p = 1  # 一个自变量
                adj_r_squared = 1 - (1 - r_squared) * (n - 1) / (n - p - 1) if n > p + 1 else r_squared

                self.regression_params = {
                    'method': 'linear',
                    'slope': slope,
                    'intercept': intercept,
                    'r_squared': r_squared,
                    'adj_r_squared': adj_r_squared,
                    'std_err': std_err,
                    'p_value': p_value,
                    'n_points': n
                }

                # 创建校准曲线函数
                def calibration_func(n_value):
                    return intercept + slope * n_value

                self.calibration_curve = calibration_func

                return True, f"线性回归校准曲线建立成功: R²={r_squared:.4f}, 斜率={slope:.4f}, 截距={intercept:.4f}"

            elif method == 'interpolation':
                # 插值方法
                x = df['n'].values
                y = df['RT'].values

                # 创建插值函数
                self.calibration_curve = interp1d(x, y, kind='linear', fill_value='extrapolate')
                self.regression_params = {'method': 'interpolation', 'n_points': len(x)}

                return True, f"插值校准曲线建立成功: 使用{len(x)}个标准品点"

            else:
                return False, f"不支持的校准方法: {method}"

        except Exception as e:
            return False, f"建立校准曲线失败: {str(e)}"

    def calculate_ppg_index(self, rt_values: Union[float, List[float], pd.Series],
                            method: str = 'interpolation') -> Union[float, List[float], pd.Series]:
        """计算PPG保留指数

        参数:
            rt_values: 保留时间值，可以是单个值、列表或pandas Series
            method: 计算方法，可选 'linear' (线性回归反算) 或 'interpolation' (线性插值)

        返回:
            对应的PPG保留指数（n*100）
        """
        if self.ppg_standards is None or len(self.ppg_standards) < 2:
            raise ValueError("PPG标准品数据未加载或不足")

        # 处理输入
        if isinstance(rt_values, (int, float)):
            is_scalar = True
            rt_list = [float(rt_values)]
        elif isinstance(rt_values, pd.Series):
            is_scalar = False
            rt_list = rt_values.tolist()
        else:
            is_scalar = False
            rt_list = list(rt_values)

        results = []

        for rt in rt_list:
            if pd.isna(rt):
                results.append(None)
                continue

            try:
                if method == 'linear' and self.regression_params and self.regression_params['method'] == 'linear':
                    # 使用线性回归反算
                    slope = self.regression_params['slope']
                    intercept = self.regression_params['intercept']

                    if abs(slope) < 1e-10:  # 避免除零
                        n = 0
                    else:
                        n = (rt - intercept) / slope

                else:
                    # 使用线性插值
                    df = self.ppg_standards.sort_values('RT')

                    # 边界处理
                    if rt < df['RT'].min():
                        # 小于最小RT，使用外推
                        rt_min = df['RT'].iloc[0]
                        rt_next = df['RT'].iloc[1]
                        n_min = df['n'].iloc[0]
                        n_next = df['n'].iloc[1]

                        if rt_next - rt_min != 0:
                            n = n_min + (rt - rt_min) / (rt_next - rt_min) * (n_next - n_min)
                        else:
                            n = n_min

                    elif rt > df['RT'].max():
                        # 大于最大RT，使用外推
                        rt_max = df['RT'].iloc[-1]
                        rt_prev = df['RT'].iloc[-2]
                        n_max = df['n'].iloc[-1]
                        n_prev = df['n'].iloc[-2]

                        if rt_max - rt_prev != 0:
                            n = n_max + (rt - rt_max) / (rt_max - rt_prev) * (n_max - n_prev)
                        else:
                            n = n_max

                    else:
                        # 找到所在的区间
                        for i in range(len(df) - 1):
                            rt_i = df['RT'].iloc[i]
                            rt_ip1 = df['RT'].iloc[i + 1]

                            if rt_i <= rt <= rt_ip1:
                                n_i = df['n'].iloc[i]
                                n_ip1 = df['n'].iloc[i + 1]

                                if rt_ip1 - rt_i != 0:
                                    # 线性插值
                                    n = n_i + (rt - rt_i) / (rt_ip1 - rt_i) * (n_ip1 - n_i)
                                else:
                                    n = n_i
                                break
                        else:
                            n = None

                if n is not None:
                    # 乘以100得到标准化的保留指数（如方案中所述）
                    n_final = n * 100
                else:
                    n_final = None

                results.append(n_final)

            except Exception as e:
                print(f"计算保留指数时出错 (RT={rt}): {str(e)}")
                results.append(None)

        # 返回格式与输入一致
        if is_scalar:
            return results[0] if results else None
        elif isinstance(rt_values, pd.Series):
            return pd.Series(results, index=rt_values.index, name='PPG_Index')
        else:
            return results

    def calculate_ppg_index_batch(self, rt_df: pd.DataFrame, rt_column: str = '保留时间(RT)') -> pd.DataFrame:
        """批量计算PPG保留指数

        参数:
            rt_df: 包含保留时间的DataFrame
            rt_column: 保留时间列名

        返回:
            添加了PPG保留指数列的DataFrame
        """
        if rt_column not in rt_df.columns:
            raise ValueError(f"DataFrame中未找到列: {rt_column}")

        result_df = rt_df.copy()
        result_df['PPG_Index'] = self.calculate_ppg_index(rt_df[rt_column])

        return result_df

    def save_calibration_report(self, output_dir: str, base_name: str = "PPG校准报告") -> Tuple[bool, str, List[str]]:
        """保存PPG校准报告"""
        try:
            if self.ppg_standards is None:
                return False, "没有PPG标准品数据可保存", []

            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)

            saved_files = []
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

            # 保存校准数据
            excel_file = output_path / f"{base_name}_{timestamp}.xlsx"

            with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
                # PPG标准品数据
                self.ppg_standards.to_excel(writer, sheet_name='PPG标准品', index=False)

                # 校准曲线参数
                if self.regression_params:
                    params_df = pd.DataFrame([self.regression_params])
                    params_df.to_excel(writer, sheet_name='校准参数', index=False)

                # 统计信息
                if self.ppg_standards is not None:
                    stats_data = {
                        '统计项目': ['标准品数量', '聚合度范围', '保留时间范围', '平均保留时间', '保留时间标准差'],
                        '数值': [
                            len(self.ppg_standards),
                            f"{self.ppg_standards['n'].min():.1f} - {self.ppg_standards['n'].max():.1f}",
                            f"{self.ppg_standards['RT'].min():.3f} - {self.ppg_standards['RT'].max():.3f}",
                            f"{self.ppg_standards['RT'].mean():.3f}",
                            f"{self.ppg_standards['RT'].std():.3f}"
                        ]
                    }

                    if self.regression_params and self.regression_params['method'] == 'linear':
                        stats_data['统计项目'].extend(['R²', '调整R²', '斜率', '截距', '标准误差'])
                        stats_data['数值'].extend([
                            f"{self.regression_params['r_squared']:.4f}",
                            f"{self.regression_params['adj_r_squared']:.4f}",
                            f"{self.regression_params['slope']:.4f}",
                            f"{self.regression_params['intercept']:.4f}",
                            f"{self.regression_params['std_err']:.4f}"
                        ])

                    stats_df = pd.DataFrame(stats_data)
                    stats_df.to_excel(writer, sheet_name='统计信息', index=False)

            saved_files.append(str(excel_file))

            return True, f"校准报告已保存到 {excel_file.name}", saved_files

        except Exception as e:
            return False, f"保存校准报告失败: {str(e)}", []


class CompoundMatcher:
    """化合物匹配器 - 从Excel加载化合物信息并进行匹配"""

    def __init__(self):
        """初始化匹配器"""
        self.compounds_df = None
        self.match_results = []
        self.mz_calculator = MzCalculator()  # m/z计算器
        self.match_settings = {
            'ppm_tolerance': 10,  # ppm误差
            'rt_window': 30,  # 保留时间窗口（秒）
            'intensity_threshold': 1000,  # 强度阈值
            'ion_mode': 'M+H',  # 离子化模式
            'calculate_mz': True  # 是否自动计算m/z
        }

    def load_compounds_from_excel(self, excel_file: str,
                                  calculate_mz: bool = True,
                                  ion_mode: str = 'M+H') -> Tuple[bool, str, pd.DataFrame]:
        """从Excel文件加载化合物信息，可自动计算m/z值"""
        try:
            file_path = Path(excel_file)

            if not file_path.exists():
                return False, "文件不存在", None

            # 尝试读取Excel文件
            try:
                # 尝试所有可能的sheet
                xls = pd.ExcelFile(file_path)
                sheet_names = xls.sheet_names

                # 尝试第一个sheet
                df = pd.read_excel(file_path, sheet_name=sheet_names[0])

                # 查找必要的列
                mz_col = self._find_column(df, ['mz', 'M/Z', 'm/z', 'mass', 'Mass', 'MZ', '质荷比'])
                name_col = self._find_column(df, ['name', 'Name', 'compound', 'Compound', 'ID', '名称', '化合物名称'])
                formula_col = self._find_column(df, ['formula', 'Formula', '分子式'])
                rt_col = self._find_column(df, ['rt', 'RT', '保留时间', 'retention_time'])
                smiles_col = self._find_column(df, ['smiles', 'SMILES', 'Smiles'])
                inchi_col = self._find_column(df, ['inchi', 'InChI', 'InChIKey'])

                # 重命名列
                rename_dict = {}
                if mz_col:
                    rename_dict[mz_col] = 'mz'
                if name_col:
                    rename_dict[name_col] = 'name'
                if formula_col:
                    rename_dict[formula_col] = 'formula'
                if rt_col:
                    rename_dict[rt_col] = 'rt_reference'
                if smiles_col:
                    rename_dict[smiles_col] = 'smiles'
                if inchi_col:
                    rename_dict[inchi_col] = 'inchi'

                if rename_dict:
                    df = df.rename(columns=rename_dict)

                # 确保必要的列存在
                if 'name' not in df.columns:
                    df['name'] = [f"化合物_{i + 1}" for i in range(len(df))]

                # 如果需要计算m/z且m/z列不存在或需要计算
                if calculate_mz and ('mz' not in df.columns or df['mz'].isna().all()):
                    # 批量计算m/z
                    df = self.mz_calculator.batch_calculate_mz(
                        df,
                        smiles_col='smiles' if 'smiles' in df.columns else None,
                        inchi_col='inchi' if 'inchi' in df.columns else None,
                        formula_col='formula' if 'formula' in df.columns else None,
                        ion_mode=ion_mode
                    )

                    # 统计计算情况
                    calculated_count = df['mz_source'].notna().sum()
                    if calculated_count > 0:
                        print(f"自动计算了 {calculated_count} 个化合物的m/z值")

                # 确保m/z列是数值
                if 'mz' in df.columns:
                    df['mz'] = pd.to_numeric(df['mz'], errors='coerce')
                    # 删除m/z为NaN的行（如果允许）
                    df_original_len = len(df)
                    df = df.dropna(subset=['mz'])
                    dropped_count = df_original_len - len(df)
                    if dropped_count > 0:
                        print(f"删除了 {dropped_count} 个没有m/z值的化合物")

                self.compounds_df = df

                return True, f"成功加载 {len(df)} 个化合物", df

            except Exception as e:
                return False, f"读取Excel文件失败: {str(e)}", None

        except Exception as e:
            return False, f"加载化合物信息失败: {str(e)}", None

    def _find_column(self, df: pd.DataFrame, possible_names: List[str]) -> Optional[str]:
        """在DataFrame中查找可能的列名"""
        for name in possible_names:
            if name in df.columns:
                return name
        return None

    def set_match_settings(self, ppm_tolerance: float = 10, rt_window: float = 30,
                           intensity_threshold: float = 1000, ion_mode: str = 'M+H',
                           calculate_mz: bool = True):
        """设置匹配参数"""
        self.match_settings = {
            'ppm_tolerance': ppm_tolerance,
            'rt_window': rt_window,
            'intensity_threshold': intensity_threshold,
            'ion_mode': ion_mode,
            'calculate_mz': calculate_mz
        }

    def calculate_mz_for_compound(self, identifier: str, identifier_type: str = 'smiles',
                                  ion_mode: str = 'M+H') -> Tuple[bool, str, float]:
        """为单个化合物计算m/z值

        参数:
            identifier: SMILES、InChI或分子式字符串
            identifier_type: 'smiles', 'inchi', 或 'formula'
            ion_mode: 离子化模式

        返回:
            (成功与否, 消息, m/z值)
        """
        try:
            if identifier_type.lower() == 'smiles':
                mz = self.mz_calculator.calculate_mz_from_smiles(identifier, ion_mode)
                return True, f"从SMILES成功计算m/z", mz
            elif identifier_type.lower() == 'inchi':
                mz = self.mz_calculator.calculate_mz_from_inchi(identifier, ion_mode)
                return True, f"从InChI成功计算m/z", mz
            elif identifier_type.lower() in ['formula', 'mf', 'molformula']:
                mz = self.mz_calculator.calculate_mz_from_formula(identifier, ion_mode)
                return True, f"从分子式成功计算m/z", mz
            else:
                return False, f"不支持的标识符类型: {identifier_type}", 0.0
        except Exception as e:
            return False, f"计算m/z失败: {str(e)}", 0.0

    def match_compounds(self, peaks_data: Union[List[Dict], pd.DataFrame, str]) -> Tuple[bool, str]:
        """匹配化合物和峰数据"""
        try:
            # 检查是否有化合物数据
            if self.compounds_df is None or len(self.compounds_df) == 0:
                return False, "未加载化合物信息"

            # 转换峰数据为DataFrame
            if isinstance(peaks_data, str):
                # 文件路径
                peaks_df = pd.read_csv(peaks_data)
            elif isinstance(peaks_data, pd.DataFrame):
                peaks_df = peaks_data
            elif isinstance(peaks_data, list):
                peaks_df = pd.DataFrame(peaks_data)
            else:
                return False, "不支持的峰数据类型"

            # 检查峰数据是否包含必要列
            required_columns = ['质荷比(m/z)', '保留时间(RT)', '强度']
            missing_cols = [col for col in required_columns if col not in peaks_df.columns]

            if missing_cols:
                # 尝试查找其他可能的列名
                column_mapping = {
                    '质荷比(m/z)': ['mz', 'M/Z', 'm/z'],
                    '保留时间(RT)': ['rt', 'RT', '保留时间'],
                    '强度': ['intensity', 'Intensity', '强度']
                }

                for req_col, possible_names in column_mapping.items():
                    if req_col not in peaks_df.columns:
                        for name in possible_names:
                            if name in peaks_df.columns:
                                peaks_df = peaks_df.rename(columns={name: req_col})
                                break

            # 再次检查
            missing_cols = [col for col in required_columns if col not in peaks_df.columns]
            if missing_cols:
                return False, f"峰数据缺少必要列: {', '.join(missing_cols)}"

            # 确保数据类型正确
            peaks_df['质荷比(m/z)'] = pd.to_numeric(peaks_df['质荷比(m/z)'], errors='coerce')
            peaks_df['保留时间(RT)'] = pd.to_numeric(peaks_df['保留时间(RT)'], errors='coerce')
            peaks_df['强度'] = pd.to_numeric(peaks_df['强度'], errors='coerce')

            # 过滤低于阈值的峰
            peaks_df = peaks_df[peaks_df['强度'] >= self.match_settings['intensity_threshold']]

            self.match_results = []

            # 对每个化合物进行匹配
            for _, compound in self.compounds_df.iterrows():
                compound_mz = compound.get('mz')
                compound_name = compound.get('name', '未知化合物')
                formula = compound.get('formula', '')
                rt_reference = compound.get('rt_reference', None)
                smiles = compound.get('smiles', '')
                inchi = compound.get('inchi', '')
                mz_source = compound.get('mz_source', '用户提供')

                if pd.isna(compound_mz):
                    continue

                # 计算ppm误差
                ppm_tolerance = self.match_settings['ppm_tolerance']
                mz_tolerance = compound_mz * ppm_tolerance / 1e6

                # 查找匹配的峰
                matching_peaks = peaks_df[
                    (peaks_df['质荷比(m/z)'] >= compound_mz - mz_tolerance) &
                    (peaks_df['质荷比(m/z)'] <= compound_mz + mz_tolerance)
                    ].copy()

                # 如果有参考保留时间，进一步筛选
                if rt_reference is not None and not pd.isna(rt_reference):
                    rt_window = self.match_settings['rt_window']
                    matching_peaks = matching_peaks[
                        (matching_peaks['保留时间(RT)'] >= rt_reference - rt_window) &
                        (matching_peaks['保留时间(RT)'] <= rt_reference + rt_window)
                        ]

                if len(matching_peaks) > 0:
                    # 选择最接近的峰
                    matching_peaks['mz_difference'] = abs(matching_peaks['质荷比(m/z)'] - compound_mz)
                    matching_peaks['mz_difference_ppm'] = matching_peaks['mz_difference'] / compound_mz * 1e6

                    # 按m/z误差排序
                    matching_peaks = matching_peaks.sort_values('mz_difference')

                    for i, (_, peak) in enumerate(matching_peaks.iterrows()):
                        if i >= 3:  # 只取前3个最接近的匹配
                            break

                        result = {
                            '化合物名称': compound_name,
                            '分子式': formula,
                            '理论m/z': compound_mz,
                            'm/z来源': mz_source,
                            '实测m/z': peak['质荷比(m/z)'],
                            'm/z误差(Da)': peak['mz_difference'],
                            'm/z误差(ppm)': peak['mz_difference_ppm'],
                            '保留时间(RT)': peak['保留时间(RT)'],
                            '强度': peak['强度'],
                            '匹配状态': '匹配成功',
                            '匹配排名': i + 1
                        }

                        # 添加结构信息
                        if smiles:
                            result['SMILES'] = smiles
                        if inchi:
                            result['InChI'] = inchi

                        # 添加参考保留时间信息
                        if rt_reference is not None and not pd.isna(rt_reference):
                            result['参考RT'] = rt_reference
                            result['RT误差'] = abs(peak['保留时间(RT)'] - rt_reference)

                        self.match_results.append(result)
                else:
                    # 未找到匹配
                    result = {
                        '化合物名称': compound_name,
                        '分子式': formula,
                        '理论m/z': compound_mz,
                        'm/z来源': mz_source,
                        '实测m/z': None,
                        'm/z误差(Da)': None,
                        'm/z误差(ppm)': None,
                        '保留时间(RT)': None,
                        '强度': None,
                        '匹配状态': '未匹配',
                        '匹配排名': None
                    }

                    # 添加结构信息
                    if smiles:
                        result['SMILES'] = smiles
                    if inchi:
                        result['InChI'] = inchi

                    if rt_reference is not None and not pd.isna(rt_reference):
                        result['参考RT'] = rt_reference

                    self.match_results.append(result)

            if len(self.match_results) == 0:
                return False, "未找到任何匹配"

            # 统计匹配情况
            successful_matches = [r for r in self.match_results if r['匹配状态'] == '匹配成功']

            return True, f"匹配完成: 共处理 {len(self.compounds_df)} 个化合物，找到 {len(successful_matches)} 个匹配"

        except Exception as e:
            return False, f"匹配过程中出错: {str(e)}"

    def add_ppg_index(self, ppg_calculator: PPGIndexCalculator) -> Tuple[bool, str]:
        """为匹配结果添加PPG保留指数"""
        try:
            if not self.match_results:
                return False, "没有匹配结果可添加PPG指数"

            if ppg_calculator.calibration_curve is None:
                return False, "PPG校准曲线未建立"

            # 将匹配结果转换为DataFrame
            results_df = pd.DataFrame(self.match_results)

            # 计算PPG保留指数
            if '保留时间(RT)' in results_df.columns:
                results_df['PPG保留指数'] = ppg_calculator.calculate_ppg_index(results_df['保留时间(RT)'])

            # 更新匹配结果
            self.match_results = results_df.to_dict('records')

            # 计算成功匹配的PPG指数统计
            successful_matches = [r for r in self.match_results if
                                  r['匹配状态'] == '匹配成功' and r.get('PPG保留指数') is not None]
            if successful_matches:
                ppg_values = [r['PPG保留指数'] for r in successful_matches]
                avg_ppg = sum(ppg_values) / len(ppg_values)
                min_ppg = min(ppg_values)
                max_ppg = max(ppg_values)

                return True, f"PPG指数添加成功: {len(successful_matches)}个匹配化合物，PPG指数范围: {min_ppg:.1f} - {max_ppg:.1f}，平均: {avg_ppg:.1f}"
            else:
                return True, "PPG指数添加完成，但无成功匹配的化合物"

        except Exception as e:
            return False, f"添加PPG指数失败: {str(e)}"

    def save_match_results(self, output_dir: str, base_name: str = "化合物匹配结果",
                           include_ppg: bool = True) -> Tuple[bool, str, List[str]]:
        """保存匹配结果到Excel文件"""
        try:
            if not self.match_results:
                return False, "没有匹配结果可保存", []

            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)

            saved_files = []

            # 创建DataFrame
            results_df = pd.DataFrame(self.match_results)

            # 保存为Excel
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

            # 主结果文件
            excel_file = output_path / f"{base_name}_{timestamp}.xlsx"

            with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
                # 所有结果
                results_df.to_excel(writer, sheet_name='所有匹配结果', index=False)

                # 成功匹配的结果
                successful_matches = results_df[results_df['匹配状态'] == '匹配成功'].copy()
                if len(successful_matches) > 0:
                    successful_matches = successful_matches.sort_values(['化合物名称', '匹配排名'])
                    successful_matches.to_excel(writer, sheet_name='成功匹配', index=False)

                # 按化合物汇总
                compound_summary = []
                for compound_name in results_df['化合物名称'].unique():
                    compound_matches = results_df[results_df['化合物名称'] == compound_name]
                    successful = compound_matches[compound_matches['匹配状态'] == '匹配成功']

                    summary = {
                        '化合物名称': compound_name,
                        '总匹配次数': len(successful),
                        '最佳m/z误差(ppm)': successful['m/z误差(ppm)'].min() if len(successful) > 0 else None,
                        '最高强度': successful['强度'].max() if len(successful) > 0 else None,
                        '平均RT': successful['保留时间(RT)'].mean() if len(successful) > 0 else None,
                        '匹配状态': '已匹配' if len(successful) > 0 else '未匹配'
                    }

                    # 添加m/z来源信息
                    if 'm/z来源' in compound_matches.columns:
                        sources = compound_matches['m/z来源'].unique()
                        summary['m/z来源'] = ', '.join([s for s in sources if pd.notna(s)])

                    # 添加PPG指数信息
                    if include_ppg and 'PPG保留指数' in successful.columns:
                        ppg_values = successful['PPG保留指数'].dropna()
                        if len(ppg_values) > 0:
                            summary['平均PPG指数'] = ppg_values.mean()
                            summary['PPG指数范围'] = f"{ppg_values.min():.1f} - {ppg_values.max():.1f}"

                    compound_summary.append(summary)

                summary_df = pd.DataFrame(compound_summary)
                summary_df.to_excel(writer, sheet_name='化合物汇总', index=False)

                # 统计信息
                stats = {
                    '统计项目': ['化合物总数', '成功匹配数', '未匹配数', '总匹配峰数',
                                 '平均m/z误差(ppm)', '最小m/z误差(ppm)', '最大m/z误差(ppm)'],
                    '数值': [
                        len(results_df['化合物名称'].unique()),
                        len(results_df[results_df['匹配状态'] == '匹配成功']['化合物名称'].unique()),
                        len(results_df[results_df['匹配状态'] == '未匹配']['化合物名称'].unique()),
                        len(results_df[results_df['匹配状态'] == '匹配成功']),
                        results_df['m/z误差(ppm)'].mean() if 'm/z误差(ppm)' in results_df.columns else 0,
                        results_df['m/z误差(ppm)'].min() if 'm/z误差(ppm)' in results_df.columns else 0,
                        results_df['m/z误差(ppm)'].max() if 'm/z误差(ppm)' in results_df.columns else 0
                    ]
                }

                # 添加m/z来源统计
                if 'm/z来源' in results_df.columns:
                    mz_sources = results_df['m/z来源'].value_counts()
                    if len(mz_sources) > 0:
                        for source, count in mz_sources.items():
                            if pd.notna(source):
                                stats['统计项目'].append(f"m/z来源: {source}")
                                stats['数值'].append(count)

                # 添加PPG指数统计
                if include_ppg and 'PPG保留指数' in results_df.columns:
                    ppg_values = results_df['PPG保留指数'].dropna()
                    if len(ppg_values) > 0:
                        stats['统计项目'].extend(['平均PPG指数', '最小PPG指数', '最大PPG指数', 'PPG指数标准差'])
                        stats['数值'].extend([
                            ppg_values.mean(),
                            ppg_values.min(),
                            ppg_values.max(),
                            ppg_values.std()
                        ])

                stats_df = pd.DataFrame(stats)
                stats_df.to_excel(writer, sheet_name='统计信息', index=False)

            saved_files.append(str(excel_file))

            # 保存为CSV（便于其他软件使用）
            csv_file = output_path / f"{base_name}_{timestamp}.csv"
            results_df.to_csv(csv_file, index=False, encoding='utf-8-sig')
            saved_files.append(str(csv_file))

            return True, f"结果已保存到 {excel_file.name}", saved_files

        except Exception as e:
            return False, f"保存结果失败: {str(e)}", []


class MzMLPeakExtractor:
    """mzML峰提取器 - 提取峰数据供化合物匹配"""

    def __init__(self, progress_callback=None, log_callback=None):
        """初始化提取器"""
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

    def extract_peaks_with_pymzml(self, mzml_path: str, intensity_threshold: float = 1000) -> Tuple[
        Optional[List[Dict]], str]:
        """使用pymzML提取峰数据"""
        try:
            import pymzml

            self.log(f"使用pymzML加载文件: {Path(mzml_path).name}")
            peaks = []

            run = pymzml.run.Reader(mzml_path)
            spectrum_count = 0

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
                if spectrum_count % 100 == 0:
                    progress = min(100, spectrum_count // 10)
                    self.update_progress(progress, f"处理谱图 {spectrum_count}")

            if peaks:
                return peaks, f"pymzML提取成功: {len(peaks)}个峰"
            else:
                return None, "pymzML未找到符合条件的峰"

        except ImportError:
            return None, "pymzML未安装"
        except Exception as e:
            return None, f"pymzML处理失败: {str(e)}"

    def extract_peaks_from_mzml(self, mzml_file: str, intensity_threshold: float = 1000) -> Tuple[
        bool, str, Optional[pd.DataFrame]]:
        """从mzML文件提取峰数据"""
        try:
            self.log("开始提取峰数据...")
            self.update_progress(10, "提取峰数据...")

            # 提取峰
            peaks, msg = self.extract_peaks_with_pymzml(mzml_file, intensity_threshold)

            if peaks is None:
                return False, msg, None

            # 转换为DataFrame
            df = pd.DataFrame(peaks)

            # 重命名列以便与化合物匹配器兼容
            column_mapping = {
                '质荷比(m/z)': '质荷比(m/z)',
                '保留时间(RT)': '保留时间(RT)',
                '强度': '强度'
            }

            for old_col, new_col in column_mapping.items():
                if old_col in df.columns:
                    df = df.rename(columns={old_col: new_col})

            self.update_progress(100, "峰提取完成")
            return True, f"成功提取 {len(df)} 个峰", df

        except Exception as e:
            return False, f"提取峰数据失败: {str(e)}", None

    def extract_peaks_from_csv(self, csv_file: str) -> Tuple[bool, str, Optional[pd.DataFrame]]:
        """从CSV文件加载峰数据"""
        try:
            self.log("从CSV文件加载峰数据...")

            df = pd.read_csv(csv_file)

            # 检查必要的列
            required_columns = ['质荷比(m/z)', '保留时间(RT)', '强度']

            # 尝试查找其他可能的列名
            column_mapping = {
                '质荷比(m/z)': ['mz', 'M/Z', 'm/z'],
                '保留时间(RT)': ['rt', 'RT', '保留时间'],
                '强度': ['intensity', 'Intensity', '强度']
            }

            for req_col, possible_names in column_mapping.items():
                if req_col not in df.columns:
                    for name in possible_names:
                        if name in df.columns:
                            df = df.rename(columns={name: req_col})
                            break

            # 再次检查
            missing_cols = [col for col in required_columns if col not in df.columns]
            if missing_cols:
                return False, f"CSV文件缺少必要列: {', '.join(missing_cols)}", None

            return True, f"成功加载 {len(df)} 个峰", df

        except Exception as e:
            return False, f"加载CSV文件失败: {str(e)}", None


class MzMLCompoundMatcherGUI:
    """mzML峰提取与化合物匹配GUI"""

    def __init__(self, root):
        """初始化GUI"""
        self.root = root
        self.root.title("mzML峰提取与化合物匹配及保留指数计算程序")
        self.root.geometry("1100x950")

        # 设置图标
        try:
            self.root.iconbitmap(default='icon.ico')
        except:
            pass

        # 处理器
        self.peak_extractor = None
        self.compound_matcher = None
        self.ppg_calculator = None
        self.processing_thread = None
        self.is_processing = False

        # 创建UI
        self.setup_ui()

        # 绑定关闭事件
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)

    def setup_ui(self):
        """设置UI界面"""
        # 创建主框架和滚动条
        main_frame = tk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True)

        # 创建Canvas和滚动条
        canvas = tk.Canvas(main_frame)
        scrollbar = ttk.Scrollbar(main_frame, orient=tk.VERTICAL, command=canvas.yview)
        scrollable_frame = tk.Frame(canvas)

        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )

        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # 主内容
        content_frame = tk.Frame(scrollable_frame, padx=20, pady=20)
        content_frame.pack(fill=tk.BOTH, expand=True)

        # 标题
        title_label = tk.Label(content_frame, text="mzML峰提取、化合物匹配与PPG保留指数计算程序",
                               font=("Arial", 18, "bold"))
        title_label.pack(pady=(0, 10))

        # 版本信息
        version_label = tk.Label(content_frame, text="版本 9.0 - 支持从SMILES/InChI计算m/z并匹配",
                                 font=("Arial", 10))
        version_label.pack(pady=(0, 20))

        # ==================== 第1步: 峰数据源选择 ====================
        step1_frame = tk.LabelFrame(content_frame, text="第1步: 选择峰数据源", font=("Arial", 12, "bold"),
                                    padx=10, pady=10)
        step1_frame.pack(fill=tk.X, pady=(0, 20))

        # 数据源选择
        self.data_source_var = tk.StringVar(value="mzml")
        ttk.Radiobutton(step1_frame, text="从mzML文件提取峰数据",
                        variable=self.data_source_var, value="mzml",
                        command=self.on_data_source_change).grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)

        ttk.Radiobutton(step1_frame, text="从现有CSV文件加载峰数据",
                        variable=self.data_source_var, value="csv",
                        command=self.on_data_source_change).grid(row=0, column=1, sticky=tk.W, padx=5, pady=5)

        # mzML文件选择
        self.mzml_frame = tk.Frame(step1_frame)
        self.mzml_frame.grid(row=1, column=0, columnspan=2, sticky=tk.W, padx=5, pady=5)

        ttk.Label(self.mzml_frame, text="mzML文件:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5))

        self.mzml_file_var = tk.StringVar()
        mzml_file_entry = ttk.Entry(self.mzml_frame, textvariable=self.mzml_file_var, width=70)
        mzml_file_entry.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=(0, 5))

        ttk.Button(self.mzml_frame, text="浏览...", command=self.browse_mzml_file).grid(row=0, column=2, padx=(0, 5))

        # CSV文件选择
        self.csv_frame = tk.Frame(step1_frame)
        self.csv_frame.grid(row=2, column=0, columnspan=2, sticky=tk.W, padx=5, pady=5)

        ttk.Label(self.csv_frame, text="CSV文件:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5))

        self.csv_file_var = tk.StringVar()
        csv_file_entry = ttk.Entry(self.csv_frame, textvariable=self.csv_file_var, width=70)
        csv_file_entry.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=(0, 5))

        ttk.Button(self.csv_frame, text="浏览...", command=self.browse_csv_file).grid(row=0, column=2, padx=(0, 5))

        # 强度阈值
        ttk.Label(step1_frame, text="强度阈值:").grid(row=3, column=0, sticky=tk.W, padx=5, pady=5)

        self.intensity_threshold_var = tk.StringVar(value="1000")
        threshold_entry = ttk.Entry(step1_frame, textvariable=self.intensity_threshold_var, width=15)
        threshold_entry.grid(row=3, column=1, sticky=tk.W, padx=5, pady=5)

        # 初始隐藏CSV框架
        self.csv_frame.grid_remove()

        # ==================== 第2步: 化合物信息加载 ====================
        step2_frame = tk.LabelFrame(content_frame, text="第2步: 加载化合物信息", font=("Arial", 12, "bold"),
                                    padx=10, pady=10)
        step2_frame.pack(fill=tk.X, pady=(0, 20))

        ttk.Label(step2_frame, text="化合物Excel文件:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)

        self.compound_file_var = tk.StringVar()
        compound_file_entry = ttk.Entry(step2_frame, textvariable=self.compound_file_var, width=70)
        compound_file_entry.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=(0, 5), pady=5)

        ttk.Button(step2_frame, text="浏览...", command=self.browse_compound_file).grid(row=0, column=2, pady=5)

        # 自动计算m/z复选框
        self.calculate_mz_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(step2_frame, text="自动计算m/z（从SMILES/InChI/MF）",
                        variable=self.calculate_mz_var).grid(row=1, column=0, sticky=tk.W, padx=5, pady=5)

        # 离子模式选择
        ttk.Label(step2_frame, text="离子模式:").grid(row=1, column=1, sticky=tk.W, padx=5, pady=5)

        self.ion_mode_var = tk.StringVar(value="M+H")
        ion_mode_combo = ttk.Combobox(step2_frame, textvariable=self.ion_mode_var,
                                      values=["M+H", "M-H", "M+Na", "M+K", "M+NH4", "M+CH3COO",
                                              "M+2H", "M+2Na", "M+FA-H"], width=15)
        ion_mode_combo.grid(row=1, column=2, sticky=tk.W, padx=5, pady=5)

        # 预览按钮
        ttk.Button(step2_frame, text="预览化合物", command=self.preview_compounds).grid(row=2, column=0, sticky=tk.W,
                                                                                        padx=5, pady=5)

        # 加载状态
        self.compound_status_var = tk.StringVar(value="未加载化合物信息")
        ttk.Label(step2_frame, textvariable=self.compound_status_var).grid(row=2, column=1, columnspan=2, sticky=tk.W,
                                                                           padx=5, pady=5)

        # ==================== 第3步: PPG保留指数计算设置 ====================
        step3_frame = tk.LabelFrame(content_frame, text="第3步: PPG保留指数计算设置", font=("Arial", 12, "bold"),
                                    padx=10, pady=10)
        step3_frame.pack(fill=tk.X, pady=(0, 20))

        # 启用PPG计算复选框
        self.enable_ppg_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(step3_frame, text="启用PPG保留指数计算",
                        variable=self.enable_ppg_var,
                        command=self.on_ppg_enable_change).grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)

        # PPG标准品文件选择
        self.ppg_frame = tk.Frame(step3_frame)
        self.ppg_frame.grid(row=1, column=0, columnspan=3, sticky=tk.W, padx=5, pady=5)

        ttk.Label(self.ppg_frame, text="PPG标准品Excel文件:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5))

        self.ppg_file_var = tk.StringVar()
        ppg_file_entry = ttk.Entry(self.ppg_frame, textvariable=self.ppg_file_var, width=60)
        ppg_file_entry.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=(0, 5))

        ttk.Button(self.ppg_frame, text="浏览...", command=self.browse_ppg_file).grid(row=0, column=2, padx=(0, 5))

        # PPG校准方法
        ttk.Label(self.ppg_frame, text="校准方法:").grid(row=1, column=0, sticky=tk.W, padx=(0, 5), pady=(10, 5))

        self.ppg_method_var = tk.StringVar(value="interpolation")
        ttk.Combobox(self.ppg_frame, textvariable=self.ppg_method_var,
                     values=["interpolation", "linear"], width=15, state="readonly").grid(row=1, column=1, sticky=tk.W,
                                                                                          padx=(0, 5), pady=(10, 5))

        # 预览PPG按钮
        ttk.Button(self.ppg_frame, text="预览PPG标准品", command=self.preview_ppg_standards).grid(row=1, column=2,
                                                                                                  pady=(10, 5))

        # 初始禁用PPG框架
        self.ppg_frame.grid_remove()

        # ==================== 第4步: 匹配参数设置 ====================
        step4_frame = tk.LabelFrame(content_frame, text="第4步: 设置匹配参数", font=("Arial", 12, "bold"),
                                    padx=10, pady=10)
        step4_frame.pack(fill=tk.X, pady=(0, 20))

        # m/z误差 (ppm)
        ttk.Label(step4_frame, text="m/z误差 (ppm):").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)

        self.ppm_tolerance_var = tk.StringVar(value="10")
        ppm_entry = ttk.Entry(step4_frame, textvariable=self.ppm_tolerance_var, width=15)
        ppm_entry.grid(row=0, column=1, sticky=tk.W, padx=5, pady=5)

        ttk.Label(step4_frame, text="ppm").grid(row=0, column=2, sticky=tk.W, padx=(0, 20), pady=5)

        # RT窗口 (秒)
        ttk.Label(step4_frame, text="RT窗口 (秒):").grid(row=0, column=3, sticky=tk.W, padx=5, pady=5)

        self.rt_window_var = tk.StringVar(value="30")
        rt_entry = ttk.Entry(step4_frame, textvariable=self.rt_window_var, width=15)
        rt_entry.grid(row=0, column=4, sticky=tk.W, padx=5, pady=5)

        ttk.Label(step4_frame, text="秒").grid(row=0, column=5, sticky=tk.W, padx=(0, 5), pady=5)

        # 输出目录
        ttk.Label(step4_frame, text="输出目录:").grid(row=1, column=0, sticky=tk.W, padx=5, pady=(10, 5))

        self.output_dir_var = tk.StringVar()
        output_dir_entry = ttk.Entry(step4_frame, textvariable=self.output_dir_var, width=70)
        output_dir_entry.grid(row=1, column=1, columnspan=4, sticky=(tk.W, tk.E), padx=(0, 5), pady=(10, 5))

        ttk.Button(step4_frame, text="浏览...", command=self.browse_output_dir).grid(row=1, column=5, pady=(10, 5))

        # ==================== 第5步: 单化合物m/z计算 ====================
        step5_frame = tk.LabelFrame(content_frame, text="第5步: 单化合物m/z计算（可选）", font=("Arial", 12, "bold"),
                                    padx=10, pady=10)
        step5_frame.pack(fill=tk.X, pady=(0, 20))

        # 输入类型选择
        ttk.Label(step5_frame, text="输入类型:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)

        self.single_input_type_var = tk.StringVar(value="smiles")
        ttk.Radiobutton(step5_frame, text="SMILES", variable=self.single_input_type_var,
                        value="smiles").grid(row=0, column=1, sticky=tk.W, padx=5, pady=5)
        ttk.Radiobutton(step5_frame, text="InChI", variable=self.single_input_type_var,
                        value="inchi").grid(row=0, column=2, sticky=tk.W, padx=5, pady=5)
        ttk.Radiobutton(step5_frame, text="分子式", variable=self.single_input_type_var,
                        value="formula").grid(row=0, column=3, sticky=tk.W, padx=5, pady=5)

        # 输入框
        ttk.Label(step5_frame, text="输入:").grid(row=1, column=0, sticky=tk.W, padx=5, pady=5)

        self.single_input_var = tk.StringVar()
        single_input_entry = ttk.Entry(step5_frame, textvariable=self.single_input_var, width=70)
        single_input_entry.grid(row=1, column=1, columnspan=3, sticky=(tk.W, tk.E), padx=(0, 5), pady=5)

        # 计算按钮
        ttk.Button(step5_frame, text="计算m/z", command=self.calculate_single_mz).grid(row=1, column=4, padx=5, pady=5)

        # 结果显示
        self.single_result_var = tk.StringVar(value="输入化合物信息以计算m/z")
        ttk.Label(step5_frame, textvariable=self.single_result_var).grid(row=2, column=0, columnspan=5,
                                                                         sticky=tk.W, padx=5, pady=5)

        # ==================== 第6步: 进度显示 ====================
        step6_frame = tk.LabelFrame(content_frame, text="第6步: 处理进度", font=("Arial", 12, "bold"),
                                    padx=10, pady=10)
        step6_frame.pack(fill=tk.X, pady=(0, 20))

        self.progress_var = tk.StringVar(value="准备就绪")
        ttk.Label(step6_frame, textvariable=self.progress_var).pack(anchor=tk.W, pady=(0, 5))

        self.progress_bar = ttk.Progressbar(step6_frame, mode='determinate', length=950)
        self.progress_bar.pack(fill=tk.X, pady=(0, 5))

        self.progress_percent = tk.StringVar(value="0%")
        ttk.Label(step6_frame, textvariable=self.progress_percent).pack(anchor=tk.E)

        # ==================== 第7步: 日志显示 ====================
        step7_frame = tk.LabelFrame(content_frame, text="第7步: 处理日志", font=("Arial", 12, "bold"),
                                    padx=10, pady=10)
        step7_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 20))

        # 日志文本框
        self.log_text = scrolledtext.ScrolledText(step7_frame, width=130, height=15,
                                                  wrap=tk.WORD, font=("Consolas", 10))
        self.log_text.pack(fill=tk.BOTH, expand=True)

        # 配置标签颜色
        self.log_text.tag_config("INFO", foreground="black")
        self.log_text.tag_config("SUCCESS", foreground="green")
        self.log_text.tag_config("WARNING", foreground="orange")
        self.log_text.tag_config("ERROR", foreground="red")

        # ==================== 按钮区域 ====================
        button_frame = tk.Frame(content_frame)
        button_frame.pack(fill=tk.X, pady=(0, 10))

        self.process_button = ttk.Button(button_frame, text="开始处理",
                                         command=self.start_processing, width=20)
        self.process_button.pack(side=tk.LEFT, padx=5)

        self.cancel_button = ttk.Button(button_frame, text="取消处理",
                                        command=self.cancel_processing, width=20,
                                        state=tk.DISABLED)
        self.cancel_button.pack(side=tk.LEFT, padx=5)

        ttk.Button(button_frame, text="清空日志", command=self.clear_log).pack(side=tk.LEFT, padx=5)

        ttk.Button(button_frame, text="打开输出目录", command=self.open_output_dir).pack(side=tk.LEFT, padx=5)

        ttk.Button(button_frame, text="退出", command=self.on_closing).pack(side=tk.LEFT, padx=5)

        # 状态栏
        self.status_var = tk.StringVar(value="就绪")
        status_bar = ttk.Label(content_frame, textvariable=self.status_var, relief=tk.SUNKEN)
        status_bar.pack(fill=tk.X, pady=(0, 10))

        # 初始化
        self.on_data_source_change()

    def on_data_source_change(self):
        """数据源改变事件"""
        source = self.data_source_var.get()

        if source == "mzml":
            self.mzml_frame.grid()
            self.csv_frame.grid_remove()
        else:
            self.mzml_frame.grid_remove()
            self.csv_frame.grid()

    def on_ppg_enable_change(self):
        """PPG计算启用状态改变事件"""
        if self.enable_ppg_var.get():
            self.ppg_frame.grid()
        else:
            self.ppg_frame.grid_remove()

    def browse_mzml_file(self):
        """浏览mzML文件"""
        file_types = [("mzML文件", "*.mzML *.mzML.gz"), ("所有文件", "*.*")]
        file_path = filedialog.askopenfilename(title="选择mzML文件", filetypes=file_types)

        if file_path:
            self.mzml_file_var.set(file_path)

    def browse_csv_file(self):
        """浏览CSV文件"""
        file_types = [("CSV文件", "*.csv"), ("所有文件", "*.*")]
        file_path = filedialog.askopenfilename(title="选择CSV文件", filetypes=file_types)

        if file_path:
            self.csv_file_var.set(file_path)

    def browse_compound_file(self):
        """浏览化合物Excel文件"""
        file_types = [("Excel文件", "*.xlsx *.xls"), ("所有文件", "*.*")]
        file_path = filedialog.askopenfilename(title="选择化合物Excel文件", filetypes=file_types)

        if file_path:
            self.compound_file_var.set(file_path)
            # 尝试预览
            self.preview_compounds()

    def browse_ppg_file(self):
        """浏览PPG标准品Excel文件"""
        file_types = [("Excel文件", "*.xlsx *.xls"), ("所有文件", "*.*")]
        file_path = filedialog.askopenfilename(title="选择PPG标准品Excel文件", filetypes=file_types)

        if file_path:
            self.ppg_file_var.set(file_path)
            # 尝试预览
            self.preview_ppg_standards()

    def browse_output_dir(self):
        """浏览输出目录"""
        dir_path = filedialog.askdirectory(title="选择输出目录")

        if dir_path:
            self.output_dir_var.set(dir_path)

    def preview_compounds(self):
        """预览化合物信息"""
        compound_file = self.compound_file_var.get()

        if not compound_file:
            messagebox.showwarning("警告", "请先选择化合物Excel文件")
            return

        try:
            matcher = CompoundMatcher()
            calculate_mz = self.calculate_mz_var.get()
            ion_mode = self.ion_mode_var.get()

            success, msg, df = matcher.load_compounds_from_excel(
                compound_file,
                calculate_mz=calculate_mz,
                ion_mode=ion_mode
            )

            if success:
                # 创建预览窗口
                preview_window = tk.Toplevel(self.root)
                preview_window.title("化合物信息预览")
                preview_window.geometry("900x600")

                # 创建Treeview显示数据
                tree_frame = tk.Frame(preview_window)
                tree_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

                # 滚动条
                tree_scroll_y = ttk.Scrollbar(tree_frame)
                tree_scroll_y.pack(side=tk.RIGHT, fill=tk.Y)

                tree_scroll_x = ttk.Scrollbar(tree_frame, orient=tk.HORIZONTAL)
                tree_scroll_x.pack(side=tk.BOTTOM, fill=tk.X)

                # Treeview
                tree = ttk.Treeview(tree_frame, yscrollcommand=tree_scroll_y.set,
                                    xscrollcommand=tree_scroll_x.set)
                tree.pack(fill=tk.BOTH, expand=True)

                tree_scroll_y.config(command=tree.yview)
                tree_scroll_x.config(command=tree.xview)

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
                info_frame = tk.Frame(preview_window)
                info_frame.pack(fill=tk.X, padx=10, pady=(0, 10))

                tk.Label(info_frame, text=f"共加载 {len(df)} 个化合物").pack(side=tk.LEFT)

                if 'mz' in df.columns:
                    mz_range = f"m/z范围: {df['mz'].min():.4f} - {df['mz'].max():.4f}"
                    tk.Label(info_frame, text=mz_range).pack(side=tk.LEFT, padx=20)

                if 'mz_source' in df.columns:
                    calculated_count = df['mz_source'].notna().sum()
                    tk.Label(info_frame, text=f"自动计算m/z: {calculated_count}个").pack(side=tk.LEFT, padx=20)

                # 关闭按钮
                tk.Button(preview_window, text="关闭", command=preview_window.destroy).pack(pady=10)

                self.compound_status_var.set(f"已加载 {len(df)} 个化合物")
                self.log_message(f"✓ 成功加载 {len(df)} 个化合物", "SUCCESS")
            else:
                messagebox.showerror("错误", msg)
                self.compound_status_var.set("加载失败")
                self.log_message(f"✗ 加载化合物失败: {msg}", "ERROR")

        except Exception as e:
            messagebox.showerror("错误", f"预览失败: {str(e)}")

    def preview_ppg_standards(self):
        """预览PPG标准品信息"""
        ppg_file = self.ppg_file_var.get()

        if not ppg_file:
            messagebox.showwarning("警告", "请先选择PPG标准品Excel文件")
            return

        try:
            calculator = PPGIndexCalculator()
            success, msg, df = calculator.load_ppg_standards(ppg_file)

            if success:
                # 创建预览窗口
                preview_window = tk.Toplevel(self.root)
                preview_window.title("PPG标准品信息预览")
                preview_window.geometry("800x600")

                # 创建Treeview显示数据
                tree_frame = tk.Frame(preview_window)
                tree_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

                # 滚动条
                tree_scroll_y = ttk.Scrollbar(tree_frame)
                tree_scroll_y.pack(side=tk.RIGHT, fill=tk.Y)

                tree_scroll_x = ttk.Scrollbar(tree_frame, orient=tk.HORIZONTAL)
                tree_scroll_x.pack(side=tk.BOTTOM, fill=tk.X)

                # Treeview
                tree = ttk.Treeview(tree_frame, yscrollcommand=tree_scroll_y.set,
                                    xscrollcommand=tree_scroll_x.set)
                tree.pack(fill=tk.BOTH, expand=True)

                tree_scroll_y.config(command=tree.yview)
                tree_scroll_x.config(command=tree.xview)

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
                info_frame = tk.Frame(preview_window)
                info_frame.pack(fill=tk.X, padx=10, pady=(0, 10))

                tk.Label(info_frame, text=f"共加载 {len(df)} 个PPG标准品").pack(side=tk.LEFT)

                if 'n' in df.columns and 'RT' in df.columns:
                    n_range = f"聚合度范围: {df['n'].min():.1f} - {df['n'].max():.1f}"
                    rt_range = f"保留时间范围: {df['RT'].min():.3f} - {df['RT'].max():.3f}"
                    tk.Label(info_frame, text=n_range).pack(side=tk.LEFT, padx=20)
                    tk.Label(info_frame, text=rt_range).pack(side=tk.LEFT, padx=20)

                # 关闭按钮
                tk.Button(preview_window, text="关闭", command=preview_window.destroy).pack(pady=10)

                self.log_message(f"✓ 成功加载 {len(df)} 个PPG标准品", "SUCCESS")
            else:
                messagebox.showerror("错误", msg)
                self.log_message(f"✗ 加载PPG标准品失败: {msg}", "ERROR")

        except Exception as e:
            messagebox.showerror("错误", f"预览失败: {str(e)}")

    def calculate_single_mz(self):
        """计算单个化合物的m/z值"""
        input_text = self.single_input_var.get().strip()
        input_type = self.single_input_type_var.get()
        ion_mode = self.ion_mode_var.get()

        if not input_text:
            messagebox.showwarning("警告", "请输入化合物信息")
            return

        if not RDKIT_AVAILABLE:
            messagebox.showerror("错误", "RDKit未安装，无法计算m/z值")
            self.log_message("✗ RDKit未安装，无法计算m/z值", "ERROR")
            return

        try:
            matcher = CompoundMatcher()
            success, msg, mz = matcher.calculate_mz_for_compound(
                input_text, input_type, ion_mode
            )

            if success:
                result_text = f"计算成功: {msg}\nm/z = {mz:.6f}"
                self.single_result_var.set(result_text)
                self.log_message(f"✓ 计算成功: {input_type} → {mz:.6f} ({ion_mode})", "SUCCESS")
            else:
                self.single_result_var.set(f"计算失败: {msg}")
                self.log_message(f"✗ 计算失败: {msg}", "ERROR")

        except Exception as e:
            error_msg = f"计算过程中出错: {str(e)}"
            self.single_result_var.set(error_msg)
            self.log_message(f"✗ {error_msg}", "ERROR")

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
        # 检查数据源
        source = self.data_source_var.get()

        if source == "mzml":
            mzml_file = self.mzml_file_var.get().strip()
            if not mzml_file:
                return False, "请选择mzML文件"
            if not os.path.exists(mzml_file):
                return False, f"mzML文件不存在: {mzml_file}"
        else:
            csv_file = self.csv_file_var.get().strip()
            if not csv_file:
                return False, "请选择CSV文件"
            if not os.path.exists(csv_file):
                return False, f"CSV文件不存在: {csv_file}"

        # 检查化合物文件
        compound_file = self.compound_file_var.get().strip()
        if not compound_file:
            return False, "请选择化合物Excel文件"
        if not os.path.exists(compound_file):
            return False, f"化合物文件不存在: {compound_file}"

        # 检查PPG设置
        if self.enable_ppg_var.get():
            ppg_file = self.ppg_file_var.get().strip()
            if not ppg_file:
                return False, "请选择PPG标准品Excel文件"
            if not os.path.exists(ppg_file):
                return False, f"PPG标准品文件不存在: {ppg_file}"

        # 检查输出目录
        output_dir = self.output_dir_var.get().strip()
        if not output_dir:
            return False, "请选择输出目录"

        # 检查参数
        try:
            ppm_tolerance = float(self.ppm_tolerance_var.get())
            if ppm_tolerance <= 0:
                return False, "m/z误差必须大于0"
        except ValueError:
            return False, "m/z误差必须是数字"

        try:
            rt_window = float(self.rt_window_var.get())
            if rt_window < 0:
                return False, "RT窗口不能为负数"
        except ValueError:
            return False, "RT窗口必须是数字"

        try:
            intensity_threshold = float(self.intensity_threshold_var.get())
            if intensity_threshold < 0:
                return False, "强度阈值不能为负数"
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
        source = self.data_source_var.get()
        compound_file = self.compound_file_var.get()
        output_dir = self.output_dir_var.get()
        enable_ppg = self.enable_ppg_var.get()
        calculate_mz = self.calculate_mz_var.get()
        ion_mode = self.ion_mode_var.get()

        ppm_tolerance = float(self.ppm_tolerance_var.get())
        rt_window = float(self.rt_window_var.get())
        intensity_threshold = float(self.intensity_threshold_var.get())

        # 创建处理器
        self.peak_extractor = MzMLPeakExtractor(
            progress_callback=self.update_progress,
            log_callback=self.log_message
        )

        self.compound_matcher = CompoundMatcher()

        if enable_ppg:
            self.ppg_calculator = PPGIndexCalculator()
            ppg_file = self.ppg_file_var.get()
            ppg_method = self.ppg_method_var.get()
        else:
            self.ppg_calculator = None
            ppg_file = None
            ppg_method = None

        # 在新线程中处理
        if source == "mzml":
            mzml_file = self.mzml_file_var.get()
            self.processing_thread = threading.Thread(
                target=self.process_mzml_match,
                args=(mzml_file, compound_file, output_dir, enable_ppg, calculate_mz,
                      ppm_tolerance, rt_window, intensity_threshold, ion_mode,
                      ppg_file, ppg_method)
            )
        else:
            csv_file = self.csv_file_var.get()
            self.processing_thread = threading.Thread(
                target=self.process_csv_match,
                args=(csv_file, compound_file, output_dir, enable_ppg, calculate_mz,
                      ppm_tolerance, rt_window, intensity_threshold, ion_mode,
                      ppg_file, ppg_method)
            )

        self.is_processing = True
        self.processing_thread.start()

        # 启动进度监视器
        self.root.after(100, self.check_processing_status)

    def process_mzml_match(self, mzml_file: str, compound_file: str, output_dir: str,
                           enable_ppg: bool, calculate_mz: bool,
                           ppm_tolerance: float, rt_window: float,
                           intensity_threshold: float, ion_mode: str,
                           ppg_file: str, ppg_method: str):
        """处理mzML文件和化合物匹配（在线程中）"""
        try:
            self.log_message("步骤1: 从mzML文件提取峰数据...")
            self.update_progress(20, "提取峰数据...")

            # 提取峰数据
            success, msg, peaks_df = self.peak_extractor.extract_peaks_from_mzml(
                mzml_file, intensity_threshold
            )

            if not success:
                self.root.after(0, lambda: messagebox.showerror("错误", msg))
                self.processing_finished()
                return

            self.log_message(f"✓ 成功提取 {len(peaks_df)} 个峰", "SUCCESS")

            # 进行化合物匹配
            self.process_compound_matching(peaks_df, compound_file, output_dir,
                                           enable_ppg, calculate_mz,
                                           ppm_tolerance, rt_window,
                                           mzml_file, ion_mode,
                                           ppg_file, ppg_method)

        except Exception as e:
            error_msg = f"处理过程中出错: {str(e)}"
            self.root.after(0, lambda: self.log_message(error_msg, "ERROR"))
            self.root.after(0, self.processing_finished)

    def process_csv_match(self, csv_file: str, compound_file: str, output_dir: str,
                          enable_ppg: bool, calculate_mz: bool,
                          ppm_tolerance: float, rt_window: float,
                          intensity_threshold: float, ion_mode: str,
                          ppg_file: str, ppg_method: str):
        """处理CSV文件和化合物匹配（在线程中）"""
        try:
            self.log_message("步骤1: 从CSV文件加载峰数据...")
            self.update_progress(20, "加载峰数据...")

            # 加载峰数据
            success, msg, peaks_df = self.peak_extractor.extract_peaks_from_csv(csv_file)

            if not success:
                self.root.after(0, lambda: messagebox.showerror("错误", msg))
                self.processing_finished()
                return

            self.log_message(f"✓ 成功加载 {len(peaks_df)} 个峰", "SUCCESS")

            # 进行化合物匹配
            self.process_compound_matching(peaks_df, compound_file, output_dir,
                                           enable_ppg, calculate_mz,
                                           ppm_tolerance, rt_window,
                                           csv_file, ion_mode,
                                           ppg_file, ppg_method)

        except Exception as e:
            error_msg = f"处理过程中出错: {str(e)}"
            self.root.after(0, lambda: self.log_message(error_msg, "ERROR"))
            self.root.after(0, self.processing_finished)

    def process_compound_matching(self, peaks_df: pd.DataFrame, compound_file: str,
                                  output_dir: str, enable_ppg: bool, calculate_mz: bool,
                                  ppm_tolerance: float, rt_window: float,
                                  source_file: str, ion_mode: str,
                                  ppg_file: str, ppg_method: str):
        """执行化合物匹配"""
        try:
            self.log_message("步骤2: 加载化合物信息...")
            self.update_progress(40, "加载化合物信息...")

            # 设置匹配参数
            self.compound_matcher.set_match_settings(
                ppm_tolerance=ppm_tolerance,
                rt_window=rt_window,
                intensity_threshold=1000,  # 已在提取时过滤
                ion_mode=ion_mode,
                calculate_mz=calculate_mz
            )

            # 加载化合物信息
            success, msg, _ = self.compound_matcher.load_compounds_from_excel(
                compound_file,
                calculate_mz=calculate_mz,
                ion_mode=ion_mode
            )
            if not success:
                self.root.after(0, lambda: messagebox.showerror("错误", msg))
                self.processing_finished()
                return

            self.log_message(f"✓ 成功加载化合物信息", "SUCCESS")

            # 如果需要，加载PPG标准品并建立校准曲线
            if enable_ppg and self.ppg_calculator and ppg_file:
                self.log_message("步骤2.5: 加载PPG标准品并建立校准曲线...")
                self.update_progress(45, "建立PPG校准曲线...")

                success, msg, _ = self.ppg_calculator.load_ppg_standards(ppg_file)
                if not success:
                    self.root.after(0, lambda: messagebox.showerror("PPG错误", msg))
                    self.processing_finished()
                    return

                self.log_message(f"✓ 成功加载PPG标准品", "SUCCESS")

                success, msg = self.ppg_calculator.build_calibration_curve(method=ppg_method)
                if not success:
                    self.root.after(0, lambda: messagebox.showerror("PPG错误", msg))
                    self.processing_finished()
                    return

                self.log_message(f"✓ {msg}", "SUCCESS")

            self.log_message("步骤3: 匹配化合物...")
            self.update_progress(60, "匹配化合物...")

            # 执行匹配
            success, msg = self.compound_matcher.match_compounds(peaks_df)
            if not success:
                self.root.after(0, lambda: messagebox.showwarning("警告", msg))

            self.log_message(f"✓ 化合物匹配完成", "SUCCESS")

            # 如果需要，添加PPG保留指数
            if enable_ppg and self.ppg_calculator:
                self.log_message("步骤3.5: 计算PPG保留指数...")
                self.update_progress(70, "计算PPG保留指数...")

                success, msg = self.compound_matcher.add_ppg_index(self.ppg_calculator)
                if not success:
                    self.log_message(f"✗ {msg}", "WARNING")
                else:
                    self.log_message(f"✓ {msg}", "SUCCESS")

                # 保存PPG校准报告
                success, msg, _ = self.ppg_calculator.save_calibration_report(output_dir)
                if success:
                    self.log_message(f"✓ PPG校准报告已保存", "SUCCESS")

            self.log_message("步骤4: 保存匹配结果...")
            self.update_progress(80, "保存匹配结果...")

            # 保存结果
            base_name = Path(source_file).stem
            success, msg, saved_files = self.compound_matcher.save_match_results(
                output_dir, f"{base_name}_化合物匹配", include_ppg=enable_ppg
            )

            if success:
                self.log_message(f"✓ 结果已保存", "SUCCESS")
                self.log_message(f"  生成的文件:", "INFO")
                for file_path in saved_files:
                    self.log_message(f"    • {Path(file_path).name}", "INFO")

                # 显示成功对话框
                self.root.after(0, lambda: messagebox.showinfo(
                    "处理成功",
                    f"化合物匹配完成！\n结果已保存到: {output_dir}"
                ))
            else:
                self.log_message(f"✗ 保存结果失败: {msg}", "ERROR")
                self.root.after(0, lambda: messagebox.showerror("错误", f"保存结果失败: {msg}"))

            self.update_progress(100, "处理完成")
            self.processing_finished()

        except Exception as e:
            error_msg = f"化合物匹配过程中出错: {str(e)}"
            self.root.after(0, lambda: self.log_message(error_msg, "ERROR"))
            self.root.after(0, self.processing_finished)

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
        self.peak_extractor = None
        self.compound_matcher = None
        self.ppg_calculator = None

        # 启用开始按钮，禁用取消按钮
        self.process_button.config(state=tk.NORMAL)
        self.cancel_button.config(state=tk.DISABLED)

        # 更新状态
        self.update_status("处理完成")

    def cancel_processing(self):
        """取消处理"""
        if self.peak_extractor and self.is_processing:
            self.peak_extractor.cancel_requested = True
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
                if self.peak_extractor:
                    self.peak_extractor.cancel_requested = True
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
        'scipy': '科学计算库（用于PPG指数计算）',
        'openpyxl': 'Excel文件支持',
        'pymzml': 'mzML文件读取 (峰提取)',
        'rdkit': '化学信息学库（用于从SMILES/InChI计算m/z）'
    }

    missing = []
    warnings = []

    for lib, desc in dependencies.items():
        try:
            __import__(lib)
            print(f"✓ {lib}: {desc}")
        except ImportError:
            if lib == 'rdkit':
                print(f"⚠ {lib}: {desc} - 未安装，m/z计算功能不可用")
                warnings.append(lib)
            else:
                print(f"✗ {lib}: {desc} - 未安装")
                missing.append(lib)

    if missing:
        print(f"\n缺少以下库:")
        for lib in missing:
            print(f"  pip install {lib}")

    if warnings:
        print(f"\n以下库未安装，相关功能不可用:")
        for lib in warnings:
            print(f"  conda install -c conda-forge {lib}")

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
    window_width = 1100
    window_height = 950
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()
    x = (screen_width - window_width) // 2
    y = (screen_height - window_height) // 2
    root.geometry(f"{window_width}x{window_height}+{x}+{y}")

    # 创建GUI
    app = MzMLCompoundMatcherGUI(root)

    # 运行主循环
    root.mainloop()


if __name__ == "__main__":
    print("mzML峰提取、化合物匹配与PPG保留指数计算程序")
    print("版本: 9.0")
    print("功能:")
    print("  1. 从SMILES/InChI/分子式计算化合物的m/z值")
    print("  2. 从mzML文件提取峰数据")
    print("  3. 从现有CSV文件加载峰数据")
    print("  4. 从Excel文件加载化合物信息")
    print("  5. 根据化合物m/z匹配峰数据")
    print("  6. 计算化合物的PPG保留指数")
    print("  7. 生成详细的匹配结果和保留指数Excel文件")
    print("=" * 70)

    main()