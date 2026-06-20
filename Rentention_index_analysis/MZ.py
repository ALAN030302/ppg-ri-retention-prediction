#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
mzML peak extraction, compound matching, and retention index calculation program
:
1. Calculate compound m/z values from SMILES/InChI
2. Extract peak data from mzML files
3. Load compound information from Excel tables
4. Match peak data by compound m/z
5. Calculate compound PPG retention indices
6. Save matching results and retention indices to a new Excel table
Author: Qianlei Yao
Version: 9.0
Date: 2025
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

# Suppress warnings
warnings.filterwarnings('ignore')

try:
    import pandas as pd
    import numpy as np
    from scipy import stats
    from scipy.interpolate import interp1d
except ImportError:
    print("Error: please install pandas, numpy, and scipy first")
    print("Installation command: pip install pandas numpy scipy")
    sys.exit(1)

# Try to import cheminformatics libraries
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    from rdkit.Chem import AllChem
    from rdkit.Chem import rdMolDescriptors

    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("Warning: RDKit is not installed, so m/z cannot be calculated from SMILES/InChI")
    print("Installation command: conda install -c conda-forge rdkit")

# Try to import GUI libraries
try:
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox, scrolledtext

    GUI_AVAILABLE = True
except ImportError:
    GUI_AVAILABLE = False
    print("Warning: tkinter is not installed, so the GUI is unavailable")


class MzCalculator:
    """m/z calculator - calculate compound m/z values from SMILES/InChI"""

    def __init__(self):
        """Initialize the m/z calculator"""
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
            'M+FA-H': ('-', 44.998203), # formate
            'M+Cl': ('-', 34.969402) # chloride ion
        }

        # Masses of common elements
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
        """Calculate monoisotopic mass from a molecular formula"""
        if not RDKIT_AVAILABLE:
            raise ImportError("RDKit is not installed, so formula mass cannot be calculated")

        try:
            # Use RDKit to calculate mass from a molecular formula
            mol = Chem.MolFromSmiles('') # Create an empty molecule
            mol = Chem.AddHs(mol) # Add hydrogens

            # Parse molecular formula
            pattern = r'([A-Z][a-z]*)(\d*)'
            elements = re.findall(pattern, formula)

            # calculate
            total_mass = 0.0
            for element, count in elements:
                element = element.capitalize()
                if count == '':
                    count = 1
                else:
                    count = int(count)

                # translated note
                atomic_num = Chem.GetPeriodicTable().GetAtomicNumber(element)
                # translated note
                isotopes = Chem.GetPeriodicTable().GetMostCommonIsotopeMass(atomic_num)
                if element in self.element_masses:
                    mass = self.element_masses[element] * count
                else:
                    # UseRDKit
                    mass = isotopes * count

                total_mass += mass

            return total_mass

        except Exception as e:
            # RDKitfailed, calculate
            try:
                return self._simple_formula_mass(formula)
            except:
                raise ValueError(f"Parse molecular formula {formula}: {str(e)}")

    def _simple_formula_mass(self, formula: str) -> float:
        """molecular_formulacalculate (UseRDKit)"""
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
                # Use
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
                    raise ValueError(f": {element}")

        return total_mass

    def calculate_mz_from_smiles(self, smiles: str, ion_mode: str = 'M+H', charge: int = 1) -> float:
        """SMILEScalculatem/z

        Parameters:
            smiles: SMILES
            ion_mode:
            charge:

        Returns:
            m/z
        """
        if not RDKIT_AVAILABLE:
            raise ImportError("RDKit, SMILEScalculatem/z")

        try:
            # SMILESmolecule
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f"SMILES: {smiles}")

            # calculate
            mass = rdMolDescriptors.CalcExactMolWt(mol)

            # translated note
            mz = self._apply_ionization(mass, ion_mode, charge)

            return mz

        except Exception as e:
            raise ValueError(f"SMILEScalculatem/zfailed: {str(e)}")

    def calculate_mz_from_inchi(self, inchi: str, ion_mode: str = 'M+H', charge: int = 1) -> float:
        """InChIcalculatem/z"""
        if not RDKIT_AVAILABLE:
            raise ImportError("RDKit, InChIcalculatem/z")

        try:
            # InChImolecule
            mol = Chem.MolFromInchi(inchi)
            if mol is None:
                raise ValueError(f"InChI: {inchi}")

            # calculate
            mass = rdMolDescriptors.CalcExactMolWt(mol)

            # translated note
            mz = self._apply_ionization(mass, ion_mode, charge)

            return mz

        except Exception as e:
            raise ValueError(f"InChIcalculatem/zfailed: {str(e)}")

    def calculate_mz_from_formula(self, formula: str, ion_mode: str = 'M+H', charge: int = 1) -> float:
        """molecular_formulacalculatem/z"""
        try:
            # calculate
            mass = self.calculate_monoisotopic_mass_from_formula(formula)

            # translated note
            mz = self._apply_ionization(mass, ion_mode, charge)

            return mz

        except Exception as e:
            raise ValueError(f"molecular_formulacalculatem/zfailed: {str(e)}")

    def _apply_ionization(self, mass: float, ion_mode: str, charge: int) -> float:
        """"""
        if ion_mode not in self.ion_modes:
            raise ValueError(f": {ion_mode}")

        polarity, offset = self.ion_modes[ion_mode]

        # calculatem/z
        if charge == 0:
            raise ValueError("0")

        mz = (mass + offset) / abs(charge)

        return mz

    def detect_ion_mode_from_mz(self, observed_mz: float, calculated_mass: float,
                                tolerance_ppm: float = 10) -> List[Tuple[str, float]]:
        """m/zcalculate

        Returns:
            column, namematch(ppm)
        """
        possible_modes = []

        for mode, (polarity, offset) in self.ion_modes.items():
            for charge in [1, 2, 3]: # 1, 2, 3
                expected_mz = (calculated_mass + offset) / charge
                error_ppm = abs((observed_mz - expected_mz) / expected_mz * 1e6)

                if error_ppm <= tolerance_ppm:
                    possible_modes.append((f"{mode} ({charge}+)" if charge > 1 else mode, error_ppm))

        # translated note
        possible_modes.sort(key=lambda x: x[1])

        return possible_modes

    def batch_calculate_mz(self, compounds_df: pd.DataFrame,
                           smiles_col: str = None,
                           inchi_col: str = None,
                           formula_col: str = None,
                           ion_mode: str = 'M+H') -> pd.DataFrame:
        """calculatem/z"""
        result_df = compounds_df.copy()

        if 'mz' not in result_df.columns:
            result_df['mz'] = np.nan

        if 'mz_source' not in result_df.columns:
            result_df['mz_source'] = ''

        for idx, row in result_df.iterrows():
            calculated = False

            # SMILEScalculate
            if smiles_col and smiles_col in row and pd.notna(row[smiles_col]):
                try:
                    mz = self.calculate_mz_from_smiles(str(row[smiles_col]), ion_mode)
                    result_df.at[idx, 'mz'] = mz
                    result_df.at[idx, 'mz_source'] = f'SMILES ({ion_mode})'
                    calculated = True
                except Exception as e:
                    pass

            # InChIcalculate
            if not calculated and inchi_col and inchi_col in row and pd.notna(row[inchi_col]):
                try:
                    mz = self.calculate_mz_from_inchi(str(row[inchi_col]), ion_mode)
                    result_df.at[idx, 'mz'] = mz
                    result_df.at[idx, 'mz_source'] = f'InChI ({ion_mode})'
                    calculated = True
                except Exception as e:
                    pass

            # molecular_formulacalculate
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
    """PPG retention index calculator"""

    def __init__(self):
        """Initialize the retention index calculator"""
        self.ppg_standards = None # PPG standard data
        self.regression_params = None # linear regression parameters
        self.calibration_curve = None # calibration curve

    def load_ppg_standards(self, excel_file: str) -> Tuple[bool, str, pd.DataFrame]:
        """Load PPG standard data from an Excel file"""
        try:
            file_path = Path(excel_file)

            if not file_path.exists():
                return False, "PPG standard file does not exist", None

            # Try to read the Excel file
            try:
                # Try all possible worksheets
                xls = pd.ExcelFile(file_path)
                sheet_names = xls.sheet_names

                # Try the first worksheet
                df = pd.read_excel(file_path, sheet_name=sheet_names[0])

                # Find required columns
                rt_col = self._find_column(df, ['RT', 'rt', 'retention_time', 'RetentionTime', 'retention_time'])
                n_col = self._find_column(df, ['n', 'N', 'degree_of_polymerization', 'degree_of_polymerizationn', 'degree', 'Degree'])
                name_col = self._find_column(df, ['name', 'Name', 'compound', 'Compound', 'PPG'])

                if rt_col is None or n_col is None:
                    return False, "Retention-time or degree-of-polymerization column was not found; ensure that the Excel file contains RT and n information", None

                # column
                rename_dict = {}
                if rt_col:
                    rename_dict[rt_col] = 'RT'
                if n_col:
                    rename_dict[n_col] = 'n'
                if name_col:
                    rename_dict[name_col] = 'name'

                if rename_dict:
                    df = df.rename(columns=rename_dict)

                # data
                if 'RT' in df.columns:
                    df['RT'] = pd.to_numeric(df['RT'], errors='coerce')
                if 'n' in df.columns:
                    df['n'] = pd.to_numeric(df['n'], errors='coerce')

                # translated note
                df = df.sort_values('n').reset_index(drop=True)

                self.ppg_standards = df

                return True, f"load {len(df)} PPGstandard", df

            except Exception as e:
                return False, f"PPGstandardfilefailed: {str(e)}", None

        except Exception as e:
            return False, f"loadPPGstandardfailed: {str(e)}", None

    def _find_column(self, df: pd.DataFrame, possible_names: List[str]) -> Optional[str]:
        """DataFramecolumn"""
        for name in possible_names:
            if name in df.columns:
                return name
        return None

    def build_calibration_curve(self, method: str = 'linear') -> Tuple[bool, str]:
        """PPGcalibration curve

        Parameters:
            method: method, 'linear' () 'interpolation' ()
        """
        if self.ppg_standards is None or len(self.ppg_standards) < 2:
            return False, "PPG standard data, 2standard"

        try:
            df = self.ppg_standards.copy()
            df = df.dropna(subset=['RT', 'n'])

            if len(df) < 2:
                return False, "PPG standard data"

            if method == 'linear':
                # method
                x = df['n'].values
                y = df['RT'].values

                # translated note
                slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

                # calculate
                y_pred = intercept + slope * x
                residuals = y - y_pred

                # calculateR²R²
                ss_res = np.sum(residuals ** 2)
                ss_tot = np.sum((y - np.mean(y)) ** 2)
                r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0
                n = len(x)
                p = 1 #
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

                # calibration curve
                def calibration_func(n_value):
                    return intercept + slope * n_value

                self.calibration_curve = calibration_func

                return True, f"calibration curve: R²={r_squared:.4f}, ={slope:.4f}, ={intercept:.4f}"

            elif method == 'interpolation':
                # method
                x = df['n'].values
                y = df['RT'].values

                # translated note
                self.calibration_curve = interp1d(x, y, kind='linear', fill_value='extrapolate')
                self.regression_params = {'method': 'interpolation', 'n_points': len(x)}

                return True, f"calibration curve: Use{len(x)}standard"

            else:
                return False, f"method: {method}"

        except Exception as e:
            return False, f"calibration curvefailed: {str(e)}"

    def calculate_ppg_index(self, rt_values: Union[float, List[float], pd.Series],
                            method: str = 'interpolation') -> Union[float, List[float], pd.Series]:
        """calculatePPG retention indices

        Parameters:
            rt_values: retention_time, , columnpandas Series
            method: calculatemethod, 'linear' () 'interpolation' ()

        Returns:
            PPG retention indices (n*100)
        """
        if self.ppg_standards is None or len(self.ppg_standards) < 2:
            raise ValueError("PPG standard dataload")

        # input
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
                    # Use
                    slope = self.regression_params['slope']
                    intercept = self.regression_params['intercept']

                    if abs(slope) < 1e-10: #
                        n = 0
                    else:
                        n = (rt - intercept) / slope

                else:
                    # Use
                    df = self.ppg_standards.sort_values('RT')

                    # translated note
                    if rt < df['RT'].min():
                        # RT, Use
                        rt_min = df['RT'].iloc[0]
                        rt_next = df['RT'].iloc[1]
                        n_min = df['n'].iloc[0]
                        n_next = df['n'].iloc[1]

                        if rt_next - rt_min != 0:
                            n = n_min + (rt - rt_min) / (rt_next - rt_min) * (n_next - n_min)
                        else:
                            n = n_min

                    elif rt > df['RT'].max():
                        # RT, Use
                        rt_max = df['RT'].iloc[-1]
                        rt_prev = df['RT'].iloc[-2]
                        n_max = df['n'].iloc[-1]
                        n_prev = df['n'].iloc[-2]

                        if rt_max - rt_prev != 0:
                            n = n_max + (rt - rt_max) / (rt_max - rt_prev) * (n_max - n_prev)
                        else:
                            n = n_max

                    else:
                        # translated note
                        for i in range(len(df) - 1):
                            rt_i = df['RT'].iloc[i]
                            rt_ip1 = df['RT'].iloc[i + 1]

                            if rt_i <= rt <= rt_ip1:
                                n_i = df['n'].iloc[i]
                                n_ip1 = df['n'].iloc[i + 1]

                                if rt_ip1 - rt_i != 0:
                                    # translated note
                                    n = n_i + (rt - rt_i) / (rt_ip1 - rt_i) * (n_ip1 - n_i)
                                else:
                                    n = n_i
                                break
                        else:
                            n = None

                if n is not None:
                    # 100retention_index ()
                    n_final = n * 100
                else:
                    n_final = None

                results.append(n_final)

            except Exception as e:
                print(f"calculateretention_index (RT={rt}): {str(e)}")
                results.append(None)

        # Returnsinput
        if is_scalar:
            return results[0] if results else None
        elif isinstance(rt_values, pd.Series):
            return pd.Series(results, index=rt_values.index, name='PPG_Index')
        else:
            return results

    def calculate_ppg_index_batch(self, rt_df: pd.DataFrame, rt_column: str = 'retention_time(RT)') -> pd.DataFrame:
        """calculatePPG retention indices

        Parameters:
            rt_df: retention_timeDataFrame
            rt_column: retention_timecolumn

        Returns:
            PPG retention indicescolumnDataFrame
        """
        if rt_column not in rt_df.columns:
            raise ValueError(f"DataFramecolumn: {rt_column}")

        result_df = rt_df.copy()
        result_df['PPG_Index'] = self.calculate_ppg_index(rt_df[rt_column])

        return result_df

    def save_calibration_report(self, output_dir: str, base_name: str = "PPG") -> Tuple[bool, str, List[str]]:
        """savePPG"""
        try:
            if self.ppg_standards is None:
                return False, "PPG standard datasave", []

            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)

            saved_files = []
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

            # savedata
            excel_file = output_path / f"{base_name}_{timestamp}.xlsx"

            with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
                # PPG standard data
                self.ppg_standards.to_excel(writer, sheet_name='PPGstandard', index=False)

                # calibration curveParameters
                if self.regression_params:
                    params_df = pd.DataFrame([self.regression_params])
                    params_df.to_excel(writer, sheet_name='Parameters', index=False)

                # translated note
                if self.ppg_standards is not None:
                    stats_data = {
                        '': ['standard', 'degree_of_polymerization', 'retention_time', 'retention_time', 'retention_time'],
                        '': [
                            len(self.ppg_standards),
                            f"{self.ppg_standards['n'].min():.1f} - {self.ppg_standards['n'].max():.1f}",
                            f"{self.ppg_standards['RT'].min():.3f} - {self.ppg_standards['RT'].max():.3f}",
                            f"{self.ppg_standards['RT'].mean():.3f}",
                            f"{self.ppg_standards['RT'].std():.3f}"
                        ]
                    }

                    if self.regression_params and self.regression_params['method'] == 'linear':
                        stats_data[''].extend(['R²', 'R²', '', '', ''])
                        stats_data[''].extend([
                            f"{self.regression_params['r_squared']:.4f}",
                            f"{self.regression_params['adj_r_squared']:.4f}",
                            f"{self.regression_params['slope']:.4f}",
                            f"{self.regression_params['intercept']:.4f}",
                            f"{self.regression_params['std_err']:.4f}"
                        ])

                    stats_df = pd.DataFrame(stats_data)
                    stats_df.to_excel(writer, sheet_name='', index=False)

            saved_files.append(str(excel_file))

            return True, f"save {excel_file.name}", saved_files

        except Exception as e:
            return False, f"savefailed: {str(e)}", []


class CompoundMatcher:
    """compound matching - Excelloadcompoundmatch"""

    def __init__(self):
        """match"""
        self.compounds_df = None
        self.match_results = []
        self.mz_calculator = MzCalculator() # m/zcalculate
        self.match_settings = {
            'ppm_tolerance': 10, # ppm
            'rt_window': 30, # retention_time ()
            'intensity_threshold': 1000, #
            'ion_mode': 'M+H', #
            'calculate_mz': True # calculatem/z
        }

    def load_compounds_from_excel(self, excel_file: str,
                                  calculate_mz: bool = True,
                                  ion_mode: str = 'M+H') -> Tuple[bool, str, pd.DataFrame]:
        """Excelfileloadcompound, calculatem/z"""
        try:
            file_path = Path(excel_file)

            if not file_path.exists():
                return False, "file does not exist", None

            # Try to read the Excel file
            try:
                # Try all possible worksheets
                xls = pd.ExcelFile(file_path)
                sheet_names = xls.sheet_names

                # Try the first worksheet
                df = pd.read_excel(file_path, sheet_name=sheet_names[0])

                # Find required columns
                mz_col = self._find_column(df, ['mz', 'M/Z', 'm/z', 'mass', 'Mass', 'MZ', ''])
                name_col = self._find_column(df, ['name', 'Name', 'compound', 'Compound', 'ID', 'name', 'compound_name'])
                formula_col = self._find_column(df, ['formula', 'Formula', 'molecular_formula'])
                rt_col = self._find_column(df, ['rt', 'RT', 'retention_time', 'retention_time'])
                smiles_col = self._find_column(df, ['smiles', 'SMILES', 'Smiles'])
                inchi_col = self._find_column(df, ['inchi', 'InChI', 'InChIKey'])

                # column
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

                # column
                if 'name' not in df.columns:
                    df['name'] = [f"compound_{i + 1}" for i in range(len(df))]

                # calculatem/zm/zcolumncalculate
                if calculate_mz and ('mz' not in df.columns or df['mz'].isna().all()):
                    # calculatem/z
                    df = self.mz_calculator.batch_calculate_mz(
                        df,
                        smiles_col='smiles' if 'smiles' in df.columns else None,
                        inchi_col='inchi' if 'inchi' in df.columns else None,
                        formula_col='formula' if 'formula' in df.columns else None,
                        ion_mode=ion_mode
                    )

                    # calculate
                    calculated_count = df['mz_source'].notna().sum()
                    if calculated_count > 0:
                        print(f"calculate {calculated_count} compoundm/z")

                # m/zcolumn
                if 'mz' in df.columns:
                    df['mz'] = pd.to_numeric(df['mz'], errors='coerce')
                    # m/zNaN ()
                    df_original_len = len(df)
                    df = df.dropna(subset=['mz'])
                    dropped_count = df_original_len - len(df)
                    if dropped_count > 0:
                        print(f" {dropped_count} m/zcompound")

                self.compounds_df = df

                return True, f"load {len(df)} compound", df

            except Exception as e:
                return False, f"Read Excel filefailed: {str(e)}", None

        except Exception as e:
            return False, f"loadcompoundfailed: {str(e)}", None

    def _find_column(self, df: pd.DataFrame, possible_names: List[str]) -> Optional[str]:
        """DataFramecolumn"""
        for name in possible_names:
            if name in df.columns:
                return name
        return None

    def set_match_settings(self, ppm_tolerance: float = 10, rt_window: float = 30,
                           intensity_threshold: float = 1000, ion_mode: str = 'M+H',
                           calculate_mz: bool = True):
        """matchParameters"""
        self.match_settings = {
            'ppm_tolerance': ppm_tolerance,
            'rt_window': rt_window,
            'intensity_threshold': intensity_threshold,
            'ion_mode': ion_mode,
            'calculate_mz': calculate_mz
        }

    def calculate_mz_for_compound(self, identifier: str, identifier_type: str = 'smiles',
                                  ion_mode: str = 'M+H') -> Tuple[bool, str, float]:
        """compoundcalculatem/z

        Parameters:
            identifier: SMILES, InChImolecular_formula
            identifier_type: 'smiles', 'inchi', 'formula'
            ion_mode:

        Returns:
            (, message, m/z)
        """
        try:
            if identifier_type.lower() == 'smiles':
                mz = self.mz_calculator.calculate_mz_from_smiles(identifier, ion_mode)
                return True, f"SMILEScalculatem/z", mz
            elif identifier_type.lower() == 'inchi':
                mz = self.mz_calculator.calculate_mz_from_inchi(identifier, ion_mode)
                return True, f"InChIcalculatem/z", mz
            elif identifier_type.lower() in ['formula', 'mf', 'molformula']:
                mz = self.mz_calculator.calculate_mz_from_formula(identifier, ion_mode)
                return True, f"molecular_formulacalculatem/z", mz
            else:
                return False, f": {identifier_type}", 0.0
        except Exception as e:
            return False, f"calculatem/zfailed: {str(e)}", 0.0

    def match_compounds(self, peaks_data: Union[List[Dict], pd.DataFrame, str]) -> Tuple[bool, str]:
        """matchcompoundpeakdata"""
        try:
            # compound data
            if self.compounds_df is None or len(self.compounds_df) == 0:
                return False, "loadcompound"

            # peakdataDataFrame
            if isinstance(peaks_data, str):
                # filepath
                peaks_df = pd.read_csv(peaks_data)
            elif isinstance(peaks_data, pd.DataFrame):
                peaks_df = peaks_data
            elif isinstance(peaks_data, list):
                peaks_df = pd.DataFrame(peaks_data)
            else:
                return False, "peakdata"

            # peakdatacolumn
            required_columns = ['(m/z)', 'retention_time(RT)', '']
            missing_cols = [col for col in required_columns if col not in peaks_df.columns]

            if missing_cols:
                # column
                column_mapping = {
                    '(m/z)': ['mz', 'M/Z', 'm/z'],
                    'retention_time(RT)': ['rt', 'RT', 'retention_time'],
                    '': ['intensity', 'Intensity', '']
                }

                for req_col, possible_names in column_mapping.items():
                    if req_col not in peaks_df.columns:
                        for name in possible_names:
                            if name in peaks_df.columns:
                                peaks_df = peaks_df.rename(columns={name: req_col})
                                break

            # translated note
            missing_cols = [col for col in required_columns if col not in peaks_df.columns]
            if missing_cols:
                return False, f"peakdatacolumn: {', '.join(missing_cols)}"

            # data
            peaks_df['(m/z)'] = pd.to_numeric(peaks_df['(m/z)'], errors='coerce')
            peaks_df['retention_time(RT)'] = pd.to_numeric(peaks_df['retention_time(RT)'], errors='coerce')
            peaks_df[''] = pd.to_numeric(peaks_df[''], errors='coerce')

            # peak
            peaks_df = peaks_df[peaks_df[''] >= self.match_settings['intensity_threshold']]

            self.match_results = []

            # compoundmatch
            for _, compound in self.compounds_df.iterrows():
                compound_mz = compound.get('mz')
                compound_name = compound.get('name', 'compound')
                formula = compound.get('formula', '')
                rt_reference = compound.get('rt_reference', None)
                smiles = compound.get('smiles', '')
                inchi = compound.get('inchi', '')
                mz_source = compound.get('mz_source', '')

                if pd.isna(compound_mz):
                    continue

                # calculateppm
                ppm_tolerance = self.match_settings['ppm_tolerance']
                mz_tolerance = compound_mz * ppm_tolerance / 1e6

                # matchpeak
                matching_peaks = peaks_df[
                    (peaks_df['(m/z)'] >= compound_mz - mz_tolerance) &
                    (peaks_df['(m/z)'] <= compound_mz + mz_tolerance)
                    ].copy()

                # retention_time,
                if rt_reference is not None and not pd.isna(rt_reference):
                    rt_window = self.match_settings['rt_window']
                    matching_peaks = matching_peaks[
                        (matching_peaks['retention_time(RT)'] >= rt_reference - rt_window) &
                        (matching_peaks['retention_time(RT)'] <= rt_reference + rt_window)
                        ]

                if len(matching_peaks) > 0:
                    # peak
                    matching_peaks['mz_difference'] = abs(matching_peaks['(m/z)'] - compound_mz)
                    matching_peaks['mz_difference_ppm'] = matching_peaks['mz_difference'] / compound_mz * 1e6

                    # m/z
                    matching_peaks = matching_peaks.sort_values('mz_difference')

                    for i, (_, peak) in enumerate(matching_peaks.iterrows()):
                        if i >= 3: # 3match
                            break

                        result = {
                            'compound_name': compound_name,
                            'molecular_formula': formula,
                            'm/z': compound_mz,
                            'm/z': mz_source,
                            'm/z': peak['(m/z)'],
                            'm/z(Da)': peak['mz_difference'],
                            'm/z(ppm)': peak['mz_difference_ppm'],
                            'retention_time(RT)': peak['retention_time(RT)'],
                            '': peak[''],
                            'match_status': 'match',
                            'match_rank': i + 1
                        }

                        # translated note
                        if smiles:
                            result['SMILES'] = smiles
                        if inchi:
                            result['InChI'] = inchi

                        # retention_time
                        if rt_reference is not None and not pd.isna(rt_reference):
                            result['RT'] = rt_reference
                            result['RT'] = abs(peak['retention_time(RT)'] - rt_reference)

                        self.match_results.append(result)
                else:
                    # match
                    result = {
                        'compound_name': compound_name,
                        'molecular_formula': formula,
                        'm/z': compound_mz,
                        'm/z': mz_source,
                        'm/z': None,
                        'm/z(Da)': None,
                        'm/z(ppm)': None,
                        'retention_time(RT)': None,
                        '': None,
                        'match_status': 'unmatched',
                        'match_rank': None
                    }

                    # translated note
                    if smiles:
                        result['SMILES'] = smiles
                    if inchi:
                        result['InChI'] = inchi

                    if rt_reference is not None and not pd.isna(rt_reference):
                        result['RT'] = rt_reference

                    self.match_results.append(result)

            if len(self.match_results) == 0:
                return False, "match"

            # match
            successful_matches = [r for r in self.match_results if r['match_status'] == 'match']

            return True, f"match: {len(self.compounds_df)} compound, {len(successful_matches)} match"

        except Exception as e:
            return False, f"match: {str(e)}"

    def add_ppg_index(self, ppg_calculator: PPGIndexCalculator) -> Tuple[bool, str]:
        """matching resultsPPG retention indices"""
        try:
            if not self.match_results:
                return False, "matching resultsPPGindex"

            if ppg_calculator.calibration_curve is None:
                return False, "PPGcalibration curve"

            # matching resultsDataFrame
            results_df = pd.DataFrame(self.match_results)

            # calculatePPG retention indices
            if 'retention_time(RT)' in results_df.columns:
                results_df['PPG retention indices'] = ppg_calculator.calculate_ppg_index(results_df['retention_time(RT)'])

            # matching results
            self.match_results = results_df.to_dict('records')

            # calculatesuccessful matchesPPGindex
            successful_matches = [r for r in self.match_results if
                                  r['match_status'] == 'match' and r.get('PPG retention indices') is not None]
            if successful_matches:
                ppg_values = [r['PPG retention indices'] for r in successful_matches]
                avg_ppg = sum(ppg_values) / len(ppg_values)
                min_ppg = min(ppg_values)
                max_ppg = max(ppg_values)

                return True, f"PPGindex: {len(successful_matches)}matchcompound, PPGindex: {min_ppg:.1f} - {max_ppg:.1f}, : {avg_ppg:.1f}"
            else:
                return True, "PPGindex, successful matchescompound"

        except Exception as e:
            return False, f"PPGindexfailed: {str(e)}"

    def save_match_results(self, output_dir: str, base_name: str = "compound matchingresults",
                           include_ppg: bool = True) -> Tuple[bool, str, List[str]]:
        """savematching resultsExcelfile"""
        try:
            if not self.match_results:
                return False, "matching resultssave", []

            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)

            saved_files = []

            # DataFrame
            results_df = pd.DataFrame(self.match_results)

            # saveExcel
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

            # resultsfile
            excel_file = output_path / f"{base_name}_{timestamp}.xlsx"

            with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
                # results
                results_df.to_excel(writer, sheet_name='all matching results', index=False)

                # successful matchesresults
                successful_matches = results_df[results_df['match_status'] == 'match'].copy()
                if len(successful_matches) > 0:
                    successful_matches = successful_matches.sort_values(['compound_name', 'match_rank'])
                    successful_matches.to_excel(writer, sheet_name='successful matches', index=False)

                # compound
                compound_summary = []
                for compound_name in results_df['compound_name'].unique():
                    compound_matches = results_df[results_df['compound_name'] == compound_name]
                    successful = compound_matches[compound_matches['match_status'] == 'match']

                    summary = {
                        'compound_name': compound_name,
                        'match': len(successful),
                        'm/z(ppm)': successful['m/z(ppm)'].min() if len(successful) > 0 else None,
                        '': successful[''].max() if len(successful) > 0 else None,
                        'RT': successful['retention_time(RT)'].mean() if len(successful) > 0 else None,
                        'match_status': 'match' if len(successful) > 0 else 'unmatched'
                    }

                    # m/z
                    if 'm/z' in compound_matches.columns:
                        sources = compound_matches['m/z'].unique()
                        summary['m/z'] = ', '.join([s for s in sources if pd.notna(s)])

                    # PPGindex
                    if include_ppg and 'PPG retention indices' in successful.columns:
                        ppg_values = successful['PPG retention indices'].dropna()
                        if len(ppg_values) > 0:
                            summary['PPGindex'] = ppg_values.mean()
                            summary['PPGindex'] = f"{ppg_values.min():.1f} - {ppg_values.max():.1f}"

                    compound_summary.append(summary)

                summary_df = pd.DataFrame(compound_summary)
                summary_df.to_excel(writer, sheet_name='compound', index=False)

                # translated note
                stats = {
                    '': ['compound', 'successful matches', 'unmatched', 'matchpeak',
                                 'm/z(ppm)', 'm/z(ppm)', 'm/z(ppm)'],
                    '': [
                        len(results_df['compound_name'].unique()),
                        len(results_df[results_df['match_status'] == 'match']['compound_name'].unique()),
                        len(results_df[results_df['match_status'] == 'unmatched']['compound_name'].unique()),
                        len(results_df[results_df['match_status'] == 'match']),
                        results_df['m/z(ppm)'].mean() if 'm/z(ppm)' in results_df.columns else 0,
                        results_df['m/z(ppm)'].min() if 'm/z(ppm)' in results_df.columns else 0,
                        results_df['m/z(ppm)'].max() if 'm/z(ppm)' in results_df.columns else 0
                    ]
                }

                # m/z
                if 'm/z' in results_df.columns:
                    mz_sources = results_df['m/z'].value_counts()
                    if len(mz_sources) > 0:
                        for source, count in mz_sources.items():
                            if pd.notna(source):
                                stats[''].append(f"m/z: {source}")
                                stats[''].append(count)

                # PPGindex
                if include_ppg and 'PPG retention indices' in results_df.columns:
                    ppg_values = results_df['PPG retention indices'].dropna()
                    if len(ppg_values) > 0:
                        stats[''].extend(['PPGindex', 'PPGindex', 'PPGindex', 'PPGindex'])
                        stats[''].extend([
                            ppg_values.mean(),
                            ppg_values.min(),
                            ppg_values.max(),
                            ppg_values.std()
                        ])

                stats_df = pd.DataFrame(stats)
                stats_df.to_excel(writer, sheet_name='', index=False)

            saved_files.append(str(excel_file))

            # Save as CSV (Use)
            csv_file = output_path / f"{base_name}_{timestamp}.csv"
            results_df.to_csv(csv_file, index=False, encoding='utf-8-sig')
            saved_files.append(str(csv_file))

            return True, f"resultssave {excel_file.name}", saved_files

        except Exception as e:
            return False, f"saveresultsfailed: {str(e)}", []


class MzMLPeakExtractor:
    """mzMLpeak - peakdatacompound matching"""

    def __init__(self, progress_callback=None, log_callback=None):
        """"""
        self.progress_callback = progress_callback
        self.log_callback = log_callback
        self.cancel_requested = False

    def log(self, message: str, level: str = "INFO"):
        """Record log messages"""
        timestamp = datetime.now().strftime("%H:%M:%S")
        log_message = f"[{timestamp}] {message}"

        if self.log_callback:
            self.log_callback(log_message, level)
        else:
            print(log_message)

    def update_progress(self, value: int, message: str = ""):
        """Update progress"""
        if self.progress_callback:
            self.progress_callback(value, message)

    def extract_peaks_with_pymzml(self, mzml_path: str, intensity_threshold: float = 1000) -> Tuple[
        Optional[List[Dict]], str]:
        """UsepymzMLpeakdata"""
        try:
            import pymzml

            self.log(f"UsepymzMLloadfile: {Path(mzml_path).name}")
            peaks = []

            run = pymzml.run.Reader(mzml_path)
            spectrum_count = 0

            for spectrum in run:
                if self.cancel_requested:
                    return None, "user cancelled the operation"

                if spectrum.ms_level == 1:
                    mz_array = spectrum.mz
                    intensity_array = spectrum.i

                    if len(mz_array) > 0:
                        for idx, (mz, intensity) in enumerate(zip(mz_array, intensity_array)):
                            if intensity > intensity_threshold:
                                peaks.append({
                                    'spectrum_index': spectrum_count,
                                    'retention_time(RT)': spectrum.scan_time[0] if hasattr(spectrum,
                                                                                     'scan_time') and spectrum.scan_time else spectrum_count * 0.5,
                                    '(m/z)': float(mz),
                                    '': float(intensity),
                                    'MS': 1,
                                    'ID': f"F{spectrum_count:05d}_{idx:04d}"
                                })

                        spectrum_count += 1

                # Update progress
                if spectrum_count % 100 == 0:
                    progress = min(100, spectrum_count // 10)
                    self.update_progress(progress, f" {spectrum_count}")

            if peaks:
                return peaks, f"pymzML: {len(peaks)}peak"
            else:
                return None, "pymzMLconditionpeak"

        except ImportError:
            return None, "pymzML"
        except Exception as e:
            return None, f"pymzMLfailed: {str(e)}"

    def extract_peaks_from_mzml(self, mzml_file: str, intensity_threshold: float = 1000) -> Tuple[
        bool, str, Optional[pd.DataFrame]]:
        """Extract peak data from mzML files"""
        try:
            self.log("peakdata...")
            self.update_progress(10, "peakdata...")

            # peak
            peaks, msg = self.extract_peaks_with_pymzml(mzml_file, intensity_threshold)

            if peaks is None:
                return False, msg, None

            # DataFrame
            df = pd.DataFrame(peaks)

            # columncompound matching
            column_mapping = {
                '(m/z)': '(m/z)',
                'retention_time(RT)': 'retention_time(RT)',
                '': ''
            }

            for old_col, new_col in column_mapping.items():
                if old_col in df.columns:
                    df = df.rename(columns={old_col: new_col})

            self.update_progress(100, "peak")
            return True, f" {len(df)} peak", df

        except Exception as e:
            return False, f"peakdatafailed: {str(e)}", None

    def extract_peaks_from_csv(self, csv_file: str) -> Tuple[bool, str, Optional[pd.DataFrame]]:
        """CSVfileloadpeakdata"""
        try:
            self.log("CSVfileloadpeakdata...")

            df = pd.read_csv(csv_file)

            # column
            required_columns = ['(m/z)', 'retention_time(RT)', '']

            # column
            column_mapping = {
                '(m/z)': ['mz', 'M/Z', 'm/z'],
                'retention_time(RT)': ['rt', 'RT', 'retention_time'],
                '': ['intensity', 'Intensity', '']
            }

            for req_col, possible_names in column_mapping.items():
                if req_col not in df.columns:
                    for name in possible_names:
                        if name in df.columns:
                            df = df.rename(columns={name: req_col})
                            break

            # translated note
            missing_cols = [col for col in required_columns if col not in df.columns]
            if missing_cols:
                return False, f"CSVfilecolumn: {', '.join(missing_cols)}", None

            return True, f"load {len(df)} peak", df

        except Exception as e:
            return False, f"loadCSVfilefailed: {str(e)}", None


class MzMLCompoundMatcherGUI:
    """mzMLpeakcompound matchingGUI"""

    def __init__(self, root):
        """GUI"""
        self.root = root
        self.root.title("mzMLpeakcompound matchingretention_indexcalculate")
        self.root.geometry("1100x950")

        # translated note
        try:
            self.root.iconbitmap(default='icon.ico')
        except:
            pass

        # translated note
        self.peak_extractor = None
        self.compound_matcher = None
        self.ppg_calculator = None
        self.processing_thread = None
        self.is_processing = False

        # UI
        self.setup_ui()

        # translated note
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)

    def setup_ui(self):
        """UI"""
        # translated note
        main_frame = tk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Canvas
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

        # translated note
        content_frame = tk.Frame(scrollable_frame, padx=20, pady=20)
        content_frame.pack(fill=tk.BOTH, expand=True)

        # translated note
        title_label = tk.Label(content_frame, text="mzMLpeak, compound matchingPPG retention indicescalculate",
                               font=("Arial", 18, "bold"))
        title_label.pack(pady=(0, 10))

        # Version
        version_label = tk.Label(content_frame, text="Version 9.0 - SMILES/InChIcalculatem/zmatch",
                                 font=("Arial", 10))
        version_label.pack(pady=(0, 20))

        # ==================== 1: peakdata ====================
        step1_frame = tk.LabelFrame(content_frame, text="1: peakdata", font=("Arial", 12, "bold"),
                                    padx=10, pady=10)
        step1_frame.pack(fill=tk.X, pady=(0, 20))

        # data
        self.data_source_var = tk.StringVar(value="mzml")
        ttk.Radiobutton(step1_frame, text="Extract peak data from mzML files",
                        variable=self.data_source_var, value="mzml",
                        command=self.on_data_source_change).grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)

        ttk.Radiobutton(step1_frame, text="CSVfileloadpeakdata",
                        variable=self.data_source_var, value="csv",
                        command=self.on_data_source_change).grid(row=0, column=1, sticky=tk.W, padx=5, pady=5)

        # mzMLfile
        self.mzml_frame = tk.Frame(step1_frame)
        self.mzml_frame.grid(row=1, column=0, columnspan=2, sticky=tk.W, padx=5, pady=5)

        ttk.Label(self.mzml_frame, text="mzMLfile:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5))

        self.mzml_file_var = tk.StringVar()
        mzml_file_entry = ttk.Entry(self.mzml_frame, textvariable=self.mzml_file_var, width=70)
        mzml_file_entry.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=(0, 5))

        ttk.Button(self.mzml_frame, text="...", command=self.browse_mzml_file).grid(row=0, column=2, padx=(0, 5))

        # CSVfile
        self.csv_frame = tk.Frame(step1_frame)
        self.csv_frame.grid(row=2, column=0, columnspan=2, sticky=tk.W, padx=5, pady=5)

        ttk.Label(self.csv_frame, text="CSVfile:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5))

        self.csv_file_var = tk.StringVar()
        csv_file_entry = ttk.Entry(self.csv_frame, textvariable=self.csv_file_var, width=70)
        csv_file_entry.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=(0, 5))

        ttk.Button(self.csv_frame, text="...", command=self.browse_csv_file).grid(row=0, column=2, padx=(0, 5))

        # translated note
        ttk.Label(step1_frame, text=":").grid(row=3, column=0, sticky=tk.W, padx=5, pady=5)

        self.intensity_threshold_var = tk.StringVar(value="1000")
        threshold_entry = ttk.Entry(step1_frame, textvariable=self.intensity_threshold_var, width=15)
        threshold_entry.grid(row=3, column=1, sticky=tk.W, padx=5, pady=5)

        # CSV
        self.csv_frame.grid_remove()

        # ==================== 2: compoundload ====================
        step2_frame = tk.LabelFrame(content_frame, text="2: loadcompound", font=("Arial", 12, "bold"),
                                    padx=10, pady=10)
        step2_frame.pack(fill=tk.X, pady=(0, 20))

        ttk.Label(step2_frame, text="compoundExcelfile:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)

        self.compound_file_var = tk.StringVar()
        compound_file_entry = ttk.Entry(step2_frame, textvariable=self.compound_file_var, width=70)
        compound_file_entry.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=(0, 5), pady=5)

        ttk.Button(step2_frame, text="...", command=self.browse_compound_file).grid(row=0, column=2, pady=5)

        # calculatem/z
        self.calculate_mz_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(step2_frame, text="calculatem/z (SMILES/InChI/MF)",
                        variable=self.calculate_mz_var).grid(row=1, column=0, sticky=tk.W, padx=5, pady=5)

        # translated note
        ttk.Label(step2_frame, text=":").grid(row=1, column=1, sticky=tk.W, padx=5, pady=5)

        self.ion_mode_var = tk.StringVar(value="M+H")
        ion_mode_combo = ttk.Combobox(step2_frame, textvariable=self.ion_mode_var,
                                      values=["M+H", "M-H", "M+Na", "M+K", "M+NH4", "M+CH3COO",
                                              "M+2H", "M+2Na", "M+FA-H"], width=15)
        ion_mode_combo.grid(row=1, column=2, sticky=tk.W, padx=5, pady=5)

        # translated note
        ttk.Button(step2_frame, text="compound", command=self.preview_compounds).grid(row=2, column=0, sticky=tk.W,
                                                                                        padx=5, pady=5)

        # load
        self.compound_status_var = tk.StringVar(value="loadcompound")
        ttk.Label(step2_frame, textvariable=self.compound_status_var).grid(row=2, column=1, columnspan=2, sticky=tk.W,
                                                                           padx=5, pady=5)

        # ==================== 3: PPG retention indicescalculate ====================
        step3_frame = tk.LabelFrame(content_frame, text="3: PPG retention indicescalculate", font=("Arial", 12, "bold"),
                                    padx=10, pady=10)
        step3_frame.pack(fill=tk.X, pady=(0, 20))

        # PPGcalculate
        self.enable_ppg_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(step3_frame, text="PPG retention indicescalculate",
                        variable=self.enable_ppg_var,
                        command=self.on_ppg_enable_change).grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)

        # PPGstandardfile
        self.ppg_frame = tk.Frame(step3_frame)
        self.ppg_frame.grid(row=1, column=0, columnspan=3, sticky=tk.W, padx=5, pady=5)

        ttk.Label(self.ppg_frame, text="PPGstandardExcelfile:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5))

        self.ppg_file_var = tk.StringVar()
        ppg_file_entry = ttk.Entry(self.ppg_frame, textvariable=self.ppg_file_var, width=60)
        ppg_file_entry.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=(0, 5))

        ttk.Button(self.ppg_frame, text="...", command=self.browse_ppg_file).grid(row=0, column=2, padx=(0, 5))

        # PPGmethod
        ttk.Label(self.ppg_frame, text="method:").grid(row=1, column=0, sticky=tk.W, padx=(0, 5), pady=(10, 5))

        self.ppg_method_var = tk.StringVar(value="interpolation")
        ttk.Combobox(self.ppg_frame, textvariable=self.ppg_method_var,
                     values=["interpolation", "linear"], width=15, state="readonly").grid(row=1, column=1, sticky=tk.W,
                                                                                          padx=(0, 5), pady=(10, 5))

        # PPG
        ttk.Button(self.ppg_frame, text="PPGstandard", command=self.preview_ppg_standards).grid(row=1, column=2,
                                                                                                  pady=(10, 5))

        # PPG
        self.ppg_frame.grid_remove()

        # ==================== 4: matchParameters ====================
        step4_frame = tk.LabelFrame(content_frame, text="4: matchParameters", font=("Arial", 12, "bold"),
                                    padx=10, pady=10)
        step4_frame.pack(fill=tk.X, pady=(0, 20))

        # m/z (ppm)
        ttk.Label(step4_frame, text="m/z (ppm):").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)

        self.ppm_tolerance_var = tk.StringVar(value="10")
        ppm_entry = ttk.Entry(step4_frame, textvariable=self.ppm_tolerance_var, width=15)
        ppm_entry.grid(row=0, column=1, sticky=tk.W, padx=5, pady=5)

        ttk.Label(step4_frame, text="ppm").grid(row=0, column=2, sticky=tk.W, padx=(0, 20), pady=5)

        # RT ()
        ttk.Label(step4_frame, text="RT ():").grid(row=0, column=3, sticky=tk.W, padx=5, pady=5)

        self.rt_window_var = tk.StringVar(value="30")
        rt_entry = ttk.Entry(step4_frame, textvariable=self.rt_window_var, width=15)
        rt_entry.grid(row=0, column=4, sticky=tk.W, padx=5, pady=5)

        ttk.Label(step4_frame, text="").grid(row=0, column=5, sticky=tk.W, padx=(0, 5), pady=5)

        # outputdirectory
        ttk.Label(step4_frame, text="outputdirectory:").grid(row=1, column=0, sticky=tk.W, padx=5, pady=(10, 5))

        self.output_dir_var = tk.StringVar()
        output_dir_entry = ttk.Entry(step4_frame, textvariable=self.output_dir_var, width=70)
        output_dir_entry.grid(row=1, column=1, columnspan=4, sticky=(tk.W, tk.E), padx=(0, 5), pady=(10, 5))

        ttk.Button(step4_frame, text="...", command=self.browse_output_dir).grid(row=1, column=5, pady=(10, 5))

        # ==================== 5: compoundm/zcalculate ====================
        step5_frame = tk.LabelFrame(content_frame, text="5: compoundm/zcalculate ()", font=("Arial", 12, "bold"),
                                    padx=10, pady=10)
        step5_frame.pack(fill=tk.X, pady=(0, 20))

        # input
        ttk.Label(step5_frame, text="input:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)

        self.single_input_type_var = tk.StringVar(value="smiles")
        ttk.Radiobutton(step5_frame, text="SMILES", variable=self.single_input_type_var,
                        value="smiles").grid(row=0, column=1, sticky=tk.W, padx=5, pady=5)
        ttk.Radiobutton(step5_frame, text="InChI", variable=self.single_input_type_var,
                        value="inchi").grid(row=0, column=2, sticky=tk.W, padx=5, pady=5)
        ttk.Radiobutton(step5_frame, text="molecular_formula", variable=self.single_input_type_var,
                        value="formula").grid(row=0, column=3, sticky=tk.W, padx=5, pady=5)

        # input
        ttk.Label(step5_frame, text="input:").grid(row=1, column=0, sticky=tk.W, padx=5, pady=5)

        self.single_input_var = tk.StringVar()
        single_input_entry = ttk.Entry(step5_frame, textvariable=self.single_input_var, width=70)
        single_input_entry.grid(row=1, column=1, columnspan=3, sticky=(tk.W, tk.E), padx=(0, 5), pady=5)

        # calculate
        ttk.Button(step5_frame, text="calculatem/z", command=self.calculate_single_mz).grid(row=1, column=4, padx=5, pady=5)

        # results
        self.single_result_var = tk.StringVar(value="inputcompoundcalculatem/z")
        ttk.Label(step5_frame, textvariable=self.single_result_var).grid(row=2, column=0, columnspan=5,
                                                                         sticky=tk.W, padx=5, pady=5)

        # ==================== 6: ====================
        step6_frame = tk.LabelFrame(content_frame, text="6: ", font=("Arial", 12, "bold"),
                                    padx=10, pady=10)
        step6_frame.pack(fill=tk.X, pady=(0, 20))

        self.progress_var = tk.StringVar(value="")
        ttk.Label(step6_frame, textvariable=self.progress_var).pack(anchor=tk.W, pady=(0, 5))

        self.progress_bar = ttk.Progressbar(step6_frame, mode='determinate', length=950)
        self.progress_bar.pack(fill=tk.X, pady=(0, 5))

        self.progress_percent = tk.StringVar(value="0%")
        ttk.Label(step6_frame, textvariable=self.progress_percent).pack(anchor=tk.E)

        # ==================== 7: ====================
        step7_frame = tk.LabelFrame(content_frame, text="7: ", font=("Arial", 12, "bold"),
                                    padx=10, pady=10)
        step7_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 20))

        # translated note
        self.log_text = scrolledtext.ScrolledText(step7_frame, width=130, height=15,
                                                  wrap=tk.WORD, font=("Consolas", 10))
        self.log_text.pack(fill=tk.BOTH, expand=True)

        # translated note
        self.log_text.tag_config("INFO", foreground="black")
        self.log_text.tag_config("SUCCESS", foreground="green")
        self.log_text.tag_config("WARNING", foreground="orange")
        self.log_text.tag_config("ERROR", foreground="red")

        # ==================== ====================
        button_frame = tk.Frame(content_frame)
        button_frame.pack(fill=tk.X, pady=(0, 10))

        self.process_button = ttk.Button(button_frame, text="",
                                         command=self.start_processing, width=20)
        self.process_button.pack(side=tk.LEFT, padx=5)

        self.cancel_button = ttk.Button(button_frame, text="",
                                        command=self.cancel_processing, width=20,
                                        state=tk.DISABLED)
        self.cancel_button.pack(side=tk.LEFT, padx=5)

        ttk.Button(button_frame, text="", command=self.clear_log).pack(side=tk.LEFT, padx=5)

        ttk.Button(button_frame, text="outputdirectory", command=self.open_output_dir).pack(side=tk.LEFT, padx=5)

        ttk.Button(button_frame, text="", command=self.on_closing).pack(side=tk.LEFT, padx=5)

        # translated note
        self.status_var = tk.StringVar(value="")
        status_bar = ttk.Label(content_frame, textvariable=self.status_var, relief=tk.SUNKEN)
        status_bar.pack(fill=tk.X, pady=(0, 10))

        # translated note
        self.on_data_source_change()

    def on_data_source_change(self):
        """data"""
        source = self.data_source_var.get()

        if source == "mzml":
            self.mzml_frame.grid()
            self.csv_frame.grid_remove()
        else:
            self.mzml_frame.grid_remove()
            self.csv_frame.grid()

    def on_ppg_enable_change(self):
        """PPGcalculate"""
        if self.enable_ppg_var.get():
            self.ppg_frame.grid()
        else:
            self.ppg_frame.grid_remove()

    def browse_mzml_file(self):
        """mzMLfile"""
        file_types = [("mzMLfile", "*.mzML *.mzML.gz"), ("file", "*.*")]
        file_path = filedialog.askopenfilename(title="mzMLfile", filetypes=file_types)

        if file_path:
            self.mzml_file_var.set(file_path)

    def browse_csv_file(self):
        """CSVfile"""
        file_types = [("CSVfile", "*.csv"), ("file", "*.*")]
        file_path = filedialog.askopenfilename(title="CSVfile", filetypes=file_types)

        if file_path:
            self.csv_file_var.set(file_path)

    def browse_compound_file(self):
        """compoundExcelfile"""
        file_types = [("Excelfile", "*.xlsx *.xls"), ("file", "*.*")]
        file_path = filedialog.askopenfilename(title="compoundExcelfile", filetypes=file_types)

        if file_path:
            self.compound_file_var.set(file_path)
            # translated note
            self.preview_compounds()

    def browse_ppg_file(self):
        """PPGstandardExcelfile"""
        file_types = [("Excelfile", "*.xlsx *.xls"), ("file", "*.*")]
        file_path = filedialog.askopenfilename(title="PPGstandardExcelfile", filetypes=file_types)

        if file_path:
            self.ppg_file_var.set(file_path)
            # translated note
            self.preview_ppg_standards()

    def browse_output_dir(self):
        """outputdirectory"""
        dir_path = filedialog.askdirectory(title="outputdirectory")

        if dir_path:
            self.output_dir_var.set(dir_path)

    def preview_compounds(self):
        """compound"""
        compound_file = self.compound_file_var.get()

        if not compound_file:
            messagebox.showwarning("Warning", "compoundExcelfile")
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
                # translated note
                preview_window = tk.Toplevel(self.root)
                preview_window.title("compound")
                preview_window.geometry("900x600")

                # Treeviewdata
                tree_frame = tk.Frame(preview_window)
                tree_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

                # translated note
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

                # column
                tree["columns"] = list(df.columns)
                tree["show"] = "headings"

                # column
                for col in df.columns:
                    tree.heading(col, text=col)
                    tree.column(col, width=100)

                # data
                for _, row in df.iterrows():
                    tree.insert("", tk.END, values=list(row))

                # translated note
                info_frame = tk.Frame(preview_window)
                info_frame.pack(fill=tk.X, padx=10, pady=(0, 10))

                tk.Label(info_frame, text=f"load {len(df)} compound").pack(side=tk.LEFT)

                if 'mz' in df.columns:
                    mz_range = f"m/z: {df['mz'].min():.4f} - {df['mz'].max():.4f}"
                    tk.Label(info_frame, text=mz_range).pack(side=tk.LEFT, padx=20)

                if 'mz_source' in df.columns:
                    calculated_count = df['mz_source'].notna().sum()
                    tk.Label(info_frame, text=f"calculatem/z: {calculated_count}").pack(side=tk.LEFT, padx=20)

                # translated note
                tk.Button(preview_window, text="", command=preview_window.destroy).pack(pady=10)

                self.compound_status_var.set(f"load {len(df)} compound")
                self.log_message(f"✓ load {len(df)} compound", "SUCCESS")
            else:
                messagebox.showerror("", msg)
                self.compound_status_var.set("loadfailed")
                self.log_message(f"✗ loadcompoundfailed: {msg}", "ERROR")

        except Exception as e:
            messagebox.showerror("", f"failed: {str(e)}")

    def preview_ppg_standards(self):
        """PPGstandard"""
        ppg_file = self.ppg_file_var.get()

        if not ppg_file:
            messagebox.showwarning("Warning", "PPGstandardExcelfile")
            return

        try:
            calculator = PPGIndexCalculator()
            success, msg, df = calculator.load_ppg_standards(ppg_file)

            if success:
                # translated note
                preview_window = tk.Toplevel(self.root)
                preview_window.title("PPGstandard")
                preview_window.geometry("800x600")

                # Treeviewdata
                tree_frame = tk.Frame(preview_window)
                tree_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

                # translated note
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

                # column
                tree["columns"] = list(df.columns)
                tree["show"] = "headings"

                # column
                for col in df.columns:
                    tree.heading(col, text=col)
                    tree.column(col, width=100)

                # data
                for _, row in df.iterrows():
                    tree.insert("", tk.END, values=list(row))

                # translated note
                info_frame = tk.Frame(preview_window)
                info_frame.pack(fill=tk.X, padx=10, pady=(0, 10))

                tk.Label(info_frame, text=f"load {len(df)} PPGstandard").pack(side=tk.LEFT)

                if 'n' in df.columns and 'RT' in df.columns:
                    n_range = f"degree_of_polymerization: {df['n'].min():.1f} - {df['n'].max():.1f}"
                    rt_range = f"retention_time: {df['RT'].min():.3f} - {df['RT'].max():.3f}"
                    tk.Label(info_frame, text=n_range).pack(side=tk.LEFT, padx=20)
                    tk.Label(info_frame, text=rt_range).pack(side=tk.LEFT, padx=20)

                # translated note
                tk.Button(preview_window, text="", command=preview_window.destroy).pack(pady=10)

                self.log_message(f"✓ load {len(df)} PPGstandard", "SUCCESS")
            else:
                messagebox.showerror("", msg)
                self.log_message(f"✗ loadPPGstandardfailed: {msg}", "ERROR")

        except Exception as e:
            messagebox.showerror("", f"failed: {str(e)}")

    def calculate_single_mz(self):
        """calculatecompoundm/z"""
        input_text = self.single_input_var.get().strip()
        input_type = self.single_input_type_var.get()
        ion_mode = self.ion_mode_var.get()

        if not input_text:
            messagebox.showwarning("Warning", "inputcompound")
            return

        if not RDKIT_AVAILABLE:
            messagebox.showerror("", "RDKit, calculatem/z")
            self.log_message("✗ RDKit, calculatem/z", "ERROR")
            return

        try:
            matcher = CompoundMatcher()
            success, msg, mz = matcher.calculate_mz_for_compound(
                input_text, input_type, ion_mode
            )

            if success:
                result_text = f"calculate: {msg}\nm/z = {mz:.6f}"
                self.single_result_var.set(result_text)
                self.log_message(f"✓ calculate: {input_type} → {mz:.6f} ({ion_mode})", "SUCCESS")
            else:
                self.single_result_var.set(f"calculatefailed: {msg}")
                self.log_message(f"✗ calculatefailed: {msg}", "ERROR")

        except Exception as e:
            error_msg = f"calculate: {str(e)}"
            self.single_result_var.set(error_msg)
            self.log_message(f"✗ {error_msg}", "ERROR")

    def log_message(self, message: str, level: str = "INFO"):
        """message"""
        self.log_text.insert(tk.END, message + "\n", level)
        self.log_text.see(tk.END)
        self.root.update()

    def update_progress(self, value: int, message: str = ""):
        """Update progress"""
        self.progress_bar['value'] = value
        if message:
            self.progress_var.set(message)
        self.progress_percent.set(f"{value}%")
        self.root.update()

    def update_status(self, message: str):
        """"""
        self.status_var.set(message)
        self.root.update()

    def clear_log(self):
        """"""
        self.log_text.delete(1.0, tk.END)

    def validate_inputs(self) -> Tuple[bool, str]:
        """inputParameters"""
        # data
        source = self.data_source_var.get()

        if source == "mzml":
            mzml_file = self.mzml_file_var.get().strip()
            if not mzml_file:
                return False, "mzMLfile"
            if not os.path.exists(mzml_file):
                return False, f"mzMLfile does not exist: {mzml_file}"
        else:
            csv_file = self.csv_file_var.get().strip()
            if not csv_file:
                return False, "CSVfile"
            if not os.path.exists(csv_file):
                return False, f"CSVfile does not exist: {csv_file}"

        # compoundfile
        compound_file = self.compound_file_var.get().strip()
        if not compound_file:
            return False, "compoundExcelfile"
        if not os.path.exists(compound_file):
            return False, f"compoundfile does not exist: {compound_file}"

        # PPG
        if self.enable_ppg_var.get():
            ppg_file = self.ppg_file_var.get().strip()
            if not ppg_file:
                return False, "PPGstandardExcelfile"
            if not os.path.exists(ppg_file):
                return False, f"PPG standard file does not exist: {ppg_file}"

        # outputdirectory
        output_dir = self.output_dir_var.get().strip()
        if not output_dir:
            return False, "outputdirectory"

        # Parameters
        try:
            ppm_tolerance = float(self.ppm_tolerance_var.get())
            if ppm_tolerance <= 0:
                return False, "m/z0"
        except ValueError:
            return False, "m/z"

        try:
            rt_window = float(self.rt_window_var.get())
            if rt_window < 0:
                return False, "RT"
        except ValueError:
            return False, "RT"

        try:
            intensity_threshold = float(self.intensity_threshold_var.get())
            if intensity_threshold < 0:
                return False, ""
        except ValueError:
            return False, ""

        return True, ""

    def start_processing(self):
        """"""
        # input
        is_valid, message = self.validate_inputs()
        if not is_valid:
            messagebox.showerror("input", message)
            return

        # ,
        self.process_button.config(state=tk.DISABLED)
        self.cancel_button.config(state=tk.NORMAL)

        # translated note
        self.progress_bar['value'] = 0
        self.progress_var.set("...")
        self.progress_percent.set("0%")
        self.clear_log()

        # Parameters
        source = self.data_source_var.get()
        compound_file = self.compound_file_var.get()
        output_dir = self.output_dir_var.get()
        enable_ppg = self.enable_ppg_var.get()
        calculate_mz = self.calculate_mz_var.get()
        ion_mode = self.ion_mode_var.get()

        ppm_tolerance = float(self.ppm_tolerance_var.get())
        rt_window = float(self.rt_window_var.get())
        intensity_threshold = float(self.intensity_threshold_var.get())

        # translated note
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

        # translated note
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

        # translated note
        self.root.after(100, self.check_processing_status)

    def process_mzml_match(self, mzml_file: str, compound_file: str, output_dir: str,
                           enable_ppg: bool, calculate_mz: bool,
                           ppm_tolerance: float, rt_window: float,
                           intensity_threshold: float, ion_mode: str,
                           ppg_file: str, ppg_method: str):
        """mzMLfilecompound matching ()"""
        try:
            self.log_message("1: Extract peak data from mzML files...")
            self.update_progress(20, "peakdata...")

            # peakdata
            success, msg, peaks_df = self.peak_extractor.extract_peaks_from_mzml(
                mzml_file, intensity_threshold
            )

            if not success:
                self.root.after(0, lambda: messagebox.showerror("", msg))
                self.processing_finished()
                return

            self.log_message(f"✓ {len(peaks_df)} peak", "SUCCESS")

            # compound matching
            self.process_compound_matching(peaks_df, compound_file, output_dir,
                                           enable_ppg, calculate_mz,
                                           ppm_tolerance, rt_window,
                                           mzml_file, ion_mode,
                                           ppg_file, ppg_method)

        except Exception as e:
            error_msg = f": {str(e)}"
            self.root.after(0, lambda: self.log_message(error_msg, "ERROR"))
            self.root.after(0, self.processing_finished)

    def process_csv_match(self, csv_file: str, compound_file: str, output_dir: str,
                          enable_ppg: bool, calculate_mz: bool,
                          ppm_tolerance: float, rt_window: float,
                          intensity_threshold: float, ion_mode: str,
                          ppg_file: str, ppg_method: str):
        """CSVfilecompound matching ()"""
        try:
            self.log_message("1: CSVfileloadpeakdata...")
            self.update_progress(20, "loadpeakdata...")

            # loadpeakdata
            success, msg, peaks_df = self.peak_extractor.extract_peaks_from_csv(csv_file)

            if not success:
                self.root.after(0, lambda: messagebox.showerror("", msg))
                self.processing_finished()
                return

            self.log_message(f"✓ load {len(peaks_df)} peak", "SUCCESS")

            # compound matching
            self.process_compound_matching(peaks_df, compound_file, output_dir,
                                           enable_ppg, calculate_mz,
                                           ppm_tolerance, rt_window,
                                           csv_file, ion_mode,
                                           ppg_file, ppg_method)

        except Exception as e:
            error_msg = f": {str(e)}"
            self.root.after(0, lambda: self.log_message(error_msg, "ERROR"))
            self.root.after(0, self.processing_finished)

    def process_compound_matching(self, peaks_df: pd.DataFrame, compound_file: str,
                                  output_dir: str, enable_ppg: bool, calculate_mz: bool,
                                  ppm_tolerance: float, rt_window: float,
                                  source_file: str, ion_mode: str,
                                  ppg_file: str, ppg_method: str):
        """compound matching"""
        try:
            self.log_message("2: loadcompound...")
            self.update_progress(40, "loadcompound...")

            # matchParameters
            self.compound_matcher.set_match_settings(
                ppm_tolerance=ppm_tolerance,
                rt_window=rt_window,
                intensity_threshold=1000, #
                ion_mode=ion_mode,
                calculate_mz=calculate_mz
            )

            # loadcompound
            success, msg, _ = self.compound_matcher.load_compounds_from_excel(
                compound_file,
                calculate_mz=calculate_mz,
                ion_mode=ion_mode
            )
            if not success:
                self.root.after(0, lambda: messagebox.showerror("", msg))
                self.processing_finished()
                return

            self.log_message(f"✓ loadcompound", "SUCCESS")

            # , loadPPGstandardcalibration curve
            if enable_ppg and self.ppg_calculator and ppg_file:
                self.log_message("2.5: loadPPGstandardcalibration curve...")
                self.update_progress(45, "PPGcalibration curve...")

                success, msg, _ = self.ppg_calculator.load_ppg_standards(ppg_file)
                if not success:
                    self.root.after(0, lambda: messagebox.showerror("PPG", msg))
                    self.processing_finished()
                    return

                self.log_message(f"✓ loadPPGstandard", "SUCCESS")

                success, msg = self.ppg_calculator.build_calibration_curve(method=ppg_method)
                if not success:
                    self.root.after(0, lambda: messagebox.showerror("PPG", msg))
                    self.processing_finished()
                    return

                self.log_message(f"✓ {msg}", "SUCCESS")

            self.log_message("3: matchcompound...")
            self.update_progress(60, "matchcompound...")

            # match
            success, msg = self.compound_matcher.match_compounds(peaks_df)
            if not success:
                self.root.after(0, lambda: messagebox.showwarning("Warning", msg))

            self.log_message(f"✓ compound matching", "SUCCESS")

            # , PPG retention indices
            if enable_ppg and self.ppg_calculator:
                self.log_message("3.5: calculatePPG retention indices...")
                self.update_progress(70, "calculatePPG retention indices...")

                success, msg = self.compound_matcher.add_ppg_index(self.ppg_calculator)
                if not success:
                    self.log_message(f"✗ {msg}", "WARNING")
                else:
                    self.log_message(f"✓ {msg}", "SUCCESS")

                # savePPG
                success, msg, _ = self.ppg_calculator.save_calibration_report(output_dir)
                if success:
                    self.log_message(f"✓ PPGsave", "SUCCESS")

            self.log_message("4: savematching results...")
            self.update_progress(80, "savematching results...")

            # saveresults
            base_name = Path(source_file).stem
            success, msg, saved_files = self.compound_matcher.save_match_results(
                output_dir, f"{base_name}_compound matching", include_ppg=enable_ppg
            )

            if success:
                self.log_message(f"✓ resultssave", "SUCCESS")
                self.log_message(f" file:", "INFO")
                for file_path in saved_files:
                    self.log_message(f"    • {Path(file_path).name}", "INFO")

                # translated note
                self.root.after(0, lambda: messagebox.showinfo(
                    "",
                    f"compound matching！\nresultssave: {output_dir}"
                ))
            else:
                self.log_message(f"✗ saveresultsfailed: {msg}", "ERROR")
                self.root.after(0, lambda: messagebox.showerror("", f"saveresultsfailed: {msg}"))

            self.update_progress(100, "")
            self.processing_finished()

        except Exception as e:
            error_msg = f"compound matching: {str(e)}"
            self.root.after(0, lambda: self.log_message(error_msg, "ERROR"))
            self.root.after(0, self.processing_finished)

    def check_processing_status(self):
        """"""
        if self.is_processing and self.processing_thread.is_alive():
            # translated note
            self.root.after(100, self.check_processing_status)
        elif self.is_processing:
            # , results
            self.root.after(100, self.check_processing_status)

    def processing_finished(self):
        """"""
        self.is_processing = False
        self.peak_extractor = None
        self.compound_matcher = None
        self.ppg_calculator = None

        # ,
        self.process_button.config(state=tk.NORMAL)
        self.cancel_button.config(state=tk.DISABLED)

        # translated note
        self.update_status("")

    def cancel_processing(self):
        """"""
        if self.peak_extractor and self.is_processing:
            self.peak_extractor.cancel_requested = True
            self.log_message("...", "WARNING")
            self.update_status("...")

    def open_output_dir(self):
        """outputdirectory"""
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
                messagebox.showerror("", f"directory: {str(e)}")
        else:
            messagebox.showwarning("Warning", "outputdirectory")

    def on_closing(self):
        """"""
        if self.is_processing:
            if messagebox.askyesno("", ", ？"):
                if self.peak_extractor:
                    self.peak_extractor.cancel_requested = True
                self.root.destroy()
        else:
            self.root.destroy()


def check_dependencies():
    """"""
    print("=" * 70)
    print("...")
    print("=" * 70)

    dependencies = {
        'pandas': 'data',
        'numpy': 'calculate',
        'scipy': 'calculate (PPGindexcalculate)',
        'openpyxl': 'Excelfile',
        'pymzml': 'mzMLfile (peak)',
        'rdkit': ' (SMILES/InChIcalculatem/z)'
    }

    missing = []
    warnings = []

    for lib, desc in dependencies.items():
        try:
            __import__(lib)
            print(f"✓ {lib}: {desc}")
        except ImportError:
            if lib == 'rdkit':
                print(f"⚠ {lib}: {desc} - , m/zcalculate")
                warnings.append(lib)
            else:
                print(f"✗ {lib}: {desc} - ")
                missing.append(lib)

    if missing:
        print(f"\n:")
        for lib in missing:
            print(f"  pip install {lib}")

    if warnings:
        print(f"\n, :")
        for lib in warnings:
            print(f"  conda install -c conda-forge {lib}")

    print("\n" + "=" * 70)
    return len(missing) == 0


def main():
    """"""

    # GUI
    if not GUI_AVAILABLE:
        print(": tkinter, GUI")
        print("tkinter:")
        print("  Ubuntu/Debian: sudo apt-get install python3-tk")
        print(" Windows/macOS: ")
        return

    # translated note
    if not check_dependencies():
        response = input("\n, ? (y/n): ")
        if response.lower() != 'y':
            return

    # translated note
    root = tk.Tk()

    # translated note
    window_width = 1100
    window_height = 950
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()
    x = (screen_width - window_width) // 2
    y = (screen_height - window_height) // 2
    root.geometry(f"{window_width}x{window_height}+{x}+{y}")

    # GUI
    app = MzMLCompoundMatcherGUI(root)

    # translated note
    root.mainloop()


if __name__ == "__main__":
    print("mzMLpeak, compound matchingPPG retention indicescalculate")
    print("Version: 9.0")
    print(":")
    print(" 1. SMILES/InChI/molecular_formulacalculatecompoundm/z")
    print(" 2. Extract peak data from mzML files")
    print(" 3. CSVfileloadpeakdata")
    print(" 4. Excelfileloadcompound")
    print(" 5. Match peak data by compound m/z")
    print(" 6. Calculate compound PPG retention indices")
    print(" 7. matching resultsretention_indexExcelfile")
    print("=" * 70)

    main()