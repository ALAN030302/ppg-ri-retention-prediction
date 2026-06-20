"""
SDF CSV/Excel -
inputfilepath
"""

import os
import sys
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import PandasTools, Descriptors, AllChem
from rdkit.Chem import inchi
from pathlib import Path
import time
from datetime import datetime
import traceback
import re
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

# translated note
class CustomProgressBar:
    def __init__(self, total, desc=""):
        self.total = total
        self.desc = desc
        self.current = 0
        self.start_time = time.time()
        self.bar_length = 40

    def update(self, n=1):
        self.current += n
        elapsed = time.time() - self.start_time
        if self.current > 0:
            avg_time = elapsed / self.current
            eta = avg_time * (self.total - self.current)
        else:
            eta = 0

        percent = self.current / self.total
        filled_length = int(self.bar_length * percent)
        bar = '█' * filled_length + '░' * (self.bar_length - filled_length)

        # translated note
        sys.stdout.write(f'\r{self.desc}: |{bar}| {self.current}/{self.total} '
                        f'[{percent:.1%}] : {elapsed:.1f}s ETA: {eta:.1f}s')
        sys.stdout.flush()

    def close(self):
        elapsed = time.time() - self.start_time
        sys.stdout.write(f'\r{self.desc}: ! : {elapsed:.1f}s\n')
        sys.stdout.flush()

class InteractiveSDFConverter:
    def __init__(self):
        """"""
        self.version = "2.0.0"
        self.author = "SDF Converter Tool"
        self.supported_formats = ['.sdf', '.sd', '.mol']
        self.output_dir = Path.cwd() / "SDF_Output"
        self.log_file = self.output_dir / "conversion_log.txt"

        # outputdirectory
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def print_banner(self):
        """"""
        print("\n" + "="*60)
        print(" SDF file v2.0")
        print("="*60)
        print(f"Version: {self.version}")
        print(f": {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print("-"*60)
        print(":")
        print(" ✓ SDF file CSV ")
        print(" ✓ SDF file Excel ()")
        print(" ✓ moleculecalculate")
        print(" ✓ ")
        print(" ✓ ")
        print(" ✓ ")
        print("="*60)
        print()

    def check_dependencies(self):
        """"""
        print("...")

        required_libs = {
            'pandas': 'data',
            'rdkit': '',
            'openpyxl': 'Excelfile',
            'tqdm': ''
        }

        missing_libs = []
        for lib, desc in required_libs.items():
            try:
                if lib == 'rdkit':
                    import rdkit
                    version = getattr(rdkit, '__version__', '')
                    print(f" ✓ {lib} ({desc}) - Version: {version}")
                elif lib == 'tqdm':
                    import tqdm
                    print(f" ✓ {lib} ({desc}) - Version: {tqdm.__version__}")
                else:
                    __import__(lib)
                    print(f"  ✓ {lib} ({desc})")
            except ImportError:
                missing_libs.append(lib)
                print(f" ✗ {lib} ({desc}) - ")

        if missing_libs:
            print(f"\n⚠ : {', '.join(missing_libs)}")
            print("? (y/n): ", end="")
            response = input().strip().lower()

            if response == 'y':
                self.install_dependencies(missing_libs)
            else:
                print("\n:")
                print("pip install pandas rdkit-pypi openpyxl tqdm")
                print("\n...")
                input()

        print("!\n")

    def install_dependencies(self, libs):
        """"""
        print("...")

        install_commands = {
            'rdkit': 'pip install rdkit-pypi',
            'pandas': 'pip install pandas',
            'openpyxl': 'pip install openpyxl',
            'tqdm': 'pip install tqdm',
            'xlsxwriter': 'pip install xlsxwriter'
        }

        for lib in libs:
            if lib in install_commands:
                print(f" {lib}...")
                os.system(install_commands[lib])
                time.sleep(1)

        print("! . ")
        input("...")
        sys.exit(0)

    def get_input_file(self):
        """inputSDFfilepath"""
        print("\n" + "="*60)
        print(" 1: inputSDFfilepath")
        print("="*60)

        while True:
            print("\ninput:")
            print(" 1. inputSDFfilepath")
            print(" 2. file")
            print(" 3. file")
            print(" 4. ")
            print("\n (1-4): ", end="")

            choice = input().strip()

            if choice == '4' or choice.lower() == 'exit':
                print("Use, !")
                sys.exit(0)

            if choice == '1':
                print("\ninputSDFfilepath:")
                print("input: ", end="")
                user_input = input().strip().strip('"').strip("'")
                return self.validate_input_file(user_input)

            elif choice == '2':
                print("\nSDFfile, :")
                print("file: ", end="")
                user_input = input().strip().strip('"').strip("'")
                return self.validate_input_file(user_input)

            elif choice == '3':
                return self.handle_batch_mode()

            else:
                print("❌ , input")

    def validate_input_file(self, user_input):
        """inputfile"""
        if not user_input or user_input.lower() == 'exit':
            return []

        path = Path(user_input)

        if not path.exists():
            print(f"❌ : path - {user_input}")
            return []

        if path.is_file():
            if path.suffix.lower() in self.supported_formats:
                print(f"✓ file: {path.name}")
                return [path]
            else:
                print(f"❌ : unsupported file format - {path.suffix}")
                print(f": {', '.join(self.supported_formats)}")
                return []
        elif path.is_dir():
            return self.find_sdf_files_in_folder(path)

        return []

    def find_sdf_files_in_folder(self, folder_path):
        """fileSDFfile"""
        print(f"fileSDFfile: {folder_path}")

        sdf_files = []
        for ext in self.supported_formats:
            sdf_files.extend(folder_path.glob(f"*{ext}"))
            sdf_files.extend(folder_path.glob(f"*{ext.upper()}"))
            sdf_files.extend(folder_path.glob(f"**/*{ext}"))
            sdf_files.extend(folder_path.glob(f"**/*{ext.upper()}"))

        # translated note
        sdf_files = sorted(list(set(sdf_files)))

        if not sdf_files:
            print(f"❌ : fileSDFfile")
            return []

        print(f"✓ {len(sdf_files)} SDFfile")

        # filecolumn
        print("\nfilecolumn:")
        for i, f in enumerate(sdf_files[:10], 1):
            rel_path = f.relative_to(folder_path) if f.is_relative_to(folder_path) else f
            print(f"  {i:3d}. {rel_path}")

        if len(sdf_files) > 10:
            print(f" ... {len(sdf_files) - 10} file")

        print("\nfile? (y/n): ", end="")
        response = input().strip().lower()

        if response == 'y':
            return sdf_files
        else:
            return []

    def handle_batch_mode(self):
        """"""
        print("\n" + "="*60)
        print("")
        print("="*60)

        print("\ninputSDFfilefilepath:")
        print("(file)")
        print("input: ", end="")

        folder_path = input().strip().strip('"').strip("'")
        folder = Path(folder_path)

        if not folder.exists() or not folder.is_dir():
            print("❌ : file")
            return []

        return self.find_sdf_files_in_folder(folder)

    def get_output_options(self, input_files):
        """output"""
        print("\n" + "="*60)
        print(" 2: output")
        print("="*60)

        options = {
            'formats': [],
            'output_dir': self.output_dir,
            'prefix': '',
            'add_calculated_props': True,
            'include_structures': False,
            'excel_engine': 'openpyxl',
            'csv_encoding': 'utf-8-sig'
        }

        # output
        print("\noutput:")
        print(" 1. CSV (.csv) - ")
        print(" 2. Excel (.xlsx) - ")
        print(" 3. ")
        print(" 4. ")
        print("\n (1-4): ", end="")

        format_choice = input().strip()

        if format_choice == '1':
            options['formats'] = ['csv']
        elif format_choice == '2':
            options['formats'] = ['excel']
        elif format_choice == '3':
            options['formats'] = ['csv', 'excel']
        elif format_choice == '4':
            print("input (csv, excel): ", end="")
            custom = input().strip().lower().split(',')
            if 'csv' in [f.strip() for f in custom]:
                options['formats'].append('csv')
            if 'excel' in [f.strip() for f in custom] or 'xlsx' in [f.strip() for f in custom]:
                options['formats'].append('excel')

        if not options['formats']:
            print("⚠ , UseCSV")
            options['formats'] = ['csv']

        # outputdirectory
        print(f"\noutputdirectory: {options['output_dir']}")
        print("outputdirectory? (y/n): ", end="")
        if input().strip().lower() == 'y':
            print("inputoutputdirectory: ", end="")
            new_dir = input().strip().strip('"').strip("'")
            options['output_dir'] = Path(new_dir)

        # outputdirectory
        options['output_dir'].mkdir(parents=True, exist_ok=True)

        # file
        if len(input_files) == 1:
            default_prefix = input_files[0].stem
        else:
            default_prefix = datetime.now().strftime("%Y%m%d_%H%M")

        print(f"\nfile (Use '{default_prefix}'): ", end="")
        prefix = input().strip()
        options['prefix'] = prefix if prefix else default_prefix

        # calculate
        print("\ncalculate (SMILES, molecule)? (y/n): ", end="")
        options['add_calculated_props'] = input().strip().lower() == 'y'

        # Excel
        if 'excel' in options['formats']:
            print("\nExcel:")
            print(" 1. openpyxl (, )")
            print(" 2. xlsxwriter ()")
            print(" (1-2, Use): ", end="")
            engine_choice = input().strip()
            if engine_choice == '2':
                options['excel_engine'] = 'xlsxwriter'

        # CSV
        if 'csv' in options['formats']:
            print("\nCSVfile:")
            print(" 1. utf-8-sig (, )")
            print("  2. utf-8")
            print(" 3. gb18030 ()")
            print(" (1-3, Useutf-8-sig): ", end="")
            encoding_choice = input().strip()
            if encoding_choice == '2':
                options['csv_encoding'] = 'utf-8'
            elif encoding_choice == '3':
                options['csv_encoding'] = 'gb18030'

        return options

    def read_sdf_with_progress(self, sdf_file):
        """UseSDFfile"""
        print(f"file: {sdf_file.name}")

        try:
            # method1: UseSDMolSupplier +
            print("method1: UseSDMolSupplier...")
            suppl = Chem.SDMolSupplier(str(sdf_file))

            # molecule ()
            try:
                # translated note
                total_mols = len(suppl)
            except:
                # ,
                print("File size...")
                with open(sdf_file, 'r', encoding='latin-1') as f:
                    content = f.read()
                    # : "$$$$" molecule
                    total_mols = content.count('$$$$')

            print(f"molecule: {total_mols}")

            # Usemolecule
            molecules = []
            properties_list = []

            progress_bar = CustomProgressBar(total_mols, "molecule")

            for i, mol in enumerate(suppl):
                if mol is not None:
                    molecules.append(mol)

                    # translated note
                    props = self.extract_mol_properties(mol)
                    properties_list.append(props)

                progress_bar.update(1)

                # 100moleculesave
                if i > 0 and i % 100 == 0:
                    progress_bar.desc = f" {i}/{total_mols}"

            progress_bar.close()

            if not molecules:
                print("❌ : molecule")
                return None

            # DataFrame
            df = pd.DataFrame(properties_list)
            df['ROMol'] = molecules

            print(f"✓ {len(molecules)} molecule")
            return df

        except Exception as e:
            print(f"❌ method1failed: {e}")

            # method2: UsePandasTools ()
            try:
                print("method2: UsePandasTools...")
                df = PandasTools.LoadSDF(str(sdf_file))
                print(f"✓ UsePandasTools: {len(df)} molecule")
                return df
            except Exception as e2:
                print(f"❌ method2failed: {e2}")
                return None

    def extract_mol_properties(self, mol):
        """molecule"""
        props = {}

        # translated note
        if mol.HasProp("_Name"):
            props['Name'] = mol.GetProp("_Name")

        # translated note
        prop_names = mol.GetPropNames()
        for prop in prop_names:
            if prop != "_Name": #
                try:
                    value = mol.GetProp(prop)
                    # (Excel)
                    clean_prop = self.clean_column_name(prop)
                    props[clean_prop] = value
                except:
                    props[prop] = ""

        return props

    def clean_column_name(self, column_name):
        """columnExcel"""
        # translated note
        cleaned = re.sub(r'[\[\]{}()<>+=!@#$%^&*|\\/~`]', '', column_name)
        # translated note
        cleaned = cleaned.replace(' ', '_').replace('-', '_')
        # translated note
        if cleaned and not cleaned[0].isalpha():
            cleaned = 'C_' + cleaned
        # translated note
        if len(cleaned) > 50:
            cleaned = cleaned[:50]

        return cleaned

    def add_calculated_properties(self, df, progress_bar=None):
        """calculateDataFrame"""
        if 'ROMol' not in df.columns:
            return df

        print("calculatemolecule...")

        # column
        new_columns = {
            'SMILES': [],
            'InChI': [],
            'MolecularWeight': [],
            'MolecularFormula': [],
            'NumHeavyAtoms': [],
            'NumAtoms': [],
            'NumBonds': [],
            'NumRotatableBonds': [],
            'NumHDonors': [],
            'NumHAcceptors': [],
            'TPSA': [], #
            'LogP': [], # /
            'RingCount': [] #
        }

        total_mols = len(df)
        if progress_bar is None:
            progress_bar = CustomProgressBar(total_mols, "calculate")

        for idx, mol in enumerate(df['ROMol']):
            if mol is not None:
                # SMILES
                try:
                    new_columns['SMILES'].append(Chem.MolToSmiles(mol))
                except:
                    new_columns['SMILES'].append('')

                # InChI
                try:
                    new_columns['InChI'].append(inchi.MolToInchi(mol))
                except:
                    new_columns['InChI'].append('')

                # molecule
                try:
                    new_columns['MolecularWeight'].append(Descriptors.MolWt(mol))
                except:
                    new_columns['MolecularWeight'].append(0)

                # molecular_formula
                try:
                    new_columns['MolecularFormula'].append(Chem.rdMolDescriptors.CalcMolFormula(mol))
                except:
                    new_columns['MolecularFormula'].append('')

                # translated note
                new_columns['NumHeavyAtoms'].append(mol.GetNumHeavyAtoms())
                new_columns['NumAtoms'].append(mol.GetNumAtoms())
                new_columns['NumBonds'].append(mol.GetNumBonds())

                # translated note
                try:
                    new_columns['NumRotatableBonds'].append(Descriptors.NumRotatableBonds(mol))
                except:
                    new_columns['NumRotatableBonds'].append(0)

                try:
                    new_columns['NumHDonors'].append(Descriptors.NumHDonors(mol))
                except:
                    new_columns['NumHDonors'].append(0)

                try:
                    new_columns['NumHAcceptors'].append(Descriptors.NumHAcceptors(mol))
                except:
                    new_columns['NumHAcceptors'].append(0)

                try:
                    new_columns['TPSA'].append(Descriptors.TPSA(mol))
                except:
                    new_columns['TPSA'].append(0)

                try:
                    new_columns['LogP'].append(Descriptors.MolLogP(mol))
                except:
                    new_columns['LogP'].append(0)

                try:
                    new_columns['RingCount'].append(Descriptors.RingCount(mol))
                except:
                    new_columns['RingCount'].append(0)
            else:
                # molecule,
                for col in new_columns:
                    if col in ['SMILES', 'InChI', 'MolecularFormula']:
                        new_columns[col].append('')
                    else:
                        new_columns[col].append(0)

            progress_bar.update(1)

        progress_bar.close()

        # columnDataFrame
        for col_name, values in new_columns.items():
            df[col_name] = values

        return df

    def fix_dataframe_for_excel(self, df):
        """DataFrameExcel"""
        print("dataExcel...")

        # translated note
        df_fixed = df.copy()

        # ROMolcolumn (Excel)
        if 'ROMol' in df_fixed.columns:
            df_fixed = df_fixed.drop(columns=['ROMol'])

        # column
        df_fixed.columns = [self.clean_column_name(str(col)) for col in df_fixed.columns]

        # data
        for col in df_fixed.columns:
            # translated note
            if df_fixed[col].dtype == 'object':
                try:
                    # translated note
                    sample = df_fixed[col].dropna().iloc[0] if not df_fixed[col].dropna().empty else None
                    if sample is not None and not isinstance(sample, (str, int, float, bool)):
                        df_fixed[col] = df_fixed[col].astype(str)
                except:
                    df_fixed[col] = df_fixed[col].astype(str)

            # NaN
            df_fixed[col] = df_fixed[col].fillna('')

        # (Excel1048576)
        if len(df_fixed) > 1000000:
            print(f"⚠ Warning: data({len(df_fixed)})Excel, ")
            df_fixed = df_fixed.head(1000000)

        # column (Excelcolumn16384)
        if len(df_fixed.columns) > 1000:
            print(f"⚠ Warning: datacolumn({len(df_fixed.columns)}), retention1000column")
            df_fixed = df_fixed.iloc[:, :1000]

        print(f"✓ data: {len(df_fixed)} x {len(df_fixed.columns)}column")
        return df_fixed

    def export_to_csv(self, df, output_path, encoding='utf-8-sig'):
        """CSVfile"""
        try:
            print(f"CSV: {output_path.name}")

            # data
            df_fixed = self.fix_dataframe_for_excel(df)

            # CSV
            df_fixed.to_csv(output_path, index=False, encoding=encoding)

            print(f"✓ CSV file saved: {output_path}")
            print(f" File size: {os.path.getsize(output_path) / 1024:.1f} KB")

            return True
        except Exception as e:
            print(f"❌ CSVfailed: {e}")
            return False

    def export_to_excel(self, df, output_path, engine='openpyxl'):
        """Excelfile ()"""
        try:
            print(f"Excel: {output_path.name}")
            print(f"Use: {engine}")

            # data
            df_fixed = self.fix_dataframe_for_excel(df)

            # Excel
            with pd.ExcelWriter(output_path, engine=engine) as writer:
                df_fixed.to_excel(writer, sheet_name='Molecules', index=False)

                # openpyxl, column
                if engine == 'openpyxl':
                    worksheet = writer.sheets['Molecules']

                    # column
                    for column in worksheet.columns:
                        max_length = 0
                        column_letter = column[0].column_letter

                        for cell in column:
                            try:
                                cell_length = len(str(cell.value))
                                if cell_length > max_length:
                                    max_length = cell_length
                            except:
                                pass

                        adjusted_width = min(max_length + 2, 50) # 50
                        worksheet.column_dimensions[column_letter].width = adjusted_width

            print(f"✓ Excelfilesave: {output_path}")
            print(f" File size: {os.path.getsize(output_path) / 1024:.1f} KB")

            return True
        except Exception as e:
            print(f"❌ Excelfailed: {e}")

            # translated note
            if engine == 'openpyxl':
                print("Usexlsxwriter...")
                return self.export_to_excel(df, output_path, engine='xlsxwriter')
            else:
                print("Useopenpyxl...")
                return self.export_to_excel(df, output_path, engine='openpyxl')

    def convert_sdf_file(self, sdf_file, options):
        """SDFfile"""
        print(f"\n" + "="*60)
        print(f"file: {sdf_file.name}")
        print("="*60)

        start_time = time.time()
        result = {
            'file': sdf_file.name,
            'success': False,
            'molecules_processed': 0,
            'molecules_successful': 0,
            'output_files': [],
            'error': None,
            'time_taken': 0
        }

        try:
            # 1. SDFfile
            df = self.read_sdf_with_progress(sdf_file)

            if df is None or len(df) == 0:
                result['error'] = "SDFfilemolecule"
                result['time_taken'] = time.time() - start_time
                return result

            result['molecules_processed'] = len(df)

            # 2. calculate
            if options['add_calculated_props']:
                progress_bar = CustomProgressBar(len(df), "calculatemolecule")
                df = self.add_calculated_properties(df, progress_bar)

            result['molecules_successful'] = len(df)

            # 3. outputfile
            if options['prefix']:
                base_name = f"{options['prefix']}_{sdf_file.stem}"
            else:
                base_name = sdf_file.stem

            # 4.
            output_files = []
            for fmt in options['formats']:
                if fmt == 'csv':
                    output_path = options['output_dir'] / f"{base_name}.csv"
                    if self.export_to_csv(df, output_path, options['csv_encoding']):
                        output_files.append(output_path)

                elif fmt == 'excel':
                    output_path = options['output_dir'] / f"{base_name}.xlsx"
                    if self.export_to_excel(df, output_path, options['excel_engine']):
                        output_files.append(output_path)

            result['output_files'] = output_files
            result['success'] = True

        except Exception as e:
            result['error'] = str(e)
            print(f"❌ failed: {e}")
            traceback.print_exc()

        result['time_taken'] = time.time() - start_time
        return result

    def generate_report(self, results, options):
        """"""
        print("\n" + "="*60)
        print("")
        print("="*60)

        total_files = len(results)
        successful_files = sum(1 for r in results if r['success'])
        total_molecules = sum(r['molecules_processed'] for r in results)
        successful_molecules = sum(r['molecules_successful'] for r in results)

        # output
        print(f"\n📊 :")
        print(f" file: {total_files}")
        print(f" file: {successful_files}")
        print(f" molecule: {total_molecules}")
        print(f" molecule: {successful_molecules}")

        if total_molecules > 0:
            success_rate = (successful_molecules / total_molecules) * 100
            print(f" : {success_rate:.1f}%")

        total_time = sum(r['time_taken'] for r in results)
        print(f" : {total_time:.1f}")

        if successful_molecules > 0:
            print(f" : {successful_molecules/total_time:.1f} molecule/")

        # outputfilecolumn
        print(f"\n📁 outputfile ({options['output_dir']}):")
        all_output_files = []
        for result in results:
            if result['success']:
                for file in result['output_files']:
                    if file.exists():
                        file_size = file.stat().st_size / 1024
                        print(f"   ✓ {file.name} ({file_size:.1f} KB)")
                        all_output_files.append(file)

        # file
        report_file = options['output_dir'] / f"{options['prefix']}_report.txt"
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("SDFfile\n")
            f.write("="*50 + "\n")
            f.write(f": {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Version: {self.version}\n")
            f.write(f"outputdirectory: {options['output_dir']}\n")
            f.write(f"output: {', '.join(options['formats'])}\n")
            f.write(f"calculate: {'' if options['add_calculated_props'] else ''}\n")
            f.write("\nresults:\n")

            for i, result in enumerate(results, 1):
                f.write(f"\n[{i}] file: {result['file']}\n")
                f.write(f" : {'' if result['success'] else 'failed'}\n")
                f.write(f" molecule: {result['molecules_processed']}\n")
                f.write(f" molecule: {result['molecules_successful']}\n")
                f.write(f" : {result.get('time_taken', 0):.2f}\n")

                if result['output_files']:
                    f.write(" outputfile:\n")
                    for out_file in result['output_files']:
                        file_size = out_file.stat().st_size / 1024 if out_file.exists() else 0
                        f.write(f"     - {out_file.name} ({file_size:.1f} KB)\n")

                if result['error']:
                    f.write(f" : {result['error']}\n")

        print(f"\n📝 save: {report_file.name}")

        # Record log messages
        self.log_conversion(results, options)

        return report_file

    def log_conversion(self, results, options):
        """"""
        try:
            with open(self.log_file, 'a', encoding='utf-8') as f:
                timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                f.write(f"\n{'='*60}\n")
                f.write(f": {timestamp}\n")
                f.write(f"outputdirectory: {options['output_dir']}\n")

                for result in results:
                    status = "" if result['success'] else "failed"
                    f.write(f"{result['file']}: {status} ({result['molecules_successful']}/{result['molecules_processed']})\n")
        except:
            pass

    def show_data_preview(self, output_dir):
        """outputfile"""
        print("\n" + "="*60)
        print("data")
        print("="*60)

        # file
        csv_files = list(output_dir.glob("*.csv"))
        excel_files = list(output_dir.glob("*.xlsx"))

        preview_files = (csv_files[:2] + excel_files[:1])[:3] # 3file

        if not preview_files:
            print("outputfile")
            return

        for file in preview_files:
            try:
                print(f"\n📋 {file.name}:")
                print("-"*40)

                if file.suffix == '.csv':
                    # 5
                    df_preview = pd.read_csv(file, nrows=5)
                else:
                    # Excel5
                    df_preview = pd.read_excel(file, nrows=5)

                # translated note
                print(f": {len(pd.read_csv(file) if file.suffix == '.csv' else pd.read_excel(file))}")
                print(f"column: {len(df_preview.columns)}")

                # translated note
                print("\n5data:")
                pd.set_option('display.max_columns', 10) # 10column
                pd.set_option('display.width', 120)
                print(df_preview)

                # column
                print(f"\n10column:")
                for i, col in enumerate(df_preview.columns[:10], 1):
                    print(f"  {i:2d}. {col}")

                if len(df_preview.columns) > 10:
                    print(f" ... {len(df_preview.columns) - 10} column")

            except Exception as e:
                print(f" {file.name}: {e}")

    def open_output_directory(self, output_dir):
        """outputdirectory"""
        try:
            if sys.platform == 'win32':
                os.startfile(output_dir)
            elif sys.platform == 'darwin':
                os.system(f'open "{output_dir}"')
            elif sys.platform == 'linux':
                os.system(f'xdg-open "{output_dir}"')
            print(f"✓ outputdirectory")
        except:
            print(f"⚠ directory, : {output_dir}")

    def main(self):
        """"""
        try:
            # 1.
            self.print_banner()

            # 2.
            self.check_dependencies()

            # 3. inputfile
            input_files = self.get_input_file()

            if not input_files:
                print("file, . ")
                return

            # 4. output
            options = self.get_output_options(input_files)

            # 5. file
            print("\n" + "="*60)
            print(" 3: file")
            print("="*60)

            results = []
            total_start_time = time.time()

            for i, sdf_file in enumerate(input_files, 1):
                print(f"\n[{i}/{len(input_files)}] ", end="")
                result = self.convert_sdf_file(sdf_file, options)
                results.append(result)

                # translated note
                progress = (i / len(input_files)) * 100
                elapsed = time.time() - total_start_time
                print(f": {progress:.1f}% | : {elapsed:.1f}")

            # 6.
            report_file = self.generate_report(results, options)

            # 7.
            self.show_data_preview(options['output_dir'])

            # 8.
            print("\n" + "="*60)
            print("🎉 !")
            print("="*60)

            total_time = time.time() - total_start_time
            print(f"\n: {total_time:.1f}")
            print(f"outputdirectory: {options['output_dir']}")
            print(f": {report_file.name}")

            # outputdirectory
            print("\noutputdirectory? (y/n): ", end="")
            if input().strip().lower() == 'y':
                self.open_output_directory(options['output_dir'])

            print("\nUse SDF file!")
            print("...")
            input()

        except KeyboardInterrupt:
            print("\n\n. ")
        except Exception as e:
            print(f"\n❌ : {e}")
            print("\n:")
            traceback.print_exc()
            print("\n...")
            input()

# translated note
def quick_convert_mode():
    """ - """
    print(" SDF ")
    print("=" * 50)

    converter = InteractiveSDFConverter()

    # file
    print("\nSDFfileinputfilepath:")
    file_path = input("file: ").strip().strip('"').strip("'")

    if not file_path:
        print("inputfilepath")
        return

    path = Path(file_path)
    if not path.exists():
        print(f"file does not exist: {file_path}")
        return

    # translated note
    print("\noutput:")
    print("1. CSV")
    print("2. Excel")
    print("3. ")
    choice = input(" (1-3): ").strip()

    formats = []
    if choice in ['1', '3']:
        formats.append('csv')
    if choice in ['2', '3']:
        formats.append('excel')

    if not formats:
        formats = ['csv']

    # translated note
    result = converter.convert_sdf_file(
        path,
        {
            'formats': formats,
            'output_dir': converter.output_dir,
            'prefix': path.stem,
            'add_calculated_props': True,
            'include_structures': False,
            'excel_engine': 'openpyxl',
            'csv_encoding': 'utf-8-sig'
        }
    )

    if result['success']:
        print(f"\n✅ !")
        for out_file in result['output_files']:
            print(f"   ✓ {out_file.name}")

        # directory
        print("\noutputdirectory? (y/n): ", end="")
        if input().strip().lower() == 'y':
            converter.open_output_directory(converter.output_dir)
    else:
        print(f"\n❌ failed: {result['error']}")

if __name__ == "__main__":
    # Parameters
    if len(sys.argv) > 1:
        if sys.argv[1] in ['--quick', '-q']:
            quick_convert_mode()
        elif sys.argv[1] in ['--help', '-h']:
            print("SDF file v2.0")
            print("\nUsemethod:")
            print(" python sdf_converter.py # ")
            print(" python sdf_converter.py --quick # ")
            print(" python sdf_converter.py --help # ")
            print("\n:")
            print(" - SDFCSVExcel")
            print(" - calculatemolecule")
            print(" - ")
            print(" - Excel")
        else:
            # filepath,
            converter = InteractiveSDFConverter()
            path = Path(sys.argv[1])
            if path.exists():
                result = converter.convert_sdf_file(
                    path,
                    {
                        'formats': ['csv', 'excel'],
                        'output_dir': converter.output_dir,
                        'prefix': path.stem,
                        'add_calculated_props': True,
                        'include_structures': False,
                        'excel_engine': 'openpyxl',
                        'csv_encoding': 'utf-8-sig'
                    }
                )

                if result['success']:
                    print(f"\n✅ !")
                    for out_file in result['output_files']:
                        print(f"   ✓ {out_file.name}")
                else:
                    print(f"\n❌ failed: {result['error']}")
            else:
                print(f"file does not exist: {sys.argv[1]}")
    else:
        # translated note
        converter = InteractiveSDFConverter()
        converter.main()