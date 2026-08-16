"""
SDF 转 CSV/Excel 交互式工具 - 带进度条和格式修正
运行后按照提示输入文件路径和选项即可完成转换
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

# 自定义进度条样式
class CustomProgressBar:
    def __init__(self, total, desc="处理中"):
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

        # 进度显示
        sys.stdout.write(f'\r{self.desc}: |{bar}| {self.current}/{self.total} '
                        f'[{percent:.1%}] 耗时: {elapsed:.1f}s ETA: {eta:.1f}s')
        sys.stdout.flush()

    def close(self):
        elapsed = time.time() - self.start_time
        sys.stdout.write(f'\r{self.desc}: 完成! 共耗时: {elapsed:.1f}s\n')
        sys.stdout.flush()

class InteractiveSDFConverter:
    def __init__(self):
        """初始化交互式转换器"""
        self.version = "2.0.0"
        self.author = "SDF Converter Tool"
        self.supported_formats = ['.sdf', '.sd', '.mol']
        self.output_dir = Path.cwd() / "SDF_Output"
        self.log_file = self.output_dir / "conversion_log.txt"

        # 创建输出目录
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def print_banner(self):
        """打印程序横幅"""
        print("\n" + "="*60)
        print("         SDF 文件转换工具 v2.0")
        print("="*60)
        print(f"版本: {self.version}")
        print(f"时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print("-"*60)
        print("功能:")
        print("  ✓ 将 SDF 文件转换为 CSV 格式")
        print("  ✓ 将 SDF 文件转换为 Excel 格式（已修正格式）")
        print("  ✓ 自动提取分子属性和计算描述符")
        print("  ✓ 批量处理支持")
        print("  ✓ 实时进度条显示")
        print("  ✓ 智能格式修正")
        print("="*60)
        print()

    def check_dependencies(self):
        """检查必要的依赖库"""
        print("正在检查依赖库...")

        required_libs = {
            'pandas': '数据处理',
            'rdkit': '化学信息学处理',
            'openpyxl': 'Excel文件写入',
            'tqdm': '进度条显示'
        }

        missing_libs = []
        for lib, desc in required_libs.items():
            try:
                if lib == 'rdkit':
                    import rdkit
                    version = getattr(rdkit, '__version__', '未知')
                    print(f"  ✓ {lib} ({desc}) - 版本: {version}")
                elif lib == 'tqdm':
                    import tqdm
                    print(f"  ✓ {lib} ({desc}) - 版本: {tqdm.__version__}")
                else:
                    __import__(lib)
                    print(f"  ✓ {lib} ({desc})")
            except ImportError:
                missing_libs.append(lib)
                print(f"  ✗ {lib} ({desc}) - 未安装")

        if missing_libs:
            print(f"\n⚠ 缺少必要的库: {', '.join(missing_libs)}")
            print("是否尝试自动安装? (y/n): ", end="")
            response = input().strip().lower()

            if response == 'y':
                self.install_dependencies(missing_libs)
            else:
                print("\n请手动安装缺少的库:")
                print("pip install pandas rdkit-pypi openpyxl tqdm")
                print("\n按回车键继续尝试运行...")
                input()

        print("依赖检查完成!\n")

    def install_dependencies(self, libs):
        """尝试安装依赖库"""
        print("正在安装依赖库...")

        install_commands = {
            'rdkit': 'pip install rdkit-pypi',
            'pandas': 'pip install pandas',
            'openpyxl': 'pip install openpyxl',
            'tqdm': 'pip install tqdm',
            'xlsxwriter': 'pip install xlsxwriter'
        }

        for lib in libs:
            if lib in install_commands:
                print(f"正在安装 {lib}...")
                os.system(install_commands[lib])
                time.sleep(1)

        print("安装完成! 请重启程序。")
        input("按回车键退出...")
        sys.exit(0)

    def get_input_file(self):
        """获取用户输入的SDF文件路径"""
        print("\n" + "="*60)
        print("步骤 1: 输入SDF文件路径")
        print("="*60)

        while True:
            print("\n请选择输入方式:")
            print("  1. 输入单个SDF文件路径")
            print("  2. 拖拽文件到此处")
            print("  3. 批量处理文件夹")
            print("  4. 退出程序")
            print("\n选择 (1-4): ", end="")

            choice = input().strip()

            if choice == '4' or choice.lower() == 'exit':
                print("感谢使用，再见!")
                sys.exit(0)

            if choice == '1':
                print("\n请输入SDF文件完整路径:")
                print("输入: ", end="")
                user_input = input().strip().strip('"').strip("'")
                return self.validate_input_file(user_input)

            elif choice == '2':
                print("\n请拖拽SDF文件到此处，然后按回车:")
                print("拖拽文件: ", end="")
                user_input = input().strip().strip('"').strip("'")
                return self.validate_input_file(user_input)

            elif choice == '3':
                return self.handle_batch_mode()

            else:
                print("❌ 无效选择，请重新输入")

    def validate_input_file(self, user_input):
        """验证输入文件"""
        if not user_input or user_input.lower() == 'exit':
            return []

        path = Path(user_input)

        if not path.exists():
            print(f"❌ 错误: 路径不存在 - {user_input}")
            return []

        if path.is_file():
            if path.suffix.lower() in self.supported_formats:
                print(f"✓ 找到文件: {path.name}")
                return [path]
            else:
                print(f"❌ 错误: 不支持的文件格式 - {path.suffix}")
                print(f"支持的格式: {', '.join(self.supported_formats)}")
                return []
        elif path.is_dir():
            return self.find_sdf_files_in_folder(path)

        return []

    def find_sdf_files_in_folder(self, folder_path):
        """在文件夹中查找SDF文件"""
        print(f"在文件夹中查找SDF文件: {folder_path}")

        sdf_files = []
        for ext in self.supported_formats:
            sdf_files.extend(folder_path.glob(f"*{ext}"))
            sdf_files.extend(folder_path.glob(f"*{ext.upper()}"))
            sdf_files.extend(folder_path.glob(f"**/*{ext}"))
            sdf_files.extend(folder_path.glob(f"**/*{ext.upper()}"))

        # 去重并排序
        sdf_files = sorted(list(set(sdf_files)))

        if not sdf_files:
            print(f"❌ 错误: 文件夹中没有找到SDF文件")
            return []

        print(f"✓ 找到 {len(sdf_files)} 个SDF文件")

        # 显示文件列表
        print("\n文件列表:")
        for i, f in enumerate(sdf_files[:10], 1):
            rel_path = f.relative_to(folder_path) if f.is_relative_to(folder_path) else f
            print(f"  {i:3d}. {rel_path}")

        if len(sdf_files) > 10:
            print(f"  ... 还有 {len(sdf_files) - 10} 个文件")

        print("\n是否处理所有这些文件? (y/n): ", end="")
        response = input().strip().lower()

        if response == 'y':
            return sdf_files
        else:
            return []

    def handle_batch_mode(self):
        """处理批量转换模式"""
        print("\n" + "="*60)
        print("批量转换模式")
        print("="*60)

        print("\n请输入包含SDF文件的文件夹路径:")
        print("(支持拖拽文件夹到此处)")
        print("输入: ", end="")

        folder_path = input().strip().strip('"').strip("'")
        folder = Path(folder_path)

        if not folder.exists() or not folder.is_dir():
            print("❌ 错误: 文件夹不存在或无效")
            return []

        return self.find_sdf_files_in_folder(folder)

    def get_output_options(self, input_files):
        """获取输出选项"""
        print("\n" + "="*60)
        print("步骤 2: 设置输出选项")
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

        # 选择输出格式
        print("\n请选择输出格式:")
        print("  1. CSV 格式 (.csv) - 推荐")
        print("  2. Excel 格式 (.xlsx) - 已修正格式")
        print("  3. 两种格式都要")
        print("  4. 自定义选择")
        print("\n选择 (1-4): ", end="")

        format_choice = input().strip()

        if format_choice == '1':
            options['formats'] = ['csv']
        elif format_choice == '2':
            options['formats'] = ['excel']
        elif format_choice == '3':
            options['formats'] = ['csv', 'excel']
        elif format_choice == '4':
            print("请输入格式 (csv, excel): ", end="")
            custom = input().strip().lower().split(',')
            if 'csv' in [f.strip() for f in custom]:
                options['formats'].append('csv')
            if 'excel' in [f.strip() for f in custom] or 'xlsx' in [f.strip() for f in custom]:
                options['formats'].append('excel')

        if not options['formats']:
            print("⚠ 未选择格式，默认使用CSV格式")
            options['formats'] = ['csv']

        # 选择输出目录
        print(f"\n当前输出目录: {options['output_dir']}")
        print("是否更改输出目录? (y/n): ", end="")
        if input().strip().lower() == 'y':
            print("请输入新的输出目录: ", end="")
            new_dir = input().strip().strip('"').strip("'")
            options['output_dir'] = Path(new_dir)

        # 创建输出目录
        options['output_dir'].mkdir(parents=True, exist_ok=True)

        # 设置文件名前缀
        if len(input_files) == 1:
            default_prefix = input_files[0].stem
        else:
            default_prefix = datetime.now().strftime("%Y%m%d_%H%M")

        print(f"\n文件名前缀 (留空使用 '{default_prefix}'): ", end="")
        prefix = input().strip()
        options['prefix'] = prefix if prefix else default_prefix

        # 选择是否添加计算属性
        print("\n是否添加计算属性 (SMILES, 分子量等)? (y/n): ", end="")
        options['add_calculated_props'] = input().strip().lower() == 'y'

        # 选择Excel引擎
        if 'excel' in options['formats']:
            print("\n选择Excel写入引擎:")
            print("  1. openpyxl (默认，功能完整)")
            print("  2. xlsxwriter (速度快)")
            print("选择 (1-2, 留空使用默认): ", end="")
            engine_choice = input().strip()
            if engine_choice == '2':
                options['excel_engine'] = 'xlsxwriter'

        # 选择CSV编码
        if 'csv' in options['formats']:
            print("\n选择CSV文件编码:")
            print("  1. utf-8-sig (推荐，支持中文)")
            print("  2. utf-8")
            print("  3. gb18030 (中文编码)")
            print("选择 (1-3, 留空使用utf-8-sig): ", end="")
            encoding_choice = input().strip()
            if encoding_choice == '2':
                options['csv_encoding'] = 'utf-8'
            elif encoding_choice == '3':
                options['csv_encoding'] = 'gb18030'

        return options

    def read_sdf_with_progress(self, sdf_file):
        """使用进度条读取SDF文件"""
        print(f"读取文件: {sdf_file.name}")

        try:
            # 方法1: 使用SDMolSupplier + 进度条
            print("方法1: 使用SDMolSupplier读取...")
            suppl = Chem.SDMolSupplier(str(sdf_file))

            # 首先获取总分子数（如果可能）
            try:
                # 尝试获取总数
                total_mols = len(suppl)
            except:
                # 如果无法直接获取，估计总数
                print("正在估算文件大小...")
                with open(sdf_file, 'r', encoding='latin-1') as f:
                    content = f.read()
                    # 粗略估算：每出现一次 "$$$$" 表示一个分子结束
                    total_mols = content.count('$$$$')

            print(f"估计分子数: {total_mols}")

            # 使用进度条读取分子
            molecules = []
            properties_list = []

            progress_bar = CustomProgressBar(total_mols, "读取分子")

            for i, mol in enumerate(suppl):
                if mol is not None:
                    molecules.append(mol)

                    # 提取属性
                    props = self.extract_mol_properties(mol)
                    properties_list.append(props)

                progress_bar.update(1)

                # 每100个分子保存一次进度
                if i > 0 and i % 100 == 0:
                    progress_bar.desc = f"已读取 {i}/{total_mols}"

            progress_bar.close()

            if not molecules:
                print("❌ 错误: 未能读取任何有效分子")
                return None

            # 创建DataFrame
            df = pd.DataFrame(properties_list)
            df['ROMol'] = molecules

            print(f"✓ 成功读取 {len(molecules)} 个分子")
            return df

        except Exception as e:
            print(f"❌ 方法1失败: {e}")

            # 方法2: 使用PandasTools（备用）
            try:
                print("尝试方法2: 使用PandasTools读取...")
                df = PandasTools.LoadSDF(str(sdf_file))
                print(f"✓ 使用PandasTools读取成功: {len(df)} 个分子")
                return df
            except Exception as e2:
                print(f"❌ 方法2也失败: {e2}")
                return None

    def extract_mol_properties(self, mol):
        """从分子对象中提取所有属性"""
        props = {}

        # 提取内置属性
        if mol.HasProp("_Name"):
            props['Name'] = mol.GetProp("_Name")

        # 提取所有自定义属性
        prop_names = mol.GetPropNames()
        for prop in prop_names:
            if prop != "_Name":  # 已经单独处理
                try:
                    value = mol.GetProp(prop)
                    # 清理属性名（Excel兼容）
                    clean_prop = self.clean_column_name(prop)
                    props[clean_prop] = value
                except:
                    props[prop] = ""

        return props

    def clean_column_name(self, column_name):
        """清理列名以兼容Excel"""
        # 移除非法字符
        cleaned = re.sub(r'[\[\]{}()<>+=!@#$%^&*|\\/~`]', '', column_name)
        # 替换空格和下划线
        cleaned = cleaned.replace(' ', '_').replace('-', '_')
        # 确保以字母开头
        if cleaned and not cleaned[0].isalpha():
            cleaned = 'C_' + cleaned
        # 限制长度
        if len(cleaned) > 50:
            cleaned = cleaned[:50]

        return cleaned

    def add_calculated_properties(self, df, progress_bar=None):
        """添加计算属性到DataFrame"""
        if 'ROMol' not in df.columns:
            return df

        print("正在计算分子属性...")

        # 初始化新列
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
            'TPSA': [],  # 拓扑极性表面积
            'LogP': [],  # 辛醇/水分配系数
            'RingCount': []  # 环计数
        }

        total_mols = len(df)
        if progress_bar is None:
            progress_bar = CustomProgressBar(total_mols, "计算属性")

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

                # 分子量
                try:
                    new_columns['MolecularWeight'].append(Descriptors.MolWt(mol))
                except:
                    new_columns['MolecularWeight'].append(0)

                # 分子式
                try:
                    new_columns['MolecularFormula'].append(Chem.rdMolDescriptors.CalcMolFormula(mol))
                except:
                    new_columns['MolecularFormula'].append('')

                # 原子数等
                new_columns['NumHeavyAtoms'].append(mol.GetNumHeavyAtoms())
                new_columns['NumAtoms'].append(mol.GetNumAtoms())
                new_columns['NumBonds'].append(mol.GetNumBonds())

                # 描述符
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
                # 如果分子无效，填充默认值
                for col in new_columns:
                    if col in ['SMILES', 'InChI', 'MolecularFormula']:
                        new_columns[col].append('')
                    else:
                        new_columns[col].append(0)

            progress_bar.update(1)

        progress_bar.close()

        # 添加新列到DataFrame
        for col_name, values in new_columns.items():
            df[col_name] = values

        return df

    def fix_dataframe_for_excel(self, df):
        """修复DataFrame以兼容Excel格式"""
        print("正在修复数据格式以兼容Excel...")

        # 创建副本
        df_fixed = df.copy()

        # 移除ROMol列（Excel无法存储）
        if 'ROMol' in df_fixed.columns:
            df_fixed = df_fixed.drop(columns=['ROMol'])

        # 确保所有列名都是字符串且合法
        df_fixed.columns = [self.clean_column_name(str(col)) for col in df_fixed.columns]

        # 处理数据格式
        for col in df_fixed.columns:
            # 转换非标量类型为字符串
            if df_fixed[col].dtype == 'object':
                try:
                    # 检查是否包含非标量对象
                    sample = df_fixed[col].dropna().iloc[0] if not df_fixed[col].dropna().empty else None
                    if sample is not None and not isinstance(sample, (str, int, float, bool)):
                        df_fixed[col] = df_fixed[col].astype(str)
                except:
                    df_fixed[col] = df_fixed[col].astype(str)

            # 处理NaN值
            df_fixed[col] = df_fixed[col].fillna('')

        # 限制行数（Excel最大行数为1048576）
        if len(df_fixed) > 1000000:
            print(f"⚠ 警告: 数据行数({len(df_fixed)})超过Excel推荐限制，将截断")
            df_fixed = df_fixed.head(1000000)

        # 限制列数（Excel最大列数为16384）
        if len(df_fixed.columns) > 1000:
            print(f"⚠ 警告: 数据列数({len(df_fixed.columns)})较多，保留前1000列")
            df_fixed = df_fixed.iloc[:, :1000]

        print(f"✓ 数据格式修复完成: {len(df_fixed)}行 × {len(df_fixed.columns)}列")
        return df_fixed

    def export_to_csv(self, df, output_path, encoding='utf-8-sig'):
        """导出为CSV文件"""
        try:
            print(f"导出为CSV: {output_path.name}")

            # 修复数据格式
            df_fixed = self.fix_dataframe_for_excel(df)

            # 导出CSV
            df_fixed.to_csv(output_path, index=False, encoding=encoding)

            print(f"✓ CSV文件已保存: {output_path}")
            print(f"  文件大小: {os.path.getsize(output_path) / 1024:.1f} KB")

            return True
        except Exception as e:
            print(f"❌ 导出CSV失败: {e}")
            return False

    def export_to_excel(self, df, output_path, engine='openpyxl'):
        """导出为Excel文件（已修正格式）"""
        try:
            print(f"导出为Excel: {output_path.name}")
            print(f"使用引擎: {engine}")

            # 修复数据格式
            df_fixed = self.fix_dataframe_for_excel(df)

            # 导出Excel
            with pd.ExcelWriter(output_path, engine=engine) as writer:
                df_fixed.to_excel(writer, sheet_name='Molecules', index=False)

                # 如果是openpyxl引擎，可以设置列宽
                if engine == 'openpyxl':
                    worksheet = writer.sheets['Molecules']

                    # 自动调整列宽
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

                        adjusted_width = min(max_length + 2, 50)  # 最大宽度50
                        worksheet.column_dimensions[column_letter].width = adjusted_width

            print(f"✓ Excel文件已保存: {output_path}")
            print(f"  文件大小: {os.path.getsize(output_path) / 1024:.1f} KB")

            return True
        except Exception as e:
            print(f"❌ 导出Excel失败: {e}")

            # 尝试备用引擎
            if engine == 'openpyxl':
                print("尝试使用xlsxwriter引擎...")
                return self.export_to_excel(df, output_path, engine='xlsxwriter')
            else:
                print("尝试使用openpyxl引擎...")
                return self.export_to_excel(df, output_path, engine='openpyxl')

    def convert_sdf_file(self, sdf_file, options):
        """转换单个SDF文件"""
        print(f"\n" + "="*60)
        print(f"处理文件: {sdf_file.name}")
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
            # 1. 读取SDF文件
            df = self.read_sdf_with_progress(sdf_file)

            if df is None or len(df) == 0:
                result['error'] = "无法读取SDF文件或无有效分子"
                result['time_taken'] = time.time() - start_time
                return result

            result['molecules_processed'] = len(df)

            # 2. 添加计算属性
            if options['add_calculated_props']:
                progress_bar = CustomProgressBar(len(df), "计算分子属性")
                df = self.add_calculated_properties(df, progress_bar)

            result['molecules_successful'] = len(df)

            # 3. 准备输出文件名
            if options['prefix']:
                base_name = f"{options['prefix']}_{sdf_file.stem}"
            else:
                base_name = sdf_file.stem

            # 4. 导出为指定格式
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
            print(f"❌ 转换失败: {e}")
            traceback.print_exc()

        result['time_taken'] = time.time() - start_time
        return result

    def generate_report(self, results, options):
        """生成转换报告"""
        print("\n" + "="*60)
        print("转换报告")
        print("="*60)

        total_files = len(results)
        successful_files = sum(1 for r in results if r['success'])
        total_molecules = sum(r['molecules_processed'] for r in results)
        successful_molecules = sum(r['molecules_successful'] for r in results)

        # 输出摘要
        print(f"\n📊 转换摘要:")
        print(f"   处理文件数: {total_files}")
        print(f"   成功文件数: {successful_files}")
        print(f"   处理分子总数: {total_molecules}")
        print(f"   成功分子数: {successful_molecules}")

        if total_molecules > 0:
            success_rate = (successful_molecules / total_molecules) * 100
            print(f"   成功率: {success_rate:.1f}%")

        total_time = sum(r['time_taken'] for r in results)
        print(f"   总耗时: {total_time:.1f}秒")

        if successful_molecules > 0:
            print(f"   平均速度: {successful_molecules/total_time:.1f} 分子/秒")

        # 输出文件列表
        print(f"\n📁 输出文件 ({options['output_dir']}):")
        all_output_files = []
        for result in results:
            if result['success']:
                for file in result['output_files']:
                    if file.exists():
                        file_size = file.stat().st_size / 1024
                        print(f"   ✓ {file.name} ({file_size:.1f} KB)")
                        all_output_files.append(file)

        # 生成详细报告文件
        report_file = options['output_dir'] / f"{options['prefix']}_report.txt"
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("SDF文件转换报告\n")
            f.write("="*50 + "\n")
            f.write(f"生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"程序版本: {self.version}\n")
            f.write(f"输出目录: {options['output_dir']}\n")
            f.write(f"输出格式: {', '.join(options['formats'])}\n")
            f.write(f"添加计算属性: {'是' if options['add_calculated_props'] else '否'}\n")
            f.write("\n转换结果:\n")

            for i, result in enumerate(results, 1):
                f.write(f"\n[{i}] 文件: {result['file']}\n")
                f.write(f"   状态: {'成功' if result['success'] else '失败'}\n")
                f.write(f"   处理分子数: {result['molecules_processed']}\n")
                f.write(f"   成功分子数: {result['molecules_successful']}\n")
                f.write(f"   耗时: {result.get('time_taken', 0):.2f}秒\n")

                if result['output_files']:
                    f.write("   输出文件:\n")
                    for out_file in result['output_files']:
                        file_size = out_file.stat().st_size / 1024 if out_file.exists() else 0
                        f.write(f"     - {out_file.name} ({file_size:.1f} KB)\n")

                if result['error']:
                    f.write(f"   错误: {result['error']}\n")

        print(f"\n📝 详细报告已保存: {report_file.name}")

        # 记录日志
        self.log_conversion(results, options)

        return report_file

    def log_conversion(self, results, options):
        """记录转换日志"""
        try:
            with open(self.log_file, 'a', encoding='utf-8') as f:
                timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                f.write(f"\n{'='*60}\n")
                f.write(f"转换时间: {timestamp}\n")
                f.write(f"输出目录: {options['output_dir']}\n")

                for result in results:
                    status = "成功" if result['success'] else "失败"
                    f.write(f"{result['file']}: {status} ({result['molecules_successful']}/{result['molecules_processed']})\n")
        except:
            pass

    def show_data_preview(self, output_dir):
        """显示输出文件的预览"""
        print("\n" + "="*60)
        print("数据预览")
        print("="*60)

        # 查找最近生成的文件
        csv_files = list(output_dir.glob("*.csv"))
        excel_files = list(output_dir.glob("*.xlsx"))

        preview_files = (csv_files[:2] + excel_files[:1])[:3]  # 最多预览3个文件

        if not preview_files:
            print("没有找到输出文件进行预览")
            return

        for file in preview_files:
            try:
                print(f"\n📋 {file.name}:")
                print("-"*40)

                if file.suffix == '.csv':
                    # 读取前5行
                    df_preview = pd.read_csv(file, nrows=5)
                else:
                    # 读取Excel前5行
                    df_preview = pd.read_excel(file, nrows=5)

                # 显示基本信息
                print(f"总行数: {len(pd.read_csv(file) if file.suffix == '.csv' else pd.read_excel(file))}")
                print(f"列数: {len(df_preview.columns)}")

                # 显示前几行
                print("\n前5行数据:")
                pd.set_option('display.max_columns', 10)  # 最多显示10列
                pd.set_option('display.width', 120)
                print(df_preview)

                # 显示列名
                print(f"\n前10个列名:")
                for i, col in enumerate(df_preview.columns[:10], 1):
                    print(f"  {i:2d}. {col}")

                if len(df_preview.columns) > 10:
                    print(f"  ... 还有 {len(df_preview.columns) - 10} 列")

            except Exception as e:
                print(f"无法预览 {file.name}: {e}")

    def open_output_directory(self, output_dir):
        """打开输出目录"""
        try:
            if sys.platform == 'win32':
                os.startfile(output_dir)
            elif sys.platform == 'darwin':
                os.system(f'open "{output_dir}"')
            elif sys.platform == 'linux':
                os.system(f'xdg-open "{output_dir}"')
            print(f"✓ 已打开输出目录")
        except:
            print(f"⚠ 无法自动打开目录，请手动访问: {output_dir}")

    def main(self):
        """主程序"""
        try:
            # 1. 显示横幅
            self.print_banner()

            # 2. 检查依赖
            self.check_dependencies()

            # 3. 获取输入文件
            input_files = self.get_input_file()

            if not input_files:
                print("没有选择文件，程序退出。")
                return

            # 4. 获取输出选项
            options = self.get_output_options(input_files)

            # 5. 处理文件
            print("\n" + "="*60)
            print("步骤 3: 正在转换文件")
            print("="*60)

            results = []
            total_start_time = time.time()

            for i, sdf_file in enumerate(input_files, 1):
                print(f"\n[{i}/{len(input_files)}] ", end="")
                result = self.convert_sdf_file(sdf_file, options)
                results.append(result)

                # 显示进度
                progress = (i / len(input_files)) * 100
                elapsed = time.time() - total_start_time
                print(f"总体进度: {progress:.1f}% | 已用时间: {elapsed:.1f}秒")

            # 6. 生成报告
            report_file = self.generate_report(results, options)

            # 7. 显示预览
            self.show_data_preview(options['output_dir'])

            # 8. 完成提示
            print("\n" + "="*60)
            print("🎉 转换完成!")
            print("="*60)

            total_time = time.time() - total_start_time
            print(f"\n总耗时: {total_time:.1f}秒")
            print(f"输出目录: {options['output_dir']}")
            print(f"详细报告: {report_file.name}")

            # 询问是否打开输出目录
            print("\n是否打开输出目录? (y/n): ", end="")
            if input().strip().lower() == 'y':
                self.open_output_directory(options['output_dir'])

            print("\n感谢使用 SDF 文件转换工具!")
            print("按回车键退出...")
            input()

        except KeyboardInterrupt:
            print("\n\n程序被用户中断。")
        except Exception as e:
            print(f"\n❌ 程序运行时发生错误: {e}")
            print("\n错误详情:")
            traceback.print_exc()
            print("\n按回车键退出...")
            input()

# 快速转换模式
def quick_convert_mode():
    """快速转换模式 - 简化流程"""
    print("快速 SDF 转换模式")
    print("=" * 50)

    converter = InteractiveSDFConverter()

    # 获取文件
    print("\n请拖拽SDF文件到此处或输入文件路径:")
    file_path = input("文件: ").strip().strip('"').strip("'")

    if not file_path:
        print("未输入文件路径")
        return

    path = Path(file_path)
    if not path.exists():
        print(f"文件不存在: {file_path}")
        return

    # 选择格式
    print("\n选择输出格式:")
    print("1. CSV")
    print("2. Excel")
    print("3. 两种都要")
    choice = input("选择 (1-3): ").strip()

    formats = []
    if choice in ['1', '3']:
        formats.append('csv')
    if choice in ['2', '3']:
        formats.append('excel')

    if not formats:
        formats = ['csv']

    # 执行转换
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
        print(f"\n✅ 转换成功!")
        for out_file in result['output_files']:
            print(f"   ✓ {out_file.name}")

        # 打开目录
        print("\n是否打开输出目录? (y/n): ", end="")
        if input().strip().lower() == 'y':
            converter.open_output_directory(converter.output_dir)
    else:
        print(f"\n❌ 转换失败: {result['error']}")

if __name__ == "__main__":
    # 解析命令行参数
    if len(sys.argv) > 1:
        if sys.argv[1] in ['--quick', '-q']:
            quick_convert_mode()
        elif sys.argv[1] in ['--help', '-h']:
            print("SDF 文件转换工具 v2.0")
            print("\n使用方法:")
            print("  python sdf_converter.py          # 启动交互式界面")
            print("  python sdf_converter.py --quick  # 快速转换模式")
            print("  python sdf_converter.py --help   # 显示帮助")
            print("\n功能:")
            print("  - 支持SDF转CSV和Excel格式")
            print("  - 自动计算分子属性和描述符")
            print("  - 实时进度条显示")
            print("  - 智能格式修正和Excel兼容性处理")
        else:
            # 如果提供了文件路径，直接处理
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
                    print(f"\n✅ 转换成功!")
                    for out_file in result['output_files']:
                        print(f"   ✓ {out_file.name}")
                else:
                    print(f"\n❌ 转换失败: {result['error']}")
            else:
                print(f"文件不存在: {sys.argv[1]}")
    else:
        # 默认启动交互式界面
        converter = InteractiveSDFConverter()
        converter.main()