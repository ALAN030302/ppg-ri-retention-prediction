import pandas as pd
import os
from pathlib import Path


def excel_to_csv(input_file, output_file=None, sheet_name=0, encoding='utf-8'):
    """
    将Excel文件转换为CSV格式

    参数:
    input_file: 输入的Excel文件路径
    output_file: 输出的CSV文件路径（可选，默认与输入文件同名）
    sheet_name: 要转换的工作表名称或索引（默认为第一个工作表）
    encoding: CSV文件的编码格式（默认为utf-8）

    返回:
    转换是否成功
    """
    try:
        # 检查输入文件是否存在
        if not os.path.exists(input_file):
            print(f"错误：输入文件不存在 - {input_file}")
            return False

        # 确定输出文件名
        if output_file is None:
            # 使用输入文件的相同名称，但扩展名改为.csv
            output_file = str(Path(input_file).with_suffix('.csv'))

        # 确保输出目录存在
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        print(f"正在读取Excel文件: {input_file}")

        # 读取Excel文件
        # 尝试不同的引擎以支持.xls和.xlsx格式
        try:
            df = pd.read_excel(input_file, sheet_name=sheet_name, engine='openpyxl')
        except:
            try:
                df = pd.read_excel(input_file, sheet_name=sheet_name, engine='xlrd')
            except:
                df = pd.read_excel(input_file, sheet_name=sheet_name)

        # 获取工作表名称（如果使用索引）
        if isinstance(sheet_name, int):
            if sheet_name == 0:
                sheet_name = list(pd.ExcelFile(input_file).sheet_names)[0]

        print(f"工作表 '{sheet_name}' 已读取，共 {len(df)} 行，{len(df.columns)} 列")

        # 保存为CSV
        df.to_csv(output_file, index=False, encoding=encoding)

        print(f"CSV文件已保存: {output_file}")
        print(f"文件大小: {os.path.getsize(output_file)} 字节")

        # 显示前几行数据
        print("\n前5行数据预览:")
        print(df.head().to_string())

        return True

    except Exception as e:
        print(f"转换过程中发生错误: {e}")
        return False


def batch_excel_to_csv(input_dir, output_dir=None, pattern="*.xls*", recursive=False):
    """
    批量转换Excel文件为CSV格式

    参数:
    input_dir: 输入目录路径
    output_dir: 输出目录路径（可选，默认与输入目录相同）
    pattern: 文件匹配模式（默认为所有Excel文件）
    recursive: 是否递归处理子目录
    """
    if not os.path.exists(input_dir):
        print(f"错误：输入目录不存在 - {input_dir}")
        return

    if output_dir is None:
        output_dir = input_dir

    # 确保输出目录存在
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 查找Excel文件
    import glob
    if recursive:
        excel_files = glob.glob(os.path.join(input_dir, "**", pattern), recursive=True)
    else:
        excel_files = glob.glob(os.path.join(input_dir, pattern))

    if not excel_files:
        print(f"在 {input_dir} 中未找到匹配 {pattern} 的Excel文件")
        return

    print(f"找到 {len(excel_files)} 个Excel文件")

    success_count = 0
    for excel_file in excel_files:
        print(f"\n处理文件: {excel_file}")

        # 构建输出路径
        rel_path = os.path.relpath(excel_file, input_dir) if recursive else os.path.basename(excel_file)
        csv_file = os.path.join(output_dir, os.path.splitext(rel_path)[0] + ".csv")

        # 确保输出子目录存在
        csv_dir = os.path.dirname(csv_file)
        if csv_dir and not os.path.exists(csv_dir):
            os.makedirs(csv_dir)

        # 转换文件
        if excel_to_csv(excel_file, csv_file):
            success_count += 1

    print(f"\n批量转换完成！成功转换 {success_count}/{len(excel_files)} 个文件")


def main():
    """主函数：提供交互式界面"""
    print("=" * 50)
    print("Excel转CSV工具")
    print("=" * 50)
    print("请选择操作模式:")
    print("1. 转换单个文件")
    print("2. 批量转换目录中的文件")
    print("3. 退出")

    choice = input("\n请输入选项 (1-3): ").strip()

    if choice == "1":
        # 单个文件转换
        input_file = input("请输入Excel文件路径: ").strip()

        # 检查文件是否存在
        if not os.path.exists(input_file):
            print(f"文件不存在: {input_file}")
            return

        # 获取默认输出文件名
        default_output = str(Path(input_file).with_suffix('.csv'))
        output_file = input(f"请输入输出CSV文件路径 (回车使用默认值 '{default_output}'): ").strip()
        if not output_file:
            output_file = default_output

        # 选择工作表
        sheet_input = input("请输入工作表名称或索引 (回车使用第一个工作表): ").strip()
        if sheet_input:
            try:
                # 尝试转换为整数（索引）
                sheet_name = int(sheet_input)
            except ValueError:
                # 否则作为字符串处理（工作表名称）
                sheet_name = sheet_input
        else:
            sheet_name = 0

        # 选择编码
        encoding = input("请输入编码格式 (回车使用默认utf-8): ").strip()
        if not encoding:
            encoding = 'utf-8'

        # 执行转换
        excel_to_csv(input_file, output_file, sheet_name, encoding)

    elif choice == "2":
        # 批量转换
        input_dir = input("请输入输入目录路径: ").strip()

        if not os.path.exists(input_dir):
            print(f"目录不存在: {input_dir}")
            return

        output_dir = input("请输入输出目录路径 (回车使用输入目录): ").strip()
        if not output_dir:
            output_dir = input_dir

        pattern = input("请输入文件匹配模式 (回车使用默认 '*.xls*'): ").strip()
        if not pattern:
            pattern = "*.xls*"

        recursive_input = input("是否递归处理子目录? (y/n, 回车默认否): ").strip().lower()
        recursive = recursive_input in ['y', 'yes', '是']

        batch_excel_to_csv(input_dir, output_dir, pattern, recursive)

    elif choice == "3":
        print("退出程序")
        return
    else:
        print("无效的选项")


# 快速转换函数（无需交互）
def quick_convert(input_file, output_file=None):
    """
    快速转换单个Excel文件为CSV
    """
    return excel_to_csv(input_file, output_file)


if __name__ == "__main__":
    # 如果您想直接转换特定文件，可以取消注释下面的代码并修改文件路径
    # quick_convert("验证集_条件1.xlsx", "验证集_条件1.csv")

    # 运行交互式界面
    main()