import pandas as pd
import os

# 设置文件路径
input_file = "C:/Users/姚钱磊/Desktop/补充实验的数据/output-1/20251205-YanZhengJi-TiaoJian10/PPG指数_20251205-YanZhengJi-TiaoJian10_20260124_235733.xlsx"

# 1. 首先查看Excel文件中有哪些工作表
try:
    # 使用ExcelFile对象查看工作表名称
    excel_file = pd.ExcelFile(input_file)
    print("Excel文件中的工作表名称:")
    for sheet in excel_file.sheet_names:
        print(f"  - {sheet}")

    # 2. 猜测正确的工作表名称
    # 根据常见的名称猜测
    possible_sheet_names = [
        "所有匹配结果", "成功匹配", "Sheet1", "Sheet2", "Sheet3",
        "匹配结果", "结果", "数据", "原始数据", "化合物匹配"
    ]

    print("\n尝试查找匹配的工作表...")
    for sheet_name in possible_sheet_names:
        if sheet_name in excel_file.sheet_names:
            print(f"找到匹配的工作表: {sheet_name}")
            actual_sheet_name = sheet_name
            break
    else:
        # 如果没有找到匹配的，使用第一个工作表
        actual_sheet_name = excel_file.sheet_names[0]
        print(f"未找到完全匹配的工作表，将使用第一个工作表: {actual_sheet_name}")

    # 3. 读取数据
    df = pd.read_excel(input_file, sheet_name=actual_sheet_name)
    print(f"\n成功读取工作表 '{actual_sheet_name}'")
    print(f"数据形状: {df.shape} (行×列)")
    print(f"数据列名: {list(df.columns[:10])}...")  # 显示前10列

    # 4. 提取每个化合物的第一个匹配结果
    # 先检查是否有"匹配排名"列
    if "匹配排名" not in df.columns:
        print("\n警告: 数据中没有'匹配排名'列，尝试查找相关列...")
        # 尝试查找包含"匹配"或"排名"的列
        rank_columns = [col for col in df.columns if "匹配" in col or "排名" in col or "rank" in col.lower()]
        if rank_columns:
            rank_col = rank_columns[0]
            print(f"使用 '{rank_col}' 作为排名列")
        else:
            # 如果没有排名列，假设所有行都是第一个匹配
            print("未找到排名列，将提取所有行")
            rank_col = None
    else:
        rank_col = "匹配排名"

    # 5. 提取数据
    if rank_col:
        df_first = df[df[rank_col] == 1]
    else:
        df_first = df

    print(f"\n提取到 {len(df_first)} 个化合物的第一个匹配结果")

    # 6. 确定要提取的列
    # 查找包含"m/z"的列
    mz_columns = [col for col in df_first.columns if "m/z" in col]
    rt_columns = [col for col in df_first.columns if "保留" in col or "RT" in col or "retention" in col.lower()]
    ppg_columns = [col for col in df_first.columns if "PPG" in col or "指数" in col or "index" in col.lower()]

    # 选择需要的列
    selected_columns = []

    # 添加化合物名称列（如果有）
    name_cols = [col for col in df_first.columns if "化合物" in col or "名称" in col or "name" in col.lower()]
    if name_cols:
        selected_columns.append(name_cols[0])
    else:
        print("警告: 未找到化合物名称列")

    # 添加理论m/z列
    theoretical_mz = [col for col in mz_columns if "理论" in col or "theor" in col.lower()]
    if theoretical_mz:
        selected_columns.append(theoretical_mz[0])
    elif mz_columns:
        selected_columns.append(mz_columns[0])

    # 添加实测m/z列
    measured_mz = [col for col in mz_columns if "实测" in col or "meas" in col.lower() or "检测" in col]
    if measured_mz:
        selected_columns.append(measured_mz[0])
    elif len(mz_columns) > 1:
        selected_columns.append(mz_columns[1])

    # 添加保留时间列
    if rt_columns:
        selected_columns.append(rt_columns[0])

    # 添加PPG指数列
    if ppg_columns:
        selected_columns.append(ppg_columns[0])

    # 确保列名不重复
    selected_columns = list(dict.fromkeys(selected_columns))

    print(f"\n将提取以下列: {selected_columns}")

    # 7. 提取数据
    df_result = df_first[selected_columns].copy()

    # 8. 重命名列使其更易懂
    column_mapping = {}
    for col in df_result.columns:
        if "理论" in col and "m/z" in col:
            column_mapping[col] = "理论mz"
        elif "实测" in col and "m/z" in col:
            column_mapping[col] = "实测mz"
        elif any(x in col for x in ["保留时间", "RT", "retention"]):
            column_mapping[col] = "保留时间"
        elif any(x in col for x in ["PPG", "指数", "index"]):
            column_mapping[col] = "保留指数"
        elif any(x in col for x in ["化合物", "名称", "name"]):
            column_mapping[col] = "化合物名称"

    if column_mapping:
        df_result = df_result.rename(columns=column_mapping)

    # 9. 保存到文件
    output_file =  "C:/Users/姚钱磊/Desktop/补充实验预测/训练/验证集_条件11.xlsx"
    df_result.to_excel(output_file, index=False)

    print(f"\n处理完成！共提取 {len(df_result)} 个化合物")
    print(f"结果保存到: {os.path.abspath(output_file)}")
    print("\n前10个结果:")
    print(df_result.head(10).to_string(index=False))

except FileNotFoundError:
    print(f"错误: 找不到文件 {input_file}")
    print(f"当前工作目录: {os.getcwd()}")
    print("请检查文件路径是否正确")
except Exception as e:
    print(f"处理文件时发生错误: {e}")
    import traceback

    traceback.print_exc()