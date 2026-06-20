import pandas as pd
import os

# Set file paths. Replace this example input with the exported PPG-index workbook
# generated for the chromatographic condition being checked.
input_file = "../training_data/example_PPG_index_input.xlsx"

# 1. First inspect the worksheet names in the Excel file
try:
    # Use an ExcelFile object to inspect worksheet names
    excel_file = pd.ExcelFile(input_file)
    print("Worksheet names in the Excel file:")
    for sheet in excel_file.sheet_names:
        print(f"  - {sheet}")

    # 2. Guess the correct worksheet name
    # Guess based on common worksheet names
    possible_sheet_names = [
        "all matching results", "successful matches", "Sheet1", "Sheet2", "Sheet3",
        "matching results", "results", "data", "raw data", "compound matching"
    ]

    print("\nTry to find a matching worksheet...")
    for sheet_name in possible_sheet_names:
        if sheet_name in excel_file.sheet_names:
            print(f"Found a matching worksheet: {sheet_name}")
            actual_sheet_name = sheet_name
            break
    else:
        # match, Use
        actual_sheet_name = excel_file.sheet_names[0]
        print(f"No exact worksheet match was found; the first worksheet will be used: {actual_sheet_name}")

    # 3. Read data
    df = pd.read_excel(input_file, sheet_name=actual_sheet_name)
    print(f"\nSuccessfully read worksheet '{actual_sheet_name}'")
    print(f"Data shape: {df.shape} (rows x columns)")
    print(f"Data column names: {list(df.columns[:10])}...") # show the first 10 columns

    # 4. Extract the first matching result for each compound
    # First check whether there is"match_rank"column
    if "match_rank" not in df.columns:
        print("\nWarning: data'match_rank'column, try to find related columns...")
        # try to find columns containing"match""rank"column
        rank_columns = [col for col in df.columns if "match" in col or "rank" in col or "rank" in col.lower()]
        if rank_columns:
            rank_col = rank_columns[0]
            print(f"Use '{rank_col}' as the ranking column")
        else:
            # If no ranking column is present, assume all rows are first matches
            print("No ranking column was found; all rows will be extracted")
            rank_col = None
    else:
        rank_col = "match_rank"

    # 5. Extract data
    if rank_col:
        df_first = df[df[rank_col] == 1]
    else:
        df_first = df

    print(f"\nExtracted {len(df_first)} first matching results for compounds")

    # 6. Determine the columns to extract
    # "m/z"column
    mz_columns = [col for col in df_first.columns if "m/z" in col]
    rt_columns = [col for col in df_first.columns if "retention" in col or "RT" in col or "retention" in col.lower()]
    ppg_columns = [col for col in df_first.columns if "PPG" in col or "index" in col or "index" in col.lower()]

    # column
    selected_columns = []

    # compound_namecolumn ()
    name_cols = [col for col in df_first.columns if "compound" in col or "name" in col or "name" in col.lower()]
    if name_cols:
        selected_columns.append(name_cols[0])
    else:
        print("Warning: compound_namecolumn")

    # m/zcolumn
    theoretical_mz = [col for col in mz_columns if "" in col or "theor" in col.lower()]
    if theoretical_mz:
        selected_columns.append(theoretical_mz[0])
    elif mz_columns:
        selected_columns.append(mz_columns[0])

    # m/zcolumn
    measured_mz = [col for col in mz_columns if "" in col or "meas" in col.lower() or "" in col]
    if measured_mz:
        selected_columns.append(measured_mz[0])
    elif len(mz_columns) > 1:
        selected_columns.append(mz_columns[1])

    # retention_timecolumn
    if rt_columns:
        selected_columns.append(rt_columns[0])

    # PPGindexcolumn
    if ppg_columns:
        selected_columns.append(ppg_columns[0])

    # column
    selected_columns = list(dict.fromkeys(selected_columns))

    print(f"\ncolumn: {selected_columns}")

    # 7. Extract data
    df_result = df_first[selected_columns].copy()

    # 8. column
    column_mapping = {}
    for col in df_result.columns:
        if "" in col and "m/z" in col:
            column_mapping[col] = "theoretical_mz"
        elif "" in col and "m/z" in col:
            column_mapping[col] = "measured_mz"
        elif any(x in col for x in ["retention_time", "RT", "retention"]):
            column_mapping[col] = "retention_time"
        elif any(x in col for x in ["PPG", "index", "index"]):
            column_mapping[col] = "retention_index"
        elif any(x in col for x in ["compound", "name", "name"]):
            column_mapping[col] = "compound_name"

    if column_mapping:
        df_result = df_result.rename(columns=column_mapping)

    # 9. Save the extracted validation set.
    output_file = "../training_data/validation_set_condition11.xlsx"
    df_result.to_excel(output_file, index=False)

    print(f"\n！ {len(df_result)} compound")
    print(f"resultssave: {os.path.abspath(output_file)}")
    print("\n10results:")
    print(df_result.head(10).to_string(index=False))

except FileNotFoundError:
    print(f": file {input_file}")
    print(f"directory: {os.getcwd()}")
    print("filepath")
except Exception as e:
    print(f"file: {e}")
    import traceback

    traceback.print_exc()
