import pandas as pd
import os
from pathlib import Path


def excel_to_csv(input_file, output_file=None, sheet_name=0, encoding='utf-8'):
    """
    Convert an Excel file to CSV format

    Parameters:
    input_file: Input Excel file path
    output_file: Output CSV file path (optional; defaults to the input file name)
    sheet_name: Worksheet name or index to convert (defaults to the first worksheet)
    encoding: CSV file encoding (defaults to utf-8)

    Returns:
    Whether conversion was successful
    """
    try:
        # Check whether the input file exists
        if not os.path.exists(input_file):
            print(f"Error: input file does not exist - {input_file}")
            return False

        # Determine the output file name
        if output_file is None:
            # Use the input file name with the extension changed to .csv
            output_file = str(Path(input_file).with_suffix('.csv'))

        # Ensure that the output directory exists
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        print(f"Reading Excel file: {input_file}")

        # Read Excel file
        # Try different engines to support .xls and .xlsx formats
        try:
            df = pd.read_excel(input_file, sheet_name=sheet_name, engine='openpyxl')
        except:
            try:
                df = pd.read_excel(input_file, sheet_name=sheet_name, engine='xlrd')
            except:
                df = pd.read_excel(input_file, sheet_name=sheet_name)

        # Get the worksheet name when an index is used
        if isinstance(sheet_name, int):
            if sheet_name == 0:
                sheet_name = list(pd.ExcelFile(input_file).sheet_names)[0]

        print(f" '{sheet_name}' has been read, with {len(df)} , {len(df.columns)} column")

        # Save as CSV
        df.to_csv(output_file, index=False, encoding=encoding)

        print(f"CSV file saved: {output_file}")
        print(f"File size: {os.path.getsize(output_file)} bytes")

        # data
        print("\nPreview of the first five rows:")
        print(df.head().to_string())

        return True

    except Exception as e:
        print(f"An error occurred during conversion: {e}")
        return False


def batch_excel_to_csv(input_dir, output_dir=None, pattern="*.xls*", recursive=False):
    """
    Batch-convert Excel files to CSV format

    Parameters:
    input_dir: Input directory path
    output_dir: Output directory path (, defaults to the input directory)
    pattern: File matching pattern (defaults to all Excel files)
    recursive: Whether to process subdirectories recursively
    """
    if not os.path.exists(input_dir):
        print(f": inputdirectory - {input_dir}")
        return

    if output_dir is None:
        output_dir = input_dir

    # Ensure that the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Excelfile
    import glob
    if recursive:
        excel_files = glob.glob(os.path.join(input_dir, "**", pattern), recursive=True)
    else:
        excel_files = glob.glob(os.path.join(input_dir, pattern))

    if not excel_files:
        print(f" {input_dir} match {pattern} Excelfile")
        return

    print(f" {len(excel_files)} Excelfile")

    success_count = 0
    for excel_file in excel_files:
        print(f"\nfile: {excel_file}")

        # outputpath
        rel_path = os.path.relpath(excel_file, input_dir) if recursive else os.path.basename(excel_file)
        csv_file = os.path.join(output_dir, os.path.splitext(rel_path)[0] + ".csv")

        # outputdirectory
        csv_dir = os.path.dirname(csv_file)
        if csv_dir and not os.path.exists(csv_dir):
            os.makedirs(csv_dir)

        # file
        if excel_to_csv(excel_file, csv_file):
            success_count += 1

    print(f"\n！ {success_count}/{len(excel_files)} file")


def main():
    """: """
    print("=" * 50)
    print("ExcelCSV")
    print("=" * 50)
    print(":")
    print("1. file")
    print("2. directoryfile")
    print("3. ")

    choice = input("\ninput (1-3): ").strip()

    if choice == "1":
        # file
        input_file = input("inputExcelfilepath: ").strip()

        # file
        if not os.path.exists(input_file):
            print(f"file does not exist: {input_file}")
            return

        # outputfile
        default_output = str(Path(input_file).with_suffix('.csv'))
        output_file = input(f"inputoutputCSVfilepath (Use '{default_output}'): ").strip()
        if not output_file:
            output_file = default_output

        # translated note
        sheet_input = input("inputname (Use): ").strip()
        if sheet_input:
            try:
                # ()
                sheet_name = int(sheet_input)
            except ValueError:
                # (name)
                sheet_name = sheet_input
        else:
            sheet_name = 0

        # translated note
        encoding = input("input (Useutf-8): ").strip()
        if not encoding:
            encoding = 'utf-8'

        # translated note
        excel_to_csv(input_file, output_file, sheet_name, encoding)

    elif choice == "2":
        # translated note
        input_dir = input("inputInput directory path: ").strip()

        if not os.path.exists(input_dir):
            print(f"directory: {input_dir}")
            return

        output_dir = input("inputOutput directory path (Useinputdirectory): ").strip()
        if not output_dir:
            output_dir = input_dir

        pattern = input("inputFile matching pattern (Use '*.xls*'): ").strip()
        if not pattern:
            pattern = "*.xls*"

        recursive_input = input("Whether to process subdirectories recursively? (y/n, ): ").strip().lower()
        recursive = recursive_input in ['y', 'yes', '']

        batch_excel_to_csv(input_dir, output_dir, pattern, recursive)

    elif choice == "3":
        print("")
        return
    else:
        print("")


# ()
def quick_convert(input_file, output_file=None):
    """
    ExcelfileCSV
    """
    return excel_to_csv(input_file, output_file)


if __name__ == "__main__":
    # file, filepath
    # quick_convert("validation set_condition1.xlsx", "validation set_condition1.csv")

    # translated note
    main()