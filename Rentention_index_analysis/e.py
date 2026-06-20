import pandas as pd
import numpy as np
from pubchempy import get_compounds, Compound
import requests
from time import sleep
import warnings

warnings.filterwarnings('ignore')


def load_known_compounds(known_file_path):
    """
    Load a known-compound list containing CAS numbers and SMILES
    Assume a CSV file with columns: cas, smiles, and optional mz
    """
    try:
        # Select the reading method according to file extension
        if known_file_path.endswith('.csv'):
            known_df = pd.read_csv(known_file_path)
        elif known_file_path.endswith('.xlsx') or known_file_path.endswith('.xls'):
            known_df = pd.read_excel(known_file_path)
        else:
            # Try CSV format
            known_df = pd.read_csv(known_file_path, sep='\t') # tab

        print(f"Successfully loaded known-compound list with {len(known_df)} records")
        return known_df
    except Exception as e:
        print(f"Failed to load known-compound list: {e}")
        return pd.DataFrame()


def search_by_mz(target_df, known_df, tolerance=0.01):
    """
    Match compounds by m/z
    """
    results = []
    matched_count = 0

    for idx, row in target_df.iterrows():
        theory_mz = row['theoretical_mz']
        found = False
        smiles = None
        cas = None

        # Find matching m/z values in the known list
        if 'mz' in known_df.columns:
            matches = known_df[np.abs(known_df['mz'] - theory_mz) <= tolerance]
            if len(matches) > 0:
                # Take the first match
                matched_row = matches.iloc[0]
                smiles = matched_row.get('smiles', matched_row.get('SMILES', None))
                cas = matched_row.get('cas', matched_row.get('CAS', None))
                found = True
                matched_count += 1

        results.append({
            'compound_name': row['compound_name'],
            'theoretical_mz': theory_mz,
            'measured_mz': row['measured_mz'],
            'retention_time': row['retention_time'],
            'retention_index': row['retention_index'],
            'match_status': 'matched_by_mz' if found else 'unmatched',
            'SMILES': smiles,
            'CAS': cas
        })

    print(f"Matched by m/z {matched_count} compound")
    return pd.DataFrame(results)


def search_by_cas(cas_number):
    """
    Query SMILES from PubChem by CAS number
    """
    if pd.isna(cas_number) or cas_number == '':
        return None

    try:
        # Remove spaces and special characters from the CAS number
        cas_str = str(cas_number).strip()

        # Try querying by CAS number
        compounds = get_compounds(cas_str, namespace='name')

        if compounds:
            # Take the first compound
            compound = compounds[0]
            return compound.canonical_smiles
        else:
            # If direct query fails, try the PubChem API
            try:
                url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{cas_str}/property/CanonicalSMILES/JSON"
                response = requests.get(url, timeout=10)
                if response.status_code == 200:
                    data = response.json()
                    smiles = data['PropertyTable']['Properties'][0].get('CanonicalSMILES', None)
                    return smiles
            except:
                pass

    except Exception as e:
        print(f"Query CAS {cas_number} failed: {e}")

    return None


def batch_search_cas(result_df):
    """
    Batch-query SMILES by CAS number
    """
    total = len(result_df[result_df['match_status'] == 'unmatched'])
    print(f"Start querying by CAS number {total} compound...")

    queried_count = 0
    for idx, row in result_df.iterrows():
        if row['match_status'] == 'unmatched' and pd.notna(row['CAS']):
            smiles = search_by_cas(row['CAS'])
            if smiles:
                result_df.at[idx, 'SMILES'] = smiles
                result_df.at[idx, 'match_status'] = 'matched_by_CAS'
                queried_count += 1

            # Avoid making requests too frequently
            sleep(0.5)

    print(f"CASmatch {queried_count} compound")
    return result_df


def main():
    # 1. Excel file with validation compounds.
    target_file = "../training_data/validation_set_condition1.xlsx"
    target_df = pd.read_excel(target_file)
    print(f"loadfile, {len(target_df)} compound")

    # 2. Known compound-column file. Replace this path with the local file used
    # for CAS/name matching when rerunning the PubChem helper workflow.
    known_file_path = "../training_data/known_compounds.xlsx"

    known_df = load_known_compounds(known_file_path)

    if known_df.empty:
        print("Warning: compoundcolumn, mzmatch")
        # resultsDataFrame
        result_df = pd.DataFrame({
            'compound_name': target_df['compound_name'],
            'theoretical_mz': target_df['theoretical_mz'],
            'measured_mz': target_df['measured_mz'],
            'retention_time': target_df['retention_time'],
            'retention_index': target_df['retention_index'],
            'match_status': '',
            'SMILES': None,
            'CAS': None
        })
    else:
        # 3. matched_by_mz
        result_df = search_by_mz(target_df, known_df, tolerance=0.01)

    # 4. unmatchedcompound, Try querying by CAS number
    if 'CAS' in result_df.columns:
        result_df = batch_search_cas(result_df)
    else:
        print("Warning: resultsCAScolumn, CAS")

    # 5. matching results
    status_counts = result_df['match_status'].value_counts()
    print("\nmatching results:")
    for status, count in status_counts.items():
        print(f"{status}: {count}")

    # 6. saveresultsExcelfile
    output_file = 'validation set_condition1_SMILES.xlsx'

    # column
    columns_order = [
        'compound_name', 'theoretical_mz', 'measured_mz', 'retention_time', 'retention_index',
        'match_status', 'SMILES', 'CAS'
    ]

    # retentioncolumn
    columns_order = [col for col in columns_order if col in result_df.columns]

    result_df = result_df[columns_order]
    result_df.to_excel(output_file, index=False)
    print(f"\nresultssave: {output_file}")

    # 7. results
    print("\n10compoundmatching results:")
    print(result_df.head(10).to_string())

    return result_df


# : compoundcolumnfile, UseCAS
def add_cas_numbers_manually():
    """
    CAS ()
    Returnscompound_nameCASDataFrame
    """
    # ,
    cas_data = {
        'compound_name': [],
        'CAS': []
    }

    # : CAS
    # cas_data['compound_name'].append('compound_1')
    # cas_data['CAS'].append('123-45-6')

    return pd.DataFrame(cas_data)


if __name__ == "__main__":
    # translated note
    result = main()

    # compoundcolumn, CAS
    # translated note
    print("\n:")
    print("1. compoundcolumnfilecolumn:")
    print(" - mz: theoretical_mz")
    print(" - smiles SMILES: SMILES")
    print(" - cas CAS: CAS ()")
    print("\n2. CASfailed, :")
    print(" - ")
    print(" - CAS")
    print(" - PubChem: https://pubchem.ncbi.nlm.nih.gov/")
