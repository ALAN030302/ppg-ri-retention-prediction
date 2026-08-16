import pandas as pd
import numpy as np
from pubchempy import get_compounds, Compound
import requests
from time import sleep
import warnings

warnings.filterwarnings('ignore')


def load_known_compounds(known_file_path):
    """
    加载已知化合物列表（包含CAS号和SMILES）
    假设文件格式为CSV，包含列：cas, smiles, mz (可选)
    """
    try:
        # 根据文件扩展名选择读取方式
        if known_file_path.endswith('.csv'):
            known_df = pd.read_csv(known_file_path)
        elif known_file_path.endswith('.xlsx') or known_file_path.endswith('.xls'):
            known_df = pd.read_excel(known_file_path)
        else:
            # 尝试CSV格式
            known_df = pd.read_csv(known_file_path, sep='\t')  # 尝试tab分隔

        print(f"成功加载已知化合物列表，共 {len(known_df)} 条记录")
        return known_df
    except Exception as e:
        print(f"加载已知化合物列表失败: {e}")
        return pd.DataFrame()


def search_by_mz(target_df, known_df, tolerance=0.01):
    """
    根据mz匹配化合物
    """
    results = []
    matched_count = 0

    for idx, row in target_df.iterrows():
        theory_mz = row['理论mz']
        found = False
        smiles = None
        cas = None

        # 在已知列表中查找匹配的mz
        if 'mz' in known_df.columns:
            matches = known_df[np.abs(known_df['mz'] - theory_mz) <= tolerance]
            if len(matches) > 0:
                # 取第一个匹配项
                matched_row = matches.iloc[0]
                smiles = matched_row.get('smiles', matched_row.get('SMILES', None))
                cas = matched_row.get('cas', matched_row.get('CAS', None))
                found = True
                matched_count += 1

        results.append({
            '化合物名称': row['化合物名称'],
            '理论mz': theory_mz,
            '实测mz': row['实测mz'],
            '保留时间': row['保留时间'],
            '保留指数': row['保留指数'],
            '匹配状态': '通过mz匹配' if found else '未匹配',
            'SMILES': smiles,
            'CAS': cas
        })

    print(f"通过mz匹配到 {matched_count} 个化合物")
    return pd.DataFrame(results)


def search_by_cas(cas_number):
    """
    通过CAS号在PubChem中查询SMILES
    """
    if pd.isna(cas_number) or cas_number == '':
        return None

    try:
        # 去除CAS号中的空格和特殊字符
        cas_str = str(cas_number).strip()

        # 尝试通过CAS号查询
        compounds = get_compounds(cas_str, namespace='name')

        if compounds:
            # 取第一个化合物
            compound = compounds[0]
            return compound.canonical_smiles
        else:
            # 如果直接查询失败，尝试通过PubChem API
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
        print(f"查询CAS {cas_number} 失败: {e}")

    return None


def batch_search_cas(result_df):
    """
    批量通过CAS号查询SMILES
    """
    total = len(result_df[result_df['匹配状态'] == '未匹配'])
    print(f"开始通过CAS号查询 {total} 个化合物...")

    queried_count = 0
    for idx, row in result_df.iterrows():
        if row['匹配状态'] == '未匹配' and pd.notna(row['CAS']):
            smiles = search_by_cas(row['CAS'])
            if smiles:
                result_df.at[idx, 'SMILES'] = smiles
                result_df.at[idx, '匹配状态'] = '通过CAS匹配'
                queried_count += 1

            # 避免请求过于频繁
            sleep(0.5)

    print(f"通过CAS号匹配到 {queried_count} 个化合物")
    return result_df


def main():
    # 1. 读取目标Excel文件
    target_file = "C:/Users/姚钱磊/Desktop/补充实验预测/训练/验证集_条件1.xlsx"
    target_df = pd.read_excel(target_file)
    print(f"加载目标文件成功，共 {len(target_df)} 个化合物")

    # 2. 加载已知化合物列表（需要您提供文件路径）
    # 请将下面的路径替换为您的已知化合物列表文件路径
    known_file_path = "C:/Users/姚钱磊/Desktop/补充实验的数据/验证化合物列表20251209.xlsx"  # 或 'known_compounds.xlsx'

    known_df = load_known_compounds(known_file_path)

    if known_df.empty:
        print("警告：未找到已知化合物列表，将跳过mz匹配步骤")
        # 创建空的结果DataFrame
        result_df = pd.DataFrame({
            '化合物名称': target_df['化合物名称'],
            '理论mz': target_df['理论mz'],
            '实测mz': target_df['实测mz'],
            '保留时间': target_df['保留时间'],
            '保留指数': target_df['保留指数'],
            '匹配状态': '等待查询',
            'SMILES': None,
            'CAS': None
        })
    else:
        # 3. 首先通过mz匹配
        result_df = search_by_mz(target_df, known_df, tolerance=0.01)

    # 4. 对于未匹配的化合物，尝试通过CAS号查询
    if 'CAS' in result_df.columns:
        result_df = batch_search_cas(result_df)
    else:
        print("警告：结果中未找到CAS列，无法通过CAS号查询")

    # 5. 统计匹配结果
    status_counts = result_df['匹配状态'].value_counts()
    print("\n匹配结果统计:")
    for status, count in status_counts.items():
        print(f"{status}: {count}")

    # 6. 保存结果到新Excel文件
    output_file = '验证集_条件1_带SMILES.xlsx'

    # 确保列的顺序
    columns_order = [
        '化合物名称', '理论mz', '实测mz', '保留时间', '保留指数',
        '匹配状态', 'SMILES', 'CAS'
    ]

    # 只保留实际存在的列
    columns_order = [col for col in columns_order if col in result_df.columns]

    result_df = result_df[columns_order]
    result_df.to_excel(output_file, index=False)
    print(f"\n结果已保存到: {output_file}")

    # 7. 显示部分结果
    print("\n前10个化合物的匹配结果:")
    print(result_df.head(10).to_string())

    return result_df


# 辅助函数：如果没有已知化合物列表文件，可以使用这个函数手动添加CAS号
def add_cas_numbers_manually():
    """
    手动添加CAS号（如果需要）
    返回一个包含化合物名称和CAS号的DataFrame
    """
    # 这里只是一个示例，您需要根据实际情况添加
    cas_data = {
        '化合物名称': [],
        'CAS': []
    }

    # 示例：添加一些CAS号
    # cas_data['化合物名称'].append('化合物_1')
    # cas_data['CAS'].append('123-45-6')

    return pd.DataFrame(cas_data)


if __name__ == "__main__":
    # 运行主程序
    result = main()

    # 如果没有已知化合物列表，可以尝试手动添加CAS号
    # 然后重新运行查询
    print("\n提示:")
    print("1. 请确保已知化合物列表文件包含以下列:")
    print("   - mz: 理论mz值")
    print("   - smiles 或 SMILES: SMILES字符串")
    print("   - cas 或 CAS: CAS号（可选）")
    print("\n2. 如果通过CAS号查询失败，可以:")
    print("   - 检查网络连接")
    print("   - 确认CAS号是否正确")
    print("   - 手动在PubChem网站查询: https://pubchem.ncbi.nlm.nih.gov/")