import os
import pandas as pd
import numpy as np
from pathlib import Path
from statsmodels.stats.multitest import multipletests
from SUtils import *


def process_smr_qtl_recursive(folder='SMR_QTL'):
    """
    递归遍历 folder 及其所有子目录下的 .csv 文件，
    按 qtl_name 分组对 p_SMR 做 Benjamini-Hochberg FDR 校正，
    结果写入 p_SMR_FDR 列，并输出为同目录下的 .xlsx 文件。
    """
    root = Path(folder)
    if not root.exists():
        raise FileNotFoundError(f"文件夹不存在: {root.resolve()}")

    # 递归查找所有 csv 文件（不区分大小写）
    csv_files = [p for p in root.rglob('*')
                 if p.is_file() and p.suffix.lower() == '.csv']

    if not csv_files:
        print(f"在 {root.resolve()} 及其子目录中未找到 CSV 文件")
        return

    for csv_path in csv_files:
        print(f"正在处理: {csv_path}")

        # 1. 读取 CSV，支持 Tab 或空格分隔
        try:
            df = pd.read_csv(csv_path, sep=r'\s+', engine='python')
        except Exception as e:
            print(f"  [错误] 读取失败: {e}")
            continue

        # 2. 检查必要列
        if 'qtl_name' not in df.columns or 'p_SMR' not in df.columns:
            print(f"  [跳过] 文件缺少 'qtl_name' 或 'p_SMR' 列")
            continue

        # 3. 按 qtl_name 分组，每组单独做 FDR (Benjamini-Hochberg)
        def fdr_correct(group):
            pvals = group['p_SMR'].to_numpy(dtype=float)
            mask = ~np.isnan(pvals)  # 排除缺失值
            fdr_arr = np.full_like(pvals, np.nan, dtype=float)

            if mask.any():
                _, fdr_adj, _, _ = multipletests(pvals[mask], method='fdr_bh')
                fdr_arr[mask] = fdr_adj

            group['p_SMR_FDR'] = fdr_arr
            return group

        df = df.groupby('qtl_name', group_keys=False).apply(fdr_correct)

        # 4. 写入同目录的 xlsx 文件（仅替换扩展名）
        xlsx_path = csv_path.with_suffix('.xlsx')

        try:
            df.to_excel(xlsx_path, index=False)
            print(f"  [完成] 已保存: {xlsx_path}")
        except Exception as e:
            print(f"  [错误] 保存失败: {e}")


if __name__ == '__main__':
    # 如果 SMR_QTL 文件夹在当前工作目录下，直接运行即可
    process_smr_qtl_recursive(path+'/SMR_QTL/')