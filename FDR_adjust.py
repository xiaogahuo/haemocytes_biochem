import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests

# 读取Excel文件，使用前两行作为多级列索引（第一层：Depression/Subtype1/Subtype2等，第二层：指标名）
df = pd.read_excel('遗传相关性分析.xlsx', header=[0, 1])

# 识别所有P值列（第二层包含'P-value'且不含'FDR'）和对应的FDR列
pval_cols = [col for col in df.columns if 'P-value' in col[1] and 'FDR' not in col[1]]
fdr_cols = [col for col in df.columns if 'FDR P-value' in col[1]]

# 建立P值列到对应FDR列的映射（基于第一层名称匹配）
p_to_f = {}
for pcol in pval_cols:
    first_level = pcol[0]
    match = [fcol for fcol in fdr_cols if fcol[0] == first_level]
    if match:
        p_to_f[pcol] = match[0]

# 收集所有非空P值及其位置（列名 + 行索引）
p_values = []
positions = []  # 元素为 (pcol, idx)
for pcol in pval_cols:
    series = df[pcol]
    mask = series.notna()
    for idx, val in series[mask].items():
        p_values.append(val)
        positions.append((pcol, idx))

# 若存在有效P值，则进行整体FDR校正（控制错误发现率）
if p_values:
    _, p_corrected, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
    # 将校正后的P值写回对应的FDR列
    for (pcol, idx), corr_val in zip(positions, p_corrected):
        fcol = p_to_f[pcol]
        df.loc[idx, fcol] = corr_val

# 保存校正后的Excel文件
df.to_excel('遗传相关性分析_FDR校正.xlsx', index=True)