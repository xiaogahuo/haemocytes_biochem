from SUtils import *
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests

# 读取Excel文件，使用前两行作为多级列索引
df = pd.read_excel(multiview_result_path+'遗传相关性分析.xlsx', header=[0, 1])

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

# ========== 修改部分：各P-value列独立做FDR校正 ==========
for pcol, fcol in p_to_f.items():
    series = df[pcol]
    mask = series.notna()
    p_values = series[mask].values
    print(len(p_values))

    if len(p_values) > 0:
        _, p_corrected, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
        # 将校正后的P值写回对应的FDR列的对应位置
        df.loc[mask, fcol] = p_corrected

# 保存校正后的Excel文件
df.to_excel(multiview_result_path+'遗传相关性分析_FDR校正.xlsx', index=True)