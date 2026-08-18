import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from SUtils import *

file_path = '/mnt/agents/upload/ld.xls'

# 读取数据
df_raw = pd.read_excel(multiview_result_path+"ld.xlsx", sheet_name='Figure', header=None)

# 提取 trait 名称、rg 和 p-value
traits = df_raw.iloc[2:, 0].values
rg_data = df_raw.iloc[2:, [1, 3, 5]].values.astype(float)
pval_data = df_raw.iloc[2:, [2, 4, 6]].values.astype(float)

df_rg = pd.DataFrame(rg_data, index=traits, columns=['Depression', 'Subtype 1', 'Subtype 2'])
df_pval = pd.DataFrame(pval_data, index=traits, columns=['Depression', 'Subtype 1', 'Subtype 2'])

# 转置：3 行 × 53 列
df_rg_T = df_rg.T
df_pval_T = df_pval.T

# p-value 转星号
def pval_to_stars(p):
    if p <= 0.001:  return '***'
    elif p <= 0.01: return '**'
    elif p <= 0.05: return '*'
    else: return ''

df_annot = df_pval_T.map(pval_to_stars)

# 绘图
fig, ax = plt.subplots(figsize=(27, 3))
sns.heatmap(
    df_rg_T,
    annot=df_annot,
    fmt='',
    annot_kws={"size": 12, "color": "grey"},  # "color": "#FFD700"
    cmap="vlag",
    center=0,
    vmin=-1.0,
    vmax=1.0,
    linewidth=.5,
    ax=ax)
ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha='right', fontsize=11)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=11)
plt.tight_layout()
plt.savefig(multiview_result_path+"ld.png", dpi=720, bbox_inches='tight')
plt.show()