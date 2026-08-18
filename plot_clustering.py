from SUtils import *
import math


# ============================================================
# 配置文件路径 —— 6 个 Subtype 的聚类结果 CSV
# ============================================================
files = [
    multiview_result_path+'/Hae_Bio_0.0_Subtype_2_clustering_results.csv',
    multiview_result_path+'/Hae_Bio_0.0_Subtype_3_clustering_results.csv',
    multiview_result_path+'/Hae_Bio_0.0_Subtype_4_clustering_results.csv',
    multiview_result_path+'/Hae_Bio_0.0_Subtype_5_clustering_results.csv',
    multiview_result_path+'/Hae_Bio_0.0_Subtype_6_clustering_results.csv',
    multiview_result_path+'/Hae_Bio_0.0_Subtype_7_clustering_results.csv'
]

# ============================================================
# 读取数据并标注 Subtype
# ============================================================
dfs = []
for f in files:
    subtype = int(f.split('Subtype_')[1].split('_')[0])
    df = pd.read_csv(f)
    df['Subtype'] = subtype
    dfs.append(df)

# ============================================================
# 格式化函数：整数显示为整数，小数保留原样
# ============================================================
def fmt_val(x):
    if x == int(x):
        return str(int(x))
    else:
        return str(x)

# ============================================================
# 绘图：6×1 纵向子图布局
#   - 每个 Subtype 独占一行，保证 80 个横轴标签全部显示且不重叠
#   - 总宽度 36 英寸、总高度 30 英寸，DPI 300
# ============================================================
def plot(y):
    global_max = max(df[y].max() for df in dfs)
    y_max = math.ceil(global_max)

    fig, axes = plt.subplots(6, 1, figsize=(36, 30), sharey=True)

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']

    for idx, (df, ax) in enumerate(zip(dfs, axes)):
        subtype = df['Subtype'].iloc[0]

        # 按用户指定格式拼接横轴标签
        x_labels = df.apply(
            lambda r: f"Anchor:{int(r['Anchor'])}-β:{fmt_val(r['Beta'])}-λ:{fmt_val(r['Lambda'])}",
            axis=1
        ).tolist()

        y_values = df[y].values
        x_pos = np.arange(len(x_labels))

        # 绘制折线 + 散点
        ax.plot(x_pos, y_values,
                marker='o', markersize=4, linewidth=1.2,
                color=colors[idx], alpha=0.85)

        # 横轴刻度与标签：全部显示，90° 旋转，底部对齐
        ax.set_xticks(x_pos)
        ax.set_xticklabels(x_labels, rotation=90, ha='center', fontsize=7)

        # 纵轴与标题
        ax.set_ylabel(y, fontsize=12)
        ax.set_title(f'Subtype {subtype}', fontsize=14, fontweight='bold', loc='left')
        ax.set_ylim(0, y_max)

        # 网格线与边框美化
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.suptitle(f'{y} with Multi-view Clustering Parameters (Anchor, β, λ) on different Clustering Subtype Number: 2-7',
                 fontsize=18, fontweight='bold', y=0.998)
    plt.tight_layout()
    plt.savefig(multiview_result_path+f'/{y}.png',
                dpi=300, bbox_inches='tight')
    plt.show()


if __name__ == '__main__':
    plot(y='Initial Silhouette')
    plot(y='Initial CH')
    plot(y='Initial DB')