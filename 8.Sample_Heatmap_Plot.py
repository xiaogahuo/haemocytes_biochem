from SUtils import *


def heatmap_plot_bio_chem(cohort_stage):
    df_MDD = pd.read_csv(data_path + f"Diag_MDD-{cohort_stage}.0.csv", low_memory=False)
    #df_Y = pd.read_csv(result_path + "{}/{}_components_cluster_rst_{}.csv".format("NMF", 2, 3), low_memory=False,
                       #usecols=['eid', 'final_cluster'])
    #df_Y = pd.read_csv("/data4/myuan/Hae_Bio/F7_results.csv", low_memory=False, usecols=['eid', 'final_cluster'])
    df_Y = load_final_cluster(2)
    df_MDD = pd.merge(df_MDD, df_Y, how='inner', on='eid', sort=True, suffixes=('_x', '_y'), copy=True, indicator=False)
    # bio_chem_lst = [item for item in cate17518_new_cols if item not in ['RF', 'OEST']]

    df = df_MDD[['final_cluster']+Hae_Bio]
    print(df.shape[0])
    df.dropna(thresh=20, inplace=True)
    print(df.shape[0])
    null_ratios = df.isnull().sum() / df.shape[0]
    # 筛选出空值比例小于等于20%的列
    columns_to_keep = null_ratios[null_ratios <= 0.1].index
    keep_cols = columns_to_keep.tolist()
    print("{}, {}".format(len(keep_cols), keep_cols))
    del_col = [item for item in Hae_Bio if item not in keep_cols]
    print("del_col: " + str(del_col))

    df = df[keep_cols].sort_values(by='final_cluster')
    for col in keep_cols:
        df[col] = (df[col] - df[col].mean()) / df[col].std()

    # 设置索引并选择值
    heatmap_data = df.T
    plt.figure(figsize=(12, 12))  # 设置宽度为30，高度为12
    sns.heatmap(heatmap_data, annot=False, cmap='RdBu_r', center=0, vmin=-3, vmax=3)
    plt.xticks([], [])
    plt.title('Heatmap by Feature Order')
    plt.ylabel('Features')
    plt.savefig(multiview_result_path + "Heatmap_bio_chem.png", dpi=720)
    plt.show()
    plt.clf()

if __name__ == '__main__':
    heatmap_plot_bio_chem(cohort_stage)

