import pandas as pd
from SUtils import *
import gwaslab as gl
import re
from statsmodels.stats.multitest import multipletests


def Sample_QC(l1, l2, f_name):
    df_Y = load_final_cluster(2)
    df_NC = pd.read_csv(data_path + "Diag_NC-0.0.csv", low_memory=False, usecols=['eid'])
    df_NC['final_cluster'] = 0
    df_C = pd.concat([df_NC, df_Y])

    # 1 vs 2  ----- df_1为对照组，df_2为病例组
    df_1, df_2 = df_C[df_C['final_cluster'].isin(l1)], df_C[df_C['final_cluster'].isin(l2)]
    df_C = pd.concat([df_1, df_2])

    df_bed_sample = pd.read_csv(bed_path + "SNP_QC_raw_data/c1_b0_v3.fam", sep=r'[ \t]', header=None, engine='python')
    df_bed_sample.rename(columns={0: 'fid', 1: 'eid', 2: 'pid', 3: 'mid', 4: 'sid'}, inplace=True)

    df = pd.read_csv(data_path + "Diag-0.0.csv", usecols=['eid', 'Sex_cate'], low_memory=False)  # 'L1_Female', 'L2_Male'
    df_MDS = pd.read_csv(root_path + "Sample_QC.csv")
    df_MDS = reduce(lambda left, right: pd.merge(left, right, on=['eid'], how='inner'), [df_MDS, df, df_C, df_bed_sample])
    df_MDS['Sex_cate'] = df_MDS['Sex_cate'].replace({'L1_Female': 0, 'L2_Male': 1})

    print("df: {}".format(df_MDS.shape[0]))
    df_MDS = df_MDS[df_MDS["22006-0.0"] == 1]  # 22006-0.0 代表高加索人群
    print("Sample_MDS_QC after filter: {}".format(df_MDS.shape[0]))
    df_MDS = df_MDS[df_MDS["22021-0.0"] == 0]  # Genetic relatedness to other participants
    print("filter after 22021: {}".format(df_MDS.shape[0]))
    # 比较两列22001和Sex_cate并赋值给新列 'Sex_match', 相等的值被赋值为1，不相等的值被赋值为0
    df_MDS['Sex_match'] = (df_MDS['22001-0.0'] == df_MDS['Sex_cate']).astype(int)
    df_MDS = df_MDS[df_MDS["Sex_match"] == 1]
    print("filter after Sex Match Compare: {}".format(df_MDS.shape[0]))
    print(df_MDS.groupby('final_cluster').size())

    df_MDS = df_MDS.sort_values(by='eid')
    if not os.path.exists(bed_path + "{}/{}/".format(exp_name, f_name)):
        os.makedirs(bed_path + "{}/{}/".format(exp_name, f_name))
    df_MDS.to_csv(bed_path + "{}/{}/{}.csv".format(exp_name, f_name, f_name), columns=['fid', 'eid'], sep='\t', header=None, index=False)


def trans_genotyping_array(x):
    if -11 <= x <= -1:
        return 1   # 'BiLEVE'
    elif 1 <= x <= 95:
        return 2  # 'Axiom'
    else:
        return None  # 或者根据需要返回其他值


# .fam文件是储存样本信息的文件，为没有标题行的文本文件，每个样本占一行，共有以下6列：
# 1. 家系ID, Family ID (FID)
# 2. 个体ID, Within-family ID (IID;不能为0)
# 3. 父本个体ID(如果父本不在样本中则为0)
# 4. 母本个体ID(如果母本不在样本中则为0)
# 5. 性别(1为雄性,2为雌性, 0表示未知)
# 6. 表型值 (1为对照组,2为病例组,-9、0和非数值在二分类性状时表示缺失值)。如果存在除-9，0，1，2以外的任何数字表型值，则表型为数量性状，而不是病例/对照的二分类性状。在这种情况下，-9仍表示表型缺失值。
covar_num = 6  # <=41
def split_cluster(l1, l2, f_name):
    df_Y = load_final_cluster(2) # subtype == 3
    df_NC = pd.read_csv(data_path + "Diag_NC-0.0.csv", low_memory=False, usecols=['eid'])
    df_NC['final_cluster'] = 0
    df_C = pd.concat([df_NC, df_Y])

    # 0 vs 1  ----- df_1为对照组，df_2为病例组
    df_1, df_2 = df_C[df_C['final_cluster'].isin(l1)], df_C[df_C['final_cluster'].isin(l2)]
    df_1['final_cluster'], df_2['final_cluster'] = 1, 2 # 如果值为 1, 2 后续 GWAS分析无法进行 --linear 分析,可以进行--logistic 分析
    df_C = pd.concat([df_1, df_2])
    df = pd.read_csv(data_path + "Diag-0.0.csv", usecols=['eid', 'Age'], low_memory=False)

    df_MDS = pd.read_csv(root_path + "Sample_QC.csv")
    PCA_original_cols, PCA_new_cols = ["22009-0." + str(i) for i in range(1, covar_num)], ["PCA" + str(i) for i in range(1, covar_num)]
    df_MDS.rename(columns=dict(zip(PCA_original_cols, PCA_new_cols)), inplace=True)
    df_MDS['genotyping_array'] = df_MDS['22000-0.0'].apply(trans_genotyping_array)

    for chr in range(1, 23):
        df_bed_fam = pd.read_csv(bed_path + "{}/{}/step2_chr{}.fam.backup".format(exp_name, f_name, chr), sep=r'[ \t]', header=None, engine='python')
        df_bed_fam.rename(columns={0: 'fid', 1: 'eid', 2: 'pid', 3: 'mid', 4: 'sid'}, inplace=True)
        df_Temp = reduce(lambda left, right: pd.merge(left, right, on=['eid'], how='inner'), [df_bed_fam, df_C, df, df_MDS])
        print("df_Temp after merge: {}".format(df_Temp.shape))
        print(df_Temp.columns.tolist())
        df_Temp.to_csv(bed_path + "{}/{}/step2_chr{}.fam".format(exp_name, f_name, chr),
                  columns=['fid', 'eid', 'pid', 'mid', 'sid', 'final_cluster'], sep='\t', header=None, index=False)
        df_Temp.to_csv(bed_path + "{}/{}/covr_chr{}.txt".format(exp_name, f_name, chr),
                  columns=['fid', 'eid', 'Age', 'genotyping_array'] + PCA_new_cols, sep='\t', header=None, index=False)
        print(df_Temp.groupby('final_cluster').size())


def plot_manhattan_QQ(file_name, cut=20, sig_level=5e-8):
    df = pd.read_csv(bed_path + "{}/{}/{}_gwas_summary.txt.gz".format(exp_name, file_name, file_name), sep=r'\s+', engine='python')
    print(df.columns.tolist()) # 'SNP','CHR','BP','A1','A2','P','OR','BETA','SE'
    print("file_name: {}, snp number: {}".format(file_name, df.shape[0]))
    df_m = gl.Sumstats(df, snpid="SNP", chrom="CHR", pos="BP", ea="A1", nea="A2", eaf=None,
             beta="BETA", se="SE", p="P")
    df_m.plot_mqq(skip=0, cut=cut, anno=None, sig_level=sig_level, sig_level_lead=sig_level,
        colors=["#ff0000", "#8cc43c", "#fc4444", "#029658", "#fc6404", "#fcd444", "#1abc9c", "#5bc0de", "#6454ac", "#fc8c84"],
        fontsize=8, anno_fontsize=15, title_fontsize=13, marker_size=(5, 25), figargs={"figsize": (15, 5), "dpi":300},
        save=bed_path + "{}/{}/{}manhattan_QQ_plots.png".format(exp_name, file_name, file_name))

def cnt_sig_pvalue():
    # f_s = ["0_VS_123", "0_VS_1", "0_VS_2", "0_VS_3", "1_VS_2", "1_VS_3", "2_VS_3"]
    # f_s = ["0_VS_1", "0_VS_2"]
    f_s = ["0_VS_2", "0_VS_1"]
    for file_name in f_s:
        df = pd.read_csv(bed_path + "{}/{}/logistic_results.assoc.logistic2".format(exp_name, file_name), sep=r'\s+', engine='python')
        df.rename(columns={'STAT': 'Z', 'NMISS': 'n', 'P': 'p'}, inplace=True)
        df_freq = pd.read_csv(bed_path + "{}/{}/MAF.frq".format(exp_name, file_name), sep=r'\s+', engine='python')
        df_freq.rename(columns={'MAF': 'freq'}, inplace=True)
        df_freq = df_freq[["SNP", "A2", "freq"]]
        print(df.head())
        print(df_freq.head())
        df = pd.merge(df, df_freq, on="SNP")
        # df = df.sort_values(by="P")
        df['b'] = df['OR'].apply(lambda x: np.log(x))
        df['se'] = df['b'] / df['Z']
        df.rename(columns={'SNP': 'snp'}, inplace=True)
        print(df.columns.tolist())
        # SNP    A1  A2  freq    b   se  p   n , 'OR', 'Z', 'CHR', 'BP'
        df[['snp', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'n']].to_csv(bed_path + "{}/{}/{}_gwas_summary.txt.gz".format(exp_name, file_name, file_name), sep='\t', header=True, index=False, compression='gzip')  # sep='\t'
        df.rename(columns={'snp': 'snpid', 'CHR': 'hg18chr', 'BP': 'bp', 'A1': 'a1', 'A2': 'a2', 'OR': 'or', 'p': 'pval'}, inplace=True)
        df[['snpid', 'hg18chr', 'bp', 'a1', 'a2', 'or', 'pval']].to_csv(
            bed_path + "{}/{}/logistic_results.assoc.logistic3".format(exp_name, file_name), sep='\t', header=True,
            index=False)
        df = df[df['pval'] <= 1e-5]
        df.to_csv(bed_path + "{}/{}/logistic_results.1e-5.txt".format(exp_name, file_name), sep='\t', index=False)


def format_LDGC_rst():
    # f_s = ["0_VS_1", "0_VS_2", "1_VS_2"]
    f_s = ["0_VS_2"]
    for file_name in f_s:
        result_path = bed_path + "{}/{}/LDGC/".format(exp_name, file_name)
        f = open(bed_path + "{}/{}/{}_LDGC.csv".format(exp_name, file_name, file_name), "w")
        header_str = "Name,rg,se,z,p,h2_obs,h2_obs_se,h2_int,h2_int_se,gcov_int,gcov_int_se"
        print(header_str, file=f)
        for item in os.listdir(result_path):
            if item.endswith('.log'):
                file_path = os.path.join(result_path, item)
                item = item.replace('.log', '')
                print(item)
                with open(file_path, 'r') as file:
                    lines = file.readlines()
                    extract_line = lines[-4]  # 读取倒数第4行
                    match_lst = extract_line.split()
                    match_lst = [item] + match_lst[2:]
                    print(match_lst)
                    print(','.join(map(str, match_lst)), file=f)


if __name__ == '__main__':
    # # before GWAS  -  1.1
    # Sample_QC([0], [1, 2], '0_VS_12')
    # Sample_QC([0], [1], '0_VS_1')
    # Sample_QC([0], [2], '0_VS_2')
    # Sample_QC([1], [2], '1_VS_2')

    # # before GWAS  - 2.1
    # split_cluster([0], [1, 2], '0_VS_12')
    # split_cluster([0], [1], '0_VS_1')
    # split_cluster([0], [2], '0_VS_2')
    # split_cluster([1], [2], '1_VS_2')

    # # after GWAS  - 2.3
    cnt_sig_pvalue()
    # plot_manhattan_QQ('0_VS_12', cut=9)
    # plot_manhattan_QQ('0_VS_1', cut=9)
    # plot_manhattan_QQ('0_VS_2', cut=9)
    # plot_manhattan_QQ('1_VS_2', cut=9)


    # # after LDGC  - 3.2
    # format_LDGC_rst()
    pass


