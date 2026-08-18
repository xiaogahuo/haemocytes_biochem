import pandas as pd
import os
import numpy as np
from scipy.stats import chi2_contingency, zscore
from scipy.stats.mstats import kruskalwallis
from sklearn.feature_selection import SelectKBest, f_regression
from scipy.stats import ttest_ind
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm
from functools import reduce
from sklearn.metrics import roc_curve


exp_name = "Hae_Biochem_Clusterings_NUM2"
bed_path = "/data/UKBioBank/bed_Data/"
root_path = "/data/UKBioBank/FileHandler_And_Dataset/"
path = root_path +"Hae_Bio/"
data_path = path + "dataset/"
PRS_path = root_path + "PRS_MDD_Scores/"
# bgen_path = "/data/UKBioBank/bgen_Data/raw_data/"

cohort_stage = 0  # 要替换 {} 的数字
final_path = path + f"final_cluster_{cohort_stage}/"
result_path = final_path
multiview_result_path = path + "multiview_cluster_0/"

def load_final_cluster(subtype_num):
    csv_path = multiview_result_path+f"Hae_Bio_0.0_Subtype_{subtype_num}_anchor_2_beta_0.2500_lambda_0.0010_results.csv"
    df = pd.read_csv(csv_path, low_memory=False, usecols=['eid', 'final_cluster'])
    print(f"成功加载文件: {csv_path}")
    return df

# 自动创建目录
# os.makedirs(img_path, exist_ok=True)
# os.makedirs(result_path, exist_ok=True)

Haemocytes = ['WBC', 'RBC', 'HGB', 'HCT', 'MCV', 'MCH', 'MCHC', 'RDW', 'PLT', 'PCT', 'MPV', 'PDW', 'LYcnt', 'MOcnt', 'NEcnt', 'EOcnt', 'BAnt', 'NRBCcnt', 'Lyp', 'Mop', 'Nep', 'Eop', 'Bap', 'NRBCp', 'RETp', 'RETcnt', 'MRV', 'MSCV', 'IRF', 'HLRp','HLRc']   #31
Biochem = ['ALB', 'ALP', 'ALT', 'APOA.bio', 'APOB.bio', 'AST', 'DBIL', 'UREA', 'CALC', 'CHOL', 'CREA', 'CRP', 'CYSC', 'GGT', 'GLU', 'HbA1c', 'HDL', 'IGF1', 'LDL', 'LPa', 'OEST', 'PHOS', 'RF', 'SHBG.bio', 'TBIL', 'TEST', 'TP', 'TRIG', 'UA', 'VITD']   #30
Hae_Bio = Haemocytes + Biochem
Top17 =['APOA.bio', 'CRP', 'CYSC', 'GGT', 'HDL', 'HLRc', 'NEcnt', 'RBC', 'RETcnt', 'SHBG.bio', 'TBIL', 'TEST', 'TRIG', 'UA', 'UREA', 'VITD', 'WBC']

cat_cols = ['Sex_cate', 'Townsend_index_cate', 'Qualifications_cate', 'Smoking_status_cate', 'Alcohol_status_cate', 'BMI_cate', 'IPAQ_cate', 'family_history']
conti_cols = ['Age', 'chronic_num']
covr10 = cat_cols + conti_cols
covr3 = ['Age', 'BMI_cate', 'Sex_cate']


NC, MDD, UNNORMAL = 0, 1, -1
cate100042_cols = ['2188-0.0','2178-0.0','2296-0.0','2306-0.0']
cate100057_cols = ['1160-0.0','1170-0.0','1180-0.0','1190-0.0','1200-0.0','1210-0.0','1220-0.0']
cate100060_cols = [c+'-0.0' for c in ['20123','20124','20125','4609','4620','5375','5386','1920','1930','1940','1950','1960','1970','1980','1990','2000','2010','2020','2030','2040','2090','2100','4598','4631','4642','4653','4526','4548','4559','4570','4581','4537','2050','2060','2070','2080','6156','5663','5674','6145','10721','20122','20126','20127']]
cate100061_cols = [c+'-0.0' for c in ['1031','6160','2110','10740']]
cate100065_cols = [c+'-0.0' for c in ['21000']]
cate145_cols = [c+'-0.0' for c in ['20487','20488','20489','20490','20491','20521','20522','20523','20524','20525','20526','20527','20528','20529','20530','20531','20494','20495','20496','20497','20498']]


def get_best_cutoff(y_true, y_scores):
    fpr, tpr, thresholds = roc_curve(y_true, y_scores, pos_label=1)
    optimal_idx = np.argmax(tpr - fpr)
    optimal_threshold = thresholds[optimal_idx]
    return optimal_threshold


# 将labels转换为numpy数组，并使用广播和向量化操作来计算相似度矩阵。这样可以避免使用双重循环，提高代码性能。
def compute_similarity_matrix(labels):
    labels = np.array(labels)
    similarity_matrix = (labels[:, None] == labels).astype(int)
    return similarity_matrix


def annot_trans(pval):
    if pval < 0.001:
        return '***'
    elif pval < 0.01:
        return '**'
    elif pval < 0.05:
        return '*'
    else:
        return ' '


def social_stats(exps, cmp_col):  # 'MDD'
    for exp in exps:
        df = pd.read_csv(data_path + exp + '.csv')
        print('---------------exp: {}, number: {}------------------------'.format(exp, df.shape[0]))
        for col in cat_cols+[cmp_col]:
            count_1feature(col, df)
        for col in Haemocytes + Biochem + conti_cols:
            count_mean_std(col, df)

        for col in cat_cols:
            count_2features(col, cmp_col, df)
            count_Chi2(col, cmp_col, df)
        for col in Haemocytes + Biochem + conti_cols:
            count_mean_std_group_by_div_feature(col, cmp_col, df)
            kwtest(col, cmp_col, df)


def count_1feature(feature_name, df):
    print('----------打印字段: <' + feature_name + "> 的人数分组统计情况----------")
    target_feature_value_list = df.drop_duplicates(subset=[feature_name], keep='first')[feature_name].tolist()
    target_feature_value_list.sort()
    for value in target_feature_value_list:
        print("字段< %s > 类别 < %s > 人数统计: %d" % (feature_name, str(value), df[df[feature_name] == value].shape[0]))


def count_2features(feature_name1, feature_name2, df):
    print('----------打印字段: <' + feature_name1 + "> 和 字段< " + feature_name2 + " > 的分组统计情况----------")
    target_feature_value_list1 = df.drop_duplicates(subset=[feature_name1], keep='first')[feature_name1].tolist()
    target_feature_value_list1.sort()
    target_feature_value_list2 = df.drop_duplicates(subset=[feature_name2], keep='first')[feature_name2].tolist()
    target_feature_value_list2.sort()

    for value1 in target_feature_value_list1:
        for value2 in target_feature_value_list2:
            print("字段< %s > 类别 < %s > 字段< %s > 类别 < %s > 人数统计: %d" % (
            feature_name1, str(value1), feature_name2, str(value2),
            df[(df[feature_name1] == value1) & (df[feature_name2] == value2)].shape[0]))


def count_Chi2(cnt_feature, div_feature, df):
    print('卡方检验结果: <' + cnt_feature + "> on <" + div_feature + ">")
    # print(df.drop_duplicates(subset=[target_feature], keep='first')[target_feature])
    div_feature_value_list = sorted(df.drop_duplicates(subset=[div_feature], keep='first')[div_feature].tolist())
    cnt_feature_value_list = sorted(df.drop_duplicates(subset=[cnt_feature], keep='first')[cnt_feature].tolist())
    print("{}: {}".format(div_feature, div_feature_value_list))
    print("{}: {}".format(cnt_feature, cnt_feature_value_list))

    kf_data = []
    for div in div_feature_value_list:
        c_g =[]
        for value in cnt_feature_value_list:
            count = df[(df[div_feature] == div) & (df[cnt_feature] == value)].shape[0]
            c_g.append(count)
        kf_data.append(c_g)
    print(kf_data)
    kf = chi2_contingency(kf_data)
    print('chisq-statistic=%.4f, p-value=%.4f, df=%i expected_frep=%s' % kf)
    print('-'*25)
    return kf, kf_data


def kwtest(feature_name, target_feature, df):
    print('kwtest结果: <' + feature_name + "> on <" + target_feature + ">")
    # print(df.drop_duplicates(subset=[target_feature], keep='first')[target_feature])
    target_feature_value_list = df.drop_duplicates(subset=[target_feature], keep='first')[target_feature].tolist()

    kw_data = []
    for value in target_feature_value_list:
        feature_name_value_list = df[df[target_feature] == value][feature_name].tolist()
        kw_data.append(feature_name_value_list)
    statistic, pvalue = kruskalwallis(*kw_data) #  将kw_data中的值作为参数传入kruskalwallis中, 第1个值是第1个参数， 第2个值是第2个参数
    print("column: {}, statistic: {}, pvalue: {}".format(feature_name, statistic, pvalue))
    print('-'*25)


def count_mean_std_group_by_div_feature(feature_name, div_feature_name, df):
    print(
        '----------打印字段: <' + feature_name + "> 根据字段< " + div_feature_name + " > 的分组 均值(Mean) 标准差(Std) 统计情况----------")
    print('字段< %s > 总体Mean: %f' % (feature_name, df[feature_name].mean()))
    print('字段< %s > 总体Std: %f' % (feature_name, df[feature_name].std(ddof=1)))
    print('-' * 25)
    div_feature_value_list = df.drop_duplicates(subset=[div_feature_name], keep='first')[div_feature_name].tolist()
    div_feature_value_list.sort()

    for div_value in div_feature_value_list:
        print("字段< %s > 按字段< %s > 类别 < %s > 统计的均值 (Mean): %f" % (
            feature_name, div_feature_name, str(div_value), df[df[div_feature_name] == div_value][feature_name].mean()))
        print("字段< %s > 按字段< %s > 类别 < %s > 统计的标准差 (Std): %f" % (
            feature_name, div_feature_name, str(div_value),
            df[df[div_feature_name] == div_value][feature_name].std(ddof=1)))


def count_mean_std(feature_name, df):
    print('----------打印字段: <' + feature_name + "> 的均值(Mean) 标准差(Std) 统计情况----------")
    print('字段< %s > 总体Mean: %f' % (feature_name, df[feature_name].mean()))
    print('字段< %s > 总体Std: %f' % (feature_name, df[feature_name].std(ddof=1)))
    print('-' * 25)


def count_T(feature_name, target_feature, df):
    count_mean_std_group_by_div_feature(feature_name, target_feature, df)
    print('T检验结果: <' + feature_name + "> on <" + target_feature + ">")
    # print(df.drop_duplicates(subset=[target_feature], keep='first')[target_feature])
    target_feature_value_list = df.drop_duplicates(subset=[target_feature], keep='first')[target_feature].tolist()
    c1 = target_feature_value_list[0]
    c2 = target_feature_value_list[1]
    v1_list = df[df[target_feature] == c1][feature_name].tolist()
    v2_list = df[df[target_feature] == c2][feature_name].tolist()
    v1_list = [v for v in v1_list if not pd.isnull(v)]
    v2_list = [v for v in v2_list if not pd.isnull(v)]
    print("size:{}".format(len(v1_list)))
    print("size:{}".format(len(v2_list)))

    # Run a two sample t-test to compare the two samples
    tstat, pval = ttest_ind(a=v1_list, b=v2_list, alternative="two-sided")

    # Display results
    print("t-stat: {:.2f}   pval: {:.4f}".format(tstat, pval))
    return tstat, pval


#palette = {"HC": "#009432", "Subtpye1": "#F79F1F", "Subtpye2": "#EA2027", "Subtpye3": "#112027","Subtpye4": "#1E90FF", "MDD": "#EAB543"}
palette = {
    "HC": "#009432",        # Healthy Control
    "Subtpye1": "#F79F1F",  # Subtype 1
    "Subtpye2": "#EA2027",  # Subtype 2
    "Subtpye3": "#112027",  # Subtype 3
    "Subtpye4": "#1E90FF",  # Subtype 4
    "Subtpye5": "#8A2BE2",  # Subtype 5 - 新增
    "Subtpye6": "#FF6347",  # Subtype 6 - 新增
    "Subtpye7": "#4682B4",  # Subtype 7 - 新增
    "MDD": "#EAB543"        # Major Depressive Disorder
}
def plot_bar(df, file_name, feature_lst, width, high, hue_order, y_label):
    df_bar = []
    for feature in feature_lst:
        df_sub = df[[feature, 'final_cluster']].copy()
        df_sub['Features'] = feature
        df_sub.rename(columns={feature: y_label}, inplace=True)
        # print(df_sub.head())
        df_bar.append(df_sub)
    df_bar = pd.concat(df_bar, ignore_index=True)
    # print("before: {}".format(df_bar.shape))
    print(df_bar.head())
    # df_bar = df_bar[(df_bar['z-score']<=5) & (df_bar['z-score']>=-5)]
    print(" after: {}".format(df_bar.shape))

    plt.figure(figsize=(width, high))
    ax = sns.barplot(hue="final_cluster", hue_order=hue_order, y=y_label, x="Features", data=df_bar,
                     capsize=.2, errorbar=("ci", 95), err_kws={"color": ".9", "linewidth": 0.6},
                palette=palette, width=.5) # errorbar=None   ,gap=0.1   width=0.5
    ax.tick_params(axis='x', labelrotation=0, labelsize=6)
    ax.tick_params(axis='y', labelrotation=0, labelsize=6)
    plt.legend(loc='upper center', bbox_to_anchor=(1.07, 1.02))  # 使用这行代码把图例放到图外。
    plt.savefig(result_path + "Bar_{}.png".format(file_name), dpi=720)
    plt.clf()


def plot_violin(df, file_name, feature_lst, width, high, hue_order, y_label):
    df_violin = []
    for feature in feature_lst:
        df_sub = df[[feature, 'final_cluster']].copy()
        df_sub['Features'] = feature
        df_sub.rename(columns={feature: y_label}, inplace=True)
        # print(df_sub.head())
        df_violin.append(df_sub)
    df_violin = pd.concat(df_violin, ignore_index=True)
    # print("before: {}".format(df_violin.shape))
    print(df_violin.head())

    plt.figure(figsize=(width, high))
    ax = sns.violinplot(hue="final_cluster", hue_order=hue_order, y=y_label, x="Features", data=df_violin,
                palette=palette, width=.9, gap=0.1, legend=False) # errorbar=None
    ax.tick_params(axis='x', labelrotation=0, labelsize=8)
    ax.tick_params(axis='y', labelrotation=0, labelsize=8)
    # plt.legend(loc='upper center', bbox_to_anchor=(1.07, 1.02))  # 使用这行代码把图例放到图外。
    sns.despine(ax=ax, top=True, right=True) # 去掉顶部和右侧框线
    plt.savefig(multiview_result_path + "Violin_{}.png".format(file_name), dpi=720)
    plt.show()
    plt.clf()


def cnt_OR_CI95(f_list):
    for f_name in f_list:
        z = norm.ppf(0.975)  # 对应于95%置信水平的z值
        df = pd.read_csv(result_path+'{}.csv'.format(f_name))
        # df['clust1_OR'] = df.apply(lambda row: np.exp(row['clust1_Estimate']), axis=1)
        df['clust1_lower95'] = df.apply(lambda row: row['clust1_Estimate']-z*row['clust1_StdError'], axis=1)
        df['clust1_upper95'] = df.apply(lambda row: row['clust1_Estimate']+z*row['clust1_StdError'], axis=1)
        # df['clust1'] = 'clust1'

        # df['clust2_OR'] = df.apply(lambda row: np.exp(row['clust2_Estimate']), axis=1)
        df['clust2_lower95'] = df.apply(lambda row: row['clust2_Estimate']-z*row['clust2_StdError'], axis=1)
        df['clust2_upper95'] = df.apply(lambda row: row['clust2_Estimate']+z*row['clust2_StdError'], axis=1)
        # df['clust2'] = 'clust2'

        # df_rst = pd.DataFrame({'Feature': df['Col'].tolist()+df['Col'].tolist(),
        #                        'OR': df['clust1_OR'].tolist()+df['clust2_OR'].tolist(),
        #                        'CI_lower': df['clust1_lower95'].tolist()+df['clust2_lower95'].tolist(),
        #                        'CI_upper': df['clust1_upper95'].tolist()+df['clust2_upper95'].tolist(),
        #                        'Clust': df['clust1'].tolist()+df['clust2'].tolist()})
        df.to_csv(result_path+'{}.csv'.format(f_name), index=False)


if __name__ == '__main__':
    # 创建一个示例dataframe
    # rst = set(pd.read_csv(path + "proteomics.csv").columns.tolist())-set(Prote_cols_all)
    # print(rst)

    # df1 = pd.read_excel(path+'protein2923.xlsx', sheet_name='Sheet1', usecols=['coding', 'meaning'])
    # df2 = pd.read_excel(path+'protein2923.xlsx', sheet_name='Sheet2', usecols=['Assay Target', 'Protein panel'])
    # df = pd.merge(df1, df2, left_on='coding', right_on='Assay Target', how='inner')
    # df.to_csv(path + "protein2923.csv", index=None)
    # print(df.shape[0])
    load_final_cluster()
    pass
