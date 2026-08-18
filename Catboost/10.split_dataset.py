from functools import reduce
import pandas as pd
from sklearn.model_selection import train_test_split, StratifiedKFold
import sys
import os
sys.path.append(os.path.abspath(".."))
from SUtils import *


def cross_valid_exp(df, dataset_Name, clu_list, hc_list):
    print(df.columns.tolist())
    df_new_MDD, df_new_NC = df[df['final_cluster'].isin(clu_list)], df[df['final_cluster'].isin(hc_list)]
    df_new_MDD['final_cluster'], df_new_NC['final_cluster'] = 1, 0

    df_ALL = pd.concat([df_new_MDD, df_new_NC], ignore_index=True)
    print(df_ALL.groupby('final_cluster').size())
    # 重置索引
    df_ALL.reset_index(drop=True, inplace=True)
    # # 使用train_test_split分割DataFrame，确保'target'字段的比例在训练集和验证集中相似
    # train_df, test_df = train_test_split(df, test_size=0.3, stratify=df['final_cluster'], random_state=0)

    # 创建 StratifiedKFold 对象，设置折数为 5
    skf = StratifiedKFold(n_splits=5, random_state=42, shuffle=True)

    all_cols=Hae_Bio+covr10+['PRS_Score_5e-02']
    X, y = df_ALL[all_cols], df_ALL['final_cluster']
    # 使用 StratifiedKFold 对象的 split 方法进行数据集划分
    for fold, (train_index, test_index) in enumerate(skf.split(X, y)):
        if not os.path.exists(data_path + "cross_valid/{}/{}/".format(fold, dataset_Name)):
            os.makedirs(data_path + "cross_valid/{}/{}/".format(fold, dataset_Name))
        # print("TRAIN:", train_index, "TEST:", test_index)
        X_train, X_test = X.loc[train_index], X.loc[test_index]
        y_train, y_test = y.loc[train_index], y.loc[test_index]
        X_train.to_csv(data_path + "cross_valid/{}/{}/X_train.csv".format(fold, dataset_Name), index=None)
        y_train.to_csv(data_path + "cross_valid/{}/{}/y_train.csv".format(fold, dataset_Name), index=None)
        X_test.to_csv(data_path + "cross_valid/{}/{}/X_test.csv".format(fold, dataset_Name), index=None)
        y_test.to_csv(data_path + "cross_valid/{}/{}/y_test.csv".format(fold, dataset_Name), index=None)


if __name__ == '__main__':
    subtype_num = 2
    df_Y = load_final_cluster(subtype_num)
    df_MDD = pd.read_csv(data_path + "Diag_MDD-0.0.csv", low_memory=False)
    df_MDD = pd.merge(df_MDD, df_Y, on='eid', how='inner')
    print(df_MDD.columns.tolist())

    df_NC = pd.read_csv(data_path + "Diag_NC-0.0.csv", low_memory=False)
    df_NC['final_cluster'] = 0
    df = pd.concat([df_MDD, df_NC], ignore_index=True)

    cross_valid_exp(df.copy(), 'ALL_cluster12_vs_0', [1, 2], [0])
    cross_valid_exp(df.copy(), 'ALL_cluster1_vs_0', [1], [0])
    cross_valid_exp(df.copy(), 'ALL_cluster2_vs_0', [2], [0])

    df_Male = df[df['Sex_cate'] == 'L2_Male']
    cross_valid_exp(df_Male.copy(), 'Male_cluster12_vs_0', [1, 2], [0])
    cross_valid_exp(df_Male.copy(), 'Male_cluster1_vs_0', [1], [0])
    cross_valid_exp(df_Male.copy(), 'Male_cluster2_vs_0', [2], [0])

    df_Female = df[df['Sex_cate'] == 'L1_Female']
    cross_valid_exp(df_Female.copy(), 'Female_cluster12_vs_0', [1, 2], [0])
    cross_valid_exp(df_Female.copy(), 'Female_cluster1_vs_0', [1], [0])
    cross_valid_exp(df_Female.copy(), 'Female_cluster2_vs_0', [2], [0])

    df_BMI = df[df['BMI_cate'].isin(['BMI1', 'BMI3'])]
    cross_valid_exp(df_BMI.copy(), 'BMI_cluster12_vs_0', [1, 2], [0])
    cross_valid_exp(df_BMI.copy(), 'BMI_cluster1_vs_0', [1], [0])
    cross_valid_exp(df_BMI.copy(), 'BMI_cluster2_vs_0', [2], [0])

    df_BMI_Male = df[(df['BMI_cate'].isin(['BMI1', 'BMI3'])) & (df['Sex_cate'] == 'L2_Male')]
    cross_valid_exp(df_BMI_Male.copy(), 'BMI_Male_cluster12_vs_0', [1, 2], [0])
    cross_valid_exp(df_BMI_Male.copy(), 'BMI_Male_cluster1_vs_0', [1], [0])
    cross_valid_exp(df_BMI_Male.copy(), 'BMI_Male_cluster2_vs_0', [2], [0])

    df_BMI_Female = df[(df['BMI_cate'].isin(['BMI1', 'BMI3'])) & (df['Sex_cate'] == 'L1_Female')]
    cross_valid_exp(df_BMI_Female.copy(), 'BMI_Female_cluster12_vs_0', [1, 2], [0])
    cross_valid_exp(df_BMI_Female.copy(), 'BMI_Female_cluster1_vs_0', [1], [0])
    cross_valid_exp(df_BMI_Female.copy(), 'BMI_Female_cluster2_vs_0', [2], [0])

    df_Chronic = df[df['chronic_num'] == 0]
    cross_valid_exp(df_Chronic.copy(), 'Chronic_cluster12_vs_0', [1, 2], [0])
    cross_valid_exp(df_Chronic.copy(), 'Chronic_cluster1_vs_0', [1], [0])
    cross_valid_exp(df_Chronic.copy(), 'Chronic_cluster2_vs_0', [2], [0])

    df_Chronic_Male = df[(df['chronic_num'] == 0) & (df['Sex_cate'] == 'L2_Male')]
    cross_valid_exp(df_Chronic_Male.copy(), 'Chronic_Male_cluster12_vs_0', [1, 2], [0])
    cross_valid_exp(df_Chronic_Male.copy(), 'Chronic_Male_cluster1_vs_0', [1], [0])
    cross_valid_exp(df_Chronic_Male.copy(), 'Chronic_Male_cluster2_vs_0', [2], [0])

    df_Chronic_Female = df[(df['chronic_num'] == 0) & (df['Sex_cate'] == 'L1_Female')]
    cross_valid_exp(df_Chronic_Female.copy(), 'Chronic_Female_cluster12_vs_0', [1, 2], [0])
    cross_valid_exp(df_Chronic_Female.copy(), 'Chronic_Female_cluster1_vs_0', [1], [0])
    cross_valid_exp(df_Chronic_Female.copy(), 'Chronic_Female_cluster2_vs_0', [2], [0])
    


 