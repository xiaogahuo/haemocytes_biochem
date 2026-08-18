import shap
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import os
sys.path.append(os.path.abspath(".."))
from SUtils import *


# 设置 matplotlib 使用无显示后端
import matplotlib
matplotlib.use('Agg')

def get_topN_reason(exp_Name, feature_importance, cols, model_name, top_num=20):
    feature_importance_dict = {f: i for i, f in zip(feature_importance, cols)}
    new_dict = dict(sorted(feature_importance_dict.items(), key=lambda e: e[1], reverse=True))
    
    # 确保保存路径存在
    save_dir = os.path.join(multiview_result_path, exp_Name)
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
        
    with open(os.path.join(save_dir, f"{model_name}_SHAP.csv"), "w") as f:
        print("key,value", file=f)
        for k, v in new_dict.items():
            if top_num > 0:
                print(f"{k},{v}", file=f)
                top_num -= 1
            else:
                break

def plot_shap(exp_Name, model_name, cols):
    # 读取数据并检查文件是否存在
    shap_values, X_test = [], []
    for fold in range(5):
        test = pd.read_csv(data_path + "cross_valid/{}/{}/X_test.csv".format(fold, exp_Name), usecols=cols)
        X_test.append(test)
        shap_result = np.load(multiview_result_path + "{}/{}/{}_test_shap.npy".format(fold, exp_Name, model_name))
        shap_values.append(shap_result)

    X_test = pd.concat(X_test, ignore_index=True, axis=0)
    shap_values = np.concatenate(shap_values, axis=0)

    # 绘制 SHAP 图像
    plt.figure(figsize=(20, 8), dpi=720)
    plt.subplot(1, 2, 1)

    shap.plots.violin(shap_values, features=X_test, feature_names=cols, plot_type="layered_violin", show=False, max_display=20)
    plt.xlabel('SHAP value')
    plt.subplot(1, 2, 2)
    shap.summary_plot(shap_values, features=X_test, feature_names=cols, plot_type="bar", show=False, max_display=20)
    plt.xlabel('Mean(|SHAP|)')
    
    # 确保图片保存路径存在
    save_dir = os.path.join(multiview_result_path, exp_Name)
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    plt.savefig(os.path.join(save_dir, f"{model_name}_SHAP.png"), dpi=720)
    plt.clf()  # 清除图形

    # 计算特征重要性并保存
    feature_importance = np.absolute(shap_values).mean(axis=0)
    get_topN_reason(exp_Name, feature_importance, cols, model_name, top_num=len(cols))


if __name__ == '__main__':
    pass
    # exps = ['ALL_cluster12_vs_0', 'ALL_cluster1_vs_0', 'ALL_cluster2_vs_0']
    # for exp_Name in exps:
    #     plot_shap(exp_Name, 'Hae_Bio', Hae_Bio)
    #     plot_shap(exp_Name, 'covr10', covr10)
    #     plot_shap(exp_Name, 'Hae_Bio_covr10', Hae_Bio + covr10)
    #     plot_shap(exp_Name, 'Top17', Top17)
    #     plot_shap(exp_Name, 'Top17_covr10', Top17 + covr10)
    #     plot_shap(exp_Name, 'Top17_covr3', Top17 + covr3)
    #
    #     plot_shap(exp_Name, 'Hae_Bio_PRS', Hae_Bio + ['PRS_Score_5e-02'])
    #     plot_shap(exp_Name, 'covr10_PRS', covr10 + ['PRS_Score_5e-02'])
    #     plot_shap(exp_Name, 'Hae_Bio_covr10_PRS', Hae_Bio + covr10 + ['PRS_Score_5e-02'])
    #     plot_shap(exp_Name, 'Top17_PRS', Top17 + ['PRS_Score_5e-02'])
    #     plot_shap(exp_Name, 'Top17_covr10_PRS', Top17 + covr10 + ['PRS_Score_5e-02'])
    #     plot_shap(exp_Name, 'Top17_covr3_PRS', Top17 + covr3 + ['PRS_Score_5e-02'])
    #
    #
    # exps = ['Male_cluster12_vs_0', 'Male_cluster1_vs_0', 'Male_cluster2_vs_0',
    #             'Female_cluster12_vs_0', 'Female_cluster1_vs_0', 'Female_cluster2_vs_0']
    # covr10_no_Gender, cat_no_Gender, covr3_no_Gender = covr10.copy(), cat_cols.copy(), covr3.copy()
    # covr10_no_Gender.remove('Sex_cate'), cat_no_Gender.remove('Sex_cate'), covr3_no_Gender.remove('Sex_cate')
    # for exp_Name in exps:
    #     plot_shap(exp_Name, 'Hae_Bio', Hae_Bio)
    #     plot_shap(exp_Name, 'covr10', covr10_no_Gender)
    #     plot_shap(exp_Name, 'Hae_Bio_covr10', Hae_Bio + covr10_no_Gender)
    #     plot_shap(exp_Name, 'Top17', Top17)
    #     plot_shap(exp_Name, 'Top17_covr10', Top17 + covr10_no_Gender)
    #     plot_shap(exp_Name, 'Top17_covr3', Top17 + covr3_no_Gender)
    #
    #     plot_shap(exp_Name, 'Hae_Bio_PRS', Hae_Bio + ['PRS_Score_5e-02'])
    #     plot_shap(exp_Name, 'covr10_PRS', covr10_no_Gender + ['PRS_Score_5e-02'])
    #     plot_shap(exp_Name, 'Hae_Bio_covr10_PRS', Hae_Bio + covr10_no_Gender + ['PRS_Score_5e-02'])
    #     plot_shap(exp_Name, 'Top17_PRS', Top17 + ['PRS_Score_5e-02'])
    #     plot_shap(exp_Name, 'Top17_covr10_PRS', Top17 + covr10_no_Gender + ['PRS_Score_5e-02'])
    #     plot_shap(exp_Name, 'Top17_covr3_PRS', Top17 + covr3_no_Gender + ['PRS_Score_5e-02'])
    #
    #
    # exps = ['BMI_cluster12_vs_0', 'BMI_cluster1_vs_0', 'BMI_cluster2_vs_0']
    # covr10_no_BMI, cat_no_BMI, covr3_no_BMI = covr10.copy(), cat_cols.copy(), covr3.copy()
    # covr10_no_BMI.remove('BMI_cate'), cat_no_BMI.remove('BMI_cate'), covr3_no_BMI.remove('BMI_cate')
    # for exp_Name in exps:
    #     plot_shap(exp_Name, 'Hae_Bio', Hae_Bio)
    #     plot_shap(exp_Name, 'covr10', covr10_no_BMI)
    #     plot_shap(exp_Name, 'Hae_Bio_covr10', Hae_Bio + covr10_no_BMI)
    #     plot_shap(exp_Name, 'Top17', Top17)
    #     plot_shap(exp_Name, 'Top17_covr10', Top17 + covr10_no_BMI)
    #     plot_shap(exp_Name, 'Top17_covr3', Top17 + covr3_no_BMI)
    #
    #     plot_shap(exp_Name, 'Hae_Bio_PRS', Hae_Bio + ['PRS_Score_5e-02'])
    #     plot_shap(exp_Name, 'covr10_PRS', covr10_no_BMI + ['PRS_Score_5e-02'])
    #     plot_shap(exp_Name, 'Hae_Bio_covr10_PRS', Hae_Bio + covr10_no_BMI + ['PRS_Score_5e-02'])
    #     plot_shap(exp_Name, 'Top17_PRS', Top17 + ['PRS_Score_5e-02'])
    #     plot_shap(exp_Name, 'Top17_covr10_PRS', Top17 + covr10_no_BMI + ['PRS_Score_5e-02'])
    #     plot_shap(exp_Name, 'Top17_covr3_PRS', Top17 + covr3_no_BMI + ['PRS_Score_5e-02'])
    #
    # exps = ['BMI_Male_cluster12_vs_0', 'BMI_Male_cluster1_vs_0', 'BMI_Male_cluster2_vs_0',
    #         'BMI_Female_cluster12_vs_0', 'BMI_Female_cluster1_vs_0', 'BMI_Female_cluster2_vs_0']
    # covr10_no_BMI_Sex, cat_no_BMI_Sex, covr3_no_BMI_Sex = covr10.copy(), cat_cols.copy(), covr3.copy()
    # covr10_no_BMI_Sex.remove('BMI_cate'), cat_no_BMI_Sex.remove('BMI_cate'), covr3_no_BMI_Sex.remove('BMI_cate')
    # covr10_no_BMI_Sex.remove('Sex_cate'), cat_no_BMI_Sex.remove('Sex_cate'), covr3_no_BMI_Sex.remove('Sex_cate')
    # for exp_Name in exps:
    #     plot_shap(exp_Name, 'Hae_Bio', Hae_Bio)
    #     plot_shap(exp_Name, 'covr10', covr10_no_BMI_Sex)
    #     plot_shap(exp_Name, 'Hae_Bio_covr10', Hae_Bio + covr10_no_BMI_Sex)
    #     plot_shap(exp_Name, 'Top17', Top17)
    #     plot_shap(exp_Name, 'Top17_covr10', Top17 + covr10_no_BMI_Sex)
    #     plot_shap(exp_Name, 'Top17_covr3', Top17 + covr3_no_BMI_Sex)
    #
    #     plot_shap(exp_Name, 'Hae_Bio_PRS', Hae_Bio + ['PRS_Score_5e-02'])
    #     plot_shap(exp_Name, 'covr10_PRS', covr10_no_BMI_Sex + ['PRS_Score_5e-02'])
    #     plot_shap(exp_Name, 'Hae_Bio_covr10_PRS', Hae_Bio + covr10_no_BMI_Sex + ['PRS_Score_5e-02'])
    #     plot_shap(exp_Name, 'Top17_PRS', Top17 + ['PRS_Score_5e-02'])
    #     plot_shap(exp_Name, 'Top17_covr10_PRS', Top17 + covr10_no_BMI_Sex + ['PRS_Score_5e-02'])
    #     plot_shap(exp_Name, 'Top17_covr3_PRS', Top17 + covr3_no_BMI_Sex + ['PRS_Score_5e-02'])
    #
    #
    # exps = ['Chronic_cluster12_vs_0', 'Chronic_cluster1_vs_0', 'Chronic_cluster2_vs_0']
    # covr10_no_Chronic = covr10.copy()
    # covr10_no_Chronic.remove('chronic_num')
    # for exp_Name in exps:
    #     plot_shap(exp_Name, 'Hae_Bio', Hae_Bio)
    #     plot_shap(exp_Name, 'covr10', covr10_no_Chronic)
    #     plot_shap(exp_Name, 'Hae_Bio_covr10', Hae_Bio + covr10_no_Chronic)
    #     plot_shap(exp_Name, 'Top17', Top17)
    #     plot_shap(exp_Name, 'Top17_covr10', Top17 + covr10_no_Chronic)
    #     plot_shap(exp_Name, 'Top17_covr3', Top17 + covr3)
    #
    #     plot_shap(exp_Name, 'Hae_Bio_PRS', Hae_Bio + ['PRS_Score_5e-02'])
    #     plot_shap(exp_Name, 'covr10_PRS', covr10_no_Chronic + ['PRS_Score_5e-02'])
    #     plot_shap(exp_Name, 'Hae_Bio_covr10_PRS', Hae_Bio + covr10_no_Chronic + ['PRS_Score_5e-02'])
    #     plot_shap(exp_Name, 'Top17_PRS', Top17 + ['PRS_Score_5e-02'])
    #     plot_shap(exp_Name, 'Top17_covr10_PRS', Top17 + covr10_no_Chronic + ['PRS_Score_5e-02'])
    #     plot_shap(exp_Name, 'Top17_covr3_PRS', Top17 + covr3 + ['PRS_Score_5e-02'])
    #
    #
    # exps = ['Chronic_Male_cluster12_vs_0', 'Chronic_Male_cluster1_vs_0', 'Chronic_Male_cluster2_vs_0',
    #         'Chronic_Female_cluster12_vs_0', 'Chronic_Female_cluster1_vs_0', 'Chronic_Female_cluster2_vs_0']
    # covr10_no_Chronic_Sex, cat_no_Sex, covr3_no_Sex = covr10.copy(), cat_cols.copy(), covr3.copy()
    # covr10_no_Chronic_Sex.remove('Sex_cate'), cat_no_Sex.remove('Sex_cate'), covr3_no_Sex.remove('Sex_cate')
    # covr10_no_Chronic_Sex.remove('chronic_num')
    # for exp_Name in exps:
    #     plot_shap(exp_Name, 'Hae_Bio', Hae_Bio)
    #     plot_shap(exp_Name, 'covr10', covr10_no_Chronic_Sex)
    #     plot_shap(exp_Name, 'Hae_Bio_covr10', Hae_Bio + covr10_no_Chronic_Sex)
    #     plot_shap(exp_Name, 'Top17', Top17)
    #     plot_shap(exp_Name, 'Top17_covr10', Top17 + covr10_no_Chronic_Sex)
    #     plot_shap(exp_Name, 'Top17_covr3', Top17 + covr3_no_Sex)
    #
    #     plot_shap(exp_Name, 'Hae_Bio_PRS', Hae_Bio + ['PRS_Score_5e-02'])
    #     plot_shap(exp_Name, 'covr10_PRS', covr10_no_Chronic_Sex + ['PRS_Score_5e-02'])
    #     plot_shap(exp_Name, 'Hae_Bio_covr10_PRS', Hae_Bio + covr10_no_Chronic_Sex + ['PRS_Score_5e-02'])
    #     plot_shap(exp_Name, 'Top17_PRS', Top17 + ['PRS_Score_5e-02'])
    #     plot_shap(exp_Name, 'Top17_covr10_PRS', Top17 + covr10_no_Chronic_Sex + ['PRS_Score_5e-02'])
    #     plot_shap(exp_Name, 'Top17_covr3_PRS', Top17 + covr3_no_Sex + ['PRS_Score_5e-02'])
