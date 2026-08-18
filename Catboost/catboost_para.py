from catboost import CatBoostClassifier
import csv
import shap
from sklearn.model_selection import GridSearchCV
import pandas as pd
import sys
import os
sys.path.append(os.path.abspath(".."))
from SUtils import *


def run_grid_search(exp_Name, model_name, cols, n_estimators, class_weights):
  
    best_params_list = []
    
    # 定义参数网格
    param_grid = {
        'n_estimators': [500, 600, 700],  # 尝试不同的树的数量
        'learning_rate': [0.01, 0.02, 0.05],  # 不同的学习率
        'max_depth': [5],  # 决策树的深度
        'subsample': [0.5,0.55,0.6],  # 不同的子采样比例
       
    }

    for fold in range(5):
        # 定义每一折的结果路径
        fold_result_path = os.path.join(result_path, str(fold), exp_Name)
        if not os.path.exists(fold_result_path):
            os.makedirs(fold_result_path)

        # 读取训练和测试数据
        X_train = pd.read_csv(data_path + "cross_valid/{}/{}/X_train.csv".format(fold, exp_Name))[cols]
        y_train = pd.read_csv(data_path + "cross_valid/{}/{}/y_train.csv".format(fold, exp_Name))['final_cluster']

        # 初始化 CatBoostClassifier
        model = CatBoostClassifier(loss_function='Logloss', eval_metric='F1', cat_features=None, verbose=0)

        # 使用 GridSearchCV 进行参数搜索
        grid_search = GridSearchCV(estimator=model, param_grid=param_grid, cv=5, scoring='f1', verbose=1)
        grid_search.fit(X_train, y_train)

        # 保存最佳参数
        best_params = grid_search.best_params_
        best_params_list.append(best_params)
        print("Fold {}: Best parameters found: {}".format(fold, best_params))
        print("Fold {}: Best F1 score: {}".format(fold, grid_search.best_score_))

    # 保存所有折的最佳参数到 CSV 文件
    final_result_path = os.path.join(result_path, exp_Name)
    if not os.path.exists(final_result_path):
        os.makedirs(final_result_path)
    best_params_df = pd.DataFrame(best_params_list)
    best_params_df.to_csv(os.path.join(final_result_path, "{}_best_params.csv".format(model_name)), index=False)

if __name__ == '__main__':
    all_cols = Haemocytes_cols + Biochem_cols
    # 更新 `cols_l` 为包含所有特征的列表
    cols_l = [all_cols]
    # 现在 `filtered_feature_lists` 中的特征列表都已过滤掉缺失列
    # 分配给 `cols_l`，与 `m_names` 对应
    # 对 `cols_l` 中的每个特征列表进行去重
    cols_l = [list(set(cols)) for cols in cols_l]

    m_names = ["catboost"]

    for model_name, cols in zip(m_names, cols_l):

        run_grid_search("cluster12_vs_0", model_name, cols, 700, {NC: 1, MDD: 20})
        run_grid_search("cluster1_vs_other", model_name, cols, 500, {NC: 1, MDD: 45})
        run_grid_search("cluster2_vs_other", model_name, cols, 500, {NC: 1, MDD: 45})
