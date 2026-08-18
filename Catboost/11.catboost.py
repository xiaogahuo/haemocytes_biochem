from catboost import CatBoostClassifier
import sys
import os
sys.path.append(os.path.abspath(".."))
from SUtils import *
import csv
import shap
from plot_shap import plot_shap


def run_cat(exp_Name, model_name, cols, cat_features, n_estimators):
    for fold in range(5):
        X_train = pd.read_csv(data_path + "cross_valid/{}/{}/X_train.csv".format(fold, exp_Name))[cols]
        X_test = pd.read_csv(data_path + "cross_valid/{}/{}/X_test.csv".format(fold, exp_Name))[cols]
        y_train = pd.read_csv(data_path + "cross_valid/{}/{}/y_train.csv".format(fold, exp_Name))['final_cluster']
        y_test = pd.read_csv(data_path + "cross_valid/{}/{}/y_test.csv".format(fold, exp_Name))['final_cluster']
        if not os.path.exists(multiview_result_path + "{}/{}/".format(fold, exp_Name)):
            os.makedirs(multiview_result_path + "{}/{}/".format(fold, exp_Name))

        count_0 = (y_train == 0).sum()
        count_1 = (y_train == 1).sum()
        ratio = int(count_0 / count_1)
        print(f"class_weights: {ratio}")
        class_weights = {NC: 1, MDD: ratio}
        model = CatBoostClassifier(subsample=0.55, learning_rate=0.02, max_depth=5, eval_metric='F1',
                                   loss_function='Logloss', cat_features=cat_features, n_estimators=n_estimators, class_weights=class_weights)
        #将性别特征 改成 类别选择  不是浮点数
        cat_features = [col for col in cols if col in cat_cols]
        model.fit(X_train, y_train,cat_features=cat_features, eval_set=[(X_test, y_test)], use_best_model=True, verbose=False)
        # Extract the loss values for training and validation datasets
        results = model.get_evals_result()
        train_loss = results['learn']['Logloss']
        test_loss = results['validation']['Logloss']
        pd.DataFrame({'loss': train_loss}).to_csv(multiview_result_path + "{}/{}/{}_train_loss.csv".format(fold, exp_Name, model_name), index=None)
        pd.DataFrame({'loss': test_loss}).to_csv(multiview_result_path + "{}/{}/{}_test_loss.csv".format(fold, exp_Name, model_name), index=None)

        predict_test_prob = model.predict_proba(X_test)  # 返回预测概率
        df_test = pd.DataFrame()
        df_test['label'] = y_test
        df_test['pred_prob'] = predict_test_prob[:, 1] 
        df_test['pred_label'] = model.predict(X_test)
        df_test.to_csv(multiview_result_path + "{}/{}/{}_test_eval.csv".format(fold, exp_Name, model_name), index=None)

        predict_train_prob = model.predict_proba(X_train)  # 返回预测概率 
        df_train = pd.DataFrame()
        df_train['label'] = y_train
        df_train['pred_prob'] = predict_train_prob[:, 1]
        df_train['pred_label'] = model.predict(X_train)
        df_train.to_csv(multiview_result_path + "{}/{}/{}_train_eval.csv".format(fold, exp_Name, model_name), index=None)
  
        # save SHAP_values
        explainer = shap.TreeExplainer(model)  # #这里的model在准备工作中已经完成建模，模型名称就是model
        shap_values = explainer.shap_values(X_test)  # 传入特征矩阵X，计算SHAP值
        # print("shap_values shape: {}".format(shap_values.shape))
        np.save(multiview_result_path + "{}/{}/{}_test_shap.npy".format(fold, exp_Name, model_name), shap_values)


def run(exp_Name, model_name, cols, cat_features, n_estimators):
    print(f"exp_Name: {exp_Name}, model_name: {model_name}, cols: {cols}, cat_features: {cat_features}, n_estimators: {n_estimators}")
    run_cat(exp_Name, model_name, cols, cat_features, n_estimators)
    plot_shap(exp_Name, model_name, cols)


if __name__ == '__main__':
    exps = ['ALL_cluster12_vs_0', 'ALL_cluster1_vs_0', 'ALL_cluster2_vs_0']
    for exp_Name in exps:
        run(exp_Name, 'Hae_Bio', Hae_Bio, None, 500)
        run(exp_Name, 'covr10', covr10, cat_cols, 500)
        run(exp_Name, 'Hae_Bio_covr10', Hae_Bio + covr10, cat_cols, 500)
        run(exp_Name, 'Top17', Top17, None, 500)
        run(exp_Name, 'Top17_covr10', Top17 + covr10, cat_cols, 500)
        run(exp_Name, 'Top17_covr3', Top17 + covr3, ['BMI_cate', 'Sex_cate'], 500)

        run(exp_Name, 'Hae_Bio_PRS', Hae_Bio + ['PRS_Score_5e-02'], None, 500)
        run(exp_Name, 'covr10_PRS', covr10 + ['PRS_Score_5e-02'], cat_cols, 500)
        run(exp_Name, 'Hae_Bio_covr10_PRS', Hae_Bio + covr10 + ['PRS_Score_5e-02'], cat_cols, 500)
        run(exp_Name, 'Top17_PRS', Top17 + ['PRS_Score_5e-02'], None, 500)
        run(exp_Name, 'Top17_covr10_PRS', Top17 + covr10 + ['PRS_Score_5e-02'], cat_cols, 500)
        run(exp_Name, 'Top17_covr3_PRS', Top17 + covr3 + ['PRS_Score_5e-02'], ['BMI_cate', 'Sex_cate'], 500)


    exps = ['Male_cluster12_vs_0', 'Male_cluster1_vs_0', 'Male_cluster2_vs_0',
                'Female_cluster12_vs_0', 'Female_cluster1_vs_0', 'Female_cluster2_vs_0']
    covr10_no_Gender, cat_no_Gender, covr3_no_Gender = covr10.copy(), cat_cols.copy(), covr3.copy()
    covr10_no_Gender.remove('Sex_cate'), cat_no_Gender.remove('Sex_cate'), covr3_no_Gender.remove('Sex_cate')
    for exp_Name in exps:
        run(exp_Name, 'Hae_Bio', Hae_Bio, None, 500)
        run(exp_Name, 'covr10', covr10_no_Gender, cat_no_Gender, 500)
        run(exp_Name, 'Hae_Bio_covr10', Hae_Bio + covr10_no_Gender, cat_no_Gender, 500)
        run(exp_Name, 'Top17', Top17, None, 500)
        run(exp_Name, 'Top17_covr10', Top17 + covr10_no_Gender, cat_no_Gender, 500)
        run(exp_Name, 'Top17_covr3', Top17 + covr3_no_Gender, ['BMI_cate'], 500)

        run(exp_Name, 'Hae_Bio_PRS', Hae_Bio + ['PRS_Score_5e-02'], None, 500)
        run(exp_Name, 'covr10_PRS', covr10_no_Gender + ['PRS_Score_5e-02'], cat_no_Gender, 500)
        run(exp_Name, 'Hae_Bio_covr10_PRS', Hae_Bio + covr10_no_Gender + ['PRS_Score_5e-02'], cat_no_Gender, 500)
        run(exp_Name, 'Top17_PRS', Top17 + ['PRS_Score_5e-02'], None, 500)
        run(exp_Name, 'Top17_covr10_PRS', Top17 + covr10_no_Gender + ['PRS_Score_5e-02'], cat_no_Gender, 500)
        run(exp_Name, 'Top17_covr3_PRS', Top17 + covr3_no_Gender + ['PRS_Score_5e-02'], ['BMI_cate'], 500)


    exps = ['BMI_cluster12_vs_0', 'BMI_cluster1_vs_0', 'BMI_cluster2_vs_0']
    covr10_no_BMI, cat_no_BMI, covr3_no_BMI = covr10.copy(), cat_cols.copy(), covr3.copy()
    covr10_no_BMI.remove('BMI_cate'), cat_no_BMI.remove('BMI_cate'), covr3_no_BMI.remove('BMI_cate')
    for exp_Name in exps:
        run(exp_Name, 'Hae_Bio', Hae_Bio, None, 500)
        run(exp_Name, 'covr10', covr10_no_BMI, cat_no_BMI, 500)
        run(exp_Name, 'Hae_Bio_covr10', Hae_Bio + covr10_no_BMI, cat_no_BMI, 500)
        run(exp_Name, 'Top17', Top17, None, 500)
        run(exp_Name, 'Top17_covr10', Top17 + covr10_no_BMI, cat_no_BMI, 500)
        run(exp_Name, 'Top17_covr3', Top17 + covr3_no_BMI, ['Sex_cate'], 500)

        run(exp_Name, 'Hae_Bio_PRS', Hae_Bio + ['PRS_Score_5e-02'], None, 500)
        run(exp_Name, 'covr10_PRS', covr10_no_BMI + ['PRS_Score_5e-02'], cat_no_BMI, 500)
        run(exp_Name, 'Hae_Bio_covr10_PRS', Hae_Bio + covr10_no_BMI + ['PRS_Score_5e-02'], cat_no_BMI, 500)
        run(exp_Name, 'Top17_PRS', Top17 + ['PRS_Score_5e-02'], None, 500)
        run(exp_Name, 'Top17_covr10_PRS', Top17 + covr10_no_BMI + ['PRS_Score_5e-02'], cat_no_BMI, 500)
        run(exp_Name, 'Top17_covr3_PRS', Top17 + covr3_no_BMI + ['PRS_Score_5e-02'], ['Sex_cate'], 500)

    exps = ['BMI_Male_cluster12_vs_0', 'BMI_Male_cluster1_vs_0', 'BMI_Male_cluster2_vs_0',
            'BMI_Female_cluster12_vs_0', 'BMI_Female_cluster1_vs_0', 'BMI_Female_cluster2_vs_0']
    covr10_no_BMI_Sex, cat_no_BMI_Sex, covr3_no_BMI_Sex = covr10.copy(), cat_cols.copy(), covr3.copy()
    covr10_no_BMI_Sex.remove('BMI_cate'), cat_no_BMI_Sex.remove('BMI_cate'), covr3_no_BMI_Sex.remove('BMI_cate')
    covr10_no_BMI_Sex.remove('Sex_cate'), cat_no_BMI_Sex.remove('Sex_cate'), covr3_no_BMI_Sex.remove('Sex_cate')
    for exp_Name in exps:
        run(exp_Name, 'Hae_Bio', Hae_Bio, None, 500)
        run(exp_Name, 'covr10', covr10_no_BMI_Sex, cat_no_BMI_Sex, 500)
        run(exp_Name, 'Hae_Bio_covr10', Hae_Bio + covr10_no_BMI_Sex, cat_no_BMI_Sex, 500)
        run(exp_Name, 'Top17', Top17, None, 500)
        run(exp_Name, 'Top17_covr10', Top17 + covr10_no_BMI_Sex, cat_no_BMI_Sex, 500)
        run(exp_Name, 'Top17_covr3', Top17 + covr3_no_BMI_Sex, None, 500)

        run(exp_Name, 'Hae_Bio_PRS', Hae_Bio + ['PRS_Score_5e-02'], None, 500)
        run(exp_Name, 'covr10_PRS', covr10_no_BMI_Sex + ['PRS_Score_5e-02'], cat_no_BMI_Sex, 500)
        run(exp_Name, 'Hae_Bio_covr10_PRS', Hae_Bio + covr10_no_BMI_Sex + ['PRS_Score_5e-02'], cat_no_BMI_Sex, 500)
        run(exp_Name, 'Top17_PRS', Top17 + ['PRS_Score_5e-02'], None, 500)
        run(exp_Name, 'Top17_covr10_PRS', Top17 + covr10_no_BMI_Sex + ['PRS_Score_5e-02'], cat_no_BMI_Sex, 500)
        run(exp_Name, 'Top17_covr3_PRS', Top17 + covr3_no_BMI_Sex + ['PRS_Score_5e-02'], None, 500)


    exps = ['Chronic_cluster12_vs_0', 'Chronic_cluster1_vs_0', 'Chronic_cluster2_vs_0']
    covr10_no_Chronic = covr10.copy()
    covr10_no_Chronic.remove('chronic_num')
    for exp_Name in exps:
        run(exp_Name, 'Hae_Bio', Hae_Bio, None, 500)
        run(exp_Name, 'covr10', covr10_no_Chronic, cat_cols, 500)
        run(exp_Name, 'Hae_Bio_covr10', Hae_Bio + covr10_no_Chronic, cat_cols, 500)
        run(exp_Name, 'Top17', Top17, None, 500)
        run(exp_Name, 'Top17_covr10', Top17 + covr10_no_Chronic, cat_cols, 500)
        run(exp_Name, 'Top17_covr3', Top17 + covr3, ['BMI_cate', 'Sex_cate'], 500)

        run(exp_Name, 'Hae_Bio_PRS', Hae_Bio + ['PRS_Score_5e-02'], None, 500)
        run(exp_Name, 'covr10_PRS', covr10_no_Chronic + ['PRS_Score_5e-02'], cat_cols, 500)
        run(exp_Name, 'Hae_Bio_covr10_PRS', Hae_Bio + covr10_no_Chronic + ['PRS_Score_5e-02'], cat_cols, 500)
        run(exp_Name, 'Top17_PRS', Top17 + ['PRS_Score_5e-02'], None, 500)
        run(exp_Name, 'Top17_covr10_PRS', Top17 + covr10_no_Chronic + ['PRS_Score_5e-02'], cat_cols, 500)
        run(exp_Name, 'Top17_covr3_PRS', Top17 + covr3 + ['PRS_Score_5e-02'], ['BMI_cate', 'Sex_cate'], 500)


    exps = ['Chronic_Male_cluster12_vs_0', 'Chronic_Male_cluster1_vs_0', 'Chronic_Male_cluster2_vs_0',
            'Chronic_Female_cluster12_vs_0', 'Chronic_Female_cluster1_vs_0', 'Chronic_Female_cluster2_vs_0']
    covr10_no_Chronic_Sex, cat_no_Sex, covr3_no_Sex = covr10.copy(), cat_cols.copy(), covr3.copy()
    covr10_no_Chronic_Sex.remove('Sex_cate'), cat_no_Sex.remove('Sex_cate'), covr3_no_Sex.remove('Sex_cate')
    covr10_no_Chronic_Sex.remove('chronic_num')
    for exp_Name in exps:
        run(exp_Name, 'Hae_Bio', Hae_Bio, None, 500)
        run(exp_Name, 'covr10', covr10_no_Chronic_Sex, cat_no_Sex, 500)
        run(exp_Name, 'Hae_Bio_covr10', Hae_Bio + covr10_no_Chronic_Sex, cat_no_Sex, 500)
        run(exp_Name, 'Top17', Top17, None, 500)
        run(exp_Name, 'Top17_covr10', Top17 + covr10_no_Chronic_Sex, cat_no_Sex, 500)
        run(exp_Name, 'Top17_covr3', Top17 + covr3_no_Sex, ['BMI_cate'], 500)

        run(exp_Name, 'Hae_Bio_PRS', Hae_Bio + ['PRS_Score_5e-02'], None, 500)
        run(exp_Name, 'covr10_PRS', covr10_no_Chronic_Sex + ['PRS_Score_5e-02'], cat_no_Sex, 500)
        run(exp_Name, 'Hae_Bio_covr10_PRS', Hae_Bio + covr10_no_Chronic_Sex + ['PRS_Score_5e-02'], cat_no_Sex, 500)
        run(exp_Name, 'Top17_PRS', Top17 + ['PRS_Score_5e-02'], None, 500)
        run(exp_Name, 'Top17_covr10_PRS', Top17 + covr10_no_Chronic_Sex + ['PRS_Score_5e-02'], cat_no_Sex, 500)
        run(exp_Name, 'Top17_covr3_PRS', Top17 + covr3_no_Sex + ['PRS_Score_5e-02'], ['BMI_cate'], 500)
