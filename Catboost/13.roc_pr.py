# -*- coding: utf-8 -*-
import collections
import sys
import os
sys.path.append(os.path.abspath(".."))
from SUtils import *
from sklearn.metrics import auc, precision_recall_curve, average_precision_score
import matplotlib.pyplot as plt


def get_sensitivity_specificity(y_true, y_scores, threshold):
    # print('{} {}:'.format(type(y_true), type(y_scores)))
    y_pred = (y_scores >= threshold).astype(int)
    tp = np.sum((y_true == 1) & (y_pred == 1))
    tn = np.sum((y_true == 0) & (y_pred == 0))
    fp = np.sum((y_true == 0) & (y_pred == 1))
    fn = np.sum((y_true == 1) & (y_pred == 0))
    print('tp={}, tn={}, fp={}, fn={}'.format(tp, tn, fp, fn))

    precision = tp / (tp + fp)
    sensitivity = tp / (tp + fn)
    F1_score = 2*(precision * sensitivity) / (precision + sensitivity)
    specificity = tn / (tn + fp)
    Acc = (tp + tn) / (tp + tn + fp + fn)
    return precision, sensitivity, specificity, Acc, F1_score

def get_score_label(fold, exp_Name, model_name, d_type):
    # get the raw scores and labels from the test csvfile for the ROC PR curves
    df = pd.read_csv(multiview_result_path + "{}/{}/{}_{}_eval.csv".format(fold, exp_Name, model_name, d_type))
    return df['pred_prob'], df['label']


def generate_roc_pr(exp_Name, model_name, colors):
    lw = 1
    text_size = 18
    plt.figure(figsize=(14, 6), dpi=720)
    plt.subplot(1, 2, 1)
    types = ['train', 'test']
    mean_fpr = np.linspace(0, 1, 100)
    for i, d_type in enumerate(types):
        color = colors[i]
        tprs, aucs, cutoffs, sensitivities, specificities = [], [], [], [], []
        for fold in range(5):
            scores, labels = get_score_label(fold, exp_Name, model_name, d_type)
            fpr, tpr, thres = roc_curve(labels, scores, pos_label=1)
            AUC = auc(fpr, tpr)
            plt.plot(fpr, tpr, lw=lw / 2, alpha=0.15)
            interp_tpr = np.interp(mean_fpr, fpr, tpr)
            interp_tpr[0] = 0.0
            tprs.append(interp_tpr)
            aucs.append(AUC)
            cutoff = get_best_cutoff(labels, scores)
            precision, sensitivity, specificity, Acc, F1_score = get_sensitivity_specificity(labels, scores, cutoff)
            print('{}: fold{}: {}: best cutoff={:.3f}, sensitivity={:.3f}, specificity={:.3f}'.format(exp_Name, fold, d_type, cutoff, sensitivity, specificity))
            cutoffs.append(cutoff)
            sensitivities.append(sensitivity)
            specificities.append(specificity)
        mean_tpr = np.mean(tprs, axis=0)
        mean_tpr[-1] = 1.0
        mean_auc = auc(mean_fpr, mean_tpr)
        std_auc = np.std(aucs)
        plt.plot(mean_fpr, mean_tpr, color=color,
                 label='{} AUC={:.3f}$\pm${:.3f}'.format(d_type, mean_auc, std_auc),
                 lw=2, alpha=.8)
        std_tpr = np.std(tprs, axis=0)
        tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
        tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
        plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color=color, alpha=.2)
        print('{}: {}: mean best cutoff={:.3f}$\pm${:.3f}, mean sensitivity={:.3f}$\pm${:.3f}, mean specificity={:.3f}$\pm${:.3f}'.format(exp_Name, d_type, np.mean(cutoffs, axis=0), np.std(cutoffs), np.mean(sensitivities, axis=0), np.std(sensitivities), np.mean(specificities, axis=0), np.std(specificities)))
    plt.plot([0, 1], [0, 1], 'k--', lw=lw)
    plt.xlim=(-0.05, 1.05)
    plt.ylim=(-0.05, 1.05)
    legend_properties = {'weight': 'bold', 'size': 12}
    plt.xlabel('False Positive Rate', fontsize=text_size, fontweight='bold')
    plt.ylabel('True Positive Rate', fontsize=text_size, fontweight='bold')
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.legend(loc="lower right", prop=legend_properties)

    plt.subplot(1, 2, 2)
    mean_rc = np.linspace(0, 1, 100)
    for i, d_type in enumerate(types):
        color = colors[i]
        prs, aucs,  labels= [], [], []
        for fold in range(5):
            scores, labels = get_score_label(fold, exp_Name, model_name, d_type)
            pr, rc, thres = precision_recall_curve(labels, scores, pos_label=1)
            pr, rc = pr[::-1], rc[::-1]
            AUC = average_precision_score(labels, scores, pos_label=1)
            plt.plot(rc, pr, lw=lw / 2, alpha=0.15)
            interp_pr = np.interp(mean_rc, rc, pr)
            prs.append(interp_pr)
            aucs.append(AUC)
        mean_pr = np.mean(prs, axis=0)
        mean_auc = np.mean(aucs)
        std_auc = np.std(aucs)
        plt.plot(mean_rc, mean_pr, color=color,
                 label='{} AP={:.3f}$\pm${:.3f}'.format(d_type, mean_auc, std_auc),
                 lw=2, alpha=.8)
        count = collections.Counter(labels)
        ratio = count[1] / (count[1] + count[0])
        plt.plot([0, 1], [ratio, ratio], 'k--', lw=lw)
        std_pr = np.std(prs, axis=0)
        prs_upper = np.minimum(mean_pr + std_pr, 1)
        prs_lower = np.maximum(mean_pr - std_pr, 0)
        plt.fill_between(mean_rc, prs_lower, prs_upper, color=color, alpha=.2)
    legend_properties = {'weight': 'bold', 'size': 12}
    plt.yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    plt.xlim=(-0.05, 1.05)
    plt.ylim=(-0.05, 1.05)
    plt.xlabel('Recall', fontsize=text_size, fontweight='bold')
    plt.ylabel('Precision', fontsize=text_size, fontweight='bold')
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.legend(loc="upper right", prop=legend_properties)

    if not os.path.exists(multiview_result_path+"{}/".format(exp_Name)):
        os.makedirs(multiview_result_path+"{}/".format(exp_Name))
    plt.savefig(multiview_result_path+"{}/{}_ROC_PR.png".format(exp_Name, model_name), bbox_inches='tight')
    plt.clf()
    plt.close()
    return


def generate_loss_curve(exp_Name, model_name, colors):
    lw = 1
    text_size = 18
    plt.figure(figsize=(14, 6), dpi=720)
    plt.subplot(1, 2, 1)
    tpyes = ['train', 'test']
    for i, d_type in enumerate(tpyes):
        color = colors[i]
        loss_5folds = []
        for fold in range(5):
            loss = pd.read_csv(multiview_result_path + "{}/{}/{}_{}_loss.csv".format(fold, exp_Name, model_name, d_type))['loss']
            loss_5folds.append(loss)
        mean_loss = np.mean(loss_5folds, axis=0)
        plt.plot(np.arange(0, len(mean_loss)), mean_loss, color=color, label= d_type, lw=lw, alpha=.8)
    legend_properties = {'weight': 'bold', 'size': 12}
    plt.xlabel('Epoch', fontsize=text_size, fontweight='bold')
    plt.ylabel('Logloss', fontsize=text_size, fontweight='bold')
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.legend(loc="upper right", prop=legend_properties)
    plt.savefig(multiview_result_path + "{}/{}_Loss_Curve.png".format(exp_Name, model_name), bbox_inches='tight')
    plt.clf()
    plt.close()
    return


def generate_sens_speci(exp_Name, model_name, f):
    types = ['train', 'test']
    for i, d_type in enumerate(types):
        all_scores, all_labels = [], []
        for fold in range(5):
            scores, labels = get_score_label(fold, exp_Name, model_name, d_type)
            all_scores.append(scores)
            all_labels.append(labels)
        all_scores = pd.concat(all_scores)
        all_labels = pd.concat(all_labels)
        fpr, tpr, thres = roc_curve(all_labels, all_scores, pos_label=1)
        AUC = auc(fpr, tpr)
        cutoff = get_best_cutoff(all_labels, all_scores)
        precision, sensitivity, specificity, Acc, F1_score = get_sensitivity_specificity(all_labels, all_scores, cutoff)
        print('{},{},{},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f}'.format(exp_Name, model_name, d_type,
                    AUC, precision, sensitivity, specificity, Acc, F1_score), file=f)

def generate_dataset(exp_Name, model_name, f):
    generate_roc_pr(exp_Name, model_name, ['#222222', 'r'])
    generate_sens_speci(exp_Name, model_name, f)
    generate_loss_curve(exp_Name, model_name, ['#222222', 'r'])


if __name__ == "__main__":
    f = open(multiview_result_path + "Sensi_Speci_ALL_PRS.csv", "w")
    header_str = "exp_Name,model_name,d_type,AUC,precision,sensitivity,specificity,Acc,F1_score"
    print(header_str, file=f)

    exp_names = ['ALL_cluster12_vs_0', 'ALL_cluster1_vs_0', 'ALL_cluster2_vs_0',
                 'Male_cluster12_vs_0', 'Male_cluster1_vs_0', 'Male_cluster2_vs_0',
                 'Female_cluster12_vs_0', 'Female_cluster1_vs_0', 'Female_cluster2_vs_0',
                 'BMI_cluster12_vs_0', 'BMI_cluster1_vs_0', 'BMI_cluster2_vs_0',
                 'BMI_Male_cluster12_vs_0', 'BMI_Male_cluster1_vs_0', 'BMI_Male_cluster2_vs_0',
                 'BMI_Female_cluster12_vs_0', 'BMI_Female_cluster1_vs_0', 'BMI_Female_cluster2_vs_0',
                 'Chronic_cluster12_vs_0', 'Chronic_cluster1_vs_0', 'Chronic_cluster2_vs_0',
                 'Chronic_Male_cluster12_vs_0', 'Chronic_Male_cluster1_vs_0', 'Chronic_Male_cluster2_vs_0',
                 'Chronic_Female_cluster12_vs_0', 'Chronic_Female_cluster1_vs_0', 'Chronic_Female_cluster2_vs_0']

    model_names =['Hae_Bio', 'covr10', 'Hae_Bio_covr10', 'Top17', 'Top17_covr10', 'Top17_covr3',
                  'Hae_Bio_PRS', 'covr10_PRS', 'Hae_Bio_covr10_PRS', 'Top17_PRS', 'Top17_covr10_PRS', 'Top17_covr3_PRS']
    for exp in exp_names:
        for model in model_names:
            generate_dataset(exp, model, f)




