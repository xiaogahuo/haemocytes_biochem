from SUtils import *
import matplotlib.pyplot as plt
import seaborn as sns


def plot_heatmap(file_name, sheet_name):
    df = pd.read_excel(multiview_result_path + "{}.xlsx".format(file_name), sheet_name=sheet_name, index_col=0)  # Estimate
    print("max: {}, min: {}".format(df.max().max(), df.min().min()))
    plt.figure(figsize=(6, 20))
    myPlot = sns.heatmap(df, annot=df, fmt='.2f', annot_kws={"size": 9, "color": "white", "weight": "bold"}, yticklabels=True, cmap="vlag", center=0.75, vmin=0.5, vmax=1.0, linewidth=.5)
    myPlot.xaxis.tick_top()
    myPlot.set_xticklabels(myPlot.get_xticklabels(), rotation=45, ha='left')
    plt.tight_layout()
    plt.savefig(multiview_result_path + "{}_{}_Heatmap.png".format(file_name, sheet_name))
    plt.show()
    plt.clf()


if __name__ == '__main__':
    # for file in ['UKB', 'DPUK']:
    for file in ['AUC_Sensi_Speci']:
        for sheet in ['No_PRS', 'PRS']:
            plot_heatmap(file, sheet)