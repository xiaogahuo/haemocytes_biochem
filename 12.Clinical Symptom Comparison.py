from scipy.stats.mstats import kruskal
import pandas as pd
from SUtils import *
import math






def PHQ_plot(df):
    df['final_cluster'] = df['final_cluster'].replace({0: 'HC', 1: 'Subtpye1', 2: 'Subtpye2'})
    df_PHQ = df[['20127-0.0', '2050-0.0', '2060-0.0', '2070-0.0', '2080-0.0', 'PHQ', 'PHQ_Total', 'final_cluster']]
    for col in ['20127-0.0', '2050-0.0', '2060-0.0', '2070-0.0', '2080-0.0', 'PHQ', 'PHQ_Total']:
        plot_violin(df_PHQ, col, [col], 2.8, 4, ['HC', 'Subtpye1', 'Subtpye2'], 'Score')

    df_1160 = df[['1160-0.0', 'final_cluster']]
    # df_1160 = df_1160[(df_1160['1160-0.0'] > -4) & (df_1160['1160-0.0'] < 14)]
    plot_violin(df_1160, '1160-0.0', ['1160-0.0'], 2.8, 4, ['HC', 'Subtpye1', 'Subtpye2'], 'Score')


if __name__ == '__main__':
    df_MDD, df_NC = pd.read_csv(data_path + "Diag_MDD-0.0.csv", low_memory=False), pd.read_csv(data_path + "Diag_NC-0.0.csv", low_memory=False)
    df_NC['final_cluster'] = 0
    df_Y = load_final_cluster(2)
    df_MDD = pd.merge(df_MDD, df_Y, how='inner', on='eid', sort=True, suffixes=('_x', '_y'), copy=True, indicator=False)
    df = pd.concat([df_MDD, df_NC])
    print(df.columns.tolist())
    PHQ_plot(df)





