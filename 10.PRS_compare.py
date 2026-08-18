# -*- coding: utf-8 -*-
from SUtils import *
import sys


#       #FID      IID  ALLELE_CT  NAMED_ALLELE_DOSAGE_SUM  SCORE1_AVG
# 0  2046642  2046642       1490                  768.988    0.000274
# 1  4270668  4270668       1490                  905.463   -0.002423
# 2  5053312  5053312       1490                  704.180   -0.000345
# 3  1780405  1780405       1490                  721.000   -0.001148
# 4  3030933  3030933       1490                  734.004   -0.000045
def PRS_score_kwtest(df, f_name):   #
    PRS_cols = ['SCORE_5e08', 'SCORE_1e07', 'SCORE_5e07', 'SCORE_1e06', 'SCORE_5e06', 'SCORE_1e05', 'SCORE_5e05',
                'SCORE_1e04', 'SCORE_5e04', 'SCORE_1e03', 'SCORE_5e03', 'SCORE_1e02', 'SCORE_5e02']
    new_cols = ['5e-08', '1e-07', '5e-07', '1e-06', '5e-06', '1e-05', '5e-05', '1e-04', '5e-04', '1e-03', '5e-03',
                '0.01', '0.05']
    for col in PRS_cols:
        print("*"*25+f_name+"*"*25)
        kwtest(col, 'final_cluster', df)


def plot_bar_PRS():
    df = pd.read_csv(data_path + "PRS_SCORE_Dataset.csv", low_memory=False)
    
    # df['final_cluster'] = df['final_cluster'].replace({0: 'HC', 1: 'Subtpye1', 2: 'Subtpye2', 3: 'Subtpye3'})
    df['final_cluster'] = df['final_cluster'].replace({0: 'HC', 1: 'Subtpye1', 2: 'Subtpye2'})
    PRS_cols = ['SCORE_5e08', 'SCORE_1e07', 'SCORE_5e07', 'SCORE_1e06', 'SCORE_5e06', 'SCORE_1e05', 'SCORE_5e05', 'SCORE_1e04', 'SCORE_5e04', 'SCORE_1e03', 'SCORE_5e03', 'SCORE_1e02', 'SCORE_5e02']
    new_cols = ['5e-08', '1e-07', '5e-07', '1e-06', '5e-06', '1e-05', '5e-05', '1e-04', '5e-04', '1e-03', '5e-03', '0.01', '0.05']
    df.rename(columns=dict(zip(PRS_cols, new_cols)), inplace=True)
    print(df.head())
    for col in new_cols:
        df[col] = (df[col] - df[col].mean())/df[col].std()
          
    # plot_bar(df, 'PRS', new_cols, 13, 3, ['HC', 'Subtpye1', 'Subtpye2'], 'PRS_Zscore')
    # plot_bar(df, 'PRS', new_cols, 13, 4, ['HC', 'Subtpye1', 'Subtpye2','Subtpye3'], 'PRS_Zscore')
    plot_bar(df, 'PRS', new_cols, 13, 4, ['HC', 'Subtpye1', 'Subtpye2'], 'PRS_Zscore')


if __name__ == '__main__':
    df_NC = pd.read_csv(data_path + "Diag_NC-0.0.csv", low_memory=False, usecols=['eid'])  #
    df_NC['final_cluster'] = 0  #
    #df_Y = pd.read_csv("/data4/myuan/Hae_Bio/F6_results.csv", low_memory=False, usecols=['eid', 'final_cluster'])   #
    df_Y = load_final_cluster(2)
    df = pd.concat([df_Y, df_NC])
    df = [df]
    for f in ['5e08', '1e07', '5e07', '1e06', '5e06', '1e05', '5e05', '1e04', '5e04', '1e03', '5e03', '1e02', '5e02']:
        df_PRS = pd.read_csv(PRS_path + "{}_PRS.sscore".format(f), sep=r'[ \t]', engine='python',
                             usecols=['#FID', 'SCORE1_AVG'])
        df_PRS.rename(columns={'#FID': 'eid', 'SCORE1_AVG': 'SCORE_' + f}, inplace=True)
        df.append(df_PRS)
    df = reduce(lambda left, right: pd.merge(left, right, on=['eid'], how='inner'), df)
    df.to_csv(data_path + "PRS_SCORE_Dataset.csv", index=False)

    PRS_score_kwtest(df[(df['final_cluster'] == 0) | (df['final_cluster'] == 1)], '0_VS_1')
    PRS_score_kwtest(df[(df['final_cluster'] == 0) | (df['final_cluster'] == 2)], '0_VS_2')
    # PRS_score_kwtest(df[(df['final_cluster'] == 0) | (df['final_cluster'] == 3)], '0_VS_3')
    PRS_score_kwtest(df[(df['final_cluster'] == 1) | (df['final_cluster'] == 2)], '1_VS_2')
    # PRS_score_kwtest(df[(df['final_cluster'] == 1) | (df['final_cluster'] == 3)], '1_VS_3')
    # PRS_score_kwtest(df[(df['final_cluster'] == 2) | (df['final_cluster'] == 3)], '2_VS_3')
    plot_bar_PRS()
