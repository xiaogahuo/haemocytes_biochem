import pandas as pd
from SUtils import *


##############################################
#
# 执行数据预处理，包括计算随访时长、家族史
# "189-0.0": "Townsend_index"  Townsend deprivation index at recruitment Q1 (≤ -3.63) ;Q2 (-3.63 - -2.12);Q3 (-2.12–0.58) ;Q4 (>0.58);Unknown
# "1558-0.0": "Alcohol_freq",
# "6138-0.0": "Qualifications",    Qualifications  Non-university ;University ;Unknown
#    1	College or University degree
#    2	A levels/AS levels or equivalent
#    3	O levels/GCSEs or equivalent
#    4	CSEs or equivalent
#    5	NVQ or HND or HNC or equivalent
#    6	Other professional qualifications eg: nursing, teaching
#    -7	None of the above
#    -3	Prefer not to answer
# "20116-0.0": "Smoking_status",
# "20117-0.0": "Alcohol drinker status",
# "21000-0.0": "Ethnic_background",   Non-White; White; Unknown   # https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=1001
# "21001-0.0": "BMI",  <18.5 kg/m2 ;18.5–24.9 kg/m2;25.0–29.9 kg/m2 ;≥30 kg/m2 ;Unknown
# "21003-0.0": "Age",  <45  45-59  >=60
# "22032-0.0": "IPAQ"  Low Moderate High  Unknown
#
##############################################
# 0 Low 1 Moderate 2 High  Unknown
def sex(x):
    val = x['31-0.0']
    if pd.isnull(val):
        return 'L3_Unknown'
    elif val == 0:
        return 'L1_Female'
    elif val == 1:
        return 'L2_Male'

# 0 Low 1 Moderate 2 High  Unknown
def ipaq(x):
    val = x['22032-0.0'] # have only -0.0 data
    if pd.isnull(val):
        return 'L4_Unknown'
    elif val == 0:
        return 'L1_Low'
    elif val == 1:
        return 'L2_Moderate'
    else:
        return 'L3_High'


# <60 >=60
def age2(x):
    val = x['21003-0.0']
    if val < 60:
        return 'Age1'
    else:
        return 'Age2'


# <18.5 kg/m2 ;18.5–24.9 kg/m2;25.0–29.9 kg/m2 ;≥30 kg/m2 ;Unknown
def bmi(x, suffix):
    val = x['21001'+suffix]
    if pd.isnull(val):
        return 'BMI5_Unknown'
    elif val < 18.5:
        return 'BMI2'
    elif 18.5 <= val < 25.0:
        return 'BMI1'
    elif 25.0 <= val < 30:
        return 'BMI3'
    else:
        return 'BMI4'

def ethnic_background(x, suffix): # Non-White; White; Unknown
    val = x['21000'+suffix]
    if pd.isnull(val) or val in [-1, -3]:
        return 'L3_Unknown'
    elif val in [1, 1001, 2001, 3001, 4001]:
        return 'L2_White'
    else:
        return 'L1_Non-White'


def qualifications(x, suffix): # Non-university ;University ;Unknown
    val = x['6138'+suffix]
    if pd.isnull(val) or val in [-7, -3]:
        return 'L3_Unknown'
    elif val == 1:
        return 'L2_University'
    else:
        return 'L1_Non-university'


def smoking(x, suffix): # -3	Prefer not to answer ;  0 Never ; 1	Previous; 2	Current;
    val = x['20116'+suffix]
    if pd.isnull(val) or val == -3:
        return 'L4_Unknown'
    elif val == 0:
        return 'L1_Never'
    elif val == 1:
        return 'L2_Previous'
    else:
        return 'L3_Current'


def alcohol(x, suffix): # -3	Prefer not to answer ;  0 Never ; 1	Previous; 2	Current;
    val = x['20117'+suffix]
    if pd.isnull(val) or val == -3:
        return 'L4_Unknown'
    elif val == 0:
        return 'L1_Never'
    elif val == 1:
        return 'L2_Previous'
    else:
        return 'L3_Current'


def townsend_index(x):
    val = x['189-0.0'] # have only -0.0 data
    if pd.isnull(val):
        return 'Q5_unknown'
    elif val <= -3.63:
        return 'Q1'
    elif -3.63 < val <= -2.12:
        return 'Q2'
    elif -2.12 < val <= 0.58:
        return 'Q3'
    else:
        return 'Q4'

field_20107 = ['20107-0.0', '20107-0.1', '20107-0.2', '20107-0.3', '20107-0.4', '20107-0.5', '20107-0.6', '20107-0.7', '20107-0.8', '20107-0.9', '20107-1.0', '20107-1.1', '20107-1.2', '20107-1.3', '20107-1.4', '20107-1.5', '20107-1.6', '20107-1.7', '20107-1.8', '20107-1.9', '20107-2.0', '20107-2.1', '20107-2.2', '20107-2.3', '20107-2.4', '20107-2.5', '20107-2.6', '20107-2.7', '20107-2.8', '20107-2.9', '20107-3.0', '20107-3.1', '20107-3.2', '20107-3.3', '20107-3.4', '20107-3.5', '20107-3.6', '20107-3.7', '20107-3.8', '20107-3.9']
field_20110 = ['20110-0.0', '20110-0.1', '20110-0.2', '20110-0.3', '20110-0.4', '20110-0.5', '20110-0.6', '20110-0.7', '20110-0.8', '20110-0.9', '20110-0.10', '20110-1.0', '20110-1.1', '20110-1.2', '20110-1.3', '20110-1.4', '20110-1.5', '20110-1.6', '20110-1.7', '20110-1.8', '20110-1.9', '20110-1.10', '20110-2.0', '20110-2.1', '20110-2.2', '20110-2.3', '20110-2.4', '20110-2.5', '20110-2.6', '20110-2.7', '20110-2.8', '20110-2.9', '20110-2.10', '20110-3.0', '20110-3.1', '20110-3.2', '20110-3.3', '20110-3.4', '20110-3.5', '20110-3.6', '20110-3.7', '20110-3.8', '20110-3.9', '20110-3.10']
field_20111 = ['20111-0.0', '20111-0.1', '20111-0.2', '20111-0.3', '20111-0.4', '20111-0.5', '20111-0.6', '20111-0.7', '20111-0.8', '20111-0.9', '20111-0.10', '20111-0.11', '20111-1.0', '20111-1.1', '20111-1.2', '20111-1.3', '20111-1.4', '20111-1.5', '20111-1.6', '20111-1.7', '20111-1.8', '20111-1.9', '20111-1.10', '20111-1.11', '20111-2.0', '20111-2.1', '20111-2.2', '20111-2.3', '20111-2.4', '20111-2.5', '20111-2.6', '20111-2.7', '20111-2.8', '20111-2.9', '20111-2.10', '20111-2.11', '20111-3.0', '20111-3.1', '20111-3.2', '20111-3.3', '20111-3.4', '20111-3.5', '20111-3.6', '20111-3.7', '20111-3.8', '20111-3.9', '20111-3.10', '20111-3.11']
family_history_column = field_20107 + field_20110 + field_20111
# coverate_family_history = ['20107-', '20110-', '20111-']  # if contain 12  set 1   others   set 0
def family_history(x):
    if 12 in x[family_history_column].tolist():
        return 'YES' # 有家族史
    else:
        return 'NO'  # 无家族史

# 算随访时间
# 死亡['40000-0.0', '40000-1.0']  失访 ['191-0.0']  Affective-Disorder 分类4  取最早时间 算随访时长
# Affective-Disorder 分类1 2 3 拿对应的诊断时间算随访时长
# 其他受试拿当前时间算随访时长
# 入组时间 ['53-0.0']
def follow_Duration(row):
    cur_date = '2023-11-23'
    if row['MDD'] == NC:
        temp = [str(t) for t in row[['40000-0.0', '40000-1.0', '191-0.0']].tolist()]
        temp.append(cur_date)
        target_date = min([t for t in temp if not pd.isnull(t)])
    else:  # MDD
        target_date = str(row['130894-0.0'])
    target_date = datetime.strptime(target_date, "%Y-%m-%d")
    base_date = datetime.strptime(str(row['53-0.0']), "%Y-%m-%d")
    return (target_date-base_date).days/30


def get_preprocess_data(suffix):
    df = pd.read_csv(data_path + "Social_MDD{}.csv".format(suffix), low_memory=False)
    df['follow_Duration'] = df.apply(lambda x: follow_Duration(x), axis=1)
    df['Townsend_index_cate'] = df.apply(lambda x: townsend_index(x), axis=1)
    df['Qualifications_cate'] = df.apply(lambda x: qualifications(x, suffix), axis=1)
    df['Smoking_status_cate'] = df.apply(lambda x: smoking(x, suffix), axis=1)
    df['Alcohol_status_cate'] = df.apply(lambda x: alcohol(x, suffix), axis=1)
    # df['Ethnic_background_cate'] = df.apply(lambda x: ethnic_background(x, suffix), axis=1)
    df['BMI_cate'] = df.apply(lambda x: bmi(x, suffix), axis=1)
    # df['Age_cate2'] = df.apply(lambda x: age2(x), axis=1)
    df = df.rename(columns={'21003'+suffix: 'Age'})
    df['IPAQ_cate'] = df.apply(lambda x: ipaq(x), axis=1)
    df['family_history'] = df.apply(lambda x: family_history(x), axis=1)
    df['Sex_cate'] = df.apply(lambda x: sex(x), axis=1)

    df_haemocytes = pd.read_csv(data_path + "Haemocytes{}.csv".format(suffix), usecols=['eid']+Haemocytes)
    df_biochem = pd.read_csv(data_path + "Biochem{}.csv".format(suffix), usecols=['eid']+Biochem)

    df = reduce(lambda left, right: pd.merge(left, right, on=['eid'], how='inner'),
                    [df, df_haemocytes, df_biochem])
    df_PRS = pd.read_csv(PRS_path + "5e02_PRS.sscore", sep=r'[ \t]', engine='python', usecols=['#FID', 'SCORE1_AVG'])
    df_PRS.rename(columns={'#FID': 'eid', 'SCORE1_AVG': 'PRS_Score_5e-02'}, inplace=True)
    df = reduce(lambda left, right: pd.merge(left, right, on=['eid'], how='left'), [df, df_PRS])
    # df.to_csv(data_path + "ALL{}.csv".format(suffix), index=None)

    df_diag = df[((df['Btime'] == 1) & (df['PHQ'] >= 2)) | ((df['Btime'] == 0) & (df['PHQ'] == 0))]
    df_diag.to_csv(data_path + f"Diag{suffix}.csv", index=None)
    print("df_diag: {},".format(df_diag.groupby('MDD').size()))
    df_MDD, df_NC = df_diag[df_diag['MDD'] == MDD], df_diag[df_diag['MDD'] == NC]
    df_MDD.to_csv(data_path + f"Diag_MDD{suffix}.csv", index=None)
    df_NC.to_csv(data_path + f"Diag_NC{suffix}.csv", index=None)

    # df_after_diag = df[df['Btime'] == 2]
    # df_after_diag.to_csv(data_path + "{}_MDD_after_{}.csv".format("Diag", suffix), index=None)


if __name__ == '__main__':
    get_preprocess_data("-0.0")  # "Baseline"
    # get_preprocess_data("-1.0")  # "FollowUP"

    # # # 4. 社会人口学数据统计
    # exps = ['Diag-0.0']
    # social_stats(exps, 'MDD')

    # df_NC = pd.read_csv(data_path + "Diag_NC-0.0.csv", low_memory=False)
    # print(df_NC.shape[0])
    # df = df_NC[df_NC['PHQ'] == 1]
    # print(df.shape[0])



