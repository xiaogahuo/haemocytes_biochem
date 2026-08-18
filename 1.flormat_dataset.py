import pandas as pd
from SUtils import *
from datetime import datetime
from functools import reduce


# ICD10 Disease Columns - -41270   from '41270-0.0' to '41270-0.242'
diagnosis_ICD = ['41270-0.'+str(i) for i in range(243)]
# print(diagnosis_ICD)
# ICD10-Affective-Disorder  Disease Code ---- 精分双相，成瘾，脑血管相关疾病
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=41270
code_ICD10_F = ['F0', 'F112', 'F122', 'F132', 'F142', 'F162', 'F2', 'F30', 'F31', 'F34', 'F38', 'F39', 'F7']
code_ICD10_G = ['G45', 'G46', 'G30', 'G31', 'G35', 'G20']
code_ICD10_I = ['I60', 'I61', 'I62', 'I63', 'I64', 'I67']
ICD_exclude_code = code_ICD10_F + code_ICD10_G + code_ICD10_I
chronic_code_ICD10 = ['I10', 'I11', 'I12', 'I13', 'I14', 'I15', 'I20', 'I21', 'I22', 'I23', 'I24', 'I25', 'E10', 'E11']


def ICD_exclude(row):
    diag_lst = row[diagnosis_ICD].tolist()
    diag_lst = [x for x in diag_lst if not pd.isna(x)]
    exclude_num = startCnt(ICD_exclude_code, diag_lst)
    return exclude_num


def cnt_chronic_diseases(x):
    # print(x[diagnosis_ICD].tolist())
    result = [str(elem) if not pd.isna(elem) else None for elem in x[diagnosis_ICD].tolist()]
    result = list(filter(None, result))
    chronic_num = startCnt(chronic_code_ICD10, result)
    # print("ICD_disease: {}, chronic_num: {}".format(result, chronic_num))
    # return result, chronic_num
    return chronic_num


codes_mdd = [20, 21, 30, 31, 40, 41, 50, 51]
def mdd(v):
    if pd.isna(v):
        return NC  # NC
    elif int(v) in codes_mdd:
        return MDD  # MDD
    else:
        return UNNORMAL  # Un-normal data


def base_time(row):
    mdd_tab, diag_time, base_time = row['MDD'], row['130894-0.0'], row['53-0.0']
    if mdd_tab == NC:
        return 0   # NC
    elif pd.isna(diag_time):  # mdd == 1, diag_time == NA, exclude
        return -1  # drop
    elif datetime.strptime(diag_time, "%Y-%m-%d") <= datetime.strptime(base_time, "%Y-%m-%d"):  # 基线之前诊断为抑郁症
        return 1
    else:  # 基线之后诊断为抑郁症
        return 2


# 判断字符串result是否以某个子串tag开头, 统计总数
def startCnt(listTag, listResult):
    num = 0
    for result in listResult:
        for tag in listTag:
            if result.startswith(tag):
                num += 1
    return num


# 1-3 Drug_Yes    4, 5, -7 Drug_No     others   Drug_Unknown
def drug_intake(x):
    result = [x['6177-0.0'], x['6177-0.1'], x['6177-0.2'], x['6153-0.0'], x['6153-0.1'], x['6153-0.2'], x['6153-0.3']]
    result = [x for x in result if x]
    if set([1, 2, 3]).intersection(set(result)):
        return 'Drug_Yes'  # 'Drug_Yes'
    elif set([4, 5, -7]).intersection(set(result)):
        return 'Drug_No'  # 'Drug_No'
    else:
        return 'Drug_Unknown'  # 'Drug_Unknown'


def format2Dataset(suffix):
    #血常规
    df_haemocytes = pd.read_csv(path + "Cate_100081_filtered{}.csv".format(suffix), low_memory=False)
    df_codings = pd.read_excel(path + '多组学字典.xlsx', sheet_name="haemocytes")
    id_cols, names_cols = df_codings['FieldID'].tolist(), df_codings['FieldName'].tolist()
    id_cols = [str(code)+suffix for code in id_cols]
    df_haemocytes.rename(columns=dict(zip(id_cols, names_cols)), inplace=True)
    print("Haemocytes_cols: {}".format(names_cols))
    for col in names_cols:
        df_haemocytes[col] = (df_haemocytes[col] - df_haemocytes[col].mean()) / df_haemocytes[col].std(ddof=0)  # z变换
        df_haemocytes[col] = df_haemocytes[col].fillna(0)  # fill NA with 0
    df_haemocytes.to_csv(data_path+"Haemocytes{}.csv".format(suffix), index=None)
    #生化
    df_biochem = pd.read_csv(path + "Cate_17518_filtered{}.csv".format(suffix), low_memory=False)
    df_codings = pd.read_excel(path + '多组学字典.xlsx', sheet_name="biochem")
    id_cols, names_cols = df_codings['FieldID'].tolist(), df_codings['FieldName'].tolist()
    id_cols = [str(code) + suffix for code in id_cols]
    df_biochem.rename(columns=dict(zip(id_cols, names_cols)), inplace=True)
    print("Biochem_cols: {}".format(names_cols))
    for col in names_cols:
        df_biochem[col] = (df_biochem[col] - df_biochem[col].mean()) / df_biochem[col].std(ddof=0)  # z变换
        df_biochem[col] = df_biochem[col].fillna(0)  # fill NA with 0
    df_biochem.to_csv(data_path+"Biochem{}.csv".format(suffix), index=None)
    print("haemocytes: {}, biochem: {}".format(df_haemocytes.shape[0],df_biochem.shape[0]))


def format_MDD(suffix):
    df_cancer = pd.read_csv(data_path + "Cancer_gotten/Cancer_gotten.csv",
                            usecols=['eid', 'cancer_list', 'cancer_gotten'])  # 排掉 cancer_gotten == 1的 人
    print("ALL cancer: {}".format(df_cancer.groupby('cancer_gotten').size()))

    df_haemocytes = pd.read_csv(data_path+"Haemocytes{}.csv".format(suffix), usecols=['eid'])
    df_biochem = pd.read_csv(data_path+"Biochem{}.csv".format(suffix), usecols=['eid'])

    df_MDD = pd.read_csv(path + "Cate_2405_not_filtered-0.0.csv", usecols=['eid', '130894-0.0', '130895-0.0'],
                     dtype={'130894-0.0': str})
    print("df_MDD before merge: {}".format(df_MDD.shape[0]))
    # df_MDD = pd.merge(df_MDD, df_prote, how='inner', on='eid', sort=True, suffixes=('_x', '_y'), copy=True, indicator=False)
    df_MDD = reduce(lambda left, right: pd.merge(left, right, on=['eid'], how='inner'), [df_MDD, df_haemocytes, df_biochem, df_cancer])
    print("df_MDD after merge with Hae_Bio: {}".format(df_MDD.shape[0]))
    values_to_exclude = ['1900-01-01', '1901-01-01', '1902-02-02', '1903-03-03', '1909-09-09', '2037-07-07']
    df_MDD = df_MDD[~df_MDD['130894-0.0'].isin(values_to_exclude)]
    print("df_MDD diag_time excluded with ['1900-01-01', '1901-01-01', '1902-02-02', '1903-03-03', '1909-09-09', '2037-07-07'] : {}".format(df_MDD.shape[0]))
    df_MDD['MDD'] = df_MDD['130895-0.0'].map(mdd)
    print('df_MDD: {}'.format(df_MDD.groupby('MDD').size()))

    df_chunk = pd.read_csv(root_path+"1.social_demography.csv", chunksize=50000, low_memory=False)
    lst = []
    for i, chunk in enumerate(df_chunk):
        chunk = pd.merge(chunk, df_MDD, how='inner', on='eid', sort=True, suffixes=('_x', '_y'), copy=True, indicator=False)
        lst.append(chunk)
    df = pd.concat(lst, ignore_index=True)
    print('df : {}'.format(df.groupby('MDD').size()))
    df['ICD_exclude'] = df.apply(lambda row: ICD_exclude(row), axis=1)
    df = df[(df['ICD_exclude'] == 0) & (df['cancer_gotten'] == 0)]
    print("df after exclude ICD & cancer:  {}".format(df.groupby('MDD').size()))
    df["Btime"] = df.apply(lambda row: base_time(row), axis=1)
    print("df Btime: {}".format(df.groupby("Btime").size()))
    chronic_num = df.apply(lambda x: cnt_chronic_diseases(x), axis=1)
    # chronic_diseases, chronic_num = zip(*rst)
    # df['chronic_diseases'] = pd.Series(chronic_diseases)
    df['chronic_num'] = pd.Series(chronic_num)

    for col in [c+suffix for c in ['2050', '2060', '2070', '2080']]:
        df[col] = df[col] - 1  # 1, 2, 3, 4 -> 0, 1, 2, 3  item -> Score
        df[col] = df[col].fillna(0)  # fill NAN with 0
        df.loc[df[col] < 0, col] = 0  # fill neg with 0
    df['PHQ'] = df.apply(lambda x: x['2050'+suffix] + x['2060'+suffix], axis=1)
    df['PHQ_Total'] = df.apply(lambda x: x['2050'+suffix] + x['2060'+suffix] + x['2070'+suffix] + x['2080'+suffix], axis=1)
    print("df: {}".format(df.groupby('MDD').size()))
    df.to_csv(data_path + "Social_MDD{}.csv".format(suffix), index=None)


if __name__ == '__main__':
    format2Dataset("-0.0") # "Baseline"
    # format2Dataset("-1.0") # "FollowUP"
    # format2Dataset("-2.0") # "FollowUP2"

    format_MDD("-0.0")
    # format_MDD("-1.0")
    # format_MDD("-2.0")

    # # 检查每个元素是否为空值
    # df = pd.read_csv(path + "Social_MDD.csv", low_memory=False)
    # null_counts = df[['chronic_num']].isnull().sum()
    # # 打印每列的空值数量
    # print(null_counts)
