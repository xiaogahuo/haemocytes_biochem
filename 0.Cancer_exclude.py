# 统计患癌症疾病的受试者， 用于后续排除不满足条件的受试者
from SUtils import *

cancer_cols = ['20001-0.0', '20001-0.1', '20001-0.2', '20001-0.3', '20001-0.4', '20001-0.5', '20001-1.0', '20001-1.1', '20001-1.2', '20001-1.3', '20001-1.4', '20001-1.5', '20001-2.0', '20001-2.1', '20001-2.2', '20001-2.3', '20001-2.4', '20001-2.5', '20001-3.0', '20001-3.1', '20001-3.2', '20001-3.3', '20001-3.4', '20001-3.5']


def single_row(x):
    rst = [elem if (not np.isnan(elem) and elem != 99999) else None for elem in x[cancer_cols]] # 去除掉 99999 == No Cancer
    rst = list(filter(None, rst))  # 去除列表中的None元素,  去除掉 99999
    return rst, 1 if len(rst) > 0 else 0


def c_gotten():
    df = pd.read_csv(root_path + "20001.csv")
    # print(df.columns.tolist())
    row = df.apply(lambda x: single_row(x), axis=1)
    rst, cancer_gotten = zip(*row)
    df['cancer_list'] = pd.Series(rst)
    df['cancer_gotten'] = pd.Series(cancer_gotten)
    if not os.path.exists(data_path +"Cancer_gotten/"):
        os.makedirs(data_path +"Cancer_gotten/")
    df.to_csv(data_path + "Cancer_gotten/" +"Cancer_gotten.csv", columns=['eid', 'cancer_list', 'cancer_gotten'], index=None)


if __name__ == '__main__':
    c_gotten()