from SUtils import *


# 抽取特定后缀(suffix)的列
def getSubColumns(df, field_sheets, suffix):
    data_cols = df.columns.tolist()
    target_cols = []

    for sheet in field_sheets:
        cols = pd.read_excel(open(path + "fields_ID.xlsx", 'rb'), sheet_name=sheet)["Field ID"].tolist()
        cols = [str(col)+suffix for col in cols]
        target_cols = target_cols + cols
    rst_cols = ['eid']+list(set(data_cols).intersection(set(target_cols)))
    return rst_cols


def filter_dataset(filename, field_sheets, suffix):
    print(filename)
    df = pd.read_csv(root_path + filename + ".csv", low_memory=False)
    sub_cols = getSubColumns(df, field_sheets, suffix)
    df = df[sub_cols]
    df = df.dropna(axis=0, thresh=2)  # Keep those rows with at least 0.2 non-NA values.
    df.to_csv(path + filename + "_filtered{}.csv".format(suffix), index=None)


def not_filter_dataset(filename, field_sheets, suffix):
    print(filename)
    df = pd.read_csv(root_path + filename + ".csv", low_memory=False)
    sub_cols = getSubColumns(df, field_sheets, suffix)
    df = df[sub_cols]
    df.to_csv(path + filename + "_not_filtered{}.csv".format(suffix), index=None)


if __name__ == '__main__':
    not_filter_dataset("Cate_2405", ["C2405"], "-0.0")  # only have 0.0 data

    filter_dataset("Cate_17518", ["C17518"], "-0.0")
    filter_dataset("Cate_100081", ["C100081"], "-0.0")
    filter_dataset("Cate_17518", ["C17518"], "-1.0")
    filter_dataset("Cate_100081", ["C100081"], "-1.0")

