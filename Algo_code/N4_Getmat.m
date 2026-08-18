root_path = "/data4/myuan/"
path = root_path +"Hae_Bio/"
data_path = path + "dataset/"
data = readtable(data_path + 'Diag_MDD-0.0.csv', 'VariableNamingRule', 'preserve');

features_Haemocytes = ["WBC", "RBC", "HGB", "HCT", "MCV", "MCH", "MCHC", "RDW", "PLT", "PCT", "MPV", "PDW", "LYcnt", "MOcnt", "NEcnt", "EOcnt", "BAnt", "NRBCcnt", "Lyp", "Mop", "Nep", "Eop", "Bap", "NRBCp", "RETp", "RETcnt", "MRV", "MSCV", "IRF", "HLRp", "HLRc"];   % 31个特征
features_Biochem = ["ALB", "ALP", "ALT", "APOA.bio", "APOB.bio", "AST", "DBIL", "UREA", "CALC", "CHOL", "CREA", "CRP", "CYSC", "GGT", "GLU", "HbA1c", "HDL", "IGF1", "LDL", "LPa", "OEST", "PHOS", "RF", "SHBG.bio", "TBIL", "TEST", "TP", "TRIG", "UA", "VITD"];   % 30个特征
data_Haemocytes = data{:, features_Haemocytes};
data_Biochem = data{:, features_Biochem};

% 提取标签数据（MDD）
Y = data.MDD;
eid = data.eid;

% 将特征数据放入 X 中
X = cell(2, 1); % 2x1的单元数组，表示2个视图
X{1} = data_Haemocytes; % 第一视图的数据
X{2} = data_Biochem; % 第二视图的数据

% 保存为 .mat 文件
save(data_path+'Hae_Bio_0.0.mat', 'X', 'Y', 'eid');
fprintf('已生成 Hae.mat 文件，其中包含 X, Y, eid 三个变量。\nX 表示特征数据，Y 表示标签数据，eid 表示样本编号。\n');
