##### 1. 过滤 Sample
    1.1 执行GWAS_data_preproc.py 的 Sample_QC()  生成用于过滤 sample 的csv文件
    conda activate GWAS
    cd /data/UKBioBank/bed_Data
    1.2 执行 nohup run_Filter_Individual_and_QC.sh {exp_name} {file_name} > run_Filter_Individual_and_QC_{file_name}.log & 生成GWAS QC后的数据集到 /data/UKBioBank/bed_Data/{exp_name}/{file_name} 子目录
        tail -f run_Filter_Individual_and_QC_{file_name}.log
        ps: 脚本有execute权限仍然报error 则执行file run_Filter_Individual_and_QC.sh 看是否有'with CRLF line terminators'日志信息， 有则执行 dos2unix run_Filter_Individual_and_QC.sh命令对sh脚本格式进行转换 

##### 2. 执行 GWAS 关联分析
    2.1 执行GWAS_data_preproc.py 的 split_cluster() 生成用于执行 GWAS 关联分析的 fam文件及covr协变量文件 
    conda activate GWAS
    cd /data/UKBioBank/bed_Data/
    2.2 执行 nohup ./run_Association_GWAS.sh {exp_name} {file_name} > GWAS_Association.log & 执行 GWAS 关联分析
        tail -f GWAS_Association.log
    2.3 执行GWAS_data_preproc.py 的 cnt_sig_pvalue() 将 A2 等位基因信息合并到 GWAS summary 的 logistic_results.assoc.logistic3 文件中, 用于后续的LD分析

##### 3. 基于第2步的GWAS summary结果计算LD SCore、基因多效性、遗传度、遗传相关性
[//]: # (https://www.cnblogs.com/chenwenyan/p/11321272.html)
    虚拟环境名称为ldsc (python 2.7 环境): 
    conda activate ldsc 
    cd /data/UKBioBank/bed_Data
    3.1 执行 nohup ./run_LDSC.sh ${exp_name} ${file_name} > /data/UKBioBank/bed_Data/${exp_name}/${file_name}/ldsc.log & 用于生成LD 分析结果
    
    Hae_Biochem_Clusterings_NUM2
    nohup ./run_LDSC.sh Hae_Biochem_Clusterings_NUM2 0_VS_1 > /data/UKBioBank/bed_Data/Hae_Biochem_Clusterings_NUM2/0_VS_1/ldsc.log &
    nohup ./run_LDSC.sh Hae_Biochem_Clusterings_NUM2 0_VS_2 > /data/UKBioBank/bed_Data/Hae_Biochem_Clusterings_NUM2/0_VS_2/ldsc.log &
    nohup ./run_LDSC.sh Hae_Biochem_Clusterings_NUM2 0_VS_12 > /data/UKBioBank/bed_Data/Hae_Biochem_Clusterings_NUM2/0_VS_12/ldsc.log &
    nohup ./run_LDSC.sh Hae_Biochem_Clusterings_NUM2 1_VS_2 > /data/UKBioBank/bed_Data/Hae_Biochem_Clusterings_NUM2/1_VS_2/ldsc.log &
    
    3.2 执行GWAS_data_preproc.py 的 format_LDGC_rst() 将LDGC 结果写入LDGC.csv 文件中

##### 4. 绘制Locuszoom
      4.1 执行Plot_Locuszoom.py 绘制Locuszoom 图

##### 5. 执行 TWAS & PWAS
      下载Weights只需要执行一次
      nohup ./gtex_download.sh /data/UKBioBank/TWAS_PWAS_Data/Weights > gtex.out 2>&1 &
      
      修改summary_twas_dataset_lst.txt 后替换服务器/data/UKBioBank/TWAS_PWAS_Data/summary_twas_dataset_lst.txt对应文件
      cd /data/UKBioBank/bed_Data
      执行nohup ./run_PWAS.sh Hae_Biochem_Clusterings_NUM2 0_VS_2  1 5 > run_PWAS_0VS2.log &
      要提高并行度，参考 run_PWAS_ALL.sh

##### 6. 执行 SMR
      工具链接网址： https://yanglab.westlake.edu.cn/software/smr
      *QTL数据下载链接： https://yanglab.westlake.edu.cn/software/smr/#DataResource
      gwas.txt.gz解压缩命令: gunzip 0_VS_2_gwas_summary.txt.gz      

      
      --beqtl-summary /data/UKBioBank/GTEx_v10_eQTL/Brain_Frontal_Cortex_BA9/Brain_Frontal_Cortex_BA9.v10.eQTL.cis_qtl_pairs.17 \
      --beqtl-summary /data/UKBioBank/GTEx_v10_eQTL/Brain_Caudate_basal_ganglia/Brain_Caudate_basal_ganglia.v10.eQTL.cis_qtl_pairs.17 \
      --beqtl-summary /data/UKBioBank/GTEx_v10_eQTL/BrainMeta_cis_eqtl_summary/BrainMeta_cis_eQTL_chr17 \
      
      cd /data/UKBioBank/bed_Data
      ./smr --bfile /data/UKBioBank/TWAS_PWAS_Data/Weights/LDREF/1000G.EUR.17 \
            --gwas-summary /data/UKBioBank/bed_Data/Hae_Biochem_Clusterings_NUM2/0_VS_1/0_VS_1_gwas_summary.txt \
            --beqtl-summary /data/UKBioBank/GTEx_v10_eQTL/BrainMeta_cis_eqtl_summary/BrainMeta_cis_eQTL_chr17 \
            --plot --probe ENSG00000120088.15 \
            --probe-wind 100 \
            --out 0vs1_BrainMeta_Prob100 \
            --gene-list /data/UKBioBank/bed_Data/plot_smr/glist-hg38
      注意修改 --gwas-summary --beqtl-summary  --probe-wind   --out (要体现前3个参数的唯一性)  --probe ENSG00000120088.15  
      命令执行后结果生成到 /data/UKBioBank/bed_Data/plot/ 目录

      mstsc 远程登录到服务器打开 cd /data/UKBioBank/bed_Data/plot_smr 目录下的 GWAS_eQTL.R 文件
      修改 SMRData 的load 目录， 绘制SMRLocusPlot（参数不用修改） 和  SMREffectPlot（参数不用修改）
      

#### 注意： 注意执行 df -h 后查看 /data/ 目录是否已满
执行完实验将 locuszoom文件夹、pwas文件夹、 
LDGC文件夹、{file_name}manhattan_QQ_plots.png、logistic_results.1e-5.txt、 
{file_name}_gwas_summary.txt.gz、{file_name}_LDGC.csv、ldsc.log结果拷贝到本地

