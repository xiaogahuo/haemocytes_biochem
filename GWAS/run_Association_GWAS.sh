#!/bin/bash
# conda activate GWAS
# run_Association_GWAS.sh  GWAS_sample_cluster1
# nohup ./run_Association_GWAS.sh Hae_Biochem_Clustering 0_VS_123 0_VS_1 0_VS_2 0_VS_3 1_VS_2 1_VS_3 2_VS_3 > run_Association_GWAS.log &
    exp_name=$1
  # 获取第二个及以后的参数转换为数组
    file_names=("${@:2}")

    # for 循环
    for file_name in "${file_names[@]}"
    do
      echo "Start Association on each Chrome with ${file_name}"

      for chr in {1..22}; do
        # logistic
        # We will be using 10 principal components as covariates in this logistic analysis. We use the MDS components calculated from the previous tutorial: covar_mds.txt.
        plink --bfile /data/UKBioBank/bed_Data/${exp_name}/${file_name}/step2_chr${chr} --freq --out /data/UKBioBank/bed_Data/${exp_name}/${file_name}/MAF_chr${chr}
        plink --bfile /data/UKBioBank/bed_Data/${exp_name}/${file_name}/step2_chr${chr} --covar /data/UKBioBank/bed_Data/${exp_name}/${file_name}/covr_chr${chr}.txt --logistic --hide-covar --out /data/UKBioBank/bed_Data/${exp_name}/${file_name}/logistic_results_chr${chr}

        # Remove NA values, those might give problems generating plots in later steps.
        awk '!/'NA'/' /data/UKBioBank/bed_Data/${exp_name}/${file_name}/logistic_results_chr${chr}.assoc.logistic > /data/UKBioBank/bed_Data/${exp_name}/${file_name}/logistic_results_chr${chr}.assoc.logistic2
      done
      # 首先将第一个染色体文件（包含标题）复制到最终文件
      cp /data/UKBioBank/bed_Data/${exp_name}/${file_name}/MAF_chr1.frq  /data/UKBioBank/bed_Data/${exp_name}/${file_name}/MAF.frq
      cp /data/UKBioBank/bed_Data/${exp_name}/${file_name}/logistic_results_chr1.assoc.logistic2  /data/UKBioBank/bed_Data/${exp_name}/${file_name}/logistic_results.assoc.logistic2
      # 然后将其余22个文件的内容追加（跳过标题行）
      for chr in {2..22}; do
        tail -n +2 /data/UKBioBank/bed_Data/${exp_name}/${file_name}/MAF_chr${chr}.frq >> /data/UKBioBank/bed_Data/${exp_name}/${file_name}/MAF.frq
        tail -n +2 /data/UKBioBank/bed_Data/${exp_name}/${file_name}/logistic_results_chr${chr}.assoc.logistic2 >> /data/UKBioBank/bed_Data/${exp_name}/${file_name}/logistic_results.assoc.logistic2
      done
      rm -rf /data/UKBioBank/bed_Data/${exp_name}/${file_name}/*.log
      rm -rf /data/UKBioBank/bed_Data/${exp_name}/${file_name}/*.logistic
      echo "Done: Association on each Chrome with ${file_name}"
    done
    exit 0

