#!/bin/bash
# conda activate GWAS
# ./run_Filter_Individual_and_QC.sh  {exp_name} {file_names}
# nohup ./run_Filter_Individual_and_QC.sh Hae_Biochem_Clustering 1_VS_2 1_VS_3 2_VS_3 > run_Filter_Individual_and_QC_12_VS_23.log &
    exp_name=$1

    # 获取第二个及以后的参数转换为数组
    file_names=("${@:2}")

    for file_name in "${file_names[@]}"
    do

      echo "Start filter Individuals & QC on each Chrome with ${file_name}"
      for chr in {1..22}; do
        plink --bfile "/data/UKBioBank/bed_Data/SNP_QC_raw_data/c${chr}_b0_v3"  --keep "/data/UKBioBank/bed_Data/${exp_name}/${file_name}/${file_name}.csv" --make-bed --out "/data/UKBioBank/bed_Data/${exp_name}/${file_name}/selIndivi_c${chr}"
        plink --bfile "/data/UKBioBank/bed_Data/${exp_name}/${file_name}/selIndivi_c${chr}" --exclude /data/UKBioBank/bed_Data/exclude_snps.txt --make-bed --out "/data/UKBioBank/bed_Data/${exp_name}/${file_name}/chr${chr}"
        rm -rf /data/UKBioBank/bed_Data/${exp_name}/${file_name}/selIndivi_c${chr}*
        # Delete SNPs with missingness >0.02.
        plink --bfile "/data/UKBioBank/bed_Data/${exp_name}/${file_name}/chr${chr}" --geno 0.02 --make-bed --out "/data/UKBioBank/bed_Data/${exp_name}/${file_name}/step1_chr${chr}"
        rm -rf /data/UKBioBank/bed_Data/${exp_name}/${file_name}/chr${chr}*
        # Remove SNPs with a low MAF frequency.
        plink --bfile "/data/UKBioBank/bed_Data/${exp_name}/${file_name}/step1_chr${chr}" --maf 0.05 --make-bed --out "/data/UKBioBank/bed_Data/${exp_name}/${file_name}/step2_chr${chr}"
        mv -f "/data/UKBioBank/bed_Data/${exp_name}/${file_name}/step2_chr${chr}.fam" "/data/UKBioBank/bed_Data/${exp_name}/${file_name}/step2_chr${chr}.fam.backup"
        rm -rf /data/UKBioBank/bed_Data/${exp_name}/${file_name}/step1_chr${chr}*
        rm -rf /data/UKBioBank/bed_Data/${exp_name}/${file_name}/*.log
      done
      echo "Done: filter Individuals & QC on each Chrome with ${file_name}"
    done
    exit 0
