#!/bin/bash
# http://gusevlab.org/projects/fusion/
# http://nilanjanchatterjeelab.org/pwas

echo "begin"
# nohup ./run_PWAS.sh Hae_Biochem_Clusterings_NUM2 0_VS_1 > run_PWAS_0VS1.log &
exp_name=$1   #  Hae_Biochem_Clusterings_NUM2
cmp_name=$2   #  0_VS_1
GWAS_NUM=$3      #  0_VS_1: 155958   0_VS_2: 156207
start_line=$4
end_line=$5
# ROSMAP.n376.fusion.WEIGHTS_v2     # Banner.n152.fusion.WEIGHTS

if [ -z "$GWAS_NUM" ]; then
    echo "ERROR: GWAS_NUM (3rd argument) is required. Usage: $0 <exp_name> <cmp_name> <GWAS_NUM> <start_line> <end_line>"
    exit 1

fi
echo "${exp_name}, ${cmp_name}, ${GWAS_NUM}, ${start_line}, ${end_line}"

while IFS=, read -r TWAS_name pos_file pvalue enter_nan; do
    echo "列1=${TWAS_name}, 列2=${pos_file}, 列3=${pvalue}, 列4=${enter_nan} "

    out_dir="/data/UKBioBank/bed_Data/${exp_name}/${cmp_name}/pwas/${TWAS_name}/"
    if [ ! -d "${out_dir}" ]; then
      mkdir -p "${out_dir}"
    fi

    for chr in {1..22}; do

      Rscript ./fusion_twas/FUSION.assoc_test.R \
      --sumstats "/data/UKBioBank/bed_Data/${exp_name}/${cmp_name}/${cmp_name}_gwas_summary.txt.gz" \
      --weights  "/data/UKBioBank/TWAS_PWAS_Data/Weights/${TWAS_name}/${pos_file}" \
      --weights_dir "/data/UKBioBank/TWAS_PWAS_Data/Weights/${TWAS_name}"  \
      --ref_ld_chr  /data/UKBioBank/TWAS_PWAS_Data/Weights/LDREF/1000G.EUR. \
      --chr "${chr}" \
      --coloc_P 0.05 \
      --GWASN "${GWAS_NUM}" \
      --out  "${out_dir}/chr_${chr}.dat"

#      # shellcheck disable=SC2002
#      cat "${out_dir}/chr_${chr}.dat" | awk -v pv="${pvalue}" 'NR == 1 || $NF < pv' > "${out_dir}/chr_${chr}.top"
#
#      Rscript ./fusion_twas/FUSION.post_process.R \
#      --sumstats "/data/UKBioBank/bed_Data/${exp_name}/${cmp_name}/${cmp_name}_gwas_summary.txt.gz" \
#      --input "${out_dir}/chr_${chr}.top" \
#      --out "${out_dir}/chr_${chr}.top.analysis" \
#      --ref_ld_chr /data/UKBioBank/TWAS_PWAS_Data/Weights/LDREF/1000G.EUR. \
#      --chr "${chr}" \
#      --plot --locus_win 100000

      echo "Done: ${exp_name}, ${cmp_name}, ${TWAS_name}, chr${chr}"
    done

    all_chr_file="${out_dir}/all_chr.dat"
    head -n 1 "${out_dir}/chr_1.dat" > "${all_chr_file}"
    for chr in {1..22}; do
        tail -n +2 "${out_dir}/chr_${chr}.dat" >> "${all_chr_file}"
    done
    echo "Merged all chromosome dat files into: ${all_chr_file}"

    Rscript /data/UKBioBank/bed_Data/fusion_twas/calculate_FDR_noNA.R \
        --input "${all_chr_file}" \
        --out "${out_dir}/all_chr_with_FDR_noNA.dat" \
        --fdr "${pvalue}" \
        --sigout "${out_dir}/FDR_${pvalue}_genes.txt"

    Rscript /data/UKBioBank/bed_Data/fusion_twas/plot_manhattan_FDR.R \
        --input "${out_dir}/all_chr_with_FDR_noNA.dat" \
        --out "${out_dir}/ManhattanPlot_FDR_${pvalue}.pdf" \
        --fdr "${pvalue}"
    echo "Finished FDR and Manhattan plot for ${TWAS_name}"

done  < <(sed -n "${start_line},${end_line}p" /data/UKBioBank/TWAS_PWAS_Data/summary_twas_dataset_lst.txt)
echo "Done PWAS"