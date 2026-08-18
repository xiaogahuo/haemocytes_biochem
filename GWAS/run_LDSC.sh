#!/bin/bash
# https://www.cnblogs.com/chenwenyan/p/11321272.html
# https://cloud.tencent.com/developer/article/1670338
# https://zenodo.org/records/10515792

exp_name=$1
file_name=$2

line_count=$(wc -l < "/data/UKBioBank/bed_Data/${exp_name}/${file_name}/step2_chr1.fam")

echo ${line_count}
echo "${exp_name}/${file_name}"

/home/user/ldsc/munge_sumstats.py \
    --sumstats "/data/UKBioBank/bed_Data/${exp_name}/${file_name}/logistic_results.assoc.logistic3" \
    --N ${line_count} \
    --merge-alleles "/data/UKBioBank/LDScore_EUR/w_hm3.snplist" \
    --out "/data/UKBioBank/bed_Data/${exp_name}/${file_name}/scz"

/home/user/ldsc/ldsc.py \
    --h2 "/data/UKBioBank/bed_Data/${exp_name}/${file_name}/scz.sumstats.gz" \
    --ref-ld-chr "/data/UKBioBank/LDScore_EUR/LDscore/" \
    --w-ld-chr "/data/UKBioBank/LDScore_EUR/LDscore/" \
    --out "/data/UKBioBank/bed_Data/${exp_name}/${file_name}/scz_h2"


sumstats_file_name="/data/UKBioBank/LDScore_EUR/sumstats_ALL/Hae_Bio_MDD_Subtype.txt"
# 构建路径
out_path="/data/UKBioBank/bed_Data/${exp_name}/${file_name}/LDGC/"
# 判断路径是否存在，如果不存在则创建
if [ ! -d "${out_path}" ]; then
    echo "路径 ${out_path} 不存在，正在创建..."
    mkdir -p "${out_path}"
else
    echo "路径 ${out_path} 已存在。"
fi

# 按行读取文件
while IFS= read -r line; do
    echo "${line}"
    /home/user/ldsc/ldsc.py \
        --rg "/data/UKBioBank/LDScore_EUR/sumstats_ALL/${line}.sumstats.gz","/data/UKBioBank/bed_Data/${exp_name}/${file_name}/scz.sumstats.gz" \
        --ref-ld-chr "/data/UKBioBank/LDScore_EUR/LDscore/" \
        --w-ld-chr "/data/UKBioBank/LDScore_EUR/LDscore/" \
        --out "${out_path}${line}"
done < "${sumstats_file_name}"