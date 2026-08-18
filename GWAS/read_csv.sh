while IFS=, read -r TWAS_name pos_name pvalue; do
    echo "列1=${TWAS_name}, 列2=${pos_name}, 列3=${pvalue}"
done < /data/UKBioBank/TWAS_PWAS_Data/summary_twas_dataset_lst_ALL.csv