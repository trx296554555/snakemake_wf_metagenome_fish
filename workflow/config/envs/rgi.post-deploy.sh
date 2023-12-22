root_path=$(dirname `dirname \`dirname "$CONDA_PREFIX"_\``)
env_name=$(basename $CONDA_PREFIX)
config_path="${root_path}"/workflow/config
db_root=$(awk -F ': ' '/db_root:/{print $2}' "${config_path}"/config.yaml)
card_db=$(awk -F ': ' '/    rgi:/{print $2}' "${config_path}"/config.yaml)
# 使用 rgi 对预测的蛋白质进行抗药性功能注释

{
  echo "Make additional adjustments for the post-deployment of the Conda environment ${env_name} (rgi)"
  echo "Annotation the assembled Contigs using rgi."
  echo "---------------------------"
} >> "${root_path}"/logs/env.log

# Check whether the card database already exists
if [ -d "${db_root}"/"${card_db}" ];then
  echo "The card database already exists, location：${db_root}/${card_db}" >> "${root_path}"/logs/env.log
  rgi load \
    --card_json "${db_root}"/"${card_db}"/card.json \
    --debug \
    --card_annotation "${db_root}"/"${card_db}"/card_database_v3.2.8.fasta \
    --card_annotation_all_models "${db_root}"/"${card_db}"/card_database_v3.2.8_all.fasta \
    --wildcard_annotation "${db_root}"/"${card_db}"/wildcard_database_v4.0.0.fasta \
    --wildcard_annotation_all_models "${db_root}"/"${card_db}"/wildcard_database_v4.0.0_all.fasta \
    --wildcard_index "${db_root}"/"${card_db}"/index-for-model-sequences.txt \
    --wildcard_version 4.0.0 \
    --amr_kmers "${db_root}"/"${card_db}"/all_amr_61mers.txt \
    --kmer_database "${db_root}"/"${card_db}"/61_kmer_db.json \
    --kmer_size 61
else
  echo "The card database doesn't exist. Use the following command to build it yourself:" >> "${root_path}"/logs/env.log
  # 构建数据库
cat >> "${root_path}"/logs/env.log <<EOF
  ------------------------------------------------------------
  mkdir -p "${db_root}"/"${card_db}"

  # Remove any previous loads:
  rgi clean --local

  # Download CARD and WildCARD data:
  wget https://card.mcmaster.ca/latest/data
  tar -xvf data ./card.json
  wget -O wildcard_data.tar.bz2 https://card.mcmaster.ca/latest/variants
  mkdir -p wildcard
  tar -xjf wildcard_data.tar.bz2 -C wildcard
  gunzip wildcard/*.gz

  # Create annotation files (note that the parameter version_number depends upon the versions of WildCARD data downloaded, please adjust accordingly):
  rgi card_annotation -i /path/to/card.json > card_annotation.log 2>&1
  rgi wildcard_annotation -i wildcard --card_json /path/to/card.json -v version_number > wildcard_annotation.log 2>&1

  # Load all data into RGI (note that the FASTA filenames plus the parameter version_number depend on the versions of CARD and WildCARD data downloaded, please adjust accordingly):
  rgi load \
    --card_json /path/to/card.json \
    --debug --local \
    --card_annotation card_database_v3.2.4.fasta \
    --card_annotation_all_models card_database_v3.2.4_all.fasta \
    --wildcard_annotation wildcard_database_v4.0.0.fasta \
    --wildcard_annotation_all_models wildcard_database_v4.0.0_all.fasta \
    --wildcard_index /path/to/wildcard/index-for-model-sequences.txt \
    --wildcard_version 4.0.0 \
    --amr_kmers /path/to/wildcard/all_amr_61mers.txt \
    --kmer_database /path/to/wildcard/61_kmer_db.json \
    --kmer_size 61
    ------------------------------------------------------------
EOF
exit 1
fi
