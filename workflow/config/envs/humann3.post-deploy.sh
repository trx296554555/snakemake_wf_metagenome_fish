root_path=$(dirname `dirname \`dirname "$CONDA_PREFIX"_\``)
config_path="${root_path}"/workflow/config
env_name=$(basename $CONDA_PREFIX)

# Specify the version of the metaphlan3 database to match the version of humann3 (cannot use the latest version)
metaphlan3_db_version="v30_CHOCOPhlAn_201901"
echo "metaphlan3_db_version: ${metaphlan3_db_version}" >> "${root_path}"/logs/env.log
echo $metaphlan3_db_version > "${CONDA_PREFIX}"/envs/humann3/lib/python3.1/site-packages/metaphlan/metaphlan_databases/mpa_latest
# TODO
# 检测是否已经下载了metaphlan3数据库，如果有，就不再下载，并软链到humann3的数据库目录下
# 等待humann4发布后，再进行修改

# Read the configuration
db_root=$(awk -F ': ' '/db_root:/{print $2}' "${config_path}"/config.yaml)
humann_chocophlan=$(awk -F ': ' '/humann_chocophlan:/{print $2}' "${config_path}"/config.yaml)
humann_uniref=$(awk -F ': ' '/humann_uniref:/{print $2}' "${config_path}"/config.yaml)
humann_utility_mapping=$(awk -F ': ' '/humann_utility_mapping:/{print $2}' "${config_path}"/config.yaml)

{
  echo "Make additional adjustments for the post-deployment of the Conda environment ${env_name} (humann3)"
  echo db_root: "$db_root"
  echo humann_chocophlan: "$humann_chocophlan"
  echo humann_uniref: "$humann_uniref"
  echo humann_utility_mapping: "$humann_utility_mapping"
} >> "${root_path}"/logs/env.log

# Check whether the kraken2 database already exists
if [ -d "${db_root}"/"${humann_chocophlan}" ] && [ -d  "${db_root}"/"${humann_uniref}" ] && [ -d  "${db_root}"/"${humann_utility_mapping}" ];then
  echo "Humann3 database already exists, location：${db_root}/$(basename $humann_uniref)" >> "${root_path}"/logs/env.log
  humann_config --update database_folders nucleotide "${db_root}"/"${humann_chocophlan}"
  humann_config --update database_folders protein "${db_root}"/"${humann_uniref}"
  humann_config --update database_folders utility_mapping "${db_root}"/"${humann_utility_mapping}"
else
  echo "The Humann3 database doesn't exist. Use the following command to build it yourself:" >> "${root_path}"/logs/env.log
  # download the humann3 database
  echo "conda activate $CONDA_PREFIX && cd ${db_root} && humann_databases --download uniref uniref90_diamond ${humann_uniref} --update-config yes" >> "${root_path}"/logs/env.log
  echo "conda activate $CONDA_PREFIX && cd ${db_root} && humann_databases --download utility_mapping full ${humann_uniref} --update-config yes" >> "${root_path}"/logs/env.log
  echo "conda activate $CONDA_PREFIX && cd ${db_root} && humann_databases --download chocophlan full ${humann_chocophlan} --update-config yes" >> "${root_path}"/logs/env.log
fi

echo "When fisrt running humann3, metaphlan3 database will be downloaded automatically." >> "${root_path}"/logs/env.log
echo "Which will take a long time, more then" >> "${root_path}"/logs/env.log
echo "---------------------------" >> "${root_path}"/logs/env.log
