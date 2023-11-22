root_path=$(dirname `dirname \`dirname "$CONDA_PREFIX"_\``)
config_path="${root_path}"/workflow/config
env_name=$(basename $CONDA_PREFIX)

# Read the configuration
db_root=$(awk -F ': ' '/db_root:/{print $2}' "${config_path}"/config.yaml)
kraken2_db=$(awk -F ': ' '/kraken2:/{print $2}' "${config_path}"/config.yaml)
reads_length=$(awk -F ': ' '/reads_length:/{print $2}' "${config_path}"/config.yaml)

{
  echo "Make additional adjustments for the post-deployment of the Conda environment ${env_name} (kraken2)"
  echo db_root: "$db_root"
  echo kraken2_db_version: "$kraken2_db"
  echo reads_length: "$reads_length"
} >> "${root_path}"/logs/env.log

# Check whether the kraken2 database already exists
if [ -d "${db_root}"/"${kraken2_db}" ];then
  echo "Kraken 2 database already exists, location：${db_root}/${kraken2_db}" >> "${root_path}"/logs/env.log
else
  echo "The Kraken 2 database doesn't exist. Use the following command to build it yourself:" >> "${root_path}"/logs/env.log
  # 构建数据库，报错参考 https://github.com/DerrickWood/kraken2/issues/508
  echo "conda activate $CONDA_PREFIX && cd ${db_root} && kraken2-build --standard --threads 12 --db ${kraken2_db}" >> "${root_path}"/logs/env.log
fi

# Check if the bracken index already exists
if [ -f "${db_root}"/"${kraken2_db}"/database.kraken ];then
  echo "The bracken index already exists." >> "${root_path}"/logs/env.log
else
  echo "The bracken index doesn't exist. Use the following command to build it yourself:" >> "${root_path}"/logs/env.log
  # Build bracken index
  echo "conda activate $CONDA_PREFIX && cd ${db_root} && bracken-build -d ${db_root}/${kraken2_db} -t 48 -l ${reads_length}" >> "${root_path}"/logs/env.log
fi

echo "---------------------------" >> "${root_path}"/logs/env.log
