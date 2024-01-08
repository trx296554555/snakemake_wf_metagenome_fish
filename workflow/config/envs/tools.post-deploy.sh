root_path=$(dirname `dirname \`dirname "$CONDA_PREFIX"_\``)
env_name=$(basename $CONDA_PREFIX)
config_path="${root_path}"/workflow/config
db_root=$(awk -F ': ' '/db_root:/{print $2}' "${config_path}"/config.yaml)
taxonkit_db=$(awk -F ': ' '/    taxonkit:/{print $2}' "${config_path}"/config.yaml)

# 包括一些常用的基础工具
# samtools bowtie2 seqkit taxonkit

# Check whether the taxonkit database already exists
if [ -d "${db_root}"/"${taxonkit_db}" ];then
  conda env config vars set TAXONKIT_DB="${db_root}/${taxonkit_db}"
  echo "taxonkit database already exists, location：${db_root}/${taxonkit_db}" >> "${root_path}"/logs/env.log
else
  echo "The taxonkit database doesn't exist. Use the following command to build it yourself:" >> "${root_path}"/logs/env.log
  exist 1
fi

{
  echo "Make additional adjustments for the post-deployment of the Conda environment ${env_name} (tools)"
  echo "Some useful tools: samtools bowtie2 seqkit taxonkit"
  echo "---------------------------"
} >> "${root_path}"/logs/env.log