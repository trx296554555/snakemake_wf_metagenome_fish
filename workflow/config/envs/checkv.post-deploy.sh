root_path=$(dirname `dirname \`dirname "$CONDA_PREFIX"_\``)
env_name=$(basename $CONDA_PREFIX)
config_path="${root_path}"/workflow/config
db_root=$(awk -F ': ' '/db_root:/{print $2}' "${config_path}"/config.yaml)
checkv_db=$(awk -F ': ' '/    checkv:/{print $2}' "${config_path}"/config.yaml)

# 检查是否已经存在 checkv 数据库
if [ -d "${db_root}"/"${checkv_db}" ];then
  echo "checkv database already exists, location：${db_root}/${checkv_db}" >> "${root_path}"/logs/env.log
else
  echo "The vibrant database doesn't exist. Use the following command to build it yourself:" >> "${root_path}"/logs/env.log
  echo "checkv download_database ${db_root}/${checkv_db}" >> "${root_path}"/logs/env.log
  exist 1
fi

{
  echo "Make additional adjustments for the post-deployment of the Conda environment ${env_name} (checkv)"
  echo "Use checkv to detect the quality of the identified viral sequences"
  echo "---------------------------"
} >> "${root_path}"/logs/env.log