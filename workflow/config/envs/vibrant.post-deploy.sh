root_path=$(dirname `dirname \`dirname "$CONDA_PREFIX"_\``)
env_name=$(basename $CONDA_PREFIX)
config_path="${root_path}"/workflow/config
db_root=$(awk -F ': ' '/db_root:/{print $2}' "${config_path}"/config.yaml)
vibrant_db=$(awk -F ': ' '/    vibrant:/{print $2}' "${config_path}"/config.yaml)

# 检查是否已经存在 vibrant 数据库
if [ -d "${db_root}"/"${vibrant_db}" ];then
  echo "vibrant database already exists, location：${db_root}/${vibrant_db}" >> "${root_path}"/logs/env.log
else
  echo "The vibrant database doesn't exist. Use the following command to build it yourself:" >> "${root_path}"/logs/env.log
  echo "download-db.sh ${db_root}/${vibrant_db}" >> "${root_path}"/logs/env.log
  exist 1
fi

{
  echo "Make additional adjustments for the post-deployment of the Conda environment ${env_name} (vibrant)"
  echo "Virus prediction by vibrant"
  echo "---------------------------"
} >> "${root_path}"/logs/env.log