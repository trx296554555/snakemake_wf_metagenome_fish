root_path=$(dirname `dirname \`dirname "$CONDA_PREFIX"_\``)
env_name=$(basename $CONDA_PREFIX)
config_path="${root_path}"/workflow/config
db_root=$(awk -F ': ' '/db_root:/{print $2}' "${config_path}"/config.yaml)
deepvirfinder_db=$(awk -F ': ' '/    deepvirfinder:/{print $2}' "${config_path}"/config.yaml)

# 检查是否已经存在 deepvirfinder 模型文件
if [ -d "${db_root}"/"${deepvirfinder_db}/models" ];then
  echo "deepvirfinder model file already exists, location：${db_root}/${deepvirfinder_db}" >> "${root_path}"/logs/env.log
else
  echo "The deepvirfinder model file doesn't exist. Use the following command to build it yourself:" >> "${root_path}"/logs/env.log
  echo "git clone https://github.com/jessieren/DeepVirFinder" >> "${root_path}"/logs/env.log
  exist 1
fi

{
  echo "Make additional adjustments for the post-deployment of the Conda environment ${env_name} (deepvirfinder)"
  echo "Virus prediction by deepvirfinder"
  echo "---------------------------"
} >> "${root_path}"/logs/env.log