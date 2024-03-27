root_path=$(dirname `dirname \`dirname "$CONDA_PREFIX"_\``)
env_name=$(basename $CONDA_PREFIX)
config_path="${root_path}"/workflow/config
db_root=$(awk -F ': ' '/db_root:/{print $2}' "${config_path}"/config.yaml)
virsorter2_db=$(awk -F ': ' '/    virsorter2:/{print $2}' "${config_path}"/config.yaml)

# 检查是否已经存在 virsorter 2 数据库
if [ -f "${db_root}"/"${virsorter2_db}/Done_all_setup" ];then
  virsorter config --init-source --db-dir ${db_root}"/"${virsorter2_db}
  virsorter config --set HMMSEARCH_THREADS=12
  echo "virsorter2 database already exists, location：${db_root}/${virsorter2_db}" >> "${root_path}"/logs/env.log

else
  echo "The virsorter2 database doesn't exist. Use the following command to build it yourself:" >> "${root_path}"/logs/env.log
  # 构建数据库
  virsorter setup -d ${db_root}"/"${virsorter2_db} -j 24
  virsorter config --init-source --db-dir ${db_root}"/"${virsorter2_db}
  virsorter config --set HMMSEARCH_THREADS=12
fi

{
  echo "Make additional adjustments for the post-deployment of the Conda environment ${env_name} (virsorter2)"
  echo "Virus prediction by virsorter2"
  echo "---------------------------"
} >> "${root_path}"/logs/env.log