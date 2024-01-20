root_path=$(dirname `dirname \`dirname "$CONDA_PREFIX"_\``)
env_name=$(basename $CONDA_PREFIX)
config_path="${root_path}"/workflow/config
db_root=$(awk -F ': ' '/db_root:/{print $2}' "${config_path}"/config.yaml)
eggnog_db=$(awk -F ': ' '/    eggnog:/{print $2}' "${config_path}"/config.yaml)

## TODO eggnog已经更新到了6.0 但是eggnog-mapper v2 目前还不支持6.0版本的数据库 持续关注
# 使用 eggnog-mapper 对预测的蛋白质进行功能注释
{
  echo "Make additional adjustments for the post-deployment of the Conda environment ${env_name} (emapper)"
  echo "Annotation the assembled Contigs using eggnog-mapper."
  echo "---------------------------"
} >> "${root_path}"/logs/env.log


# Check whether the eggnog database already exists
if [ -d "${db_root}"/"${eggnog_db}" ];then
  echo "The eggnog database already exists, location：${db_root}/${eggnog_db}" >> "${root_path}"/logs/env.log
else
  echo "The eggnog database doesn't exist. Use the following command to build it yourself:" >> "${root_path}"/logs/env.log
  # 构建数据库
cat >> "${root_path}"/logs/env.log <<EOF
  ------------------------------------------------------------
  test -d "${db_root}"/"${eggnog_db}" || mkdir "${db_root}"/"${eggnog_db}" \
  cd "${db_root}"/"${eggnog_db}"
  conda activate $CONDA_PREFIX
  download_eggnog_data.py -y --data_dir "${db_root}"/"${eggnog_db}"
  create_dbs.py -m diamond --dbname eggnog_proteins.bact_arch --taxa Bacteria,Archaea --data_dir "${db_root}"/"${eggnog_db}"
  ------------------------------------------------------------
EOF
exit 1
fi
