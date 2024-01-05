root_path=$(dirname `dirname \`dirname "$CONDA_PREFIX"_\``)
env_name=$(basename "$CONDA_PREFIX")
config_path="${root_path}"/workflow/config
db_root=$(awk -F ': ' '/db_root:/{print $2}' "${config_path}"/config.yaml)
checkm_db=$(awk -F ': ' '/    checkm:/{print $2}' "${config_path}"/config.yaml)

# 使用 drep 对分好的bin进行合并去冗余
# 使用checkM评估从宏基因组中组装的基因组的质量。它提供了基因组完整性和污染的可靠估计。
{
  echo "Make additional adjustments for the post-deployment of the Conda environment ${env_name} (drep)"
  echo "Bins de-replication using drep."
  echo "---------------------------"
} >> "${root_path}"/logs/env.log

# mamba install -c bioconda hmmer pplacer -y
ln -s "${db_root}"/software/ANIcalculator_v1/nsimscan "$CONDA_PREFIX"/bin/
ln -s "${db_root}"/software/ANIcalculator_v1/ANIcalculator "$CONDA_PREFIX"/bin/
ln -s "$CONDA_PREFIX"/lib/libgsl.so "$CONDA_PREFIX"/lib/libgsl.so.25
ln -s "$CONDA_PREFIX"/lib/libnsl.so "$CONDA_PREFIX"/lib/libnsl.so.1
checkm data setRoot "${db_root}"/"${checkm_db}"
