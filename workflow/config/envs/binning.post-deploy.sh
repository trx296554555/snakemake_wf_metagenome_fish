root_path=$(dirname `dirname \`dirname "$CONDA_PREFIX"_\``)
env_name=$(basename $CONDA_PREFIX)

# 使用 metabat2 concoct maxbin2 对组装的contigs进行binning
# 使用 DAS tool 对binning结果refine

{
  echo "Make additional adjustments for the post-deployment of the Conda environment ${env_name} (binning)"
  echo "Binning the assembled Contigs using Metabat2 Concoct MaxBin2."
  echo "Use DAS tool to refine the binning results."
  echo "---------------------------"
} >> "${root_path}"/logs/env.log