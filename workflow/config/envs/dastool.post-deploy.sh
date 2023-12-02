root_path=$(dirname `dirname \`dirname "$CONDA_PREFIX"_\``)
env_name=$(basename $CONDA_PREFIX)

# 使用 DAS tool 对不同软件binning结果合并

{
  echo "Make additional adjustments for the post-deployment of the Conda environment ${env_name} (dastool)"
  echo "Binning the assembled Contigs using Metabat2 Concoct MaxBin2."
  echo "Use DAS tool to refine the above binning results."
  echo "---------------------------"
} >> "${root_path}"/logs/env.log