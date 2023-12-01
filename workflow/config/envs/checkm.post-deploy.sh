root_path=$(dirname `dirname \`dirname "$CONDA_PREFIX"_\``)
env_name=$(basename $CONDA_PREFIX)

# 使用checkM评估从宏基因组中组装的基因组的质量。它提供了基因组完整性和污染的可靠估计。

{
  echo "Make additional adjustments for the post-deployment of the Conda environment ${env_name} (checkM)"
  echo "Assess the quality of genomes recovered from metagenomes using checkM."
  echo "It offers robust estimates of genome completeness and contamination."
  echo "And prodigal in this environment was also used to predict protein-coding genes in contigs"
  echo "---------------------------"
} >> "${root_path}"/logs/env.log