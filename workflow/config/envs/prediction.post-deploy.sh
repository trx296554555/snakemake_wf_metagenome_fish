root_path=$(dirname `dirname \`dirname "$CONDA_PREFIX"_\``)
env_name=$(basename $CONDA_PREFIX)

# 使用prodigal预测contigs中的蛋白质编码基因,使用cd-hit进行去冗余

{
  echo "Make additional adjustments for the post-deployment of the Conda environment ${env_name} (prediction)"
  echo "Use prodigal to predict protein-coding genes in contigs and use cd-hit to remove redundancy"
  echo "Use cd-hit-est to remove redundancy in contigs"
  echo "---------------------------"
} >> "${root_path}"/logs/env.log