root_path=$(dirname `dirname \`dirname "$CONDA_PREFIX"_\``)
env_name=$(basename $CONDA_PREFIX)

# 使用 salmon 不进行比对对contig进行定量

{
  echo "Make additional adjustments for the post-deployment of the Conda environment ${env_name} (salmon)"
  echo "---------------------------"
} >> "${root_path}"/logs/env.log