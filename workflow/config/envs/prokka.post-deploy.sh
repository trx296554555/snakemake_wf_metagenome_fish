root_path=$(dirname `dirname \`dirname "$CONDA_PREFIX"_\``)
env_name=$(basename $CONDA_PREFIX)

# 使用prokkka进行注释
{
  echo "Make additional adjustments for the post-deployment of the Conda environment ${env_name} (prokka)"
  echo "Use prokka to annotate bins"
  echo "---------------------------"
} >> "${root_path}"/logs/env.log