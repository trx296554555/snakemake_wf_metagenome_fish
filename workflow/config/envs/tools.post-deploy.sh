root_path=$(dirname `dirname \`dirname "$CONDA_PREFIX"_\``)
env_name=$(basename $CONDA_PREFIX)

# 包括一些常用的基础工具
# samtools bowtie2 seqkit taxonkit

{
  echo "Make additional adjustments for the post-deployment of the Conda environment ${env_name} (tools)"
  echo "Some useful tools: samtools bowtie2 seqkit taxonkit"
  echo "---------------------------"
} >> "${root_path}"/logs/env.log