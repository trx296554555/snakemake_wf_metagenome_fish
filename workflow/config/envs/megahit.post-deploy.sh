root_path=$(dirname `dirname \`dirname "$CONDA_PREFIX"_\``)
env_name=$(basename $CONDA_PREFIX)

# MEGAHIT 的主要作用是将来自微生物群体的DNA片段组装成连续的序列，称为contigs

{
  echo "Make additional adjustments for the post-deployment of the Conda environment ${env_name} (megahit)"
  echo "MEGAHIT is used to assemble DNA reads from microbial communities into contiguous sequences (contigs)."
  echo "---------------------------"
} >> "${root_path}"/logs/env.log