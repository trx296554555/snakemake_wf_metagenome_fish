root_path=$(dirname `dirname \`dirname "$CONDA_PREFIX"_\``)
env_name=$(basename $CONDA_PREFIX)

# 使用 metabat2 concoct maxbin2 对组装的contigs进行binning
# concoct 需要一个纯净的环境否则会遇到 endless OpenBLAS Warning
# 这是由于OpenBLAS的优先级比MKL高，concoct使用OpenBLAS就会出现Warning

{
  echo "Make additional adjustments for the post-deployment of the Conda environment ${env_name} (binning)"
  echo "Binning the assembled Contigs using Metabat2 Concoct MaxBin2."
  echo "---------------------------"
} >> "${root_path}"/logs/env.log