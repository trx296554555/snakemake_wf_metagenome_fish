root_path=$(dirname `dirname \`dirname "$CONDA_PREFIX"_\``)
env_name=$(basename $CONDA_PREFIX)

# 用于组装好的线粒体序列CDS推断系统发育树并构建物种树

{
  echo "Make additional adjustments for the post-deployment of the Conda environment ${env_name} (buildtree)"
  echo "CDS is used to assemble mitochondrial sequences to infer phylogenetic trees and construct species trees"
  echo "---------------------------"
} >> "${root_path}"/logs/env.log