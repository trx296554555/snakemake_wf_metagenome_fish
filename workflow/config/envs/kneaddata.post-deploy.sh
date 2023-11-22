root_path=$(dirname `dirname \`dirname "$CONDA_PREFIX"_\``)
env_name=$(basename $CONDA_PREFIX)

# 为了适应kneaddata只能通过java -Xmx500m -jar trimmomatic* 调用 trimmomatic，将jar文件链为trimmomatic
rm $CONDA_PREFIX/bin/trimmomatic
ln -s $CONDA_PREFIX/share/trimmomatic/trimmomatic.jar $CONDA_PREFIX/bin/trimmomatic

{
  echo "Make additional adjustments for the post-deployment of the Conda environment ${env_name} (kneaddata)"
  echo "Create a symbolic link named trimmomatic for the trimmomatic.jar file to avoid potential false calls."
  echo "---------------------------"
} >> "${root_path}"/logs/env.log