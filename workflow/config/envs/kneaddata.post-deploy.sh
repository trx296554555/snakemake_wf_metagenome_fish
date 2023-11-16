# 调整kneaddata配置环境
# 为了适应kneaddata只能通过java -Xmx500m -jar trimmomatic* 调用 trimmomatic，将jar文件链到 which trimmomaticls
rm $CONDA_PREFIX/bin/trimmomatic
ln -s $CONDA_PREFIX/share/trimmomatic/trimmomatic.jar $CONDA_PREFIX/bin/trimmomatic