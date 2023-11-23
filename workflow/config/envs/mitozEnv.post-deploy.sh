root_path=$(dirname `dirname \`dirname "$CONDA_PREFIX"_\``)
env_name=$(basename $CONDA_PREFIX)

# 从动物体外粪便样本或解剖的肠道中收集的粪便样本中，宏基因组测序数据通常包含大量宿主基因组序列。
# 通过使用 mitoZ 工具，我们能够从中提取出宿主的线粒体基因组序列。
# 需要注意的是，一些体外粪便样本中，由于宿主序列较少，可能导致 mitoZ 报错并输出 "task_failed.error" 结果文件。

{
  echo "Make additional adjustments for the post-deployment of the Conda environment ${env_name} (mitozEnv)"
  echo "From fecal samples collected externally from animals or from the intestines obtained through dissection, "
  echo "the obtained metagenomic data often contain a substantial amount of host genomic sequences. "
  echo "Using mitoZ, it is possible to extract the mitochondrial genome sequences of the host from these metagenomic datasets."
  echo "However, it is important to note that in some external fecal samples, the presence of a limited amount of "
  echo "host sequences may lead to an error in mitoZ, resulting in a [task_failed.error] output in the results folder."
  echo "---------------------------"
} >> "${root_path}"/logs/env.log