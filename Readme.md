# 鱼类宏基因组分析流程 workflow by snakemake
---
## Prepare Data and Database

```shell
cd /tmp_disk
cp Mambaforge-23.3.1-1-Linux-x86_64.sh ~
for i in `cat wukong_list.txt`;do cp -r 00_cleandata/$i /data/sas-2-15t/fish/00_cleandata;done
cp -r ref_db /home/user003/fish
cp -r database /data-1
chown -R user003:user003 /data/sas-2-15t/fish
```

## Init

```shell
cd YOUR_WORK_PATH
git clone https://github.com/trx296554555/snakemake_wf_metagenome_fish.git
ln -s snakemake_wf_metagenome_fish/workflow/
bash workflow/prepare.sh
source ~/.bashrc
# edit config.yaml by yourself
```

## Update

```shell
git pull
```

## Run

```shell
cd YOUR_WORK_PATH
conda activate snakemake
snakemake -c12 --use-conda --conda-create-envs-only
snakemake -c96 --use-conda -np
nohup snakemake -c96 --use-conda &
```

---
# Note:
- [kneaddata v0.12.0](https://github.com/biobakery/biobakery/wiki/kneaddata)去宿主，单个任务12线程耗时约120分钟，RAM内存不超过20G，磁盘空间临时占用不超过100G
- [kraken2 v2.1.3](https://github.com/DerrickWood/kraken2/wiki/Manual)分类，单个任务24线程耗时约1分钟，RAM内存不超过80G，磁盘空间临时占用不超过10G
- [bracken v2.9](https://github.com/jenniferlu717/Bracken)对kraken2结果进行校正，单个任务1线程耗时约0.1分钟，RAM内存不超过1G，磁盘空间临时占用不超过1G
- [mitoZ v3.6](https://github.com/linzhi2013/MitoZ/wiki)从动物肠道提取的宏基因组不可避免的含有大量宿主的基因组序列，使用mitoZ从中提取出宿主的线粒体基因组序列，单个任务24线程耗时约10分钟，RAM内存不超过10G，磁盘空间临时占用不超过10G
- [megahit v 1.2.9](https://github.com/voutcn/megahit/wiki)将来自微生物群体的DNA片段组装成contigs,单个任务24线程耗时约30分钟，RAM内存不超过10G，磁盘空间临时占用不超过10G
- [prodigal v2.6.3](https://github.com/hyattpd/Prodigal/wiki)预测contigs的基因，单个任务1线程耗时约10分钟，RAM内存不超过1G，磁盘空间临时占用不超过1G


# 文件和目录是名词+动词的形式，如：readsClassify
# rule 是动词+名词的形式，如：classify_reads