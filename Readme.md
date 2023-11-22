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
- [kneaddata](https://github.com/biobakery/biobakery/wiki/kneaddata)去宿主，单个任务12线程耗时约120分钟，RAM内存不超过20G，磁盘空间临时占用不超过100G
- 
