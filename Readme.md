# 鱼类宏基因组分析流程 workflow by snakemake
---
## Init

```shell
cd YOUR_WORK_PATH
git clone https://github.com/trx296554555/snakemake_wf_metagenome_fish.git
ln -s snakemake_wf_metagenome_fish/workflow/
bash workflow/prepare.sh
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
snakemake -c96 --use-conda -np
nohup snakemake -c96 --use-conda &
```

---
# Note:
- kneaddata去宿主，单个任务80线程耗时约20分钟，RAM内存不超过20G，磁盘空间临时占用不超过100G
