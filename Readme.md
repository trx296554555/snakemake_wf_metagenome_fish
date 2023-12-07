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

| Software                                                                              | Version | Description | Core | Ram  | Disk  | Time    |
|---------------------------------------------------------------------------------------|---------|-------------|------|------|-------|---------|
| [kneaddata](https://github.com/biobakery/biobakery/wiki/kneaddata)                    | v0.12.0 | 去宿主         | 12   | <20G | <100G | ≈120min |
| [kraken2](https://github.com/DerrickWood/kraken2/wiki/Manual)                         | v2.1.3  | 分类          | 24   | <80G | <10G  | ≈1min   |
| [krakenTools](https://github.com/jenniferlu717/KrakenTools)                           | v1.2    | 下游分析        | 1    | <1G  | <1G   | ≈1min   |
| [bracken](https://github.com/jenniferlu717/Bracken)                                   | v2.9    | 校正          | 1    | <1G  | <1G   | ≈0.1min |
| [mitoZ](https://github.com/linzhi2013/MitoZ/wiki)                                     | v3.6    | 宿主线粒体       | 24   | <10G | <10G  | ≈10min  |
| [megahit](https://github.com/voutcn/megahit/wiki)                                     | v1.2.9  | 组装          | 24   | <10G | <10G  | ≈30min  |
| [prodigal](https://github.com/hyattpd/Prodigal/wiki)                                  | v2.6.3  | 预测基因        | 12   | <1G  | <1G   | ≈10min  |
| [metabat2](https://bitbucket.org/berkeleylab/metabat/wiki/Best%20Binning%20Practices) | v2.15   | binning     | 12   | <20G | <1G   | ≈10min  |
| [CONCOCT](https://github.com/BinPro/CONCOCT)                                          | v1.1.0  | binning     | 12   | <1G  | <1G   | ≈30min  |
| [maxbin2](https://sourceforge.net/p/maxbin/code/ci/master/tree/)                      | v2.2.7  | binning     | 12   | <1G  | <1G   | ≈1min   |
| [MetaBinner](https://github.com/ziyewang/MetaBinner)                                  | v1.4.4  | binning     | 12   | <10G | <10G  | ≈45min  |
| [DAS Tool](https://github.com/cmks/DAS_Tool)                                          | v1.1.6  | refine bins | 12   | <10G | <1G   | ≈30min  |

# 文件和目录是名词+动词的形式，如：readsClassify

# rule 是动词+名词的形式，如：classify_reads