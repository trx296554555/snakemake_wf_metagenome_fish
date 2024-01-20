# 鱼类宏基因组分析流程 workflow by snakemake
---

## Init Environment

```shell
YOUR_WORK_PATH=/home/user003/fish
REAL_LOCATION=/data/sas-2-15t/fish
mkdir -p ${REAL_LOCATION} && ln -s ${REAL_LOCATION} ${YOUR_WORK_PATH}
cd ${YOUR_WORK_PATH}
# http的方式连接失败可以尝试使用ssh的方式 https://www.cnblogs.com/zxtceq/p/14037175.html
git clone https://github.com/trx296554555/snakemake_wf_metagenome_fish.git
ln -s snakemake_wf_metagenome_fish/workflow/
bash workflow/prepare.sh
source ~/.bashrc
# edit config.yaml by yourself
```

## Prepare Data and Database

```shell
cp /tmp_disk/Mambaforge-23.3.1-1-Linux-x86_64.sh ~
cp -r -n /tmp_disk/ref_db ${YOUR_WORK_PATH}
cp -r -n /tmp_disk/database $(dirname ${REAL_LOCATION}) && ln -s $(dirname ${REAL_LOCATION})/database ~
cd ${YOUR_WORK_PATH}/workflow/scripts && python prepare_sampleList.py $(hostname)
for i in $(tail -n +2 ${YOUR_WORK_PATH}/workflow/config/sampleList.csv | cut -f2 -d ,); do cp -r -n /tmp_disk/00_cleandata/$i ${YOUR_WORK_PATH}/00_cleandata/;done && chown -R user003:user003 ${YOUR_WORK_PATH}/00_cleandata/
```

## Update

```shell
git pull
```

## Run

```shell
screen -S fish
cd ${YOUR_WORK_PATH}
conda activate snakemake
snakemake -c12 --use-conda --conda-create-envs-only
snakemake -c96 --use-conda -np
ctrl + a + d
snakemake --forceall --rulegraph |sed '1d'| dot -Tpdf > fish_rules.pdf
tar -chf $(hostname).all_bins.tar all_bins/*.fa && pigz $(hostname).all_bins.tar
tar -chf $(hostname).report.tar reports/* && pigz $(hostname).report.tar && sz $(hostname).report.tar.gz && rm -f $(hostname).report.tar.gz 
cd /home/user003/fish/10_bin_classify/ && cp MAGs_classify_report.tsv $(hostname).MAGs_classify_report.tsv && sz $(hostname).MAGs_classify_report.tsv && rm -f $(hostname).MAGs_classify_report.tsv
```

---

| Software                                                                              | Version | Description       | Core | Ram  | Disk  | Time    |
|---------------------------------------------------------------------------------------|---------|-------------------|------|------|-------|---------|
| [kneaddata](https://github.com/biobakery/biobakery/wiki/kneaddata)                    | v0.12.0 | 去宿主               | 12   | <20G | <100G | ≈120min |
| [kraken2](https://github.com/DerrickWood/kraken2/wiki/Manual)                         | v2.1.3  | 分类                | 24   | <80G | <10G  | ≈1min   |
| [krakenTools](https://github.com/jenniferlu717/KrakenTools)                           | v1.2    | 下游分析              | 1    | <1G  | <1G   | ≈1min   |
| [bracken](https://github.com/jenniferlu717/Bracken)                                   | v2.9    | 校正                | 1    | <1G  | <1G   | ≈0.1min |
| [mitoZ](https://github.com/linzhi2013/MitoZ/wiki)                                     | v3.6    | 宿主线粒体             | 24   | <10G | <10G  | ≈10min  |
| [megahit](https://github.com/voutcn/megahit/wiki)                                     | v1.2.9  | 组装                | 24   | <10G | <10G  | ≈30min  |
| [prodigal](https://github.com/hyattpd/Prodigal/wiki)                                  | v2.6.3  | 预测基因              | 12   | <1G  | <1G   | ≈10min  |
| [humman3](https://github.com/biobakery/humann)                                        | v3.0.0  | 功能注释              | 12   | <10G | <10G  | ≈30min  |
| [dbcan4](https://github.com/linnabrown/run_dbcan)                                     | v4.1.0  | 碳水化合物酶            | 12   | <10G | <10G  | ≈30min  |
| [rgi](https://github.com/arpcard/rgi)                                                 | v6.0.3  | 抗性基因              | 12   | <10G | <10G  | ≈30min  |
| [vfdb](http://www.mgc.ac.cn/VFs/main.htm)                                             | v2024   | 毒力因子              | 12   | <10G | <10G  | ≈30min  |
| [metabat2](https://bitbucket.org/berkeleylab/metabat/wiki/Best%20Binning%20Practices) | v2.15   | binning           | 12   | <20G | <1G   | ≈10min  |
| [CONCOCT](https://github.com/BinPro/CONCOCT)                                          | v1.1.0  | binning           | 12   | <1G  | <1G   | ≈30min  |
| [maxbin2](https://sourceforge.net/p/maxbin/code/ci/master/tree/)                      | v2.2.7  | binning           | 12   | <1G  | <1G   | ≈1min   |
| [MetaBinner](https://github.com/ziyewang/MetaBinner)                                  | v1.4.4  | binning           | 12   | <10G | <10G  | ≈45min  |
| [DAS Tool](https://github.com/cmks/DAS_Tool)                                          | v1.1.6  | refine bins       | 12   | <10G | <1G   | ≈30min  |
| [checkm](https://github.com/Ecogenomics/CheckM/wiki)                                  | v1.2.2  | bin quality       | 64   | <10G | <1G   | ≈10min  |
| [drep](https://github.com/MrOlm/drep)                                                 | v.3.4.5 | bin dereplication | 64   | <10G | <1G   | ≈10min  |
| [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk)                                      | v2.3.2  | bin taxonomy      | 12   | <10G | <1G   | ≈10min  |
| [prokka](https://github.com/tseemann/prokka)                                          | v1.14.6 | bins功能注释          | 12   | <10G | <1G   | ≈10min  |
| [seqkit](https://bioinf.shenwei.me/seqkit/)                                           | v2.6.1  | 序列处理              | 12   | <1G  | <1G   | ≈1min   |
| [newick_utils](https://github.com/tjunier/newick_utils)                               | v1.6    | phylogenetic tree | 12   | <1G  | <1G   | ≈1min   |

# 文件和目录是名词+动词的形式，如：readsClassify

# rule 是动词+名词的形式，如：classify_reads

##   

# TODO 这里以后可能会更新，这是run_dbcan推荐的数据库更新了，但是代码没更新的原因

ln -s "${db_root}"/"${dbcan_db}"/fam-substrate-mapping-08012023.tsv "${db_root}"/"${dbcan_db}"
/fam-substrate-mapping-08252022.tsv