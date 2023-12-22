root_path=$(dirname `dirname \`dirname "$CONDA_PREFIX"_\``)
env_name=$(basename $CONDA_PREFIX)
config_path="${root_path}"/workflow/config
db_root=$(awk -F ': ' '/db_root:/{print $2}' "${config_path}"/config.yaml)
dbcan_db=$(awk -F ': ' '/    dbcan:/{print $2}' "${config_path}"/config.yaml)
vfdb_db=$(awk -F ': ' '/    vfdb:/{print $2}' "${config_path}"/config.yaml)
# 使用 dbcan rgi vfdb 对预测的蛋白质进行功能注释

# 使用 dbcan vfdb 对预测的蛋白质进行功能注释
{
  echo "Make additional adjustments for the post-deployment of the Conda environment ${env_name} (dbcan)"
  echo "Annotation the assembled Contigs using dbcan."
  echo "Annotation the assembled Contigs using vfdb by diamond."
  echo "---------------------------"
} >> "${root_path}"/logs/env.log

# 原有run_dbcan脚本的线程限制有问题，暂时使用自己修改的版本替代 TODO 后续可能会改回正式版 在她修复之后
if [ -f "${root_path}"/workflow/scripts/dbcan4_revise/run_dbcan.py ];then
  mv $CONDA_PREFIX/lib/python3.8/site-packages/dbcan_cli/run_dbcan.py $CONDA_PREFIX/lib/python3.8/site-packages/dbcan_cli/run_dbcan_bak.py
  cp "${root_path}"/workflow/scripts/dbcan4_revise/run_dbcan.py $CONDA_PREFIX/lib/python3.8/site-packages/dbcan_cli/run_dbcan.py
  cp -f "${root_path}"/workflow/scripts/dbcan4_revise/hmmer_parser.py $CONDA_PREFIX/lib/python3.8/site-packages/dbcan_cli/hmmer_parser.py
fi

# Check whether the dbcan4 database already exists
if [ -d "${db_root}"/"${dbcan_db}" ];then
  echo "The dbcan4 database already exists, location：${db_root}/${dbcan_db}" >> "${root_path}"/logs/env.log
else
  echo "The dbcan4 database doesn't exist. Use the following command to build it yourself:" >> "${root_path}"/logs/env.log
  # 构建数据库
cat >> "${root_path}"/logs/env.log <<EOF
  ------------------------------------------------------------
  test -d "${db_root}"/"${dbcan_db}" || mkdir "${db_root}"/"${dbcan_db}" \
  cd "${db_root}"/"${dbcan_db}" \
  && wget http://bcb.unl.edu/dbCAN2/download/Databases/fam-substrate-mapping-08012023.tsv \
	&& wget http://bcb.unl.edu/dbCAN2/download/Databases/PUL_12112023.faa && mv PUL_12112023.faa PUL.faa && makeblastdb -in PUL.faa -dbtype prot \
	&& wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL_12-12-2023.xlsx \
  && wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL_12-12-2023.txt \
	&& wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL.tar.gz && tar xvf dbCAN-PUL.tar.gz \
  && wget https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN_sub.hmm && hmmpress dbCAN_sub.hmm \
  && wget https://bcb.unl.edu/dbCAN2/download/Databases/V12/CAZyDB.07262023.fa && diamond makedb --in CAZyDB.07262023.fa -d CAZy \
  && wget https://bcb.unl.edu/dbCAN2/download/Databases/V12/dbCAN-HMMdb-V12.txt && mv dbCAN-HMMdb-V12.txt dbCAN.txt && hmmpress dbCAN.txt \
  && wget https://bcb.unl.edu/dbCAN2/download/Databases/V12/tcdb.fa && diamond makedb --in tcdb.fa -d tcdb \
  && wget http://bcb.unl.edu/dbCAN2/download/Databases/V12/tf-1.hmm && hmmpress tf-1.hmm \
  && wget http://bcb.unl.edu/dbCAN2/download/Databases/V12/tf-2.hmm && hmmpress tf-2.hmm \
  && wget https://bcb.unl.edu/dbCAN2/download/Databases/V12/stp.hmm && hmmpress stp.hmm \
  ------------------------------------------------------------
EOF
exit 1
fi

# Check whether the vfdb database already exists
if [ -f "${db_root}"/"${vfdb_db}" ];then
  echo "The vfdb database already exists, location：${db_root}/${vfdb_db}" >> "${root_path}"/logs/env.log
  echo "But diamond need a perfectly matched version of the database, so check version and rebuild the db if error occurs" >> "${root_path}"/logs/env.log
else
  echo "The vfdb database doesn't exist. Use the following command to build it yourself:" >> "${root_path}"/logs/env.log
  # 构建数据库
cat >> "${root_path}"/logs/env.log <<EOF
    ------------------------------------------------------------
    # download database
    mkdir -p `dirname "${db_root}"/"${vfdb_db}"` && cd `dirname "${db_root}"/"${vfdb_db}"`
    wget -c http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz
    wget -c http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz
    tar -zxvf VFDB_setA_pro.fas.gz
    tar -zxvf VFDB_setB_pro.fas.gz
    diamond makedb --in VFDB_setA_pro.fas -d VF_core_db
    diamond makedb --in VFDB_setB_pro.fas -d VF_full_db
    ------------------------------------------------------------
EOF
exit 1
fi