# 获取kraken2.post-deploy.sh脚本所在目录的父目录（项目目录）
script_dir=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")
config_path=$(dirname "$script_dir")

# 从config.yaml中读取配置，获得db_root路径和kraken2_db的名字
db_root=$(awk -F ': ' '/db_root:/{print $2}' $config_path/config.yaml)
kraken2_db=$(awk -F ': ' '/kraken2:/{print $2}' $config_path/config.yaml)

# 检查kraken2数据库是否已存在
if [ -d $db_root/$kraken2_db ];then
  echo "kraken2数据库已存在，位置：${db_root}/${kraken2_db}."
else
  echo "kraken2数据库不存在，开始下载..."
  # 下载数据库
  cd $db_root
  conda activate kraken2
  # 构建数据库，报错参考 https://github.com/DerrickWood/kraken2/issues/508
  kraken2-build --standard --threads 12 --db $kraken2_db
  echo "数据库下载构建完成，位置：${db_root}/${kraken2_db}."
fi

# 检查bracken索引是否已存在
if [ -f $db_root/${db_root}/${kraken2_db}/database.kraken ];then
  echo "bracken索引已存在."
else
  echo "bracken索引不存在，开始下载..."
  # 构建bracken索引
  cd $db_root
