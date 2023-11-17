# Desc: 设置工作路径和库路径，安装miniforge，配置conda源，安装snakemake

#### Step 0. 设置工作路径和库路径 ####
if [ "$#" -eq "2" ];then
	prj_path=$1
	db_path=$2
#elif [ "$#" -eq "1" ];then
#	prj_path=$1
  #	修改config.yaml中的root路径和db_root路径
	sed -i '1s/.*$/root: '$1'/' $prj_path/workflow/config/config.yaml
  sed -i '2s/.*$/db_root: '$2'/' $prj_path/workflow/config/config.yaml
else
	echo "Usage: `basename $0` {prj_path(full path need)} {db_path(full path need)"
	echo "You provided $# parameters,but 2 are required."
	exit 0
fi

###### Step 1. 安装miniforge，配置conda源 ######
# 检查Conda是否已安装
if command -v conda &> /dev/null; then
    echo "Conda已安装."
else
    echo "Conda未安装，开始安装..."
    # 安装miniforge
    cd ~
    # github下载可能较慢(30min左右)
    if [ ! -f Mambaforge-23.3.1-1-Linux-x86_64.sh ];then
      wget -c https://github.com/conda-forge/miniforge/releases/download/23.3.1-1/Mambaforge-23.3.1-1-Linux-x86_64.sh
    fi
    # 执行安装miniforge，-b表示所有选项选yes(不用手动干预)
    bash Mambaforge-23.3.1-1-Linux-x86_64.sh -b
    # 将Conda添加到环境变量
    echo 'export PATH="$HOME/mambaforge/bin:$PATH"' >> $HOME/.bashrc

### 注意，EOF语句不允许前面有tab或空格，所以此处不缩进
cat >> $HOME/.bashrc <<EOF
# >>> conda initialize >>>
if [ -f "~/mambaforge/etc/profile.d/conda.sh" ]; then
    . "~/mambaforge/etc/profile.d/conda.sh"
fi
if [ -f "~/mambaforge/etc/profile.d/mamba.sh" ]; then
    . "~/mambaforge/etc/profile.d/mamba.sh"
fi
# conda activate
. ~/mambaforge/bin/activate
# <<< conda initialize <<<
EOF
    echo "Conda安装完成."
fi

# 刷新当前终端的环境变量
source $HOME/.bashrc

# 设置conda使用srtict模式，适应snakemake
# conda config --set channel_priority strict

# 配置conda源
cat > $HOME/.condarc <<EOF
channels:
  - defaults
show_channel_urls: true
default_channels:
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/msys2
custom_channels:
  conda-forge: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  msys2: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  bioconda: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  menpo: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  pytorch: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  pytorch-lts: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  simpleitk: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  deepmodeling: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/
EOF

# mamba安装snakemake
# 指定要检测的环境名称
env_name="snakemake"
# 使用conda env list命令获取已安装环境列表，并使用grep命令检查环境是否存在
conda_env_list=$(conda env list | awk '{print $1}' | grep "^$env_name$")
# 检查环境是否存在
if [[ -z "$conda_env_list" ]]; then
    echo "环境 $env_name 不存在，开始创建..."
    mamba create -c conda-forge -c bioconda -n snakemake snakemake -y
else
    echo "环境 $env_name 已存在"
fi

# 配置snakemake环境下的pip源
source activate snakemake
pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple
python -m pip install --upgrade pip
pip install httpx tqdm

## TODO
# 正式发布后需要将注释替换为英文版本