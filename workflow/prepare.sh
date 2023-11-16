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
mamba create -c conda-forge -c bioconda -n snakemake snakemake -y

# 配置snakemake环境下的pip源
conda activate snakemake
pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple
python -m pip install --upgrade pip
pip install httpx tqdm

## TODO
# 正式发布后需要将注释替换为英文版本