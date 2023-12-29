import os
import pandas as pd

include: "../scripts/download_ref.py"


def change_name(wrong_dir):
    """
    这不是工作流程中的标准步骤，如果需要，可以自己修改。

    This is not standard step in the workflow, can modify by yourself if you need.
    """
    dir_path = wrong_dir
    sample_id = os.path.basename(wrong_dir)
    change_name_log = open('logs/00_change_name.log','a')
    change_name_log.write("The following file names have been modified:\n")
    for file in os.listdir(wrong_dir):
        if file.endswith(".fq1.gz"):
            os.rename(os.path.join(dir_path,file),os.path.join(dir_path,sample_id + "_1.clean.fq.gz"))
            change_name_log.write(
                f"{os.path.join(dir_path,sample_id + '_1.clean.fq.gz')} {os.path.join(dir_path,file)}\n")
        elif file.endswith(".fq2.gz"):
            os.rename(os.path.join(dir_path,file),os.path.join(dir_path,sample_id + "_2.clean.fq.gz"))
            change_name_log.write(
                f"{os.path.join(dir_path,sample_id + '_2.clean.fq.gz')} {os.path.join(dir_path,file)}\n")
        elif file.endswith('fq1.gz.md5.txt'):
            with open(os.path.join(wrong_dir,'MD5.txt'),'a') as f:
                with open(config["root"] + "/" + config["folder"]["data"] + "/" + sample_id + "/" + file,'r') as f1:
                    f.write(f1.readline().split()[0] + '  ' + sample_id + "_1.clean.fq.gz\n")
            os.remove(os.path.join(dir_path,file))
        elif file.endswith('fq2.gz.md5.txt'):
            with open(os.path.join(wrong_dir,'MD5.txt'),'a') as f:
                with open(config["root"] + "/" + config["folder"]["data"] + "/" + sample_id + "/" + file,'r') as f1:
                    f.write(f1.readline().split()[0] + '  ' + sample_id + "_2.clean.fq.gz\n")
            os.remove(os.path.join(dir_path,file))
        else:
            raise Exception(f'Unexpected file {file} in {sample_id}')
    change_name_log.close()
    return [os.path.join(dir_path,sample_id + "_1.clean.fq.gz"), os.path.join(dir_path,sample_id + "_2.clean.fq.gz")]



def check_all_sample():
    """
    单样本模式：只会处理sampleList.csv中Run_Sample列中指定的样本
    多样本合并模式：按config.yaml中指定列生成分组，
    检查样本与对应的参考基因组是否存在，
    严格检测原始数据文件夹格式，每个样本对应一个文件夹，文件夹名为样本名
    每个文件夹下有三个文件，分别为MD5.txt，样本名_1.clean.fq.gz，样本名_2.clean.fq.gz

    Only samples present in allSampleList.csv will be processed
    Strictly detect the original data folder format, each sample corresponds to a folder
    The folder name is the sample name, which are MD5.txt, sample name_1.clean.fq.gz, sample name_2.clean.fq.gz
    """
    all_samples = get_all_sample()
    data_root = config["root"] + "/" + config["folder"]["data"]
    exist_samples = os.listdir(data_root)
    if not set(all_samples).issubset(exist_samples):
        missing_samples = set(all_samples) - set(exist_samples)
        raise ValueError(f"The following samples do not exist： {missing_samples}")

    # 对文件命名不规范的样本进行处理
    file_list = []
    for dir_ in all_samples:
        this_files = os.listdir(os.path.join(data_root,dir_))
        right_files = ['MD5.txt', f'{dir_}_1.clean.fq.gz', f'{dir_}_2.clean.fq.gz']
        if len(this_files) == 3 and sorted(this_files) == sorted(right_files):
            file_list += [os.path.join(data_root,dir_,file) for file in this_files if file.endswith('.clean.fq.gz')]
        else:
            file_list += change_name(os.path.join(data_root,dir_))

    # 预先检测宿主基因组是否存在，不存在则下载
    download_df = all_sample_df.drop_duplicates(subset='Download_link')
    download_dict = dict(zip(download_df['Ref'],download_df['Download_link']))
    if not os.path.exists(config["root"] + "/" + config["folder"]["index"] + "/all_md5.pass"):
        print("Downloading reference genome...")
        # Download reference genome
        download_ref_fa_file(download_dict,config["root"] + "/" + config["folder"]["index"])
        # Check reference genome md5
        ref_fa_file_md5_check(download_dict,config["root"] + "/" + config["folder"]["index"])
    else:
        print("Reference genome has been downloaded successfully.")
    return file_list


rule check_md5:
    priority: 98
    input:
        config["root"] + "/" + config["folder"]["data"] + "/" + "{sample}/MD5.txt",
    output:
        config["root"] + "/logs/md5_info/{sample}_check_md5.log"
    message:
        "02: Verifying the MD5 checksum of FASTQ files for processing samples ------------------------------------------"
    shell:
        "cd `dirname {input}` && md5sum -c {input} > {output}"


rule check_sample:
    priority: 99
    input:
        file=check_all_sample(),
        md5=expand(config["root"] + "/logs/md5_info/{sample}_check_md5.log",sample=get_all_sample())
    output:
        file=config["root"] + "/logs/01_check_run.log"
    message:
        "01: Checking which samples need processing ------------------------------------------"
    run:
        with open(output.file,'w') as f:
            f.write(f'Run sample num: {len(get_run_sample())}\n')
            f.write(f'Run samples:\n')
            f.write("\n".join(get_run_sample()) + "\n\n")

        if config['co_assemble']['enable']:
            result_df = pd.read_csv(config['root'] + '/workflow/config/co_assemble_list.csv')
            with open(output.file,'a') as f:
                f.write(f'co_assembly {os.listdir(config["root"] + "/logs/co_assemble")} needed.\n')
        else:
            with open(output.file,'a') as f:
                f.write('co_assembly is not needed.\n')
