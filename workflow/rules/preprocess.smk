import os
import shutil
import pandas as pd
include: "../scripts/download_ref.py"

if os.path.exists(f'{config["root"]}/{config["meta"]["sampleList"]}'):
    all_sample_df = pd.read_csv(f'{config["root"]}/{config["meta"]["sampleList"]}')
    all_sample_id = all_sample_df['Sample_ID']
else:
    raise Exception(f'No {config["meta"]["sampleList"]} in {config["root"]}')
    exit(1)
if os.path.exists('logs'):
    pass
else:
    os.mkdir('logs')


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


def get_run_sample():
    root = config["root"] + "/" + config["folder"]["data"]
    run_sample = [sample for sample in os.listdir(root) if sample in all_sample_id.values]
    return run_sample


def check_run_sample():
    """
    检查样本与对应的参考基因组是否存在，只会处理sampleList.csv中存在的样本
    严格检测原始数据文件夹格式，每个样本对应一个文件夹，文件夹名为样本名
    每个文件夹下有三个文件，分别为MD5.txt，样本名_1.clean.fq.gz，样本名_2.clean.fq.gz

    Only samples present in sampleList.csv will be processed
    Strictly detect the original data folder format, each sample corresponds to a folder
    The folder name is the sample name, which are MD5.txt, sample name_1.clean.fq.gz, sample name_2.clean.fq.gz
    """
    root = config["root"] + "/" + config["folder"]["data"]
    exist_sample = os.listdir(root)
    run_sample = [sample for sample in os.listdir(root) if sample in all_sample_id.values]

    unexpect_sample = [sample for sample in exist_sample if sample not in all_sample_id.values]
    if os.path.exists(os.path.join(root,'unexpect_sample')):
        unexpect_sample.remove('unexpect_sample')
    else:
        os.mkdir(os.path.join(root,'unexpect_sample'))
    if os.path.exists(os.path.join(root,'md5_error_sample')):
        unexpect_sample.remove('md5_error_sample')
    else:
        os.mkdir(os.path.join(root,'md5_error_sample'))
    for sample in unexpect_sample:
        shutil.move(os.path.join(root,sample),os.path.join(root,'unexpect_sample',sample))

    file_list = []
    for dir_ in run_sample:
        this_files = os.listdir(os.path.join(root,dir_))
        right_files = ['MD5.txt', f'{dir_}_1.clean.fq.gz', f'{dir_}_2.clean.fq.gz']
        if len(this_files) == 3 and sorted(this_files) == sorted(right_files):
            file_list += [os.path.join(root,dir_,file) for file in this_files if file.endswith('.clean.fq.gz')]
        else:
            file_list += change_name(os.path.join(root,dir_))

    filtered_df = all_sample_df[all_sample_df['Sample_ID'].isin(run_sample)]
    download_df = filtered_df.drop_duplicates(subset='Download_link')
    download_dict = dict(zip(download_df['Ref'],download_df['Download_link']))
    # Download reference genome
    download_ref_fa_file(download_dict, config["root"] + "/" + config["folder"]["index"])
    # Check reference genome md5
    ref_fa_file_md5_check(download_dict, config["root"] + "/" + config["folder"]["index"])
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
        "cd `dirname {input}` && md5sum -c {input} > {output} || true "


rule check_run:
    priority: 99
    input:
        file=check_run_sample(),
        md5=expand(config["root"] + "/logs/md5_info/{sample}_check_md5.log",sample=get_run_sample())
    output:
        "logs/01_check_run.log"
    message:
        "01: Checking which samples need processing ------------------------------------------"
    run:
        for i in input.md5:
            dir_name = os.path.basename(i).replace('_check_md5.log','')
            with open(i,'r') as f:
                if 'FAILED' in f.read() or '失败' in f.read():
                    try:
                        shutil.move(os.path.join(config["root"],config["folder"]["data"],dir_name),
                            os.path.join(config["root"],config["folder"]["data"],'md5_error_sample',dir_name))
                    except FileNotFoundError:
                        print(f'Dir {dir_name} not found, maybe it has been moved.')
        un_exp_sample = os.listdir(os.path.join(config["root"],config["folder"]["data"],'unexpect_sample'))
        md5_err_sample = os.listdir(os.path.join(config["root"],config["folder"]["data"],'md5_error_sample'))
        with open('logs/01_check_run.log','w') as f:
            f.write(f'Run sample num: {len(get_run_sample())}\n')
            f.write(f'Run samples:\n')
            f.write("\n".join(get_run_sample()) + "\n\n")
            f.write(f'Unexpected sample num: {len(un_exp_sample)}\n')
            f.write(f'Unexpected samples: {un_exp_sample}\n\n')
            f.write(f'MD5 error sample num: {len(md5_err_sample)}\n')
            f.write(f'MD5 error samples: {md5_err_sample}\n\n')
