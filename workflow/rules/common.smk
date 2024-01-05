import os
import pandas as pd
import numpy as np


if os.path.exists(f'{config["root"]}/{config["meta"]["sampleList"]}'):
    all_sample_df = pd.read_csv(f'{config["root"]}/{config["meta"]["sampleList"]}')
else:
    raise Exception(f'No {config["meta"]["sampleList"]} in {config["root"]}')
    exit(1)
if os.path.exists('logs'):
    pass
else:
    os.mkdir('logs')


def set_co_assemble():
    # if coAssemble is needed, then prepare the input for coAssemble
    col_list = list(config['co_assemble']['col'])
    if all(col in all_sample_df.columns for col in col_list):
        pass
    else:
        raise Exception(f'co_assemble columns not all found in {config["meta"]["sampleList"]}')

    melted_df = all_sample_df.melt(id_vars=['Sample_ID'],value_vars=col_list,var_name='Cate',value_name='Groups')
    result_df = melted_df.groupby(['Cate', 'Groups'])['Sample_ID'].apply(lambda x: ','.join(x)).reset_index()
    result_df = result_df.rename(columns={'Sample_ID': 'Samples'})
    result_df.to_csv(config['root'] + '/workflow/config/co_assemble_list.csv',index=False)

    for i in result_df.index:
        if all(sample in os.listdir(f"{config['root']}/{config['folder']['data']}") for sample in
               result_df.loc[i, 'Samples'].split(',')):
            info = f'{result_df.loc[i, "Cate"]}:The samples required for grouping {result_df.loc[i, "Groups"]}' \
                   f' are all present in the current Server. The co_assemble will be performed.'
            print(info)
            with open(config["root"] + f'/logs/co_assemble/{result_df.loc[i, "Groups"]}','w') as opt:
                opt.write(result_df.loc[i, "Samples"] + '\n')
        else:
            raise Exception(f'{result_df.loc[i, "Cate"]}:The samples required for grouping {result_df.loc[i, "Groups"]}'
                            f' are not all present in the current Server. The co_assemble will not be performed.')


def get_all_sample():
    return all_sample_df['Sample_ID'].dropna().tolist()


def get_run_sample():
    return all_sample_df['Use_Sample'].dropna().tolist()


def get_run_individual():
    run_sample = get_run_sample()
    run_individual = []
    for i in run_sample:
        if (individual := i.replace('C','').replace('P','')) in run_individual:
            pass
        else:
            run_individual.append(individual)
    return run_individual


def get_co_item():
    tmp_co_assemble_path = config["root"] + f'/logs/co_assemble'
    if os.path.exists(tmp_co_assemble_path):
        return os.listdir(tmp_co_assemble_path)
    else:
        os.mkdir(tmp_co_assemble_path)
        set_co_assemble()
        return os.listdir(tmp_co_assemble_path)


def get_report():
    reads_anno_enable = config["reads_anno"]["huamnn3_enable"] or config["reads_anno"]["rgi_enable"] or \
                        config["reads_anno"]["dbcan_enable"] or config["reads_anno"]["vfdb_enable"]
    report_list = [
        config["root"] + "/" + config["folder"]["logs"] + "/01_check_run.log",
        config["root"] + "/" + config["folder"]["reports"] + "/01_rm_host.report",
        config["root"] + "/" + config["folder"]["reports"] + "/02_reads_classify.report" if config["reads_classify"][
            "enable"] else "",
        config["root"] + "/" + config["folder"]["reports"] + "/03_get_mitogenome.report" if config["get_mitogenome"][
            "enable"] else "",
        config["root"] + "/" + config["folder"]["reports"] + "/04_assemble_contigs.report" if
        config["assemble_contigs"]["enable"] else "",
        config["root"] + "/" + config["folder"]["reports"] + "/05_gene_prediction.report" if config["gene_prediction"][
            "enable"] else "",
        config["root"] + "/" + config["folder"]["reports"] + "/06_reads_annotation.report" if reads_anno_enable else "",
        config["root"] + "/" + config["folder"]["bin_refine"] + "/all_bins/gather.done" if config["contigs_binning"][
            "enable"] else "",
        config["root"] + "/" + config["folder"]["reports"] + "/07_contigs_binning.report" if config["bins_dereplicate"][
            "enable"] else "",
        # TODO 需要一个清理中间文件的rule 在工作流结束后运行
    ]
    report_list = [i for i in report_list if i != ""]
    return report_list
