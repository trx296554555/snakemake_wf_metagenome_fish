import os
import pandas as pd

def get_co_assemble_fq(wildcards):
    co_assemble_file = config['root'] + '/workflow/config/co_assemble_list.csv'
    if wildcards.sample in get_co_item():
        co_assemble_df = pd.read_csv(co_assemble_file)
        sample_df = co_assemble_df[co_assemble_df['Groups'] == wildcards.sample]
        sample_list = sample_df['Samples'].tolist()[0].strip("'").split(',')
        res_dict = {'fq1':[], 'fq2':[]}
        for sample in sample_list:
            res_dict['fq1'].append(f'{config["root"]}/{config["folder"]["rm_host"]}/{sample}/{sample}_paired_1.fq')
            res_dict['fq2'].append(f'{config["root"]}/{config["folder"]["rm_host"]}/{sample}/{sample}_paired_2.fq')
        return res_dict
    else:
        return {'fq1':f'{config["root"]}/{config["folder"]["rm_host"]}/{wildcards.sample}/{wildcards.sample}_paired_1.fq',
                'fq2':f'{config["root"]}/{config["folder"]["rm_host"]}/{wildcards.sample}/{wildcards.sample}_paired_2.fq'}


def get_assemble_len(wildcards):
    if wildcards.sample in get_co_item():
        return config["co_assemble"]["min_ctg_len"]
    else:
        return "300"


rule assemble_contigs:
    input:
        unpack(get_co_assemble_fq)
    output:
        contigs=config["root"] + "/" + config["folder"]["assemble_contigs"] + "/{sample}/{sample}.contigs.fa"
    message:
        "09: Assemble contigs using megahit -------------------------"
    threads:
        24
    params:
        opt=config["root"] + "/" + config["folder"]["assemble_contigs"] + "/{sample}/{sample}_megahit",
        min_len=get_assemble_len
    conda:
        config["root"] + "/" + config["envs"] + "/" + "megahit.yaml"
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["assemble_contigs"] + "/{sample}.log"
    log:
        config["root"] + "/" + config["folder"]["assemble_contigs"] + "/{sample}/{sample}_run.log"
    shell:
        """
        rm -rf {params.opt}
        fq1_str=$(echo {input.fq1} | sed 's/ /,/g')
        fq2_str=$(echo {input.fq2} | sed 's/ /,/g')
        megahit -1 $fq1_str -2 $fq2_str -m 0.95 --min-contig-len {params.min_len} \
        --k-list 79,99,119,139 -t {threads} \
        --out-dir {params.opt} --out-prefix {wildcards.sample} > {log} 2>&1
        # rename contigs, replace space with underscore
        sed 's/ /_/g' {params.opt}/{wildcards.sample}.contigs.fa > {output.contigs}
        """


rule report_assemble_contigs:
    input:
        megahit_logs=expand(
            config["root"] + "/" + config["folder"][
                "assemble_contigs"] + "/{sample}/{sample}_run.log",sample=get_run_sample()),
        co_megahit_logs=expand(
            config["root"] + "/" + config["folder"][
                "assemble_contigs"] + "/{item}/{item}_run.log",item=get_co_item()) if config["co_assemble"]["enable"] else []
    output:
        assemble_contigs_report=config["root"] + "/" + config["folder"]["reports"] + "/04_assemble_contigs.report"
    run:
        log_list = input.megahit_logs + input.co_megahit_logs
        res_dict = {}
        for l_file in log_list:
            sample = os.path.basename(l_file).replace('_run.log','')
            res_dict[sample] = {}
            with open(l_file) as ipt:
                lines = ipt.readlines()
            info = lines[-2].split('-')[-1]
            use_info = info.split(',')
            res_dict[sample]['contigs_num'] = use_info[0].strip().split(' ')[0]
            for i in use_info[1:]:
                res_dict[sample][i.strip().split(' ')[0]] = i.strip().split(' ')[1]
        table = pd.DataFrame(res_dict).transpose().astype(int)
        table = table.sort_index()
        table.to_csv(output.assemble_contigs_report,sep="\t")