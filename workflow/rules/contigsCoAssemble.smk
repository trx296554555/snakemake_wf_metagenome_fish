import os
import pandas as pd


co_assemble_file = config['root'] + '/workflow/config/co_assemble_list.csv'


def get_co_assemble_fq(wildcards):
    co_assemble_df = pd.read_csv(co_assemble_file)
    sample_df = co_assemble_df[co_assemble_df['Groups'] == wildcards.item]
    sample_list = sample_df['Samples'].tolist()[0].strip("'").split(',')
    res_dict = {'co_fq1':[], 'co_fq2':[]}
    for sample in sample_list:
        res_dict['co_fq1'].append(f'{config["root"]}/{config["folder"]["rm_host"]}/{sample}/{sample}_paired_1.fq')
        res_dict['co_fq2'].append(f'{config["root"]}/{config["folder"]["rm_host"]}/{sample}/{sample}_paired_2.fq')
    return res_dict


if config['co_assemble']['flag'] and os.path.exists(co_assemble_file):
    rule co_assemble_contigs:
        input:
            unpack(get_co_assemble_fq),
            co_assemble_item=config["root"] + "/logs/co_assemble/{item}",
        output:
            contigs=config["root"] + "/" + config["folder"]["assemble_contigs"] + "/{item}/co_{item}.contigs.fa"
        message:
            "09: Assemble contigs using megahit -------------------------"
        threads:
            24
        params:
            opt=config["root"] + "/" + config["folder"]["assemble_contigs"] + "/{item}/{item}_megahit"
        conda:
            config["root"] + "/" + config["envs"] + "/" + "megahit.yaml"
        benchmark:
            config["root"] + "/benchmark/" + config["folder"]["assemble_contigs"] + "/co_{item}.log"
        log:
            config["root"] + "/" + config["folder"]["assemble_contigs"] + "/{item}/co_{item}_run.log"
        shell:
            """
            rm -rf {params.opt}
            fq1_str=$(echo {input.co_fq1} | sed 's/ /,/g')
            fq2_str=$(echo {input.co_fq2} | sed 's/ /,/g')
            megahit -1 $fq1_str -2 $fq2_str -m 0.95 --min-contig-len 1000 \
            --k-list 79,99,119,139 -t {threads} \
            --out-dir {params.opt} --out-prefix {wildcards.item} > {log} 2>&1
            # rename contigs, replace space with underscore
            sed 's/ /_/g' {params.opt}/{wildcards.item}.contigs.fa > {output.contigs}
            """

