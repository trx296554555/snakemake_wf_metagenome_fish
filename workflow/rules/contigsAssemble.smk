import os
import pandas as pd

rule assemble_contigs:
    input:
        fq1=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_1.fq",
        fq2=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_2.fq"
    output:
        contigs=config["root"] + "/" + config["folder"]["assemble_contigs"] + "/{sample}/{sample}.contigs.fa"
    message:
        "09: Assemble contigs using megahit -------------------------"
    threads:
        24
    params:
        opt=config["root"] + "/" + config["folder"]["assemble_contigs"] + "/{sample}/{sample}_megahit"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "megahit.yaml"
    log:
        config["root"] + "/" + config["folder"]["assemble_contigs"] + "/{sample}/{sample}_run.log"
    shell:
        """
        rm -rf {params.opt}
        megahit -1 {input.fq1} -2 {input.fq2} -m 0.95 --min-contig-len 300 \
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
    output:
        assemble_contigs_report=config["root"] + "/" + config["folder"]["reports"] + "/05_assemble_contigs.report"
    run:
        res_dict = {}
        for l_file in input.megahit_logs:
            sample = os.path.basename(l_file).split("_")[0]
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