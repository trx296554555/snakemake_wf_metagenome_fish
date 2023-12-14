import os
import pandas as pd

rule generate_kraken2_report:
    input:
        fq1=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_1.fq",
        fq2=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_2.fq",
    output:
        kraken2_report=config["root"] + "/" + config["folder"]["reads_classify"] + "/{sample}.kraken",
        kraken2_out=temp(config["root"] + "/" + config["folder"]["reads_classify"] + "/{sample}.out")
    message:
        "05: Execute Kraken2 for assigning taxonomic labels to DNA sequences ------------------------------------------"
    threads:
        24
    params:
        db=config["db_root"] + "/" + config["db"]["kraken2"]
    conda:
        config["root"] + "/" + config["envs"] + "/" + "kraken2.yaml"
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["reads_classify"] + "/{sample}.log"
    log:
        config["root"] + "/" + config["folder"]["reads_classify"] + "/{sample}.log"
    shell:
        """
        kraken2 --paired -db {params.db} --threads {threads} \
        --report {output.kraken2_report} --out {output.kraken2_out} {input.fq1} {input.fq2} > {log} 2>&1
        """


rule generate_bracken_report:
    input:
        kraken2_report=config["root"] + "/" + config["folder"]["reads_classify"] + "/{sample}.kraken"
    output:
        bracken_report=config["root"] + "/" + config["folder"]["reads_classify"] + "/{sample}.bracken",
        bracken_out=config["root"] + "/" + config["folder"]["reads_classify"] + "/{sample}.bracken.out",
        trx_report=config["root"] + "/" + config["folder"]["reads_classify"] + "/{sample}.trx"
    message:
        "06: Reestimate Kraken2 results at the species level using Bracken ------------------------------------------"
    threads:
        1
    params:
        db=config["db_root"] + "/" + config["db"]["kraken2"],
        reads_len=config["meta"]["reads_length"],
        tmp_log=config["root"] + "/" + config["folder"]["reads_classify"] + "/{sample}.bracken.log",
        final_log=config["root"] + "/" + config["folder"]["reads_classify"] + "/{sample}.log",
        bracken2trx=config["root"] + "/workflow/scripts/bracken_to_trx.py"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "kraken2.yaml"
    shell:
        """
        sed -i '1i ------- kraken2 log --------' {params.final_log}
        bracken -d {params.db} -i {input.kraken2_report} -o {output.bracken_out} \
        -w {output.bracken_report} -r {params.reads_len} > {params.tmp_log} 2>&1
        echo '------- bracken log --------' | cat - {params.tmp_log} >> {params.final_log} && rm -f {params.tmp_log}
        python {params.bracken2trx} -r {output.bracken_report} -o {output.trx_report}
        """


rule merge_bracken_res:
    input:
        bracken_res=expand(config["root"] + "/" + config["folder"]["reads_classify"] + "/{sample}.trx", sample=get_run_sample()),
    output:
        merge_bracken=config["root"] + "/" + config["folder"]["reads_classify"] + "/all_bracken_res.tsv",
    message:
        "07: Merge Bracken reports ----------------------------------------------------------------------------------"
    run:
        merged_df = pd.DataFrame()
        for file_name in input.bracken_res:
            df = pd.read_csv(file_name, sep='\t', index_col=['TaxID', 'Taxon'])
            if merged_df.empty:
                merged_df = df
            else:
                merged_df = pd.concat([merged_df, df], axis=1, join='outer', sort=False)
        merged_df = merged_df.fillna(0)
        merged_df = merged_df.astype(int)
        merged_df = merged_df.reindex(sorted(merged_df.columns),axis=1)
        merged_df = merged_df.reset_index()
        merged_df = merged_df.sort_values(by=['Taxon'])
        merged_df.to_csv(output.merge_bracken, sep='\t', index=False)


rule report_reads_classify:
    input:
        merge_bracken = config["root"] + "/" + config["folder"]["reads_classify"] + "/all_bracken_res.tsv",
        kraken2_report=expand(config["root"] + "/" + config["folder"]["reads_classify"] + "/{sample}.kraken", sample=get_run_sample()),
    output:
        reads_classify_report=config["root"] + "/" + config["folder"]["reports"] + "/02_reads_classify.report"
    run:
        report_list = ['Eukaryota', 'Fungi', 'Bacteria', 'Viruses', 'Archaea', 'unclassified', 'Homo sapiens']
        res_dict = {}
        for k_file in input.kraken2_report:
            sample = os.path.basename(k_file).split('.')[0]
            res_dict[sample] = {}
            with open(k_file ,'r') as ipt:
                for line in ipt:
                    line = line.strip().split('\t')
                    counts = line[1]
                    taxon = line[-1].strip()
                    if taxon in report_list:
                        res_dict[sample][taxon] = counts
        table = pd.DataFrame(res_dict).transpose().astype(int)
        table['unclassified(%)'] = table['unclassified'].div(table.sum(axis=1)).multiply(100)
        table['human(%)'] = table['Homo sapiens'].div(table.sum(axis=1)).multiply(100)
        table['known_microbiome(%)'] = 100 - table['unclassified(%)'] - table['human(%)']
        table = table.sort_index()
        table.to_csv(output.reads_classify_report, sep="\t")