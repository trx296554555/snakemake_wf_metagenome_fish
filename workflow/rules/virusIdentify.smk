rule identify_virsorter2:
    input:
        fa=config["root"] + "/" + config["folder"]["assemble_contigs"] + "/{sample}/{sample}.contigs.fa",
    output:
        result_fa=config["root"] + "/" + config["folder"][
            "virus_identify"] + "/{sample}/virsorter2/{sample}.viral_sequences.fa"
    message:
        "01: Use virsorter2 to identify viruses in metagenomic contigs --------------------------"
    threads:
        8
    params:
        db=config["db_root"] + "/" + config["virus_db"]["virsorter2"],
        raw=config["root"] + "/" + config["folder"][
            "virus_identify"] + "/{sample}/virsorter2/final-viral-combined.fa",
        out_dir=config["root"] + "/" + config["folder"]["virus_identify"] + "/{sample}/virsorter2",
        min_length=config["virus_identify"]["min_virus_len"],
    log:
        config["root"] + "/" + config["folder"]["virus_identify"] + "/{sample}/virsorter2/virsorter2.log"
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["virus_identify"] + "/{sample}_virsorter2.log"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "virsorter2.yaml"
    shell:
        """
        virsorter run -j {threads} -i {input.fa} --db-dir {params.db} -w {params.out_dir} --min-length {params.min_length} all > {log} 2>&1 || touch {params.raw}
        ln -sf {params.raw} {output.result_fa}
        """


rule identify_deepvirfinder:
    input:
        fa=config["root"] + "/" + config["folder"]["assemble_contigs"] + "/{sample}/{sample}.contigs.fa"
    output:
        score_file=config["root"] + "/" + config["folder"][
            "virus_identify"] + "/{sample}/deepvirfinder/{sample}.score.txt",
        pass_list=config["root"] + "/" + config["folder"][
            "virus_identify"] + "/{sample}/deepvirfinder/{sample}.viral_contigs.pass",
        result_fa=config["root"] + "/" + config["folder"][
            "virus_identify"] + "/{sample}/deepvirfinder/{sample}.viral_sequences.fa"
    message:
        "01: Use deepvirfinder to identify viruses in metagenomic contigs --------------------------"
    threads:
        8
    params:
        out_dir=config["root"] + "/" + config["folder"]["virus_identify"] + "/{sample}/deepvirfinder",
        db=config["db_root"] + "/" + config["virus_db"]["deepvirfinder"] + '/models',
        script=config["root"] + "/workflow/scripts/dvf.py",
        min_length=config["virus_identify"]["min_virus_len"],
    log:
        config["root"] + "/" + config["folder"]["virus_identify"] + "/{sample}/deepvirfinder/deepvirfinder.log"
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["virus_identify"] + "/{sample}_deepvirfinder.log"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "deepvirfinder.yaml"
    shell:
        """
        python {params.script} -i {input} -m {params.db} -c {threads} -o {params.out_dir} -l {params.min_length} > {log} 2>&1
        mv {params.out_dir}/{wildcards.sample}.contigs.fa_gt{params.min_length}bp_dvfpred.txt {output.score_file}
        awk -F '\\t' '{{if($3 > 0.9 && $4 < 0.01)print$1}}' {output.score_file} > {output.pass_list}
        seqkit grep -f {output.pass_list} {input.fa} --quiet > {output.result_fa}
        """


rule identify_vibrant:
    input:
        fa=config["root"] + "/" + config["folder"]["assemble_contigs"] + "/{sample}/{sample}.contigs.fa"
    output:
        result_fa=config["root"] + "/" + config["folder"][
            "virus_identify"] + "/{sample}/vibrant/{sample}.viral_sequences.fa",
        done=touch(
            config["root"] + "/" + config["folder"]["virus_identify"] + "/{sample}/vibrant/{sample}.vibrant.done")
    message:
        "01: Use vibrant to identify viruses in metagenomic contigs --------------------------"
    threads:
        8
    params:
        tmp_dir=config["root"] + "/" + config["folder"]["virus_identify"] + "/{sample}/vibrant/tmp",
        useful_file=config["root"] + "/" + config["folder"][
            "virus_identify"] + "/{sample}/vibrant/tmp/VIBRANT_{sample}.contigs/VIBRANT_phages_{sample}.contigs/{sample}.contigs.phages_combined.fna",
        db=config["db_root"] + "/" + config["virus_db"]["vibrant"],
        min_length=config["virus_identify"]["min_virus_len"],
    log:
        config["root"] + "/" + config["folder"]["virus_identify"] + "/{sample}/vibrant/vibrant.log"
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["virus_identify"] + "/{sample}_vibrant.log"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "vibrant.yaml"
    shell:
        """
        rm -rf {params.tmp_dir} && mkdir -p {params.tmp_dir}
        VIBRANT_run.py -t {threads} -d {params.db}/databases -m {params.db}/files -i {input.fa} -folder {params.tmp_dir} -l {params.min_length} -no_plot > {log} 2>&1 || true
        if [ ! -f {params.useful_file} ]; then
            touch {output.result_fa}
            touch $(dirname {params.tmp_dir})/{wildcards.sample}.failed
        else
            mv {params.useful_file} {output.result_fa}
        fi
        rm -rf {params.tmp_dir}
        """


def collect_virus_from_diff_tools(wildcards):
    sample_item = str(wildcards)
    virus_samples = []
    for tool in config["virus_identify"]["tools"]:
        virus_samples.append(config["root"] + "/" + config["folder"][
            "virus_identify"] + "/" + sample_item + "/" + tool + "/" + sample_item + ".viral_sequences.fa")
    return virus_samples


rule dereplicate_virus_tools:
    input:
        collect_virus_from_diff_tools
    output:
        tmp_all_fa=temp(
            config["root"] + "/" + config["folder"]["virus_identify"] + "/{sample}/no_redundant/{sample}.tmp.fa"),
        all_fa=config["root"] + "/" + config["folder"]["virus_identify"] + "/{sample}/no_redundant/{sample}.all.fa",
        tmp_gene=temp(
            config["root"] + "/" + config["folder"]["virus_identify"] + "/{sample}/no_redundant/{sample}.gene.tmp.fa"),
        res_fa=config["root"] + "/" + config["folder"][
            "virus_identify"] + "/{sample}/no_redundant/{sample}.prune.fa",
    message:
        "02: Use cd-hit to prune redundant sequences from different tools --------------------------"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "prediction.yaml"
    threads:
        4
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["virus_identify"] + "/{sample}_cdhit.log"
    log:
        config["root"] + "/" + config["folder"]["virus_identify"] + "/{sample}/no_redundant/{sample}_cdhit.log"
    shell:
        """
        cat {input} > {output.tmp_all_fa}
        seqkit rmdup -s --quiet -w 0 {output.tmp_all_fa} > {output.all_fa}
        cd-hit-est -T {threads} -M 0 -i {output.all_fa} -o {output.tmp_gene} -c 0.99 -aS 0.9 -g 1 -sc 1 -sf 1 -d 0 > /dev/null 2>&1
        seqkit replace --quiet -w 0 -p '.*' -r {wildcards.sample}={{nr}} {output.tmp_gene} > {output.res_fa}
        seqkit stats --quiet -T $(dirname {output.res_fa})/*.fa > {log}
        """


rule gather_all_virus:
    input:
        sample_virus=expand(config["root"] + "/" + config["folder"][
            "virus_identify"] + "/{sample}/no_redundant/{sample}.prune.fa",sample=get_run_sample()),
        co_assemble_virus=expand(config["root"] + "/" + config["folder"][
            "virus_identify"] + "/{item}/no_redundant/{item}.prune.fa",item=get_co_item()) if config["co_assemble"][
            "enable"] else [],
    output:
        all_virus=directory(config["root"] + "/" + config["folder"]["virus_dereplication"] + "/all_virus"),
        gather_done=touch(config["root"] + "/" + config["folder"]["virus_dereplication"] + "/gather.done")
    message:
        "Gather all virus sequences from different samples, The purpose of adding this step is to distribute computing"
    shell:
        """
        mkdir -p {output.all_virus}
        for eachseq in {input.sample_virus}; do
            ln -sf $eachseq {output.all_virus}
        done
        for eachseq in {input.co_assemble_virus}; do
            ln -sf $eachseq {output.all_virus}
        done
        """


rule dereplicate_virus_samples:
    input:
        config["root"] + "/" + config["folder"]["virus_dereplication"] + "/all_virus"
    output:
        all_virus_fa=config["root"] + "/" + config["folder"]["virus_dereplication"] + "/all_virus_redundant.fa",
        all_no_redundant_fa=config["root"] + "/" + config["folder"][
            "virus_dereplication"] + "/all_virus_no_redundant.fa",
        done=touch(config["root"] + "/" + config["folder"]["virus_dereplication"] + "/dereplicate_virus.done")
    message:
        "02: Dereplicate all virus sequences from different samples  --------------------------"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "prediction.yaml"
    threads:
        72
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["virus_dereplication"] + "/dereplicate_virus.log"
    log:
        config["root"] + "/" + config["folder"]["virus_dereplication"] + "/dereplicate_virus.log"
    shell:
        """
        cat {input}/*.fa > {output.all_virus_fa}
        cd-hit-est -T {threads} -M 0 -i {output.all_virus_fa} -o {output.all_no_redundant_fa} -c 0.99 -aS 0.9 -g 1 -sc 1 -sf 1 -d 0 > /dev/null 2>&1
        seqkit stats --quiet -T $(dirname {output.all_no_redundant_fa})/*.fa > {log}
        """


rule checkv_virus_quality:
    input:
        all_no_redundant_fa=config["root"] + "/" + config["folder"][
            "virus_dereplication"] + "/all_virus_no_redundant.fa",
    output:
        quality_result=config["root"] + "/" + config["folder"]["virus_checkv"] + "/quality_summary.tsv",
        virus_fa=config["root"] + "/" + config["folder"]["virus_checkv"] + "/viruses.fna",
        pro_virus_fa=config["root"] + "/" + config["folder"]["virus_checkv"] + "/proviruses.fna",
        done=touch(config["root"] + "/" + config["folder"]["virus_checkv"] + "/checkv.done")
    message:
        "03: Use checkv to evaluate the quality of viral sequences --------------------------"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "checkv.yaml"
    threads:
        96
    params:
        out_dir=config["root"] + "/" + config["folder"]["virus_checkv"],
        db=config["db_root"] + "/" + config["virus_db"]["checkv"],
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["virus_checkv"] + "/checkv.log"
    log:
        config["root"] + "/" + config["folder"]["virus_checkv"] + "/checkv.log"
    shell:
        """
        checkv end_to_end --remove_tmp -t {threads} -d {params.db} {input.all_no_redundant_fa} {params.out_dir} > {log} 2>&1
        """


rule checkv_filter:
    input:
        quality_result=config["root"] + "/" + config["folder"]["virus_checkv"] + "/quality_summary.tsv",
        virus_fa=config["root"] + "/" + config["folder"]["virus_checkv"] + "/viruses.fna",
        pro_virus_fa=config["root"] + "/" + config["folder"]["virus_checkv"] + "/proviruses.fna",
    output:
        pass_list=config["root"] + "/" + config["folder"]["virus_checkv"] + "/filter_pass.txt",
        all_fa=temp(config["root"] + "/" + config["folder"]["virus_checkv"] + "/all.tmp.fa"),
        pass_fa=config["root"] + "/" + config["folder"]["virus_checkv"] + "/pass.fa",
        done=touch(config["root"] + "/" + config["folder"]["virus_checkv"] + "/filter.done"),
    message:
        """
        Filter criteria: Completeness (column 10) is not equal to NA and greater than or equal to 50,
        the number of viral genes (column 6) is greater than the number of host genes (column 7), 
        contig length (column 2) is greater than or equal to 5000, contamination level (column 12) is 
        less than or equal to 5, the occurrence count of viral sequences in contigs (kmer_freq, column 13) is
        less than 1, and the length of contigs does not exceed 1.5 times the expected genome length
        (contig > 1.5x not present in column 14).
        过滤标准：完整度（第10列）不等于NA且大于等于50，病毒gene数（第6列）大于宿主基因数（第7列），contig长度（第二列）大于等于5000，
        污染度（第12列）小于等于5，病毒序列在contig中出现次数kmer_freq（第13列）小于1，且contig长度不超过预期genome长度的1.5倍
        (contig >1.5x不在第14列中).
        """
    conda:
        config["root"] + "/" + config["envs"] + "/" + "tools.yaml"
    params:
        completeness=50,contamination=5,contig_len=config["virus_identify"]["min_virus_len"],kmer_freq=1,
    shell:
        """
        awk -F '\\t' 'NR>1 && $10 != "NA" && $6 > $7 && $2 >= {params.contig_len} && $10 >= {params.completeness} && $12 <= {params.contamination} && $13 <= {params.kmer_freq} && $14 !~ /1\.5x/ {{print $1}}' {input.quality_result} > {output.pass_list}
        cat {input.virus_fa} {input.pro_virus_fa} > {output.all_fa}
        seqkit grep -n -f {output.pass_list} {output.all_fa} --quiet -w 0 > {output.pass_fa}
        """


rule all_vs_all_blastn:
    input:
        pass_fa=config["root"] + "/" + config["folder"]["virus_checkv"] + "/pass.fa",
    output:
        blastn_db=directory(config["root"] + "/" + config["folder"]["virus_cluster"] + "/blast_db"),
        blastn_result=config["root"] + "/" + config["folder"]["virus_cluster"] + "/all_vs_all_blast.tsv",
    message:
        "04: Use blastn to perform pairwise comparisons between all viral sequences --------------------------"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "checkv.yaml"
    threads:
        96
    params:
        prefix="all_viruses",max_target_seqs=25000,perc_identity=90,
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["virus_cluster"] + "/all_vs_all_blast.log"
    log:
        config["root"] + "/" + config["folder"]["virus_cluster"] + "/makeblastdb.log"
    shell:
        """
        makeblastdb -in {input.pass_fa} -out {output.blastn_db}/{params.prefix} -dbtype nucl > {log} 2>&1
        blastn -num_threads {threads} -query {input.pass_fa} -db {output.blastn_db}/{params.prefix} -out {output.blastn_result} -outfmt '6 std qlen slen' -max_target_seqs {params.max_target_seqs} -perc_identity {params.perc_identity}
        """


rule get_votus:
    input:
        pass_fa=config["root"] + "/" + config["folder"]["virus_checkv"] + "/pass.fa",
        blastn_result=config["root"] + "/" + config["folder"]["virus_cluster"] + "/all_vs_all_blast.tsv"
    output:
        ani_result=config["root"] + "/" + config["folder"]["virus_cluster"] + "/ani.tsv",
        cluster_result=config["root"] + "/" + config["folder"]["virus_cluster"] + "/clusters.tsv",
        final_votus=config["root"] + "/" + config["folder"]["virus_cluster"] + "/final_virus_otu.fa",
    params:
        cal_script=config["root"] + "/workflow/scripts/anical.py",
        cluster_script=config["root"] + "/workflow/scripts/aniclust.py",
        get_seq_script=config["root"] + "/workflow/scripts/get_longest_from_cluster.py",
        min_ani=95,min_qcov=0,min_tcov=85,
    conda:
        config["root"] + "/" + config["envs"] + "/" + "checkv.yaml"
    log:
        config["root"] + "/" + config["folder"]["virus_cluster"] + "/votus.log"
    shell:
        """
        python {params.cal_script} -i {input.blastn_result} -o {output.ani_result} > {log} 2>&1
        python {params.cluster_script} --fna {input.pass_fa} --ani {output.ani_result} --out {output.cluster_result} --min_ani {params.min_ani} --min_qcov {params.min_qcov} --min_tcov {params.min_tcov} >> {log} 2>&1
        python {params.get_seq_script} -i {input.pass_fa} -c {output.cluster_result} -o {output.final_votus}
        """
