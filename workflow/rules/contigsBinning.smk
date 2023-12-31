rule build_contigs_index:
    input:
        contigs=config["root"] + "/" + config["folder"]["assemble_contigs"] + "/{item}/{item}.contigs.fa",
    output:
        index=directory(config["root"] + "/" + config["folder"]["contigs_binning"] + "/{item}/{item}_contig_index"),
    message:
        "13.1: Indexing and alignment of contigs to reads --------------------------"
    threads:
        24
    conda:
        config["root"] + "/" + config["envs"] + "/" + "tools.yaml"
    shell:
        """
        mkdir -p {output.index}
        bowtie2-build -f {input.contigs} {output.index}/{wildcards.item} --threads {threads} --offrate 1 > /dev/null 2>&1
        """


rule align_reads_contigs:
    input:
        index=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{item}/{item}_contig_index",
        fq1=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_1.fq",
        fq2=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_2.fq"
    output:
        sam=temp(config["root"] + "/" + config["folder"]["contigs_binning"] + "/{item}/mapping/{sample}.sam"),
        bam=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{item}/mapping/{sample}.bam",
    message:
        "13.2: Indexing and alignment of contigs to reads --------------------------"
    threads:
        12
    conda:
        config["root"] + "/" + config["envs"] + "/" + "tools.yaml"
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["contigs_binning"] + "/mapping/{item}_{sample}.log"
    log:
        config["root"] + "/" + config["folder"]["contigs_binning"] + "/{item}/mapping/{sample}.mapping.log"
    shell:
        """
        bowtie2 -x {input.index}/{wildcards.item} -1 {input.fq1} -2 {input.fq2} -S {output.sam} --threads {threads} > {log} 2>&1
        samtools view -bS -@ {threads} {output.sam} -o {output.bam} 2>&1
        samtools sort -@ {threads} {output.bam} -o {output.bam} 2>&1
        samtools index {output.bam} 2>&1
        """


def get_sorted_sample_bam(wildcards):
    if wildcards.item in get_co_item():
        co_assemble_df = pd.read_csv(config['root'] + '/workflow/config/co_assemble_list.csv')
        sample_df = co_assemble_df[co_assemble_df['Groups'] == wildcards.item]
        sample_list = sample_df['Samples'].tolist()[0].strip("'").split(',')
        return expand(config["root"] + "/" + config["folder"][
            "contigs_binning"] + "/" + wildcards.item + "/mapping/{sample}.bam",sample=sample_list)
    else:
        return config["root"] + "/" + config["folder"][
            "contigs_binning"] + f"/{wildcards.item}/mapping/{wildcards.item}.bam"


# 对于co_assemble的contigs的比对结果，合并bam文件
rule merge_bam:
    input:
        all_bam=get_sorted_sample_bam
    output:
        merge_bam=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{item}/{item}.bam"
    message:
        "13.3: Merge bam files --------------------------"
    threads:
        1
    conda:
        config["root"] + "/" + config["envs"] + "/" + "tools.yaml"
    log:
        config["root"] + "/" + config["folder"]["contigs_binning"] + "/{item}/merged_bam.log"
    shell:
        """
        samtools merge -@ {threads} -f {output.merge_bam} {input.all_bam} 2>&1 > {log}
        samtools index {output.merge_bam} 2>&1
        """


rule calculate_coverage:
    input:
        bam=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/{sample}.bam",
    output:
        depth=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/{sample}.depth",
        coverage=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/{sample}.coverage"
    message:
        "14: Generate coverage --------------------------"
    threads:
        1
    conda:
        config["root"] + "/" + config["envs"] + "/" + "binning.yaml"
    log:
        config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/coverage.log"
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam} > {log} 2>&1
        awk -v filename={wildcards.sample} 'BEGIN {{OFS="\t"; print "contig", "cov_mean_sample_" filename}} NR>1 {{print $1, $3 + $4}}' {output.depth} > {output.coverage}
        """


rule binning_metabat2:
    input:
        contigs=config["root"] + "/" + config["folder"]["assemble_contigs"] + "/{sample}/{sample}.contigs.fa",
        depth=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/{sample}.depth",
    output:
        bins=directory(config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/metabat2_bins")
    message:
        "15: Binning with metabat2 --------------------------"
    threads:
        12
    conda:
        config["root"] + "/" + config["envs"] + "/" + "binning.yaml"
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["contigs_binning"] + "/metabat2/{sample}.log"
    log:
        config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/metabat2.log"
    shell:
        """
        metabat2 -m 1500 -i {input.contigs} -a {input.depth} -o {output.bins}/{wildcards.sample}_metabat2 -t {threads} -v > {log} 2>&1
        """


rule binning_concoct:
    input:
        contigs=config["root"] + "/" + config["folder"]["assemble_contigs"] + "/{sample}/{sample}.contigs.fa",
        coverage=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/{sample}.coverage",
    output:
        opt=temp(directory(config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/concoct_opt")),
        bins=directory(config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/concoct_bins"),
    message:
        "15: Binning using concoct --------------------------"
    threads:
        12
    conda:
        config["root"] + "/" + config["envs"] + "/" + "concoct.yaml"
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["contigs_binning"] + "/concoct/{sample}.log"
    log:
        config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/concoct.log"
    shell:
        """
        mkdir -p {output.bins}
        concoct --composition_file {input.contigs} --coverage_file {input.coverage} -b {output.opt}/ --threads {threads} > {log} 2>&1
        merge_cutup_clustering.py {output.opt}/clustering_gt1000.csv > {output.opt}/clustering_merged.csv
        awk -F ',' -v string="{wildcards.sample}_concoct." 'NR==1 {{print; next}} {{OFS=","; $2 = string $2; print}}' {output.opt}/clustering_merged.csv > {output.opt}/clustering_merged_renamed.csv
        extract_fasta_bins.py {input.contigs} {output.opt}/clustering_merged_renamed.csv --output_path {output.bins} >> {log} 2>&1
        """
# cut_up_fasta.py {input.contigs} -c 10000 -o 0 --merge_last -b {output.bed} > {output.contigs_10K}
# concoct --composition_file {output.contigs_10K} --coverage_file {output.coverage} -b {output.opt}/ --threads {threads} > {log} 2>&1


rule binning_maxbin2:
    input:
        contigs=config["root"] + "/" + config["folder"]["assemble_contigs"] + "/{sample}/{sample}.contigs.fa",
        coverage=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/{sample}.coverage",
    output:
        bins=directory(config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/maxbin2_bins"),
    message:
        "15: Binning using maxbin2 --------------------------"
    threads:
        12
    conda:
        config["root"] + "/" + config["envs"] + "/" + "binning.yaml"
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["contigs_binning"] + "/maxbin2/{sample}.log"
    log:
        config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/maxbin2.log"
    shell:
        """
        mkdir -p {output.bins}
        run_MaxBin.pl -contig {input.contigs} -abund {input.coverage} -out {output.bins}/{wildcards.sample}_maxbin2 -thread {threads} > {log} 2>&1 || true
        mkdir -p {output.bins}
        """


checkpoint refine_bins_DAS_tool:
    input:
        contigs=config["root"] + "/" + config["folder"]["assemble_contigs"] + "/{sample}/{sample}.contigs.fa",
        metabat2_bins=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/metabat2_bins",
        concoct_bins=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/concoct_bins",
        proteins=config["root"] + "/" + config["folder"]["gene_prediction"] + "/{sample}/{sample}.protein.faa",
        maxbin2_bins=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/maxbin2_bins"
    output:
        bins=directory(config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/das_merged_bins")
    message:
        "16: Refine bins using DAS tool --------------------------"
    threads:
        12
    conda:
        config["root"] + "/" + config["envs"] + "/" + "dastool.yaml"
    params:
        metabat2_c2b=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/metabat2_contigs2bin.tsv",
        concoct_c2b=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/concoct_contigs2bin.tsv",
        maxbin2_c2b=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/maxbin2_contigs2bin.tsv"
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["contigs_binning"] + "/dastool/{sample}.log"
    log:
        config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/das_merged_bins.log"
    shell:
        """
        Fasta_to_Contig2Bin.sh -e fa -i {input.metabat2_bins} > {params.metabat2_c2b} 2>&1
        Fasta_to_Contig2Bin.sh -e fa -i {input.concoct_bins} > {params.concoct_c2b} 2>&1
        Fasta_to_Contig2Bin.sh -e fasta -i {input.maxbin2_bins} > {params.maxbin2_c2b} 2>&1
        mkdir -p {output.bins}
        DAS_Tool -i {params.metabat2_c2b},{params.concoct_c2b},{params.maxbin2_c2b} -l metabat2,concoct,maxbin2 \
        -c {input.contigs} -o {output.bins}/{wildcards.sample}.contigs.fa -t {threads} \
        -p {input.proteins} --write_bins --score_threshold 0 --write_bin_evals > {log} 2>&1 || true
        """


checkpoint metabinner:
    input:
        contigs=config["root"] + "/" + config["folder"]["assemble_contigs"] + "/{sample}/{sample}.contigs.fa",
        depth=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/{sample}.depth",
    output:
        bins=directory(config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/metabinner_bins")
    message:
        "15: Binning using metabinner --------------------------"
    threads:
        12
    conda:
        config["root"] + "/" + config["envs"] + "/" + "metabinner.yaml"
    params:
        tmp_dir=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/metabinner_tmp"
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["contigs_binning"] + "/metabinner/{sample}.log"
    log:
        config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/metabinner.log"
    shell:
        """
        mkdir -p {output.bins} {params.tmp_dir}
        ln -sf {input.contigs} {params.tmp_dir}/
        awk '{{if ($2>1000) print $0 }}' {input.depth} |cut -f -1,4 > {params.tmp_dir}/{wildcards.sample}.depth.gt1000
        python $CONDA_PREFIX/bin/scripts/gen_kmer.py {params.tmp_dir}/{wildcards.sample}.contigs.fa 1000 4 
        python $CONDA_PREFIX/bin/scripts/Filter_tooshort.py {params.tmp_dir}/{wildcards.sample}.contigs.fa 1000
        run_metabinner.sh -a {params.tmp_dir}/{wildcards.sample}.contigs_1000.fa -o {output.bins} -d {params.tmp_dir}/{wildcards.sample}.depth.gt1000 \
        -k {params.tmp_dir}/{wildcards.sample}.contigs_kmer_4_f1000.csv -t {threads} -p $CONDA_PREFIX/bin > {log} 2>&1 || true
        if [ -f {output.bins}/metabinner_res/metabinner_result.tsv ]; then
            mv {output.bins}/metabinner_res/metabinner_result.tsv {output.bins}/{wildcards.sample}_metabinner.tsv
            python $CONDA_PREFIX/bin/scripts/gen_bins_from_tsv.py -f {params.tmp_dir}/{wildcards.sample}.contigs_1000.fa \
            -r {output.bins}/{wildcards.sample}_metabinner.tsv -o {output.bins}
            counter=1
            for i in {output.bins}/*.fa; do
                mv $i {output.bins}/{wildcards.sample}_metabinner.$counter.fa
                ((counter++))
            done
        fi
        rm -rf {output.bins}/metabinner_res
        rm -rf {params.tmp_dir}
        """
