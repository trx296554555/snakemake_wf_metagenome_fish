rule index_and_alignment:
    input:
        contigs=config["root"] + "/" + config["folder"]["assemble_contigs"] + "/{sample}/{sample}.contigs.fa",
        fq1=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_1.fq",
        fq2=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_2.fq"
    output:
        config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/{sample}.bam"
    message:
        "000: Indexing and alignment of contigs to reads --------------------------"
    threads:
        24
    params:
        index=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/{sample}_contig_index/{sample}",
        sam=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/{sample}.sam"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "binning.yaml"
    log:
        config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/index_and_alignment.log"
    shell:
        """
        mkdir -p `dirname {params.index}`
        bowtie2-build -f {input.contigs} {params.index} --threads {threads} --offrate 1 > {log} 2>&1
        bowtie2 -x {params.index} -1 {input.fq1} -2 {input.fq2} -S {params.sam} --threads {threads} >> {log} 2>&1
        samtools view -bS -@ {threads} {params.sam} -o {output} 2>&1
        samtools sort -@ {threads} {output} -o {output} 2>&1
        samtools index {output} 2>&1
        """

rule binning_metabat2:
    input:
        bam=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/{sample}.bam",
        contigs=config["root"] + "/" + config["folder"]["assemble_contigs"] + "/{sample}/{sample}.contigs.fa"
    output:
        bins=directory(config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/metabat2_bins")
    message:
        "001: Binning with metabat2 --------------------------"
    threads:
        24
    params:
        depth=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/{sample}_depth.txt"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "binning.yaml"
    log:
        config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/metabat2.log"
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {params.depth} {input.bam} > {log} 2>&1
        mataBat2 -m 1500 -i {input.contigs} -a {params.depth} -o {output.bins}/{wildcards.sample} -t {threads} >> {log} 2>&1
        """

rule binning_concoct:
    input:
        bam=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/{sample}.bam",
        contigs=config["root"] + "/" + config["folder"]["assemble_contigs"] + "/{sample}/{sample}.contigs.fa"
    output:
        opt=temp(directory(config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/concoct_opt")),
        bins=directory(config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/concoct_bins"),
        contigs_10K=temp(config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/{sample}.contigs_10K.fa"),
        bad=temp(config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/{sample}.bad"),
        coverage=temp(config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/{sample}.coverage_table.tsv")
    message:
        "002: Binning with concoct --------------------------"
    threads:
        24
    conda:
        config["root"] + "/" + config["envs"] + "/" + "binning.yaml"
    log:
        config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/concoct.log"
    shell:
        """
        cut_up_fasta.py {input.contigs} -c 10000 -o 0 --merge_last -b {output.bad} > {output.contigs_10K} 2>&1 > {log}
        concoct_coverage_table.py {output.bad} {input.bam} > {output.coverage} 2>&1 >> {log}
        mkdir -p {output.bins}
        concoct --composition_file {output.contigs_10K} --coverage_file {output.coverage} -b {output.opt} --threads {threads} >> {log} 2>&1
        merge_cutup_clustering.py {output.opt}/clustering_gt1000.csv > {output.opt}/clustering_merged.csv 2>&1 >> {log}
        awk -F ',' -v string="Y0119C" 'NR==1 {print; next} {OFS=","; $2 = string $2; print}' {output.opt}/clustering_merged.csv > {output.opt}/clustering_merged_renamed.csv
        extract_fasta_bins.py {input.contigs} {output.opt}/clustering_merged_renamed.csv --output_path {output.bins} >> {log} 2>&1
        """


rule binning_maxbin2:
    input:
        bam=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/{sample}.bam",
        contigs=config["root"] + "/" + config["folder"]["assemble_contigs"] + "/{sample}/{sample}.contigs.fa"
    output:
        bins=directory(config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/maxbin2_bins"),
        bad=temp(config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/{sample}.bad"),
        contigs_10K=temp(config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/{sample}.contigs_10K.fa"),
        coverage=temp(config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/{sample}.coverage_table.tsv")
    message:
        "003: Binning with maxbin2 --------------------------"
    threads:
        24
    conda:
        config["root"] + "/" + config["envs"] + "/" + "binning.yaml"
    log:
        config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/maxbin2.log"
    shell:
        """
        cut_up_fasta.py {input.contigs} -c 10000 -o 0 --merge_last -b {output.bad} > {output.contigs_10K} 2>&1 > {log}
        concoct_coverage_table.py {output.bad} {input.bam} > {output.coverage} 2>&1 >> {log}
        run_MaxBin.pl -contig {input.contigs} -abund {output.coverage} -out {output.bins}/{wildcards.sample} -thread {threads} >> {log} 2>&1
        """


