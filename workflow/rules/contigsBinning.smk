rule index_and_alignment:
    input:
        contigs=config["root"] + "/" + config["folder"]["assemble_contigs"] + "/{sample}/{sample}.contigs.fa",
        fq1=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_1.fq",
        fq2=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_2.fq"
    output:
        index=temp(directory(config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/{sample}_contig_index")),
        sam=temp(config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/{sample}.sam"),
        bam=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/{sample}.bam",
        depth=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/{sample}.depth",
        coverage=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/{sample}.coverage"
    message:
        "13: Indexing and alignment of contigs to reads --------------------------"
    threads:
        24
    conda:
        config["root"] + "/" + config["envs"] + "/" + "binning.yaml"
    log:
        config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/mapping.log"
    shell:
        """
        mkdir -p {output.index}
        bowtie2-build -f {input.contigs} {output.index}/{wildcards.sample} --threads {threads} --offrate 1 > /dev/null 2>&1
        bowtie2 -x {output.index}/{wildcards.sample} -1 {input.fq1} -2 {input.fq2} -S {output.sam} --threads {threads} > {log} 2>&1
        samtools view -bS -@ {threads} {output.sam} -o {output.bam} 2>&1
        samtools sort -@ {threads} {output.bam} -o {output.bam} 2>&1
        samtools index {output.bam} 2>&1
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {output.bam} >> {log} 2>&1
        awk -v filename={wildcards.sample} 'BEGIN {{OFS="\t"; print "contig", "cov_mean_sample_" filename}} NR>1 {{print $1, $3 + $4}}' {output.depth} > {output.coverage}
        """


rule binning_metabat2:
    input:
        depth=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/{sample}.depth",
        contigs=config["root"] + "/" + config["folder"]["assemble_contigs"] + "/{sample}/{sample}.contigs.fa"
    output:
        bins=directory(config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/metabat2_bins")
    message:
        "14: Binning with metabat2 --------------------------"
    threads:
        12
    conda:
        config["root"] + "/" + config["envs"] + "/" + "binning.yaml"
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
        "14: Binning using concoct --------------------------"
    threads:
        12
    conda:
        config["root"] + "/" + config["envs"] + "/" + "concoct.yaml"
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
        "14: Binning using maxbin2 --------------------------"
    threads:
        12
    conda:
        config["root"] + "/" + config["envs"] + "/" + "binning.yaml"
    log:
        config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/maxbin2.log"
    shell:
        """
        mkdir -p {output.bins}
        run_MaxBin.pl -contig {input.contigs} -abund {input.coverage} -out {output.bins}/{wildcards.sample}_maxbin2 -thread {threads} > {log} 2>&1 || true
        mkdir -p {output.bins}
        """


rule refine_bins_DAS_tool:
    input:
        contigs=config["root"] + "/" + config["folder"]["assemble_contigs"] + "/{sample}/{sample}.contigs.fa",
        metabat2_bins=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/metabat2_bins",
        concoct_bins=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/concoct_bins",
        proteins=config["root"] + "/" + config["folder"]["gene_prediction"] + "/{sample}/{sample}.protein.faa",
        maxbin2_bins=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/maxbin2_bins"
    output:
        bins=directory(config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/merged_das_bins")
    message:
        "15: Refine bins using DAS tool --------------------------"
    threads:
        12
    conda:
        config["root"] + "/" + config["envs"] + "/" + "dastool.yaml"
    params:
        metabat2_c2b=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/metabat2_contigs2bin.tsv",
        concoct_c2b=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/concoct_contigs2bin.tsv",
        maxbin2_c2b=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/maxbin2_contigs2bin.tsv"
    log:
        config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/merged_das_bins.log"
    shell:
        """
        Fasta_to_Contig2Bin.sh -e fasta -i {input.maxbin2_bins} > {params.maxbin2_c2b}
        Fasta_to_Contig2Bin.sh -e fa -i {input.metabat2_bins} > {params.metabat2_c2b}
        Fasta_to_Contig2Bin.sh -e fa -i {input.concoct_bins} > {params.concoct_c2b}
        mkdir -p {output.bins}
        DAS_Tool -i {params.metabat2_c2b},{params.concoct_c2b},{params.maxbin2_c2b} -l metabat2,concoct,maxbin2 \
        -c {input.contigs} -o {output.bins}/{wildcards.sample}.contigs.fa -t {threads} \
        -p {input.proteins} --write_bins --score_threshold 0 --write_bin_evals > {log} 2>&1 || true
        """


rule metabinner:
    input:
        bam=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/{sample}.bam",
        contigs=config["root"] + "/" + config["folder"]["assemble_contigs"] + "/{sample}/{sample}.contigs.fa",
        depth=config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/{sample}.depth",
    output:
        bins=directory(config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/metabinner_bins")
    message:
        "15: Refine bins using metabinner.yaml --------------------------"
    threads:
        12
    conda:
        config["root"] + "/" + config["envs"] + "/" + "metabinner.yaml.yaml"
    log:
        config["root"] + "/" + config["folder"]["contigs_binning"] + "/{sample}/metabinner.log"
    shell:
        """
        mkdir -p {output.bins}
        awk '{{if ($2>1000) print $0 }}' {input.depth} |cut -f -1,4 > {input.depth}.gt1000
        gen_kmer.py {input.contigs} 1000 4 
        Filter_tooshort.py {input.contigs} 1000
        run_metabinner.sh -a {input.contigs} -o {output.bins} -d {params.depth}.gt1000 \
        -k `dirname {input.contigs}`*.contigs_kmer_4_f1000.csv -t {threads} -p $CONDA_PREFIX/bin  > {log} 2>&1
        """
