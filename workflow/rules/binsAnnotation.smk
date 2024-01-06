rule prokka_anno_bins: # rgi,dbcan,vfdb 等等，可以结合下面定量的rule
    input:
        bin=config["root"] + "/" + config["folder"]["bins_dereplication"] + "/rep_bins_user/{bin}.fa",
    output:
        gff=config["root"] + "/" + config["folder"]["bins_anno"] + "/prokka_anno_bins/{bin}.gff",
    conda:
        config["root"] + "/" + config["envs"] + "/" + "prokka.yaml"
    message:
        " Prokka annotation"
    threads:
        12
    params:
        outdir=directory(config["root"] + "/" + config["folder"]["bins_anno"] + "/prokka_anno_bins"),
    log:
        config["root"] + "/" + config["folder"]["bins_anno"] + "/prokka_anno_bins/{bin}.log"
    shell:
        """
        prokka --cpus {threads} --outdir {params.outdir} --prefix {wildcards.bin} --force --addgenes --kingdom Bacteria \
        --gcode 11 --rfam --rnammer --usegenus {input.bin} > {log} 2>&1
        mv {params.outdir}/{wildcards.bin}.gff {output.gff}
        """


rule cat_bins:
    input:
        bins_derep2=directory(config["root"] + "/" + config["folder"]["bins_dereplication"] + "/rep_bins_user"),
    output:
        cat_bins=config["root"] + "/" + config["folder"]["bins_quant"] + "/cat_bins.fasta"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "cat_bins.yaml"
    message:
        " Concatenation of bins"
    shell:
        """
        cat {input.bins_derep2}/*.fasta > {output.cat_bins}
        """

rule make_bins_index:
    input:
        cat_bins=config["root"] + "/" + config["folder"]["bins_quant"] + "/cat_bins.fasta",
    output:
        index=config["root"] + "/" + config["folder"]["bins_quant"] + "/salmon_index"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "salmon.yaml"
    message:
        " Indexing of bins"
    threads:
        12
    log:
        config["root"] + "/" + config["folder"]["logs"] + "/make_bins_index.log"
    shell:
        """
        salmon index -t {input.cat_bins} -i {output.index} -p {threads} > {log} 2>&1
        """

rule quantify_bins_expression:
    input:
        idx=config["root"] + "/" + config["folder"]["bins_quant"] + "/salmon_index",
        fq1=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_1.fq",
        fq2=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_2.fq",
    output:
        expression=config["root"] + "/" + config["folder"]["bins_quant"] + "/quant_bins/{sample}/{sample}.sf",
    conda:
        config["root"] + "/" + config["envs"] + "/" + "salmon.yaml"
    threads:
        12
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["bins_quant"] + "/quant_bins/{sample}/{sample}.quant.log"
    log:
        config["root"] + "/" + config["folder"]["bins_quant"] + "/quant_bins/{sample}/{sample}.quant.log"
    message:
        "Quantify bins expression using salmon"
    params:
        dir=config["root"] + "/" + config["folder"]["bins_quant"] + "/quant_bins/{sample}"
    shell:
        """
        salmon quant -i {input.idx} -l A --meta -1 {input.fq1} -2 {input.fq2} \
        -o {params.dir} -p {threads} --validateMappings --seqBias --quiet > {log} 2>&1
        mv {params.dir}/quant.sf {output.expression}
        """

rule merge_bins_quant:
    input:
        config["root"] + "/" + config["folder"]["bins_quant"] + "/quant_bins/done"
    output:
        config["root"] + "/" + config["folder"]["bins_quant"] + "/bin_abundance_table.tab"
    params:
        bin_folder=config["root"] + "/" + config["folder"]["bins_dereplication"] + "/rep_bins_user",
        cat_bins=config["root"] + "/" + config["folder"]["bins_quant"] + "/cat_bins.fasta",
        alignment_dir=config["root"] + "/" + config["folder"]["bins_quant"] + "/quant_bins/",
        script_dir=config["root"] + "/" + config["folder"]["scripts"],
        quant_dir=config["root"] + "/" + config["folder"]["bins_quant"] + "/quant_files"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "py27.yaml"
    # metawrap workflow
    shell:
        """
        # "summarize salmon files..."
        cd {params.alignment_dir}
        python {params.script_dir}/summarize_salmon_files.py
        mkdir -p {params.quant_dir}
        for f in $(ls {params.alignment_dir} | grep .quant.counts); do mv {params.alignment_dir}/$f {params.quant_dir}/; done
        {params.script_dir}/split_salmon_out_into_bins.py ${params.quant_dir}/ {params.bin_folder} {params.cat_bins} > {output}
        """