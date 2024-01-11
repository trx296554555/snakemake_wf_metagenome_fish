## TODO 统计MAGs 大小，N50，编码蛋白数量
# new_MAGs_report = config["root"] + "/" + config["folder"]["reports"] + "/08_new_MAGs.report",


rule get_bins_prokka_anno:
    input:
        mag=config["root"] + "/" + config["folder"]["bins_dereplication"] + "/mag_bins/{bin}.fa",
        all_classify_tsv=config["root"] + "/" + config["folder"][
            "bins_classify"] + "/gtdbtk_classify_wf/res.all.summary.tsv"
    output:
        genes=config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/{bin}/{bin}.gene.fa",
        proteins=config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/{bin}/{bin}.protein.faa",
        prokka_anno=config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/{bin}/{bin}.prokka.txt",
    conda:
        config["root"] + "/" + config["envs"] + "/" + "prokka.yaml"
    message:
        "20: Run prokka to get bins annotation"
    threads:
        12
    params:
        mags_dir=config["root"] + "/" + config["folder"]["bins_dereplication"] + "/mag_bins",
        out_dir=config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/{bin}",
    shell:
        """
        mkdir -p {params.out_dir}
        raw_id=$(basename {input.mag} .fa)
        new_id=$(printf "%04d" $(ls $(dirname {input.mag})| nl | awk '$2=="'$raw_id.fa'" {{print $1}}'))
        kingdom=$(cat {input.all_classify_tsv} | grep $raw_id$'\t' | cut -f3 | cut -d \; -f1 | sed 's/d__//')
        genus=$(cat {input.all_classify_tsv} | grep $raw_id$'\t' | cut -f3 | cut -d \; -f6 | sed 's/g__//')
        species=$(cat {input.all_classify_tsv} | grep $raw_id$'\t' | cut -f3 | cut -d \; -f7 | sed 's/s__//'|sed 's/ /_/g')
        if [ $genus ]; then genus=$genus;else genus="''"; fi
        if [ $species ]; then species=$species;else species="''"; fi
        ln -sf {input.mag} {params.out_dir}/$new_id.fa
        
        prokka --quiet --cpus {threads} --outdir {params.out_dir} --prefix {wildcards.bin} \
        --metagenome --force --gcode 11 --kingdom ${{kingdom}} --genus ${{genus}} --species ${{species}} \
        --locustag fish_bin_${{new_id}} --centre '' --compliant {input.mag} > /dev/null 2>&1
        seqkit replace --quiet -p ".*" -r ${{raw_id}}_{{nr}} -w 0 {params.out_dir}/{wildcards.bin}.ffn > {output.genes}
        seqkit replace --quiet -p ".*" -r ${{raw_id}}_{{nr}} -w 0 {params.out_dir}/{wildcards.bin}.faa > {output.proteins}
        ln -s {params.out_dir}/{wildcards.bin}.txt {output.prokka_anno}
        """


def collect_get_MAGs(wildcards):
    mag_bins_dir = checkpoints.get_MAGs.get(**wildcards).output.mags_dir
    mag_bins = glob.glob(mag_bins_dir + "/*.f*a")
    mag_names = [os.path.basename(x).replace('.fa','') for x in mag_bins]
    prokka_gene = [config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/" + bin + "/" + bin + ".gene.fa"
                  for bin in mag_names]
    prokka_prptein = [config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/" + bin + "/" + bin + ".protein.faa"
                  for bin in mag_names]
    return {"gene":prokka_gene, "protein":prokka_prptein}


rule get_bins_prokka_anno_res:
    input:
        unpack(collect_get_MAGs),
    output:
        concatenation_mag_gene=config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/concatenation_prokka/prokka_gene.fa",
        concatenation_mag_protein=config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/concatenation_prokka/prokka_protein.faa",
    shell:
        """
        cat {input.gene} > {output.concatenation_mag_gene}
        cat {input.protein} > {output.concatenation_mag_protein}
        """


rule pruning_MAG_gene_redundancy:
    input:
        concatenation_mag_gene=config["root"] + "/" + config["folder"][
            "bins_anno_prokka"] + "/concatenation_prokka/prokka_gene.fa",
        concatenation_mag_protein=config["root"] + "/" + config["folder"][
            "bins_anno_prokka"] + "/concatenation_prokka/prokka_protein.faa",
    output:
        genes=config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/concatenation_prokka/prokka_gene.prune.fa",
        proteins=config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/concatenation_prokka/prokka_protein.prune.faa",
    message:
        "21: Pruning redundancy genes from all mag bins"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "prediction.yaml"
    threads:
        80
    params:
        clean_concatenation_mag_gene=config["root"] + "/" + config["folder"][
            "bins_anno_prokka"] + "/concatenation_prokka/prokka_gene.clean.fa",
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["bins_anno_prokka"] + "/pruning_MAG_gene_redundancy.benchmark.txt"
    log:
        config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/pruning_MAG_gene_redundancy.log"
    shell:
        """
        seqkit seq --threads {threads} -n {input.concatenation_mag_protein} |seqkit grep --threads {threads} -n -f - {input.concatenation_mag_gene} > {params.clean_concatenation_mag_gene}
        cd-hit-est -T {threads} -M 0 -i {params.clean_concatenation_mag_gene} -o {output.genes} -c 0.95 -aS 0.9 -g 1 -sc 1 -sf 1 -d 0 > /dev/null 2>&1
        seqkit seq --threads {threads} -n {output.genes} |seqkit grep --threads {threads} -n -f - {input.concatenation_mag_protein} > {output.proteins}
        sed -i 's/*//g' {output.proteins}
        seqkit stats --threads {threads} -T {output.proteins} > {log}
        """
