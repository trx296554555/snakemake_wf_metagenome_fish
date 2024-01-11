## TODO 统计MAGs 大小，N50，编码蛋白数量
# new_MAGs_report = config["root"] + "/" + config["folder"]["reports"] + "/08_new_MAGs.report",


checkpoint get_bins_prokka_anno:
    input:
        mag=config["root"] + "/" + config["folder"]["bins_dereplication"] + "/mag_bins/{bin}.fa",
        all_classify_tsv=config["root"] + "/" + config["folder"]["bins_classify"] + "/gtdbtk_classify_wf/res.all.summary.tsv"
    output:
        genes=config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/{bin}/{bin}.gene.fa",
        proteins=config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/{bin}/{bin}.protein.faa",
        prokka_anno=config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/{bin}/{bin}.prokka.tsv",
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
        new_id=$(printf "%04d" $(ls $(dirname {input.mag})| nl | grep $raw_id | cut -f1))
        kingdom=$(cat {input.all_classify_tsv} | grep $raw_id | cut -f3 | cut -d \; -f1 | sed 's/d__//')
        genus=$(cat {input.all_classify_tsv} | grep $raw_id | cut -f3 | cut -d \; -f6 | sed 's/g__//')
        species=$(cat {input.all_classify_tsv} | grep $raw_id | cut -f3 | cut -d \; -f7 | sed 's/s__//')
        ln -s {input.mag} {params.out_dir}/$new_bin.fa
        prokka --quiet --cpus {threads} --outdir {params.out_dir} --prefix {wildcards.bin} \
        --metagenome --force --kingdom $kingdom --genus $genus --species $species --gcode 11 \
        --locustag fish_bin_$new_id --centre '' --compliant {input.mag} > /dev/null 2>&1
        seqkit replace --quiet -p ".*" -r ${{raw_id}}_{{nr}} -w 0 {params.out_dir}/{wildcards.bin}.ffn > {output.genes}
        seqkit replace --quiet -p ".*" -r ${{raw_id}}_{{nr}} -w 0 {params.out_dir}/{wildcards.bin}.faa > {output.proteins}
        ln -s {params.out_dir}/{wildcards.bin}.tsv {output.prokka_anno}
        """


def collect_bins_prokka_anno(wildcards):
    metabinner_bins_dir = checkpoints.get_bins_prokka_anno.get(**wildcards).output.genes
    metabinner_bins = glob.glob(metabinner_bins_dir + "/*.protein.faa")
    return metabinner_bins


rule test:
    input:
        collect_bins_prokka_anno
    output:
        done=touch(config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/prokka.done")
    shell:
        """
        echo "done"
        """