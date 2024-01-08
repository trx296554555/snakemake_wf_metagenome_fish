## TODO 统计MAGs 大小，N50，编码蛋白数量
# new_MAGs_report = config["root"] + "/" + config["folder"]["reports"] + "/08_new_MAGs.report",

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