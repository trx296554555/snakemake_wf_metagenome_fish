rule gather_all_bins:
    input:
        das_tool_bins=get_das_tool_names,
        metabinner_bins=get_metabinner_names,
    output:
        bins=directory(config["root"] + "/" + config["folder"]["bins_dereplication"] + "/all_bins")
    run:
        shell("mkdir -p {output.bins}")
        for bin in input.das_tool_bins:
            shell("ln -s {bin} {output.bins}")
        for bin in input.metabinner_bins:
            shell("ln -s {bin} {output.bins}")


rule dereplicate_bins:
    input:
        bins=directory(config["root"] + "/" + config["folder"]["bins_dereplication"] + "/all_bins"),
    output:
        bins_derep1=directory(config["root"] + "/" + config["folder"]["bins_dereplication"] + "/rep_bins_comp70_con10_sa099_nc03"),
        bins_derep2=directory(config["root"] + "/" + config["folder"]["bins_dereplication"] + "/rep_bins_comp70_con10_sa095_nc01")
    conda:
        config["root"] + "/" + config["envs"] + "/" + "drep.yaml"
    params:
        sa=0.99,nc=0.3,comp=70,con=10
    log:
        config["root"] + "/" + config["folder"]["logs"] + "/drep_dereplicate.log"
    threads:
        72
    message:
        "Dereplicate bins ----------------------"
    shell:
        """
        dRep dereplicate {output.bins_derep1} -g {input.bins} -p {threads} -sa 0.99 -nc 0.3 -comp 70 -con 10 > {log} 2>&1
        dRep dereplicate {output.bins_derep2} -g {input.bins} -p {threads} -sa 0.95 -nc 0.1 -comp 70 -con 10 > {log} 2>&1
        """