import glob
import os


def get_das_tool_names(wildcards):
    das_bins_dir = checkpoints.refine_bins_DAS_tool.get(**wildcards).output.bins
    das_bins = glob.glob(das_bins_dir + "/*.fa_DASTool_bins/*.f*a")
    return das_bins


def get_metabinner_names(wildcards):
    metabinner_bins_dir = checkpoints.metabinner.get(**wildcards).output.bins
    metabinner_bins = glob.glob(metabinner_bins_dir + "/*.f*a")
    return metabinner_bins


rule gather_need_bins:
    input:
        das_tool_bins=get_das_tool_names,
        metabinner_bins=get_metabinner_names,
    output:
        bins=config["root"] + "/" + config["folder"]["bins_refine"] + "/need_bins/{sample}_bin.list",
    run:
        shell("mkdir -p `dirname {output.bins}`")
        with open(output.bins,"w") as f:
            for bin in input.das_tool_bins:
                f.write(bin + "\n")
            for bin in input.metabinner_bins:
                f.write(bin + "\n")


rule gather_all_bins:
    input:
        sample_bins=expand(config["root"] + "/" + config["folder"][
            "bins_refine"] + "/need_bins/" + "{sample}_bin.list",sample=get_run_sample()),
        co_assemble_bins=expand(config["root"] + "/" + config["folder"][
            "bins_refine"] + "/need_bins/" + "{item}_bin.list",item=get_co_item()) if config["co_assemble"][
            "enable"] else [],
    output:
        all_bins=directory(config["root"] + "/" + config["folder"]["bins_refine"] + "/all_bins/"),
        gather_done=touch(config["root"] + "/" + config["folder"]["bins_refine"] + "/all_bins/gather.done")
    shell:
        """
        mkdir -p {output.all_bins}
        for eachbins in {input.sample_bins}; do
            cat $eachbins | xargs -I {{}} ln -s {{}} {output.all_bins}
        done
        for eachbins in {input.co_assemble_bins}; do
            cat $eachbins | xargs -I {{}} ln -s {{}} {output.all_bins}
        done
        """


rule dereplicate_bins:
    input:
        bins=config["root"] + "/" + config["folder"]["bins_refine"] + "/all_bins",
    output:
        species_bins=directory(config["root"] + "/" + config["folder"]["bins_dereplication"] + "/species_bins_comp70_con10_ANI095"),
        checkm_results=config["root"] + "/" + config["folder"]["bins_dereplication"] + "/species_bins_comp70_con10_ANI095/data_tables/genomeInfo.csv",
        strain_bins=directory(config["root"] + "/" + config["folder"]["bins_dereplication"] + "/strain_bins_comp70_con10_ANI099"),
        done=touch(config["root"] + "/" + config["folder"]["bins_dereplication"] + "/dereplicate_bins.done")
    conda:
        config["root"] + "/" + config["envs"] + "/" + "drep.yaml"
    params:
        comp=70, con=10, species_ANI=0.95, strain_ANI=0.99, nc=0.1,
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["bins_dereplication"] + "/drep_dereplicate_bins.benchmark.txt"
    log:
        config["root"] + "/" + config["folder"]["bins_dereplication"] + "/drep_dereplicate_bins.log"
    threads:
        64
    message:
        "17: Dereplicate bins at species and strain level ----------------------"
    shell:
        """
        dRep dereplicate {output.species_bins} -g {input.bins}/*.f*a -p {threads} -sa {params.species_ANI} -nc {params.nc} -comp {params.comp} -con {params.con} > {log} 2>&1
        dRep dereplicate {output.strain_bins} -g {input.bins}/*.f*a -p {threads} -sa {params.strain_ANI} -nc {params.nc} -comp {params.comp} -con {params.con} --genomeInfo {output.checkm_results} > {log} 2>&1
        """


rule get_MAGs:
    input:
        species_bins=config["root"] + "/" + config["folder"]["bins_dereplication"] + "/species_bins_comp70_con10_ANI095",
        strain_bins=config["root"] + "/" + config["folder"]["bins_dereplication"] + "/strain_bins_comp70_con10_ANI099",
    output:
        mags_dir=directory(config["root"] + "/" + config["folder"]["bins_dereplication"] + "/mag_bins"),
        mags_done=touch(config["root"] + "/" + config["folder"]["bins_dereplication"] + "/get_MAGs_bins.done")
    params:
        anno_level = config["bins_dereplicate"]["anno_level"]
    shell:
        """
        mkdir -p {output.mags_dir}
        if [ {params.anno_level} == "species" ]; then
            for eachbins in {input.species_bins}/dereplicated_genomes/*.f*a; do
                ln -sf $eachbins {output.mags_dir}
            done
        elif [ {params.anno_level} == "strain" ]; then
            for eachbins in {input.strain_bins}/dereplicated_genomes/*.f*a; do
                ln -sf $eachbins {output.mags_dir}
            done
        else
            exit 1
        fi
        """


rule report_contigs_binning:
    input:
        all_bins = config["root"] + "/" + config["folder"]["bins_refine"] + "/all_bins/",
        gather_done = config["root"] + "/" + config["folder"]["bins_refine"] + "/all_bins/gather.done",
        mags_done=config["root"] + "/" + config["folder"]["bins_dereplication"] + "/get_MAGs_bins.done",
        checkm_results=config["root"] + "/" + config["folder"]["bins_dereplication"] + "/species_bins_comp70_con10_ANI095/data_tables/genomeInfo.csv"
    output:
        contigs_binning_report=config["root"] + "/" + config["folder"]["reports"] + "/07_contigs_binning.report"
    params:
        species_bins = config["root"] + "/" + config["folder"]["bins_dereplication"] + "/species_bins_comp70_con10_ANI095/dereplicated_genomes",
        strain_bins = config["root"] + "/" + config["folder"]["bins_dereplication"] + "/strain_bins_comp70_con10_ANI099/dereplicated_genomes",
    run:
        res_dict = {}
        for i in os.listdir(input.all_bins):
            if i.endswith('fa') or i.endswith('fasta'):
                sample = '_'.join(i.split('.')[0].split('_')[0:-1])
                software = i.split('.')[0].split('_')[-1]
                if sample not in res_dict:
                    res_dict[sample] = {}
                    res_dict[sample][software] = 1
                else:
                    if software not in res_dict[sample]:
                        res_dict[sample][software] = 1
                    else:
                        res_dict[sample][software] += 1
        table = pd.DataFrame(res_dict).T
        table['total'] = table.sum(axis=1)

        checkm_table = pd.read_csv(input.checkm_results)
        col_list = []
        for ani in ['ANI99', 'ANI95']:
            if ani == 'ANI99':
                dereplicate_bins_list = os.listdir(params.strain_bins)
            else:
                dereplicate_bins_list = os.listdir(params.species_bins)
            for comp in [70,80,90,100]:
                for con in [10,5]:
                    col_name = '{}_comp{}_con{}'.format(ani, comp, con)
                    col_list.append(col_name)
                    checkm_table[col_name] = np.nan
                    for i in range(len(checkm_table)):
                        if checkm_table.loc[i, 'genome'] in dereplicate_bins_list and checkm_table.loc[i, 'completeness'] >= comp and checkm_table.loc[i, 'contamination'] <= con:
                            checkm_table.loc[i, col_name] = 1

        checkm_table = checkm_table.dropna(subset=['ANI99_comp70_con10', 'ANI95_comp70_con10'],how='all')
        # 取出checkm_table中col_list中的列
        checkm_table = checkm_table[['genome'] + col_list]
        # 将checkm_table中的nan替换为0
        checkm_table = checkm_table.fillna(0)
        checkm_table['genome'] = checkm_table['genome'].apply(lambda x: '_'.join(x.split('.')[0].split('_')[0:-1]))
        checkm_table = checkm_table.set_index('genome')
        # 将checkm_table中索引相同的行相加
        checkm_table = checkm_table.groupby(checkm_table.index).sum()
        # 将table 和 checkm_table 合并
        table = table.join(checkm_table)
        table = table.fillna(0)
        table.to_csv(output.contigs_binning_report, sep='\t')
