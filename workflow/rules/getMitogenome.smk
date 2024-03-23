import os
import re
import subprocess
import pandas as pd

# # If the reference genome contains MT sequence, can use only host reads
# rule get_mitogenome:
#     input:
#         fq1=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_1.fq",
#         fq2=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_2.fq",
#         host_fq1=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_host_1.fq",
#         host_fq2=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_host_2.fq",
#     output:
#         all_fq1=temp(config["root"] + "/" + config["folder"]["mitogenome"] + "/{sample}/{sample}_all_1.fq"),
#         all_fq2=temp(config["root"] + "/" + config["folder"]["mitogenome"] + "/{sample}/{sample}_all_2.fq"),
#         annotation=temp(directory(config["root"] + "/" + config["folder"]["mitogenome"] + "/{sample}/mt_annotation")),
#         assembly=temp(directory(config["root"] + "/" + config["folder"]["mitogenome"] + "/{sample}/mt_assembly")),
#         res_dir=directory(config["root"] + "/" + config["folder"]["mitogenome"] + "/{sample}/{sample}.result")
#     message:
#         "07: Use mitoZ to obtain the mitochondrial genome of the host from the host sequence ------------------------------------------"
#     threads:
#         24
#     params:
#         out_path=config["root"] + "/" + config["folder"]["mitogenome"] + "/{sample}",
#     conda:
#         config["root"] + "/" + config["envs"] + "/" + "mitozEnv.yaml"
#     log:
#         config["root"] + "/" + config["folder"]["mitogenome"] + "/{sample}/{sample}_mitoz.log"
#     shell:
#         """
#         cat {input.fq1} {input.host_fq1} > {output.all_fq1}
#         cat {input.fq2} {input.host_fq2} > {output.all_fq2}
#
#         mitoz all --workdir {params.out_path} --outprefix {wildcards.sample} --fq1 {output.all_fq1} --fq2 {output.all_fq2} \
#         --skip_filter --data_size_for_mt_assembly 0 --assembler megahit --clade Chordata --requiring_taxa Chordata \
#         --insert_size 350 --kmers_megahit 139 --thread_number {threads} > {log} 2>&1 || true
#
#         if [ ! -d {params.out_path}/{wildcards.sample}.result ]; then
#             mkdir {params.out_path}/{wildcards.sample}.result
#             mkdir {params.out_path}/mt_annotation
#             touch {params.out_path}/task_failed.error
#         else
#             mv {output.res_dir}/{wildcards.sample}.megahit.result {params.out_path}
#             mv {output.res_dir}/{wildcards.sample}.{wildcards.sample}.megahit.mitogenome.fa.result/* {output.res_dir}
#             rm -rf {output.res_dir}/{wildcards.sample}.{wildcards.sample}.megahit.mitogenome.fa.result
#         fi
#         """


rule get_whole_reads:
    input:
        fq1=expand(config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_1.fq", sample=get_run_sample()),
        fq2=expand(config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_2.fq", sample=get_run_sample()),
        host_fq1=expand(config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_host_1.fq", sample=get_run_sample()),
        host_fq2=expand(config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_host_2.fq", sample=get_run_sample()),
    output:
        all_fq1=expand(config["root"] + "/" + config["folder"]["mitogenome"] + "/{individual}/{individual}_all_1.fq", individual=get_run_individual()),
        all_fq2=expand(config["root"] + "/" + config["folder"]["mitogenome"] + "/{individual}/{individual}_all_2.fq", individual=get_run_individual()),
    params: individual_list=get_run_individual()
    run:
        for i in params.individual_list:
            if not os.path.exists(config["root"] + "/" + config["folder"]["mitogenome"] + "/" + i):
                os.mkdir(config["root"] + "/" + config["folder"]["mitogenome"] + "/" + i)
            else:
                if os.path.exists(config["root"] + "/" + config["folder"]["mitogenome"] + "/" + i + "/merge.complete"):
                    continue
            fq1_list = []
            fq2_list = []
            if os.path.exists(config["root"] + "/" + config["folder"]["rm_host"] + "/" + i + 'C'):
                fq1_list.append(config["root"] + "/" + config["folder"]["rm_host"] + "/" + i + "C/" + i + 'C_paired_1.fq')
                fq1_list.append(config["root"] + "/" + config["folder"]["rm_host"] + "/" + i + "C/" + i + "C_host_1.fq")
                fq2_list.append(config["root"] + "/" + config["folder"]["rm_host"] + "/" + i + "C/" + i + 'C_paired_2.fq')
                fq2_list.append(config["root"] + "/" + config["folder"]["rm_host"] + "/" + i + "C/" + i + "C_host_2.fq")
            elif os.path.exists(config["root"] + "/" + config["folder"]["rm_host"] + "/" + i + 'P'):
                fq1_list.append(config["root"] + "/" + config["folder"]["rm_host"] + "/" + i + "P/" + i + 'P_paired_1.fq')
                fq1_list.append(config["root"] + "/" + config["folder"]["rm_host"] + "/" + i + "P/" + i + "P_host_1.fq")
                fq2_list.append(config["root"] + "/" + config["folder"]["rm_host"] + "/" + i + "P/" + i + 'P_paired_2.fq')
                fq2_list.append(config["root"] + "/" + config["folder"]["rm_host"] + "/" + i + "P/" + i + "P_host_2.fq")
            else:
                raise Exception(f"Unexpect error: {i}")
            shell(f"cat {' '.join(fq1_list)} > {config['root'] + '/' + config['folder']['mitogenome'] + '/' + i + '/' + i + '_all_1.fq'}")
            shell(f"cat {' '.join(fq2_list)} > {config['root'] + '/' + config['folder']['mitogenome'] + '/' + i + '/' + i + '_all_2.fq'}")
            shell(f"touch {config['root'] + '/' + config['folder']['mitogenome'] + '/' + i + '/merge.complete'}")


# If the reference genome contains MT sequence, can use only host reads
rule get_mitogenome:
    input:
        all_fq1=config["root"] + "/" + config["folder"]["mitogenome"] + "/{individual}/{individual}_all_1.fq",
        all_fq2=config["root"] + "/" + config["folder"]["mitogenome"] + "/{individual}/{individual}_all_2.fq",
    output:
        annotation=directory(config["root"] + "/" + config["folder"]["mitogenome"] + "/{individual}/mt_annotation"),
        assembly=directory(config["root"] + "/" + config["folder"]["mitogenome"] + "/{individual}/mt_assembly"),
        res_dir=directory(config["root"] + "/" + config["folder"]["mitogenome"] + "/{individual}/{individual}.result"),
        done=touch(config["root"] + "/" + config["folder"]["mitogenome"] + "/{individual}/mitogenome.done")
    message:
        "08: Use mitoZ to obtain the mitochondrial genome of the host from the host sequence ------------------------------------------"
    threads:
        24
    params:
        out_path=config["root"] + "/" + config["folder"]["mitogenome"] + "/{individual}",
    conda:
        config["root"] + "/" + config["envs"] + "/" + "mitozEnv.yaml"
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["mitogenome"] + "/{individual}.log"
    log:
        config["root"] + "/" + config["folder"]["mitogenome"] + "/{individual}/{individual}_mitoz.log"
    shell:
        """
        mitoz all --workdir {params.out_path} --outprefix {wildcards.individual} --fq1 {input.all_fq1} --fq2 {input.all_fq2} \
        --skip_filter --data_size_for_mt_assembly 0 --assembler megahit --clade Chordata --requiring_taxa Chordata \
        --insert_size 350 --kmers_megahit 139 --thread_number {threads} > {log} 2>&1 || true

        if [ ! -d {params.out_path}/{wildcards.individual}.result ]; then
            mkdir {params.out_path}/{wildcards.individual}.result
            mkdir -p {params.out_path}/mt_annotation
            touch {params.out_path}/task_failed.error
        else
            mv {output.res_dir}/{wildcards.individual}.megahit.result {params.out_path}
            mv {output.res_dir}/{wildcards.individual}.{wildcards.individual}.megahit.mitogenome.fa.result/* {output.res_dir}
            rm -rf {output.res_dir}/{wildcards.individual}.{wildcards.individual}.megahit.mitogenome.fa.result
        fi
        """


rule report_get_mitogenome:
    input:
        done_files=expand(config["root"] + "/" + config["folder"]["mitogenome"] + "/{individual}/mitogenome.done", individual=get_run_individual()),
        res_dir=expand(config["root"] + "/" + config["folder"]["mitogenome"] + "/{individual}/{individual}.result", individual=get_run_individual()),
    output:
        get_mitogenome_report=config["root"] + "/" + config["folder"]["reports"] + "/03_get_mitogenome.report"
    run:
        res_dict = {}
        for m_dir in input.res_dir:
            sample = os.path.basename(m_dir).split(".")[0]
            res_dict[sample] = {}
            if os.path.exists(os.path.join(os.path.dirname(m_dir), "task_failed.error")):
                res_dict[sample]["Status"] = "Failed"
            else:
                res_dict[sample]["Status"] = "Success"
                depth_file = os.path.join(m_dir, "circos.depth.txt")
                command = "awk '{sum += $NF} END {print sum/NR}' " + depth_file
                avg_depth = subprocess.check_output(command, shell=True, executable="/bin/bash", text=True)
                res_dict[sample]["Average_depth"] = float(avg_depth.strip())
                with open(os.path.join(m_dir, 'summary.txt'), 'r') as ipt:
                    next(ipt)
                    info = re.split(r'\s\s+', next(ipt).strip())
                    res_dict[sample]["Seq_id"] = info[0]
                    res_dict[sample]["Length(bp)"] = info[1]
                    res_dict[sample]["Circularity"] = info[2]
                    res_dict[sample]["Closely_related_species"] = info[3]
                    for line in ipt:
                        if line.startswith("Protein coding genes totally found"):
                            res_dict[sample]["CDS"] = line.strip().split(":")[1].strip()
                        elif line.startswith("tRNA genes totally found"):
                            res_dict[sample]["tRNA"] = line.strip().split(":")[1].strip()
                        elif line.startswith("rRNA genes totally found"):
                            res_dict[sample]["rRNA"] = line.strip().split(":")[1].strip()
                        else:
                            continue
        table = pd.DataFrame(res_dict).T
        table = table.sort_index()
        table.to_csv(output.get_mitogenome_report, sep="\t", index=True, header=True)