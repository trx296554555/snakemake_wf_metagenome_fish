import os

"""
Author: Tang Ruixiang
Date: 2024-3-7 22:16:03
Version: 1.0 

This workflow is used to get the multiple sequence alignment (MSA) of mitochondrial genes.
The MSA is used to build the phylogenetic tree of the species.

# MSA
CDS and exon sequences are aligned at the codon level in two steps. 
First, the translated amino acids are aligned using Muscle 5 and gaps are reported onto the nucleotide sequences.
This alignment is then refined using MACSE v2 to obtain a final codon alignment unaffected by frameshifts, 
misassemblies, and sequencing errors. Amino acid alignments are then filtered to reduce the impact of errors
on evolutionary inferences using TrimAl.

# Phylogenetic treetar
Phylogenetic trees are inferred using IQ-TREE v2.1.3.

"""

configfile: "workflow/config/config.yaml"

MITO_CDS_LIST = ['COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6', 'ATP6', 'ATP8']
all_cds_list = glob_wildcards(config["root"] + "/" + config["folder"]["mitotree"] + "/cds/" + "{cds}.fa")
use_for_tree_cds = [cds for cds in all_cds_list.cds if cds in MITO_CDS_LIST]


rule all:
    input:
        config["root"] + "/" + config["folder"]["mitotree"] + "/species_tree/astral_hybrid_species.done"


rule get_msa:
    input:
        config["root"] + "/" + config["folder"]["mitotree"] + "/cds/" + "{cds}.fa"
    output:
        config["root"] + "/" + config["folder"]["mitotree"] + "/msa/" + "{cds}.fa"
    message:
        "01: get msa by muscle ------------------------------------------"
    threads:
        12
    params:
        mode="-super5",
    conda:
        config["root"] + "/" + config["envs"] + "/" + "buildtree.yaml"
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["mitotree"] + "/{cds}_msa.log"
    log:
        config["root"] + "/" + config["folder"]["mitotree"] + "/msa/{cds}_msa.log"
    shell:
        """
        muscle -threads {threads} {params.mode} {input} -output {output} > {log} 2>&1
        """


rule refine_msa:
    input:
        config["root"] + "/" + config["folder"]["mitotree"] + "/msa/" + "{cds}.fa"
    output:
        refine_nt=config["root"] + "/" + config["folder"]["mitotree"] + "/msa/" + "{cds}.refine.fa",
        refine_aa=config["root"] + "/" + config["folder"]["mitotree"] + "/msa/" + "{cds}.refine.faa"
    message:
        "02: refine msa by macse ------------------------------------------"
    params:
        code_num="2",
    conda:
        config["root"] + "/" + config["envs"] + "/" + "buildtree.yaml"
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["mitotree"] + "/{cds}_refine.log"
    log:
        config["root"] + "/" + config["folder"]["mitotree"] + "/msa/{cds}_refine.log"
    shell:
        """
        macse -prog refineAlignment -gc_def {params.code_num} -align {input} -out_NT {output.refine_nt} -out_AA {output.refine_aa} > {log} 2>&1
        """


rule trim_msa:
    input:
        refine_nt=config["root"] + "/" + config["folder"]["mitotree"] + "/msa/" + "{cds}.refine.fa",
        refine_aa=config["root"] + "/" + config["folder"]["mitotree"] + "/msa/" + "{cds}.refine.faa"
    output:
        tmp_seqs_nt=temp(config["root"] + "/" + config["folder"]["mitotree"] + "/msa/" + "{cds}.tmp1.fa"),
        tmp_code_nt=temp(config["root"] + "/" + config["folder"]["mitotree"] + "/msa/" + "{cds}.tmp2.fa"),
        trimmed_nt=config["root"] + "/" + config["folder"]["mitotree"] + "/msa/" + "{cds}.trim.fa",
    message:
        "03: trim msa by trimAl ------------------------------------------"
    params:
        script=config["root"] + "/workflow/scripts/subseq_from_stop_coden.py",
    conda:
        config["root"] + "/" + config["envs"] + "/" + "buildtree.yaml"
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["mitotree"] + "/{cds}_trim.log"
    shell:
        """
        seqkit replace --quiet -s -p ! -r N {input.refine_nt} | seqkit replace --quiet -p ";.*" -r "" | seqkit sort --quiet -n -w 0 > {output.tmp_seqs_nt}
        stopcode=`seqkit replace --quiet -s -p ! -r X {input.refine_aa} |seqkit locate -i -p "*" -P | tail -n +2 | cut -f5 | sort -n | uniq | tr "\\n" ';' | sed 's/;$//'`
        python {params.script} ${{stopcode}} {output.tmp_seqs_nt} {output.tmp_code_nt}
        trimal -keepheader -automated1 -in {output.tmp_code_nt} -out {output.trimmed_nt}
        """


rule ml_bs_trees:
    input:
        config["root"] + "/" + config["folder"]["mitotree"] + "/msa/" + "{cds}.trim.fa"
    output:
        ML_tree=config["root"] + "/" + config["folder"]["mitotree"] + "/tree/" + "{cds}.treefile",
        bs_trees=config["root"] + "/" + config["folder"]["mitotree"] + "/tree/" + "{cds}.ufboot",
        done=touch(config["root"] + "/" + config["folder"]["mitotree"] + "/tree/" + "{cds}.done")
    conda:
        config["root"] + "/" + config["envs"] + "/" + "buildtree.yaml"
    threads:
        24
    params:
        prefix=config["root"] + "/" + config["folder"]["mitotree"] + "/tree/" + "{cds}",
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["mitotree"] + "/{cds}_tree.log"
    shell:
        """
        iqtree -s {input} --seqtype CODON2 -pre {params.prefix} \
        -nt {threads} -m MFP -msub mitochondrial -bb 1000 -bnni -redo > /dev/null
        """


rule species_tree:
    input:
        done_file=expand(config["root"] + "/" + config["folder"]["mitotree"] + "/tree/" + "{cds}.done", cds=use_for_tree_cds),
        ML_tree=expand(config["root"] + "/" + config["folder"]["mitotree"] + "/tree/" + "{cds}.treefile",cds=use_for_tree_cds),
        map_file=config["root"] + "/" + config["folder"]["mitotree"] + "/species_tree/map.txt"
    output:
        all_gene_tree=config["root"] + "/" + config["folder"]["mitotree"] + "/species_tree/all_gene.tree",
        sample_tree=config["root"] + "/" + config["folder"]["mitotree"] + "/species_tree/sample.tree",
        species_tree=config["root"] + "/" + config["folder"]["mitotree"] + "/species_tree/species.tree",
        done=touch(config["root"] + "/" + config["folder"]["mitotree"] + "/species_tree/astral_hybrid_species.done")
    conda:
        config["root"] + "/" + config["envs"] + "/" + "buildtree.yaml"
    threads:
        24
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["mitotree"] + "/species_tree.log"
    shell:
        """
        cat {input.ML_tree} > {output.all_gene_tree}
        sed -i '/^$/d' {output.all_gene_tree}
        astral-hybrid -t {threads} -S -R -i {output.all_gene_tree} -o {output.sample_tree} 2> /dev/null
        astral-hybrid -t {threads} -S -R -i {output.all_gene_tree} -o {output.species_tree} -a {input.map_file} 2> /dev/null
        """