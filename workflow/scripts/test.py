import os
import pandas as pd

kraken2_report = ['/home/user002/fish/02_kraken2/Y0894P.kraken']
report_list = ['Eukaryota', 'Bacteria', 'Viruses', 'Archaea', 'unclassified', 'Homo sapiens']
res_dict = {}
for k_file in kraken2_report:
    sample = os.path.basename(k_file).split('.')[0]
    res_dict[sample] = {}
    with open(k_file ,'r') as ipt:
        for line in ipt:
            line = line.strip().split('\t')
            counts = line[1]
            taxon = line[-1].strip()
            if taxon in report_list:
                res_dict[sample][taxon] = counts
table = pd.DataFrame(res_dict).transpose().astype(int)
table['unclassified'] = table['unclassified'].div(table.sum(axis=1)).multiply(100)
table['human'] = table['Homo sapiens'].div(table.sum(axis=1)).multiply(100)
table['known_microbiome'] = 1-table['unclassified']-table['unclassified']
table = table.sort_index()
table.to_csv('123.csv', sep="\t")