import subprocess
import os
import argparse

"""
This script extracts the longest sequence from each cluster in the cluster file.
Arguments:
    -i: Input fasta file
    -o: Output file
    -c: Cluster file
Usage:
    python3 get_longest_from_cluster.py -i <input_fasta_file> -o <output_file> -c <cluster_file>
Author:
    Robin Tang
Date:
    2024-3-27
"""

longest_extractor = argparse.ArgumentParser(description='Usage: python3 get_longest_from_cluster.py -i <input_fasta_file> -o <output_file> -c <cluster_file>')
longest_extractor.add_argument('-i', '--input', help='Input fasta file', required=True)
longest_extractor.add_argument('-o', '--output', help='Output file', required=True)
longest_extractor.add_argument('-c', '--cluster', help='Cluster file', required=True)
args = longest_extractor.parse_args()

input_fasta = args.input
output_file = args.output
cluster_file = args.cluster

# Check seqkit command
try:
    subprocess.check_output("seqkit -h 2>/dev/null", shell=True)
except Exception as e:
    print("Error: seqkit is not installed.")
    exit()

# Use seqkit to get a list of IDs sorted in reverse order of length for each sequence in the Input Fasta file
listing = f'seqkit sort --quiet -l -r {input_fasta} | seqkit seq -n'
listing_res = subprocess.check_output(listing, shell=True)
len_sort_list = listing_res.decode().strip().split('\n')

# Read the cluster file
cluster_list = []
with open(cluster_file, 'r') as f:
    for cluster_line in f:
        each_cluster = cluster_line.strip().split()
        cluster = list(set(each_cluster))
        cluster_list.append(cluster)

# Find the sequence ID that appears first in the len sort list in each cluster
longest_seq_list = []
for cluster in cluster_list:
    if len(cluster) == 1:
        longest_seq_list.append(cluster[0])
        continue
    else:
        now_index = 999999999
        for seq in cluster:
            if len_sort_list.index(seq) < now_index:
                now_index = len_sort_list.index(seq)
                longest_seq = seq
        longest_seq_list.append(longest_seq)

longest_seq_list = list(set(longest_seq_list))
with open(os.path.join(os.path.dirname(output_file), 'longest_seq_list.txt'), 'w') as f:
    f.write('\n'.join(longest_seq_list))

# Extract the longest sequence from the input fasta file
extract_cmd = f'seqkit grep --quiet -w 0 -n -f {os.path.join(os.path.dirname(output_file), "longest_seq_list.txt")} {input_fasta} > {output_file}'
subprocess.check_output(extract_cmd, shell=True)

# Delete temporary files
os.remove(os.path.join(os.path.dirname(output_file), 'longest_seq_list.txt'))
print("Done!")