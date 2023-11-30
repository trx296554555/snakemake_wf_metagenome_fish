import sys
import os
import pandas as pd

READ_COUNT_IDENTIFIER = "READ COUNT"


def get_read_count_type(line):
    """
    Return the count number and read type from the log line
    """

    data = line.rstrip().split(":")

    count = data[-1].strip()
    type = data[-3].lstrip().rstrip()

    return count, type


def reformat_reads(reads):
    """
    Reformat the reads dictionary to a table
    """
    new_format_reads = {}
    for sample in reads:
        new_format_reads[sample] = {}
        for type in reads[sample]:
            if 'decontaminated' in type:
                ref = type.split(" ")[1]
                new_type = type.split(" ")[0] + ' ' + type.split(" ")[-1]
                new_format_reads[sample][new_type] = reads[sample][type]
                new_format_reads[sample]['reference'] = ref
            else:
                new_format_reads[sample][type] = reads[sample][type]
    table = pd.DataFrame(new_format_reads).transpose()
    table = table[['reference', 'raw pair1', 'raw pair2',
                   'trimmed pair1', 'trimmed pair2', 'final pair1', 'final pair2',
                   'trimmed orphan1', 'trimmed orphan2', 'final orphan1', 'final orphan2']]
    table = table.sort_index()
    return table


def get_reads(file, reads=None):
    """
    Get the read counts from the file
    """
    if not reads:
        reads = {}

    with open(file) as file_handle:
        get_sample = True
        for line in file_handle:
            # get the sample name from the file
            if get_sample and line.startswith("output_dir"):
                sample = line.split("/")[-1].rstrip()
                if sample not in reads:
                    reads[sample] = {}
                get_sample = False
            if READ_COUNT_IDENTIFIER in line:
                count, type = get_read_count_type(line)
                reads[sample][type] = count

    return reads


def main():
    if len(sys.argv) < 3:
        print("Usage: python statistics_kneaddata_res.py file1.log file2.log file3.log... output.file")
        exit()
    # Find the log files
    files = sys.argv[1:]
    out_file = files.pop()
    logs = filter(lambda file: file.endswith(".log") and os.path.isfile(file), files)

    # Get the reads counts for all logs for all samples
    reads = {}
    for file in logs:
        print("Reading log: " + file)
        reads = get_reads(file, reads)

    # Write the output table
    reformat_reads(reads).to_csv(out_file, sep="\t")


if __name__ == "__main__":
    main()
