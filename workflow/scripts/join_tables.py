#!/usr/bin/env python
import argparse
import sys
import os
import re
import pandas as pd


"""
Join a set of gene, pathway, or taxonomy tables into a single table
This module will join gene and pathway tables output by salmon_anno_integration.py 
To Run: 
$ ./join_tables.py -i <input_dir> -o <gene_table.tsv>
"""


GENE_TABLE_DELIMITER = "\t"  # the delimiter used in the gene table
GENE_TABLE_COMMENT_LINE = "#"  # the line which indicates the header


def join_gene_tables(gene_tables, output, verbose=None):
    """
    Join the gene tables to a single gene table
    """
    data_df = None
    for gene_table in gene_tables:
        if verbose:
            print("Reading gene table: " + gene_table)
        current_df = pd.read_table(gene_table, comment=GENE_TABLE_COMMENT_LINE, index_col=0)
        if data_df is None:
            data_df = current_df
        else:
            data_df = data_df.merge(current_df, left_index=True, right_index=True, how="outer")
    data_df = data_df.fillna(0)
    data_df.to_csv(output, sep=GENE_TABLE_DELIMITER, index_label="Sample")


def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """

    parser = argparse.ArgumentParser(
        description="Join gene, pathway, or taxonomy tables\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-v", "--verbose",
        help="additional output is printed\n",
        action="store_true",
        default=False)
    parser.add_argument(
        "-i", "--input",
        help="the directory of tables\n",
        required=True)
    parser.add_argument(
        "-o", "--output",
        help="the table to write\n",
        required=True)
    parser.add_argument(
        "--file_name",
        help="only join tables with this string included in the file name")
    parser.add_argument(
        "-s", "--search-subdirectories",
        help="search sub-directories of input folder for files\n",
        action="store_true",
        default=False)

    return parser.parse_args()


def main():
    # Parse arguments from command line
    args = parse_arguments(sys.argv)

    # check for format of the gene tables
    input_dir = os.path.abspath(args.input)

    # check the directory exists
    if not os.path.isdir(input_dir):
        sys.exit("The input directory provided can not be found." +
                 "  Please enter a new directory.")

    gene_tables = []
    file_list = os.listdir(input_dir)
    # add in files in subdirectories, if set
    if args.search_subdirectories:
        for possible_folder in os.listdir(input_dir):
            if os.path.isdir(os.path.join(input_dir, possible_folder)):
                try:
                    file_list += [os.path.join(possible_folder, file) for file in
                                  os.listdir(os.path.join(input_dir, possible_folder))]
                except EnvironmentError:
                    pass

    # filter out files which do not meet the name requirement if set
    reduced_file_list = []
    if args.file_name:
        for file in file_list:
            if re.search(args.file_name, file):
                reduced_file_list.append(file)
    else:
        for file in file_list:
            # ignore dot files, like ".DS_Store" on Apple OS X
            if file[0] == ".":
                if args.verbose:
                    print("Not including file in input folder: " + file)
            else:
                reduced_file_list.append(file)

    file_list = reduced_file_list
    args.output = os.path.abspath(args.output)

    for file in file_list:
        gene_tables.append(os.path.join(input_dir, file))
    # sort the gene tables so they are in the same order on all platforms
    gene_tables.sort()
    # split the gene table
    if gene_tables:
        if args.verbose:
            print("Joining gene table")
        join_gene_tables(gene_tables, args.output, verbose=args.verbose)
        print("Gene table created: " + args.output)
    else:
        print("Zero gene tables were found to join.")


if __name__ == "__main__":
    main()
