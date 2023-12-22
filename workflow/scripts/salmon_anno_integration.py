# -*- coding: utf-8 -*-
import os
import re
import math
import argparse


def read_quant_file(quant_file):
    dict_quant = {}
    with open(quant_file, 'r', errors='ignore') as quant:
        next(quant)
        for line in quant:
            line = line.strip().split("\t")
            count = float(line[4])
            if line[0] not in dict_quant.keys():
                dict_quant[line[0]] = count
            else:
                raise Exception('Error: The same gene name in quant file')
    return dict_quant


def read_dbcan_file(dbcan_file, dict_info, only_one=False):
    dict_data = {}
    with open(dbcan_file, 'r', errors='ignore') as ipt:
        next(ipt)
        for line in ipt:
            line = line.strip().split('\t')
            orf_id = line[0].split(' ')[0]
            dbcan_id = ''
            # res为第三列到倒数第二列的内容，去掉括号和括号内的内容
            res = [re.sub(r'\(.*?\)', '', i) for i in line[2:-1]]
            # 去掉_和_后面的数字
            res = [re.sub(r'_[a-z0-9]+', '', i) for i in res]
            # 如果只用了一个软件，那么res中只有一个元素，直接使用他的结果
            if only_one:
                # 移除res里面值为-的元素
                res = [i for i in res if i != '-']
                if res:
                    dbcan_id = res[0]
            # 如果用了多个软件，那么res中有多个元素，判断是否有至少两个软件比对结果一致的term
            elif int(line[-1]) > 1:
                # 移除res里面值为-的元素，拆分带+的term为多个term，
                res = [i for i in res if i != '-']
                all_term = []
                for i in res:
                    all_term.extend(set(i.split('+')))
                # 判断res中是否有至少两个软件比对结果一致的term
                dbcan_id = set([i for i in all_term if all_term.count(i) > 1])
            if orf_id in dict_info and dbcan_id:
                id_count = dict_info[orf_id]/len(dbcan_id)
                for id in dbcan_id:
                    if id not in dict_data:
                        count = id_count
                    else:
                        count = float(dict_data[id][1]) + id_count
                    dict_data[id] = [id, str(math.ceil(count))]
    return dict_data


def read_rgi_file(rgi_file, dict_info):
    dict_data = {}
    with open(rgi_file, 'r', errors='ignore') as ipt:
        for line in ipt:
            line = line.strip().split('\t')
            if line[5]:  # TODO 严格模式还是宽松模式  == 'Strict'
                orf_id = line[0].split(' ')[0]
                rgi_id = line[8]
                if orf_id in dict_info:
                    if rgi_id not in dict_data:
                        count = float(dict_info[orf_id])
                    else:
                        count = float(dict_info[orf_id][1]) + float(dict_data[rgi_id][1])
                    dict_data[rgi_id] = [rgi_id, str(math.ceil(count))]
    return dict_data


def read_vfdb_file(vfdb_file, dict_info):
    with open(vfdb_file, 'r', errors='ignore') as ipt:
        dict_data = {}
        for line in ipt:
            if line.startswith('Gene'): continue
            line = line.strip().split('\t')
            #### 目前提取的是全部结果，阈值在比对时设定，即：下面这个if不起作用
            if line[-1]:
                orf_id = line[0].split(' ')[0]
                if orf_id in dict_info:
                    if line[1] not in dict_data:
                        count = float(dict_info[orf_id])
                    else:
                        count = float(dict_info[orf_id][1]) + float(dict_data[line[1]][1])
                    dict_data[line[1]] = [line[1], str(math.ceil(count))]
    return dict_data


def read_anno_file(anno_file, anno_type, counts_info):
    if anno_type == 'dbcan':
        res_dict = read_dbcan_file(anno_file, counts_info)
    elif anno_type == 'rgi':
        res_dict = read_rgi_file(anno_file, counts_info)
    elif anno_type == 'vfdb':
        res_dict = read_vfdb_file(anno_file, counts_info)
    else:
        raise Exception('Error: Please input the correct type of annotation [rgi/vfdb/dbcan]')
    return res_dict


if __name__ == '__main__':
    """
    将contig定量的count和contig的注释annotation结果整合到一起，得到最终的结果
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-q', '--quant_file', dest='q_file', required=True, help='input the quant result file')
    parser.add_argument('-a', '--anno_file', dest='a_file', required=True, help='input the annotation result file')
    parser.add_argument('-t', '--type', dest='anno_type', required=True,
                        help='input the type of annotation [rgi/vfdb/dbcan]')
    parser.add_argument('-o', '--output', dest='o_file', required=True, help='input the output file')
    args = parser.parse_args()

    # read quant file counts info
    counts_info = read_quant_file(args.q_file)
    # read annotation file and join quant value
    anno_quant_dict = read_anno_file(args.a_file, args.anno_type, counts_info)

    with open(args.o_file, 'w', errors='ignore') as output:
        output.writelines('Sample\t' + os.path.basename(args.q_file).split('.')[0] + '\n')
        for value in anno_quant_dict.values():
            output.writelines('\t'.join(value) + '\n')
