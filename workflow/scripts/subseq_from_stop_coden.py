import sys
# python {params.script} ${{stopcode}} {input.refine_nt} {output.tmp_code_nt}

stop_coden_list = sys.argv[1].split(';')
stop_coden_list = [(int(i)-1)*3 for i in stop_coden_list]
# stop_coden_list = [1596, 1602, 1605, 1626, 1641]
nt_fa = sys.argv[2]
out_fa = sys.argv[3]


with open(nt_fa, 'r') as f:
    with open(out_fa, 'w') as out:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                out.write(line+'\n')
            else:
                # stop_coden 是开始的位置, 删除原序列 从stop_coden到stop_coden+3的区域
                last_stop_coden = 0
                for stop_coden in stop_coden_list:
                    out.write(line[last_stop_coden:stop_coden])
                    last_stop_coden = stop_coden+3
                out.write(line[last_stop_coden:]+ '\n')