import os.path
import sys
import pandas as pd

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
def assign_servers(df, server_dict, total_num):
    # 初始化服务器分配字典
    server_assignments = {server: [[], 0] for server in server_dict}
    # 遍历排序后的DataFrame的每一行
    for _, row in df.iterrows():
        for now_server in server_assignments.keys():
            max_server_num = int(total_num * server_dict[now_server]) + 3  # 最后的数字需要比最小的分组数量大
            if server_assignments[now_server][1] + row['counts'] <= max_server_num:
                server_assignments[now_server][0].append(row['Site_Group'])
                server_assignments[now_server][1] += row['counts']
                break
            else:
                continue
    return server_assignments


server_dict = {'mu13': 0.15, 'shuke': 0.25, 'beijita': 0.25, 'wukong': 0.15, 'huatuo': 0.1, 'mu05': 0.1}
all_df = pd.read_csv('../config/allSampleList.csv', index_col=0)
total_sample_num = all_df.shape[0]

groups = {'Species_Group': 0, 'Site_Group': 0}
# 先计算每种分组类型的分组数
for group in groups.keys():
    groups[group] = all_df.groupby(group).size().shape[0]

# 从分组数最少的分组类型开始分配
groups = sorted(groups.items(), key=lambda x: x[1])
init_df = all_df.groupby(groups[0][0]).size().reset_index(name='counts').sort_values(by='counts', ascending=False)
init_assign = assign_servers(init_df, server_dict, total_sample_num)
print(init_assign)
# 根据初始分配结果，从第二个分组类型开始进行再分配
# 计算按照初始分配结果，每个分配那些已经满足第二种分配类型的
# 子集和全集都按照第二种分配类型进行分组，子集的那些分组的数量，不等于全集的分组数量的为不满足
final_df = pd.DataFrame()
final_incomplete_df = pd.DataFrame()
now_group_cate = groups[1][0]
for server in init_assign:
    # 取出当前服务器已分配的分组
    now_sample_df = all_df[all_df['Site_Group'].isin(init_assign[server][0])]
    full_group_counts = all_df.groupby(now_group_cate).size()
    subset_group_counts = now_sample_df.groupby(now_group_cate).size()
    # 比较子集和全集的分组数量，得到满足的分组和不满足的分组
    complete_groups = subset_group_counts[subset_group_counts.eq(full_group_counts)].index.tolist()
    incomplete_groups = subset_group_counts[~subset_group_counts.eq(full_group_counts)].index.tolist()
    # 将now_sample_df中，不满足的分组，对应的列，设置为None，copy结果到res_df
    res_df = now_sample_df.copy()
    res_df.loc[res_df[now_group_cate].isin(incomplete_groups), now_group_cate] = None
    res_df['Server'] = server
    final_df = pd.concat([final_df, res_df])

    # 统计不满足的分组，在当前子集中的数量，和对应的Sample_ID
    now_sample_df_copy = now_sample_df.copy()
    incomplete_group_df = now_sample_df_copy[now_sample_df_copy[now_group_cate].isin(incomplete_groups)].reset_index()
    # 统计每组的数量和取出每组的所有 Sample_ID
    incomplete_group_info = incomplete_group_df.groupby(now_group_cate).agg({
        'Sample_ID': list,  # 使用 list 将所有 Sample_ID 放在一个列表中
    }).reset_index()
    incomplete_group_info['Group_Count'] = incomplete_group_df.groupby(now_group_cate).size().values
    incomplete_group_info['Server'] = server
    final_incomplete_df = pd.concat([final_incomplete_df, incomplete_group_info])


# 处理第二种分配类型中，不满足的分组，并将较少个样本分配到已有更多的样本的服务器
fi_df = final_incomplete_df.copy().reset_index()
# 找出每个分组中，样本数量最多的行
max_group_count_rows = fi_df.loc[fi_df.groupby('Species_Group')['Group_Count'].idxmax()]
# 找出每个分组中，样本数量不是最多的行
other_group_count_rows = fi_df.loc[~fi_df.index.isin(max_group_count_rows.index)]
# 把other_group_count_rows中的Server替换为max_group_count_rows中Species_Group相同行的Server
nother_group = other_group_count_rows.copy()
nother_group['Server'] = nother_group.apply(
    lambda x: max_group_count_rows.loc[max_group_count_rows['Species_Group'] == x['Species_Group'], 'Server'].values[0],
    axis=1)

final_df['Use_Sample'] = final_df.index
final_df.reset_index(inplace=True)
# 遍历max_group_count_rows，将其中Sample_ID对应的Species_Group填回final_df中，之前不满足全部填为None
for _, row in max_group_count_rows.iterrows():
    for sample_id in row['Sample_ID']:
        final_df.loc[final_df['Sample_ID'] == sample_id, 'Species_Group'] = row['Species_Group']

# 遍历other_group_count_rows，将其中Sample_ID列依次添加到final_df中
for _, row in nother_group.iterrows():
    for sample_id in row['Sample_ID']:
        # 取出all_df中对应的行
        sample_row = all_df.loc[all_df.index == sample_id].copy()
        # 将sample_row中的Server列设置为row['Server']
        sample_row['Server'] = row['Server']
        # 将Site_Group列设置为None
        sample_row['Site_Group'] = None
        sample_row.reset_index(inplace=True)
        # 将sample_row添加到final_df中
        final_df = pd.concat([final_df, sample_row])
# 将final_df按照 'Server', 'Site_Group', 'Species_Group', 'Use_Species', 'Sample_ID' 进行排序
final_df.sort_values(by=['Server', 'Site_Group', 'Species_Group', 'Use_Sample', 'Sample_ID'], inplace=True)
# 将final_df中的列顺序调整
final_df = final_df[['Server', 'Sample_ID', 'Use_Sample', 'Site_Group', 'Species_Group', 'Species', 'Ref','Download_link']]

# TODO 有时候在不同服务器上会产生不一样结果 之后需要修复下
final_df.to_csv('../config/finalSampleList.csv', index=False)

out_to = sys.argv[1]
if out_to in server_dict.keys():
    final_df[final_df['Server'] == out_to].to_csv('../config/sampleList.csv', index=False)
else:
    print('Wrong server name!')
    sys.exit(1)