import numpy as np
import os
import math


def prune_gmt(file_name, out_file_name, expr_matrix_file, min_gene_num, max_gene_num):
    expr_file = np.loadtxt(expr_matrix_file, delimiter='\t', dtype=str)
    gene_universe = set(expr_file[:, 0][1:])
    print(len(gene_universe), 'genes in the expression matrix')
    pathway_dict = read_gmt(file_name)
    out_file = open(out_file_name, 'w')
    pruned_pathway_dict = dict()

    for pathway in pathway_dict.keys():
        genes_shared = set(pathway_dict[pathway]).intersection(gene_universe)
        if len(genes_shared) < min_gene_num or len(genes_shared) > max_gene_num:
            print(pathway+' removed due to gene number')
        else:
            out_file.write(pathway+'\tna\t'+'\t'.join(list(genes_shared))+'\n')
            pruned_pathway_dict[pathway] = list(genes_shared)
    # return pruned_pathway_dict


def read_gmt(file_name):
    file = open(file_name, 'r')
    info = file.readlines()
    pathway_dict = dict()
    for line in info:
        line = line.strip().split('\t')
        pathway_dict[line[0]] = line[2:]
    file.close()
    return pathway_dict


def format_deseq2_script(read_count_file_path, sample_num, group1_id, group2_id, group1_num, group2_num,
                         result_file_path, script_name):
    template_file = open('deseq2_template.R', 'r')
    template_file_info = template_file.read()
    template_file.close()
    out_script = open(script_name, 'w')
    tmp_info = template_file_info.replace('data_matrix_file_name', read_count_file_path)
    tmp_info = tmp_info.replace('group1_id', group1_id)
    tmp_info = tmp_info.replace('group2_id', group2_id)
    tmp_info = tmp_info.replace('group1_num', str(group1_num))
    tmp_info = tmp_info.replace('group2_num', str(group2_num))
    tmp_info = tmp_info.replace('Sample_num', str(sample_num))
    tmp_info = tmp_info.replace('OUTPUT.csv', result_file_path)
    out_script.write(tmp_info)
    out_script.close()
    print('DESeq2 script written to '+script_name)


# generate .rnk file for GSEA prerank from DEseq2 result
def generate_rank_file(deseq_result_file, out_file, direction=1):
    file = open(deseq_result_file, 'r')
    info = file.readlines()[1:]
    up_list = []  # up in Tumor
    down_list = []  # down in Normal
    middle_list = []
    for line in info:
        line = line.strip().split(',')
        if line[2] == 'NA' or line[-2] == 'NA' or float(line[2]) == 0.0:
            middle_list.append([0, line[0].replace('\"', '')])
        elif float(line[2])*direction < 0:
            if float(line[-2]) == 0.0:
                # assign a significance for p-value=0
                down_list.append([-350, line[0].replace('\"', '')])
            else:
                down_list.append([math.log10(float(line[-2])), line[0].replace('\"', '')])
        elif float(line[2])*direction > 0:
            if float(line[-2]) == 0.0:
                # assign a significance for p-value=0
                up_list.append([350, line[0].replace('\"', '')])
            else:
                up_list.append([-math.log10(float(line[-2])), line[0].replace('\"', '')])
    up_list.sort(key=lambda x: x[0])
    up_list.reverse()
    down_list.sort(key=lambda x: x[0])
    down_list.reverse()
    final_list = []
    final_list.extend(up_list)
    final_list.extend(middle_list)
    final_list.extend(down_list)
    rank_file = open(out_file, 'w')
    for item in final_list:
        rank_file.write(item[1] + '\t' + str(item[0]) + '\n')
    rank_file.close()

    file.close()


def create_dir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)

