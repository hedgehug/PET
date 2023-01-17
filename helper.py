import numpy as np
import os


def prune_gmt(file_name, out_file_name, gene_universe, min_gene_num, max_gene_num):
    pathway_dict = read_gmt(file_name)
    out_file = open(out_file_name, 'w')
    pruned_pathway_dict = dict()
    for pathway in pathway_dict.keys():
        genes_shared = set(pathway_dict[pathway]).intersection(set(gene_universe))
        if len(genes_shared) < min_gene_num or len(genes_shared) > max_gene_num:
            print(pathway+' removed due to gene number')
        else:
            out_file.write(pathway+'\tna\t'+'\t'.join(list(genes_shared))+'\n')
            pruned_pathway_dict[pathway] = list(genes_shared)
    return pruned_pathway_dict


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


def create_dir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)

