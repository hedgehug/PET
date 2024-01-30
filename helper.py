import numpy as np
import os
import math
import pandas as pd
import pickle as pkl
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import load_example_data


def run_DEseq2(expr_file, contrast_matrix_file, script_file_name, result_file_path=None, result_dir_path=None, groups=None, method='R'):
    if method == 'R':
        if script_file_name is None:
            print('Please provide a path for writing DESeq2 script.')
            return
        run_deseq2_r(read_count_file_path=expr_file, contrast_file_path=contrast_matrix_file,
                     result_file_path=result_file_path, result_dir_path= result_dir_path,
                     script_name=script_file_name, groups=groups)
    elif method == 'Python':
        run_deseq2_python(read_count_file_path=expr_file, contrast_file_path=contrast_matrix_file,
                          result_file_path = result_file_path)
    else:
        print('Please specify a valid method for running DESeq2. Options: R, Python.')
        return
    return


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
    return pruned_pathway_dict, gene_universe


def read_gmt(file_name):
    file = open(file_name, 'r')
    info = file.readlines()
    pathway_dict = dict()
    for line in info:
        line = line.strip().split('\t')
        pathway_dict[line[0]] = line[2:]
    file.close()
    return pathway_dict


def read_contrast_file(file_name):
    contrast_df = pd.read_csv(file_name, sep='\t', header=None)
    contrast_df.columns = ['condition', 'sample_num']
    return contrast_df


def run_deseq2_r(read_count_file_path, contrast_file_path, script_name, result_file_path=None, result_dir_path=None, groups=None):
    # read contrast file
    contrast_df = read_contrast_file(contrast_file_path)
    sample_num = contrast_df['sample_num'].sum()
    # check sample number with expression file
    header = np.loadtxt(read_count_file_path, delimiter='\t', dtype=str)[0][1:]
    if len(header) != sample_num:
        raise Exception('Found ' +str(sample_num)+' samples in contrast files and '+str(len(header))+
                        ' samples in gene expression file.\n Please double check the number, stopping now!')

    template_file = open('deseq2_template.R', 'r')
    template_file_info = template_file.read()
    template_file.close()
    out_script = open(script_name, 'w')
    tmp_info = template_file_info.replace('data_matrix_file_name', read_count_file_path)
    # format condition <- factor(c(rep("group1_id", group1_num), rep("group2_id", group2_num), ...))
    cond_rep_str = []
    cond_rep_example = 'rep(\"cond_name\", cond_num)'
    for idx in contrast_df.index:
        cond_rep_str.append(cond_rep_example.replace('cond_name', contrast_df.loc[idx]['condition']).replace('cond_num', str(contrast_df.loc[idx]['sample_num'])))
    cond_rep_str = ','.join(cond_rep_str)
    tmp_info = tmp_info.replace('COND_REP', str(cond_rep_str))
    tmp_info = tmp_info.replace('Sample_num', str(sample_num))

    if groups == None:
        # save all comparison results
        # needs a dir path to save all results file
        if result_dir_path == None:
            out_script.close()
            raise Exception('Please provide a result_dir_path to save DESeq2 analysis result files.')
        tmp_info = tmp_info.replace('OUT_DIR', result_dir_path)
    else:
        if result_file_path == None:
            out_script.close()
            raise Exception('Please provide a result_file_path to save DESeq2 analysis result file.')
        # save single comparison
        tmp_info = tmp_info.replace('COND1_name', groups[0])
        tmp_info = tmp_info.replace('COND2_name', groups[1])
        tmp_info = tmp_info.replace('SAVE_ALL <- TRUE', 'SAVE_ALL <- FALSE')
        tmp_info = tmp_info.replace('OUTPUT.csv', result_file_path)

    out_script.write(tmp_info)
    out_script.close()
    print('DESeq2 script written to '+script_name)
    print('Start running R script')
    os.system('Rscript '+script_name)
    print('Analysis Done!')


def run_deseq2_python(read_count_file_path, contrast_file_path, result_file_path):
    return


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


def extract_top_gene(rank_file, num_gene, direction=1):
    rank_file = np.loadtxt(rank_file, delimiter='\t', dtype=str)
    if direction == 1:
        return rank_file[:num_gene, 0]
    elif direction == -1:
        return rank_file[-num_gene:, 0]


def create_dir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)

