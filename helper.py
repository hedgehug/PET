import numpy as np
import os
import math
import pandas as pd
import pickle as pkl
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import load_example_data
import gseapy as gp


def run_DEseq2(expr_file, contrast_matrix_file, script_file_name=None, result_file_path=None,
               result_dir_path=None, groups=None, method='R', cpu_num=5):
    if method == 'R':
        if script_file_name is None:
            print('Please provide a path for writing DESeq2 script.')
            return
        run_deseq2_r(read_count_file_path=expr_file, contrast_file_path=contrast_matrix_file,
                     result_file_path=result_file_path, result_dir_path= result_dir_path,
                     script_name=script_file_name, groups=groups)
    elif method == 'Python':
        run_deseq2_python(read_count_file_path=expr_file, contrast_file_path=contrast_matrix_file,
                          result_file_path=result_file_path, result_dir_path=result_dir_path,
                          groups=groups, cpu_num=cpu_num)
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
    print('Start running R script '+script_name)
    os.system('Rscript '+script_name)
    print('Analysis Done!')


def run_deseq2_python(read_count_file_path, contrast_file_path,
                      result_file_path, result_dir_path, groups=None, cpu_num=5):
    counts_df = pd.read_csv(read_count_file_path, sep='\t', index_col=0, header=0)
    counts_df = counts_df.dropna()
    counts_df = counts_df.T
    print(counts_df.shape[0], 'samples')
    print(counts_df.shape[1], 'genes')
    samples = counts_df.index
    contrast_df = read_contrast_file(contrast_file_path)
    meta_df = pd.DataFrame()
    meta_df.index = samples
    conds = []
    for idx in contrast_df.index:
        conds.extend([contrast_df.loc[idx]['condition']] * contrast_df.loc[idx]['sample_num'])
    meta_df['condition'] = conds
    print(len(set(meta_df['condition'])), 'conditions')
    inference = DefaultInference(n_cpus=cpu_num)
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=meta_df,
        design_factors="condition",
        refit_cooks=True,
        inference=inference,
        # n_cpus=8, # n_cpus can be specified here or in the inference object
    )
    dds.deseq2()
    all_conds = list(set(meta_df['condition']))
    all_conds.sort()

    if groups == None:
        # compute comparison for all groups
        for idx_1 in range(0, len(all_conds)-1):
            for idx_2 in range(idx_1+1, len(all_conds)):
                print('Comparing ', all_conds[idx_1], '.vs.', all_conds[idx_2])
                stat_res = DeseqStats(dds, inference=inference,
                                      contrast=["condition", all_conds[idx_1], all_conds[idx_2]])
                stat_res.summary()
                stat_res.results_df.to_csv(result_dir_path+'/'+all_conds[idx_1]+'.vs.'+
                                           all_conds[idx_2]+'.PyDESeq2_result.csv')
                print('Results written to ', result_dir_path+'/'+all_conds[idx_1]+'.vs.'+
                                           all_conds[idx_2]+'.PyDESeq2_result.csv')
    else:
        print('Comparing ',  groups[0], '.vs.', groups[1])
        stat_res = DeseqStats(dds, inference=inference, contrast=["condition", groups[0], groups[1]])
        stat_res.summary()
        stat_res.results_df.to_csv(result_file_path)
        print('Results written to ', result_file_path)
    print('Done!')


def run_GSEA(prerank_file_path, pathway_file, out_dir, gsea_cli_path=None, method='Python',
             thread_num=5, min_size=15, max_size=500, permutation_num=1000, seed=42,
             plot=False, gsea_out_label='GSEA_result'):
    if method == 'Python':
        run_GSEA_python(prerank_file_path = prerank_file_path,
                        pathway_file = pathway_file, out_dir=out_dir,
                        plot= plot, min_size=min_size, max_size=max_size,
                        thread_num = thread_num, permutation_num=permutation_num,
                        seed=seed)
    elif method == 'cli':

        if gsea_cli_path == None:
            raise Exception('Please specify file path of GSEA gsea-cli.sh. Otherwise, please keep method=\'Python\'.')
        if not os.path.exists(gsea_cli_path):
            raise Exception(gsea_cli_path+' file does not exist!')

        # call GSEA cmd script
        run_GSEA_cmd(cli_file_path=gsea_cli_path, prerank_file_path = prerank_file_path,
                     pathway_file = pathway_file, out_dir=out_dir,
                     plot= plot, min_size=min_size, max_size=max_size,
                     permutation_num=permutation_num, seed=seed, gsea_out_label=gsea_out_label)
    return


def run_GSEA_cmd(cli_file_path, prerank_file_path, pathway_file, out_dir, gsea_out_label,
                 plot=False, min_size=15, max_size=500, permutation_num=1000, seed=42):
    """
    ~/GSEA_test/GSEA_cmd/gsea-cli.sh GSEAPreranked -gmx example_new/c2.cp.kegg.v2023.1.Hs.symbols.gmt
    -collapse No_Collapse -mode Max_probe -norm meandiv -nperm 1000 -rnk example_new/cond1.vs.cond2.rnk
    -scoring_scheme weighted -rpt_label example_test   -create_svgs false -include_only_symbols true
     -make_sets true -plot_top_x 5 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false
     -out example_new/
    :return:
    """
    print('Start running GSEA!')
    print('Commandas to run GSEA: '+cli_file_path+' GSEAPreranked -gmx '+pathway_file+
              ' -collapse No_Collapse -mode Max_probe -norm meandiv -nperm '+str(permutation_num)
              +' -rnk '+prerank_file_path+' -scoring_scheme weighted -rpt_label '+gsea_out_label
              +' -create_svgs false -include_only_symbols true -make_sets true -plot_top_x 5 -rnd_seed '+str(seed)+
              ' -set_max '+str(max_size)+' -set_min '+str(min_size)+' -zip_report false -out '+out_dir)
    os.system(cli_file_path+' GSEAPreranked -gmx '+pathway_file+
              ' -collapse No_Collapse -mode Max_probe -norm meandiv -nperm '+str(permutation_num)
              +' -rnk '+prerank_file_path+' -scoring_scheme weighted -rpt_label '+gsea_out_label
              +' -create_svgs false -include_only_symbols true -make_sets true -plot_top_x 5 -rnd_seed '+str(seed)+
              ' -set_max '+str(max_size)+' -set_min '+str(min_size)+' -zip_report false -out '+out_dir)


def run_GSEA_python(prerank_file_path, pathway_file, out_dir, plot=False, thread_num=5, min_size=15,
                    max_size=500, permutation_num=1000, seed=42):
    # read rank file
    rnk = pd.read_csv(prerank_file_path, header=None, index_col=0, sep="\t")
    # # run prerank
    # # enrichr libraries are supported by prerank module. Just provide the name
    # # use 4 process to acceralate the permutation speed
    pre_res = gp.prerank(rnk=rnk,  # or rnk = rnk,
                         gene_sets=pathway_file,
                         threads=thread_num,
                         min_size=min_size,
                         max_size=max_size,
                         permutation_num=permutation_num,  # reduce number to speed up testing
                         outdir=out_dir,  # don't write to disk
                         seed=seed,
                         plot = plot,
                         verbose=True)

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
    print('Rank file written to ', out_file)
    file.close()


def extract_top_gene(deseq_result_file, num_gene=200, pval_threshold=0.05, padj_threshold=0.05,
                     fc_threshold=1, basemean_threshold=2, direction=1):
    df = pd.read_csv(deseq_result_file, sep=',', index_col=0, header=0)
    # sort by adjusted p-value
    df = df.sort_values(by='padj')
    df['fc'] = 2**df['log2FoldChange']
    deg_dict = dict()
    selected_df = df[
        (df['baseMean'] >= basemean_threshold) & (df['padj'] <= padj_threshold) & (df['pvalue'] <= pval_threshold)
        & ((df['fc'] >= fc_threshold) | (df['fc'] <= 1/fc_threshold))]
    # print(selected_df)
    if direction == 1:
        deg_dict['up'] = selected_df.loc[(df['log2FoldChange'] > 0)].index.tolist()[:min(selected_df.loc[(df['log2FoldChange'] > 0)].shape[0], num_gene)]
        deg_dict['down'] = selected_df.loc[(df['log2FoldChange'] < 0)].index.tolist()[:min(selected_df.loc[(df['log2FoldChange'] > 0)].shape[0], num_gene)]
    elif direction == -1:
        deg_dict['down'] = selected_df.loc[(df['log2FoldChange'] > 0)].index.tolist()[:min(selected_df.loc[(df['log2FoldChange'] > 0)].shape[0], num_gene)]
        deg_dict['up'] = selected_df.loc[(df['log2FoldChange'] < 0)].index.tolist()[:min(selected_df.loc[(df['log2FoldChange'] > 0)].shape[0], num_gene)]
    else:
        raise Exception(
            'Please specify direction as 1 or -1. 1 means log2FC > 0 indicates up-regulation, '
            '-1 means log2FC < 0 indicates up-regulation')
    print('Based on:')
    print('p-value <= ', pval_threshold)
    print('adjusted p-value <= ', padj_threshold)
    print('fold change threshold >= ', fc_threshold)
    print('base Mean >= ', basemean_threshold)
    print('top N = ', num_gene)
    print(len(deg_dict['up']), 'up-regulated DEGs')
    print(len(deg_dict['down']), 'down-regulated DEGs')
    return deg_dict


def create_dir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)

