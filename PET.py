import os, sys
import numpy as np
import scipy.stats as stat
import statsmodels.stats.multitest as multitest
from helper import prune_gmt

def extract_fisher_rank(fisher_result_file):
    print('Calculating fisher rank...')
    file = open(fisher_result_file, 'r')
    file_info = file.readlines()
    pathway_pval_dict = dict()
    pathway_rank_dict = dict()
    pathway_list = []
    pval_list = []
    for line in file_info[1:]:
        # Pathway	#Overlap Gene	p-value	FDR	Rank
        line = line.strip().split('\t')
        # pathway, intersection number/combined score, p-value
        pathway_pval_dict[line[0].upper()] = float(line[2])
        pathway_list.append(line[0].upper())
        pval_list.append(float(line[2]))

    # sort by p-value, ascending
    # for fisher test, there might be lots of tie of pathways
    # here we want to keep the pathways with same p-value/combined score the same rank, for later average calculation
    pval_ranks = stat.rankdata(pval_list, method='min')
    for idx in range(len(pval_ranks)):
        pathway_rank_dict[pathway_list[idx]] = pval_ranks[idx]
    return pathway_rank_dict, pathway_pval_dict


def extract_enrichr_rank(enrichr_result_file):
    print('Calculating enrichr rank...')
    file = open(enrichr_result_file, 'r')
    file_info = file.readlines()
    pathway_pval_dict = dict()
    pathway_rank_dict = dict()
    pathway_list = []
    combined_score_list = []
    for line in file_info[1:]:
        # Pathway	p-value	Combined score	Rank
        line = line.strip().split('\t')
        pathway_pval_dict[line[0].upper()] = float(line[1])
        pathway_list.append(line[0].upper())
        combined_score_list.append(float(line[2]))

    # sort by combined, descending, and keep the ties
    combined_score_ranks = stat.rankdata(combined_score_list, method='min')
    combined_score_ranks = len(combined_score_ranks)-combined_score_ranks+1
    for idx in range(len(combined_score_ranks)):
        pathway_rank_dict[pathway_list[idx]] = combined_score_ranks[idx]
    return pathway_rank_dict, pathway_pval_dict


def process_gsea_report(report_file_name, reverse, start_rank, pathway_pval_dict, pathway_rank_dict):
    file = open(report_file_name, 'r')
    info = file.readlines()[1:]
    if reverse:
        info.reverse()
    for line in info:
        line = line.strip().split('\t')
        # some GSEA p-value might be empty, correct p-value to 1
        try:
            pval = float(line[6])
        except:
            pval = 1.0
        if reverse:
            pathway_pval_dict[line[0]] = 1.0
        else:
            pathway_pval_dict[line[0]] = pval
        pathway_pval_dict[line[0]] = pval
        pathway_rank_dict[line[0]] = start_rank
        start_rank += 1
    file.close()
    return start_rank


def extract_gsea_rank(gsea_result_dir, gsea_label):
    print('Calculating GSEA rank...')
    all_file = os.listdir(gsea_result_dir)
    all_file.sort()

    # GSEA results are stored in two tsv files
    gsea_report_prefix = 'gsea_report_for_'
    gsea_report_suffix = 'tsv'

    gsea_report_list = []

    for f in all_file:
        if f.startswith(gsea_report_prefix) and f.endswith(gsea_report_suffix):
            gsea_report_list.append(gsea_result_dir+f)

    if len(gsea_report_list) != 2:
        print('Did not find the expected 2 GSEA reports. Please check GSEA run. Stopping now.')
        return 0, 0

    # for the report of expected direction, simply take the rank, for the other report, do the reverse way
    pathway_pval_dict = dict()
    pathway_rank_dict = dict()
    if gsea_report_list[0].__contains__(gsea_label):
        pos_report = gsea_report_list[0]
        neg_report = gsea_report_list[1]
    else:
        pos_report = gsea_report_list[1]
        neg_report = gsea_report_list[0]

    start_rank = process_gsea_report(pos_report, False, 1,
                                     pathway_pval_dict, pathway_rank_dict)
    _ = process_gsea_report(neg_report, True, start_rank,
                            pathway_pval_dict, pathway_rank_dict)
    return pathway_rank_dict, pathway_pval_dict


def run_PET(fisher_result_file, enrichr_result_file, gsea_result_dir, gsea_label, pathway_dict, result_file):
    # pathway: value
    rank_dict = dict()
    pval_dict = dict()

    rank_dict['Fisher exact test'], pval_dict['Fisher exact test'] = extract_fisher_rank(fisher_result_file)
    rank_dict['Enrichr'], pval_dict['Enrichr'] = extract_enrichr_rank(enrichr_result_file)
    rank_dict['GSEA'], pval_dict['GSEA'] = extract_gsea_rank(gsea_result_dir, gsea_label)
    if rank_dict['GSEA'] == 0:
        # stop executing, simply return
        return

    methods = list(rank_dict.keys())
    methods.sort()

    print('Ensembling results...')
    pathway_list = list(pathway_dict.keys())
    results = []
    for p in pathway_list:
        tmp_rank_list = []
        tmp_pval_list = []
        for method in methods:
            try:
                tmp_rank_list.append(rank_dict[method][p])
                tmp_pval_list.append(pval_dict[method][p])
            except:
                print(p, 'is missing in', method, '. Skipped.')
                break
        if len(tmp_rank_list) == len(tmp_pval_list) == len(methods):
            # compute average rank
            avg_rank = np.mean(tmp_rank_list)
            # compute combined p-value
            combined_pval = stat.combine_pvalues(tmp_pval_list, method='stouffer')[-1]
            if np.isnan(combined_pval):
                # certain p-value lists cannot be combined, e.g. [1.0, 1.0, 0.0], assign as 1.0 for FDR correction
                combined_pval = 1.0
            results.append([p, avg_rank, combined_pval])
        else:
            continue

    # sort based on average rank
    results.sort(key=lambda x: x[1])

    # calculate FDR
    pval_list = np.array(results)[:, -1].astype(float)
    fdr_list = multitest.fdrcorrection(pval_list, is_sorted=False)[-1]

    # write to results
    out_file = open(result_file, 'w')
    # write header
    out_file.write('Pathway\tPET rank\tPET p-value\tPET FDR\tAverage rank')
    for method in methods:
        out_file.write('\t'+method+' rank\t'+method+' p-value')
    out_file.write('\n')
    for idx in range(len(results)):
        out_file.write('\t'.join([results[idx][0], str(idx+1),
                                  str(results[idx][2]), str(fdr_list[idx]), str(results[idx][1])]))
        for method in methods:
            out_file.write('\t' +str(rank_dict[method][results[idx][0]])+'\t'+str(pval_dict[method][results[idx][0]]))
        out_file.write('\n')
    out_file.close()

    print('Results written to ', result_file)













if __name__ == '__main__':
    pruned_pathway_dict, gene_universe = prune_gmt(file_name='example/c2.cp.kegg.v2023.1.Hs.symbols.gmt',
                                                   out_file_name='example/c2.cp.kegg.v2023.1.Hs.symbols.cleaned.gmt',
                                                   expr_matrix_file='example/example_data.txt',
                                                   min_gene_num=15, max_gene_num=500)
    run_PET('example/fisher_up.txt', 'example/enrichr_up.txt',
            'example/gsea_out_example/example_test.GseaPreranked.1678483249803/', 'pos', pruned_pathway_dict, 'tmp.txt')