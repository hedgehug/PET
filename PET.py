import os, sys
import numpy as np
import scipy.stats as stat
import statsmodels.stats.multitest as multitest
from helper import prune_gmt
import pandas as pd


def extract_ora_rank(ora_result_file):
    print('Calculating ora rank...')
    file = open(ora_result_file, 'r')
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
    # for ora test, there might be lots of tie of pathways
    # here we want to keep the pathways with same p-value/combined score the same rank, for later average calculation
    pval_ranks = stat.rankdata(pval_list, method='average')
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
    combined_score_ranks = stat.rankdata(combined_score_list, method='average')
    combined_score_ranks = len(combined_score_ranks)-combined_score_ranks+1
    for idx in range(len(combined_score_ranks)):
        pathway_rank_dict[pathway_list[idx]] = combined_score_ranks[idx]
    return pathway_rank_dict, pathway_pval_dict


def extract_gsea_rank(gsea_result_file, gsea_label='up'):
    print('Calculating GSEA rank...')
    pathway_pval_dict = dict()
    pathway_rank_dict = dict()
    gsea_df = pd.read_csv(gsea_result_file, sep=',', header=0, index_col=0)
    if gsea_label == 'up':
        # positive NES first
        gsea_df = gsea_df.sort_values(by='NES', ascending=False)
    else:
        # negative NES first
        gsea_df = gsea_df.sort_values(by='NES', ascending=True)

    for idx in range(gsea_df.shape[0]):
        pathway_rank_dict[gsea_df.iloc[idx][0]] = idx+1
        pathway_pval_dict[gsea_df.iloc[idx][0]] = gsea_df.iloc[idx]['NOM p-val']
    return pathway_rank_dict, pathway_pval_dict


def run(ora_result_file, enrichr_result_file, gsea_result_file, out_dir, file_prefix, pathway_dict, direction):
    # pathway: value
    rank_dict = dict()
    pval_dict = dict()

    rank_dict['ora'], pval_dict['ora'] = extract_ora_rank(ora_result_file)
    rank_dict['Enrichr'], pval_dict['Enrichr'] = extract_enrichr_rank(enrichr_result_file)
    rank_dict['GSEA'], pval_dict['GSEA'] = extract_gsea_rank(gsea_result_file, gsea_label=direction)

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
    out_file = open(out_dir+'/'+file_prefix+'.'+direction+'.PET.txt', 'w')
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

    print('Results written to ', out_dir+'/'+file_prefix+'.'+direction+'.PET.txt')


def combine_gsea_cmd_results(gsea_dir, gsea_label, key):
    # GSEA results are stored in two tsv/csv files
    gsea_report_prefix = 'gsea_report_for_'
    # gsea_report_suffix = 'tsv'

    gsea_report_list = []
    all_file = os.listdir(gsea_dir)
    all_file.sort()
    for f in all_file:
        if f.startswith(gsea_report_prefix) and not f.endswith('html'):
            gsea_report_list.append(gsea_dir + f)

    if len(gsea_report_list) != 2:
        print('Did not find the expected 2 GSEA reports. Please check GSEA run. Stopping now.')
        return 0, 0

    # for the report of expected direction, simply take the rank, for the other report, do the reverse way
    if gsea_report_list[0].__contains__(gsea_label):
        pos_report = gsea_report_list[0]
        neg_report = gsea_report_list[1]
    else:
        pos_report = gsea_report_list[1]
        neg_report = gsea_report_list[0]

    df_1 = pd.read_csv(pos_report, sep='\t')
    df_2 = pd.read_csv(neg_report, sep='\t')
    # for the neg one, correct p-value to 1
    df_2['NOM p-val'] = [1] * df_2.shape[0]
    new_df = pd.concat([df_1, df_2])
    new_df.insert(0, "Name", ['Gsea']*new_df.shape[0])
    # write to a tmp file, delete later
    new_df.to_csv(gsea_dir+'/tmp_result.'+key+'.txt', sep=',', index=False)
    return gsea_dir+'/tmp_result.'+key+'.txt'


def run_PET(file_prefix, deg_dict, pathway_dict, result_dir, gsea_path, out_dir, gsea_pos_label='pos'):
    # identify files relevant
    files = dict()
    all_files = os.listdir(result_dir)
    all_files.sort()
    for key in deg_dict.keys():
        files[key] = dict()
        for f in all_files:
            if f.__contains__(key) and f.__contains__('ora_result') and f.startswith(file_prefix):
                files[key]['ora'] = result_dir+f
            elif f.__contains__(key) and f.__contains__('enrichr_result') and f.startswith(file_prefix):
                files[key]['enrichr'] = result_dir + f
            else:
                continue
    tmp_file = []
    if os.path.isdir(gsea_path):
        print('GSEA path was provided as directory, command line version used.')
        # combine GSEA results, make it GSEAPY output like
        if gsea_pos_label == 'pos':
            files['up']['gsea'] = combine_gsea_cmd_results(gsea_dir=gsea_path, gsea_label=gsea_pos_label, key='up')
            files['down']['gsea'] = combine_gsea_cmd_results(gsea_dir=gsea_path, gsea_label='neg', key='down')
        elif gsea_pos_label == 'neg':
            files['up']['gsea'] = combine_gsea_cmd_results(gsea_dir=gsea_path, gsea_label=gsea_pos_label, key='up')
            files['down']['gsea'] = combine_gsea_cmd_results(gsea_dir=gsea_path, gsea_label='pos', key='down')
        else:
            print('Please provide a valid direction key to identify GSEA report. Options: pos, neg. '
                  'The key should indicate up-regulation in the condition.')
        tmp_file.append(files['up']['gsea'])
        tmp_file.append(files['down']['gsea'])
    else:
        print('GSEA path was provided as file, GSEAPY used.')
        # correct the p-value and delete the tmp files later
        tmp_df = pd.read_csv(gsea_path, sep=',', header=0, index_col=0)
        tmp_df.loc[tmp_df['NES'] <= 0, 'NOM p-val'] = 1
        tmp_df.to_csv(gsea_path.replace('.csv', '.up.tmp.csv'), sep=',', index=True)
        # GSEA might have NA p-value, replace with 1
        tmp_df['NOM p-val'].fillna(1, inplace=True)
        files['up']['gsea'] = gsea_path.replace('.csv', '.up.tmp.csv')
        tmp_file.append(files['up']['gsea'])
        tmp_df = pd.read_csv(gsea_path, sep=',', header=0, index_col=0)
        tmp_df.loc[tmp_df['NES'] >= 0, 'NOM p-val'] = 1
        # GSEA might have NA p-value, replace with 1
        tmp_df['NOM p-val'].fillna(1, inplace=True)
        tmp_df.to_csv(gsea_path.replace('.csv', '.down.tmp.csv'), sep=',', index=True)
        files['down']['gsea'] = gsea_path.replace('.csv', '.down.tmp.csv')
        tmp_file.append(files['down']['gsea'])

    # combine the results
    for key in deg_dict.keys():
        print('For genes in ', key, 'direction, result files found: ')
        print(files[key]['ora'], '\n', files[key]['enrichr'], '\n', files[key]['gsea'])
        run(ora_result_file=files[key]['ora'], enrichr_result_file= files[key]['enrichr'],
            gsea_result_file=files[key]['gsea'], out_dir=out_dir, file_prefix = file_prefix,
            pathway_dict=pathway_dict, direction=key)

    # remove the tmp files created for GSEA
    for f in tmp_file:
        os.remove(f)

