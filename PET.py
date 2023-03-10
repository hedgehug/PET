import os, sys
import numpy as np
import scipy.stats as stat
import statsmodels.stats.multitest as multitest


def extract_fisher_rank(fisher_result_file):
    print('Calculating fisher rank...')
    file = open(fisher_result_file, 'r')
    file_info = file.readlines()
    result_list = []
    for line in file_info:
        line = line.strip().split('\t')
        # correct p-value
        if float(line[2]) == 1.0:
            pval = 0.9999999
        elif float(line[2]) == 0.0:
            pval = 0.0000001
        else:
            pval = float(line[2])
        # Pathway	#Overlap Gene	p-value	FDR
        result_list.append([line[0].upper(), float(line[col]), pval, line[4]])

    # sort by p-value, ascending
    result_list.sort(key=lambda x: x[2])

def run_PET(fisher_result_file, enrichr_result_file, gsea_result_dir, pathway_dict):
    print()