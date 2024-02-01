import scipy.stats as stats
import statsmodels.stats.multitest as multitest


def run_ora(pathway_dict, deg_dict, gene_universe_num, out_file_prefix, out_dir):
    for item in deg_dict.keys():
        print('Running ORA for genes in ', item)
        run_fisher_test(pathway_dict, deg_dict[item], gene_universe_num, out_file_prefix, out_dir, key=item)


def run_fisher_test(pathway_dict, gene_set, gene_universe_num, out_file_prefix, out_dir, key):
    gene_set = set(gene_set)
    pval_list = []
    intersection_list = []
    pathway_list = list(pathway_dict.keys())
    # print('Start performing fisher test!')
    for pathway in pathway_list:
        # print(pathway)
        intersection_num, pval = run_fisher(gene_set, pathway_dict[pathway], gene_universe_num)
        pval_list.append(pval)
        intersection_list.append(intersection_num)
    fdr_list = multitest.fdrcorrection(pval_list, is_sorted=False)[-1]
    # sort by p-value
    pval_list, fdr_list, pathway_list, intersection_list = \
        (list(t) for t in zip(*sorted(zip(pval_list, fdr_list, pathway_list, intersection_list))))

    out_file = open(out_dir+'/'+out_file_prefix+'.'+key+'.ora_result.txt', 'w')
    out_file.write('Pathway\t#Overlap Gene\tp-value\tFDR\tRank\n')
    rank = 1
    for idx in range(len(pathway_list)):
        out_file.write('\t'.join([pathway_list[idx], str(intersection_list[idx]),
                                  str(pval_list[idx]), str(fdr_list[idx]), str(rank)])+'\n')
        rank += 1
    out_file.close()
    print('Results written to ', out_dir+'/'+out_file_prefix+'.'+key+'.ora_result.txt')
    print('*'*20)


def run_fisher(geneset_1, geneset_2, gene_universe_num):
    geneset_num_1 = len(geneset_1)
    geneset_num_2 = len(geneset_2)
    intersection = len(geneset_1.intersection(geneset_2))
    _, pval = stats.fisher_exact([[intersection, geneset_num_1 - intersection], [geneset_num_2 - intersection,
                                                                                 gene_universe_num - geneset_num_2 - geneset_num_1 + intersection]],
                                 alternative='greater')
    return intersection, pval