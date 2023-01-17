import scipy.stats as stats
import sys, random, math
import numpy as np
import os
import statsmodels.stats.multitest as multitest


def extract_deg_file(file_name):
    file = open(file_name, 'r')
    file_info = file.readlines()
    header = []
    raw_data = []
    for line in file_info:
        line = line.strip().split()
        header.append(line[0])
        raw_data.append(line[2:])
    raw_data = np.array(raw_data)
    deg_dict = dict()
    for index in range(0, len(header)):
        deg_dict[header[index]] = set(raw_data[index])
    file.close()
    return header, deg_dict


def generate_random_geneset(gene_num, permutation_num, file_name):
    all_random_set = []
    for idx in range(0, permutation_num):
        tmp_geneset = random.sample(gene_universe, gene_num)
        all_random_set.append(set(tmp_geneset))
        # intersection, fisher_pval = run_fisher(tmp_geneset, pathway_geneset)

    # pathway: [all rank in all permutations]
    permutation_rank_dict = dict()
    for idx in range(0, permutation_num):
        print('Permutation '+str(idx))
        tmp_rank_list = []
        for pathway in pathway_gmt_header:
            intersection, fisher_pval = run_fisher(all_random_set[idx], pathway_gmt_dict[pathway])
            tmp_rank_list.append([pathway, intersection, fisher_pval])
        # sort by p-value
        tmp_rank_list.sort(key=lambda x: x[2])
        # tmp_rank_list.reverse()
        rank = 1
        for item in tmp_rank_list:
            try:
                permutation_rank_dict[item[0]].append(rank)
            except:
                permutation_rank_dict[item[0]] = []
                permutation_rank_dict[item[0]].append(rank)
            rank += 1

    out_file = open(file_name, 'w')
    for pathway in permutation_rank_dict.keys():
        tmp_list = permutation_rank_dict[pathway]
        tmp_list.sort()
        out_file.write(pathway+'\t'+'\t'.join([str(i) for i in tmp_list])+'\n')
    out_file.close()


def run_fisher(geneset_1, geneset_2):
    # input gene list
    geneset_num_1 = len(geneset_1)
    # pathway
    geneset_num_2 = len(geneset_2)
    intersection = len(geneset_1.intersection(geneset_2))
    _, fisher_pval = stats.fisher_exact([[intersection, geneset_num_1 - intersection], [geneset_num_2 - intersection,
                                        gene_universe_num - geneset_num_2 - geneset_num_1 + intersection]],
                                        alternative='greater')
    return intersection, fisher_pval


def read_permutation_file(permutation_file):
    file = open(permutation_file, 'r')
    file_info = file.readlines()
    pathway_mean_rank_dict = dict()
    pathway_std_rank_dict = dict()
    for line in file_info:
        line = line.strip().split()
        pathway_mean_rank_dict[line[0]] = np.mean(np.array(line[1:]).astype(np.int))
        pathway_std_rank_dict[line[0]] = np.std(np.array(line[1:]).astype(np.int))
    file.close()
    return pathway_mean_rank_dict, pathway_std_rank_dict


def run_enrichr(geneset_1):
    # generate random gene lists
    if not os.path.exists(enrichr_permutation_dir):
        os.mkdir(enrichr_permutation_dir)

    input_len = len(geneset_1)
    permutation_num = 1000
    if pathway_gmt.__contains__('.noise'):
        permutation_file = enrichr_permutation_dir + '/' + \
                           pathway_gmt.replace('/', '_').replace('.gmt', '').replace('.gmx', '') + '_size' + \
                           str(input_len) + '_num' + str(permutation_num)
    else:
        permutation_file = enrichr_permutation_dir + '/'+pathway_gmt.split('.')[0].replace('/','_')+'_size' + str(input_len)+'_num'+str(permutation_num)

    if not os.path.exists(permutation_file):
        print('Initializing permutation...')
        generate_random_geneset(input_len, permutation_num, permutation_file)

    pathway_mean_rank_dict, pathway_std_rank_dict = read_permutation_file(permutation_file)

    # run fisher exact test first
    all_fisher_result = []
    pval_dict = dict()
    for pathway in pathway_gmt_header:
        intersection, fisher_pval = run_fisher(geneset_1, pathway_gmt_dict[pathway])
        pval_dict[pathway] = fisher_pval
        all_fisher_result.append([pathway, intersection, fisher_pval])
    # sort by p-value
    all_fisher_result.sort(key=lambda x: x[2])

    rank = 1
    rank_dict = dict()
    for item in all_fisher_result:
        rank_dict[item[0]] = rank
        rank += 1

    # compute z-score
    combined_score_dict = dict()
    for pathway in pathway_gmt_header:
        z_score = (rank_dict[pathway]-pathway_mean_rank_dict[pathway])/pathway_std_rank_dict[pathway]
        combined_score_dict[pathway] = abs(z_score*np.log10(pval_dict[pathway]))
    return pval_dict, combined_score_dict


def fisher_test():
    file_idx = 0
    while file_idx < len(ranked_file_header):
        print(ranked_file_header[file_idx])
        pos_ranked_gene = ranked_file_dict[ranked_file_header[file_idx]]
        # neg_ranked_gene = ranked_file_dict[ranked_file_header[file_idx + 1]]

        pos_all_result = []
        pos_pval_list = []

        pos_pval_dict, pos_combined_score_dict = run_enrichr(pos_ranked_gene)

        for pathway in pathway_gmt_header:
            pos_all_result.append([pathway, pos_combined_score_dict[pathway], str(pos_pval_dict[pathway]), 'pos'])
            pos_pval_list.append(pos_pval_dict[pathway])

        pos_fdr_list = multitest.fdrcorrection(pos_pval_list, is_sorted=False)[-1]

        pos_file = open(out_dir + '/' + ranked_file_header[file_idx] + '_enrichr_num', 'w')
        for idx in range(0, len(pos_all_result)):
            tmp = [str(i) for i in pos_all_result[idx]]
            # pathway, combined_score, p-value, fdr, sign
            pos_file.write('\t'.join(tmp[:-1]) + '\t' + str(pos_fdr_list[idx]) + '\t' + tmp[-1] + '\n')
        pos_file.close()

        # sort by combined score
        pos_all_result.sort(key=lambda x: x[1])
        pos_all_result.reverse()
        rank = 1
        for res in pos_all_result:
            if res[0] == ranked_file_header[file_idx]:
                break
            rank += 1
        all_result.append([ranked_file_header[file_idx], str(rank)])
        file_idx += 1


def get_gene_universe():
    file = open(gene_universe_file_name, 'r')
    file_info = file.readlines()[1:]
    for line in file_info:
        line = line.strip()
        gene_universe.add(line)
    file.close()


if __name__ == '__main__':
    ranked_file = sys.argv[1]
    pathway_gmt = sys.argv[2]
    enrichr_permutation_dir = sys.argv[4]
    gene_universe_file_name = sys.argv[5]

    out_dir = sys.argv[3] + '/'
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    gene_universe = set()
    get_gene_universe()

    gene_universe_num = len(gene_universe)

    ranked_file_header, ranked_file_dict = extract_deg_file(ranked_file)
    pathway_gmt_header, pathway_gmt_dict = extract_deg_file(pathway_gmt)

    all_result = []

    fisher_test()

    out_file = open(out_dir+'/enrichr_all_result', 'w')
    for res in all_result:
        out_file.write(res[0]+'\t'+res[1]+'\n')
    out_file.close()
    print(out_dir+'/enrichr_all_result')
