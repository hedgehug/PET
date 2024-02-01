import scipy.stats as stats
import sys, random, math
import numpy as np
import os


random.seed(42)


def run_enrichr(pathway_dict, deg_dict, gene_universe, permutation_num, permutation_file_name,
                out_file_prefix, out_dir):
    for item in deg_dict.keys():
        print('Running Enrichr for genes in ', item)
        run(pathway_dict=pathway_dict, gene_set=deg_dict[item],
            gene_universe=gene_universe, permutation_num=permutation_num,
            permutation_file_name=permutation_file_name,
            out_file_prefix=out_file_prefix, out_dir=out_dir, key=item)


def run(pathway_dict, gene_set, gene_universe, permutation_num,
        permutation_file_name, out_file_prefix, out_dir, key):
    if not os.path.exists(permutation_file_name):
        print(permutation_file_name, 'does not exist, initializing permutation')
        generate_permutation(gene_universe, len(gene_set),
                             permutation_num, pathway_dict, permutation_file_name)

        print('Permutation result written to ', permutation_file_name)
        pathway_mean_rank_dict, pathway_std_rank_dict = read_permutation_file(permutation_file_name)
    else:
        print(permutation_file_name, 'already exists, reading file now...')
        pathway_mean_rank_dict, pathway_std_rank_dict = read_permutation_file(permutation_file_name)

    # run fisher exact test
    all_fisher_result = []
    pval_dict = dict()
    pathway_names = list(pathway_dict.keys())
    for pathway in pathway_names:
        intersection, fisher_pval = run_fisher(set(gene_set) , pathway_dict[pathway], len(gene_universe))
        pval_dict[pathway] = fisher_pval
        all_fisher_result.append([pathway, intersection, fisher_pval])

    # sort by p-value
    all_fisher_result.sort(key=lambda x: x[2])

    # compute the original rank
    rank = 1
    rank_dict = dict()
    for item in all_fisher_result:
        rank_dict[item[0]] = rank
        rank += 1

    # compute z-score
    combined_score_dict = dict()
    combined_score_list = []
    for pathway in pathway_names:
        z_score = (rank_dict[pathway] - pathway_mean_rank_dict[pathway]) / pathway_std_rank_dict[pathway]
        combined_score_dict[pathway] = abs(z_score * np.log10(pval_dict[pathway]))
        combined_score_list.append(combined_score_dict[pathway])

    # sort combined score, descending
    combined_score_list, pathway_names = (list(t) for t in zip(*sorted(
        zip(combined_score_list, pathway_names), reverse=True)))

    # fdr_list = multitest.fdrcorrection(pval_list, is_sorted=False)[-1]

    out_file = open(out_dir+'/'+out_file_prefix+'.'+key+'.enrichr_result.txt', 'w')
    out_file.write('Pathway\tp-value\tCombined score\tRank\n')
    rank = 1
    for idx in range(len(pathway_names)):
        out_file.write('\t'.join([pathway_names[idx], str(pval_dict[pathway_names[idx]]),
                                  str(combined_score_list[idx]), str(rank)]) + '\n')
        rank +=1
    out_file.close()
    print('Results written to ', out_dir+'/'+out_file_prefix+'.'+key+'.enrichr_result.txt')
    print('*' * 20)


def generate_permutation(gene_universe, gene_num, permutation_num, pathway_dict, permutation_file_name):
    all_random_set = []
    # sample N genes for perm_num rounds
    for idx in range(0, permutation_num):
        tmp_geneset = random.sample(gene_universe, gene_num)
        all_random_set.append(set(tmp_geneset))

    # pathway: [all rank in all permutations]
    permutation_rank_dict = dict()
    for idx in range(0, permutation_num):
        # print('Permutation ' + str(idx))
        if idx%50 == 0:
            print('Permutation ' + str(idx))
        tmp_rank_list = []
        for pathway in list(pathway_dict.keys()):
            intersection, fisher_pval = run_fisher(all_random_set[idx], pathway_dict[pathway], len(gene_universe))
            tmp_rank_list.append([pathway, intersection, fisher_pval])
        # sort by p-value
        tmp_rank_list.sort(key=lambda x: x[2])
        rank = 1
        for item in tmp_rank_list:
            try:
                permutation_rank_dict[item[0]].append(rank)
            except:
                permutation_rank_dict[item[0]] = []
                permutation_rank_dict[item[0]].append(rank)
            rank += 1

    out_file = open(permutation_file_name, 'w')
    for pathway in permutation_rank_dict.keys():
        tmp_list = permutation_rank_dict[pathway]
        tmp_list.sort()
        out_file.write(pathway + '\t' + '\t'.join([str(i) for i in tmp_list]) + '\n')
    out_file.close()


def run_fisher(geneset_1, geneset_2, gene_universe_num):
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
        line = line.strip().split('\t')
        pathway_mean_rank_dict[line[0]] = np.mean(np.array(line[1:]).astype(int))
        pathway_std_rank_dict[line[0]] = np.std(np.array(line[1:]).astype(int))
    file.close()
    return pathway_mean_rank_dict, pathway_std_rank_dict