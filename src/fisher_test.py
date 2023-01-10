import scipy.stats as stats
import sys
import numpy as np
import os


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
        deg_dict[header[index]] = raw_data[index]
    file.close()
    return header, deg_dict


def run_fisher_test():
    rank_list = []
    expected_pval_list = []
    all_result = []
    for target in deg_dict_1.keys():
        pval_list = []
        pval_value_list = []
        gene_set_1 = set(deg_dict_1[target])
        gene_set_num_1 = len(gene_set_1)
        expected_pval = 0
        intersection_list = []
        intersection_num_list = []
        if target in deg_set_2:
            # target_file = open(out_dir + '/' + target, 'w')
            print(target)
            temp_list = []
            for interest in deg_dict_2.keys():
                temp_list.append(interest)
                gene_set_2 = set(deg_dict_2[interest])
                gene_set_num_2 = len(gene_set_2)
                intersection = len(gene_set_1.intersection(gene_set_2))
                intersection_num_list.append(intersection)
                _, pval = stats.fisher_exact([[intersection, gene_set_num_1 - intersection],
                                              [gene_set_num_2 - intersection,
                                               gene_universe_num - gene_set_num_2 - gene_set_num_1 + intersection]],
                                             alternative='greater')
                intersection_list.append([interest, intersection, pval])
                pval_value_list.append(pval)
                pval_list.append([interest, pval])
            rank = 0
            if sort_key == 'p':
                # sort by p-value
                pval_value_list.sort()
                pval_value_list = np.array(pval_value_list)
                pval_list.sort(key=lambda x: x[1])
                rank = 0
                idx = 1
                for pval in pval_list:
                    if pval[0] == target:
                        all_rank = np.where(pval_value_list == pval[1])[0] + 1
                        rank = int(np.mean(all_rank))
                    idx += 1
            elif sort_key == 'n':
                # sort by intersection number
                intersection_num_list.sort()
                intersection_num_list.reverse()
                intersection_num_list = np.array(intersection_num_list)
                intersection_list.sort(key=lambda x: x[1])
                np.savetxt(out_dir+'/'+target+'_fisher_num', np.array(intersection_list), fmt='%s', delimiter='\t')
                intersection_list.reverse()
                rank = 0
                idx = 1
                for intersection in intersection_list:
                    if intersection[0] == target:
                        all_rank = np.where(intersection_num_list == intersection[1])[0] + 1
                        rank = int(np.mean(all_rank))
                    idx += 1
            all_result.append([target, rank])
            # intersection_file.write(str(rank)+'\t')
            # intersection_file.write(','.join(list(map(str, intersection_num_list)))+'\n')
            rank_list.append(rank)
            # expected_pval_list.append(expected_pval)
    all_result.sort(key=lambda x: x[0])
    return np.array(rank_list), all_result


if __name__ == '__main__':
    deg_file_1 = sys.argv[1]
    deg_file_2 = sys.argv[2]

    sort_key = sys.argv[3]
    id = sys.argv[4]

    if sort_key == 'n':
        type = '_intersection_num'
    elif sort_key == 'p':
        type = '_pval'

    out_dir = sys.argv[5]
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    # id = deg_file_1.split('/')[-1].split('_')[0]
    gene_universe_num = 22794

    deg_header_1, deg_dict_1 = extract_deg_file(deg_file_1)
    deg_set_1 = set(deg_header_1)
    deg_header_2, deg_dict_2 = extract_deg_file(deg_file_2)
    deg_set_2 = set(deg_header_2)
    # intersection_file = open('fisher_intersection_num.txt','w')

    rank_data, all_result = fisher_test()
    np.savetxt(out_dir+'/'+id+'_fisher_test_rank'+type+'.txt', np.array(all_result), fmt='%s', delimiter='\t')
    print(out_dir+'/'+id+'_fisher_test_rank'+type+'.txt')
    """
    # np.savetxt('fisher_Test/' + id + '_fisher_test_pval.txt', pval_list, fmt='%s')
    _ = plt.hist(rank_data, bins=np.arange(min(rank_data), max(rank_data) + 5, 5))
    plt.title('Mean rank')
    plt.xlabel("rank")
    plt.ylabel("count")
    plt.figtext(.5, .005, "mean rank = " + str(np.mean(rank_data)) + ', number = ' + str(len(rank_data)),
                horizontalalignment='center')
    plt.savefig(out_dir+'/'+id+'_rank_hist'+type+'.png', dpi=600)
    plt.clf()
    intersection_file.close()
    """