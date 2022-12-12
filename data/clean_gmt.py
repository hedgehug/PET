import os
import sys


if __name__ == '__main__':    
    K562_top = 'GGSEA_test/K562_top_200.gmt'
    K562_bottom = 'GGSEA_test/K562_bottom_200.gmt'
    
    HepG2_top = 'GGSEA_test/HepG2_top_200.gmt'
    HepG2_bottom = 'GGSEA_test/HepG2_bottom_200.gmt'

    K562_gene = set()
    HepG2_gene = set()

    K562_top_file = open(K562_top, 'r')
    K562_top_dict = dict()
    K562_top_info = K562_top_file.readlines()
    for line in K562_top_info:
        gene = line.split()[0]
        K562_top_dict[gene] = line
        K562_gene.add(gene)
    K562_top_file.close()

    K562_bottom_file = open(K562_bottom, 'r')
    K562_bottom_dict = dict()
    K562_bottom_info = K562_bottom_file.readlines()
    for line in K562_bottom_info:
        gene = line.split()[0]
        K562_bottom_dict[gene] = line
        # K562_gene.add(gene)
    K562_bottom_file.close()

    HepG2_top_file = open(HepG2_top, 'r')
    HepG2_top_dict = dict()
    HepG2_top_info = HepG2_top_file.readlines()
    for line in HepG2_top_info:
        gene = line.split()[0]
        HepG2_top_dict[gene] = line
        HepG2_gene.add(gene)
    HepG2_top_file.close()

    HepG2_bottom_file = open(HepG2_bottom, 'r')
    HepG2_bottom_dict = dict()
    HepG2_bottom_info = HepG2_bottom_file.readlines()
    for line in HepG2_bottom_info:
        gene = line.split()[0]
        HepG2_bottom_dict[gene] = line
        # HepG2_gene.add(gene)
    HepG2_bottom_file.close()

    shared_gene = HepG2_gene.intersection(K562_gene)

    K562_new_top_file = open(K562_top.replace('.gmt', '_filtered.gmt'), 'w')
    for gene in K562_gene:
        if gene in shared_gene:
            K562_new_top_file.write(K562_top_dict[gene])
    K562_top_file.close()

    K562_new_bottom_file = open(K562_bottom.replace('.gmt', '_filtered.gmt'), 'w')
    for gene in K562_gene:
        if gene in shared_gene:
            K562_new_bottom_file.write(K562_bottom_dict[gene])
    K562_bottom_file.close()

    HepG2_new_top_file = open(HepG2_top.replace('.gmt', '_filtered.gmt'), 'w')
    for gene in HepG2_gene:
        if gene in shared_gene:
            HepG2_new_top_file.write(HepG2_top_dict[gene])
    HepG2_top_file.close()

    HepG2_new_bottom_file = open(HepG2_bottom.replace('.gmt', '_filtered.gmt'), 'w')
    for gene in HepG2_gene:
        if gene in shared_gene:
            HepG2_new_bottom_file.write(HepG2_bottom_dict[gene])
    HepG2_bottom_file.close()