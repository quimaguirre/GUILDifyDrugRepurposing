from scipy.stats import fisher_exact, hypergeom
from statsmodels.sandbox.stats.multicomp import multipletests
import numpy as np
import pandas as pd
import os, sys, re


def main():

    type_top = 'p' # p: percent / e: enrich
    top1 = 1
    top2 = 1

    #type_top = 'e' # p: percent / e: enrich
    #top1 = 181 # asthma
    #top1 = 290 # rheumatoid arthritis
    #top2 = 290 # rheumatoid arthritis
    #top2 = 47 # breast cancer
    #top2 = 182 # breast cancer all

    working_dir = '/Users/quim/Dropbox/UPF/PhD/Projects/GUILDify/case_studies_JMB/NSCLC'
    #working_dir = '/Users/quim/Dropbox/UPF/PhD/Projects/GUILDify/case_studies_JMB/comorbidity'

    #disease1 = 'asthma'
    #disease1 = 'rheumatoid_arthritis'
    #disease1 = 'non_small_cell_lung_carcinoma'
    disease1 = 'breast_cancer_all'
    #disease1 = 'breast_cancer'

    n1_dir = os.path.join(working_dir, disease1)
    seed_file_n1 = os.path.join(n1_dir, 'seeds.txt')
    scores_file_n1 = os.path.join(n1_dir, 'guild_scores.txt')
    if type_top == 'e':
        top1_str = 'enrich_{}'.format(top1)
        top2_str = 'enrich_{}'.format(top2)
    else:
        top1_str = top1
        top2_str = top2
    functions_seeds_file_n1_bp = os.path.join(n1_dir, 'enrichment.GObp.fdr_bh.seeds.txt')
    functions_seeds_file_n1_mf = os.path.join(n1_dir, 'enrichment.GOmf.fdr_bh.seeds.txt')
    functions_file_n1_bp = os.path.join(n1_dir, 'enrichment.GObp.fdr_bh.{}.txt'.format(top1_str))
    functions_file_n1_mf = os.path.join(n1_dir, 'enrichment.GOmf.fdr_bh.{}.txt'.format(top1_str))
    output_file = os.path.join(n1_dir, 'results_{}.tsv'.format(disease1))
    output_file_corr = os.path.join(n1_dir, 'results_corr_{}.tsv'.format(disease1))
    num_functions_bp = 8520
    seeds_n1 = read_seed_file(seed_file_n1)
    ranking_n1, node_to_info_n1 = read_scores_file(scores_file_n1)
    if type_top == 'p':
        length_top_nodes_n1 = int(round(top1 * float(len(ranking_n1)) / float(100)))
        top_nodes_n1 = set(ranking_n1[0: length_top_nodes_n1])
    elif type_top == 'e':
        top_nodes_n1 = set(ranking_n1[0: top1])
    seeds_in_top_nodes_n1 = seeds_n1 & top_nodes_n1
    print('{} seeds in the {} top nodes of {}'.format(len(seeds_in_top_nodes_n1), len(top_nodes_n1), disease1))
    top_nodes_without_seeds_n1 = top_nodes_n1 - seeds_n1

    seed_functions_n1_bp, seed_functions_to_info_n1_bp = read_functions_file(functions_seeds_file_n1_bp)
    seed_functions_n1_mf, seed_functions_to_info_n1_mf = read_functions_file(functions_seeds_file_n1_mf)
    seed_functions_n1 = set(seed_functions_n1_bp) | set(seed_functions_n1_mf)
    top_functions_n1_bp, top_functions_to_info_n1_bp = read_functions_file(functions_file_n1_bp)
    top_functions_n1_mf, top_functions_to_info_n1_mf = read_functions_file(functions_file_n1_mf)
    #top_functions_n1 = set(top_functions_n1_bp) | set(top_functions_n1_mf)
    top_functions_n1 = set(top_functions_n1_bp)
    seeds_in_top_functions_n1 = seed_functions_n1 & top_functions_n1
    print('{} seed functions in the {} top functions of {}'.format(len(seeds_in_top_functions_n1), len(top_functions_n1), disease1))

    functions_file_n1_seeds_bp = os.path.join(n1_dir, 'enrichment.GObp.fdr_bh.seeds.txt')
    functions_file_n1_seeds_mf = os.path.join(n1_dir, 'enrichment.GOmf.fdr_bh.seeds.txt')
    top_functions_seeds_n1_bp, top_functions_to_info_n1_seeds_bp = read_functions_file(functions_file_n1_seeds_bp)
    top_functions_seeds_n1_mf, top_functions_to_info_n1_seeds_mf = read_functions_file(functions_file_n1_seeds_mf)
    top_functions_seeds_n1 = set(top_functions_seeds_n1_bp) | set(top_functions_seeds_n1_mf)
    top_functions_without_seeds_n1 = top_functions_n1 - top_functions_seeds_n1


    drugs = ['afatinib', 'ceritinib', 'crizotinib', 'erlotinib', 'gefitinib', 'palbociclib']
    #drugs = ['rheumatoid_arthritis']
    #drugs = ['breast_cancer']
    #drugs = ['breast_cancer_all']

    columns = ['g-n', 'g-p', 's-n', 's-p', 'f-n', 'f-p', 'sf-n', 'sf-p']
    #columns = ['g-n', 'g-p', 'gws-n', 'gws-p', 's-n', 's-p', 'f-n', 'f-p', 'fws-n', 'fws-p', 'sf-n', 'sf-p']
    results = pd.DataFrame(columns=columns)

    for drug in drugs:

        n2_dir = os.path.join(working_dir, '{}'.format(drug))
        seed_file_n2 = os.path.join(n2_dir, 'seeds.txt')
        scores_file_n2 = os.path.join(n2_dir, 'guild_scores.txt')
        functions_seeds_file_n2_bp = os.path.join(n2_dir, 'enrichment.GObp.fdr_bh.seeds.txt')
        functions_seeds_file_n2_mf = os.path.join(n2_dir, 'enrichment.GOmf.fdr_bh.seeds.txt')
        functions_file_n2_bp = os.path.join(n2_dir, 'enrichment.GObp.fdr_bh.{}.txt'.format(top2_str))
        functions_file_n2_mf = os.path.join(n2_dir, 'enrichment.GOmf.fdr_bh.{}.txt'.format(top2_str))
        seeds_n2 = read_seed_file(seed_file_n2)
        ranking_n2, node_to_info_n2 = read_scores_file(scores_file_n2)
        if type_top == 'p':
            length_top_nodes_n2 = int(round(top2 * float(len(ranking_n2)) / float(100)))
            top_nodes_n2 = set(ranking_n2[0: length_top_nodes_n2])
        elif type_top == 'e':
            top_nodes_n2 = set(ranking_n2[0: top2])
        seeds_in_top_nodes_n2 = seeds_n2 & top_nodes_n2
        print('{} seeds in the {} top nodes of {}'.format(len(seeds_in_top_nodes_n2), len(top_nodes_n2), drug))
        top_nodes_without_seeds_n2 = top_nodes_n2 - seeds_n2

        seed_functions_n2_bp, seed_functions_to_info_n2_bp = read_functions_file(functions_seeds_file_n2_bp)
        seed_functions_n2_mf, seed_functions_to_info_n2_mf = read_functions_file(functions_seeds_file_n2_mf)
        seed_functions_n2 = set(seed_functions_n2_bp) | set(seed_functions_n2_mf)
        top_functions_n2_bp, top_functions_to_info_n2_bp = read_functions_file(functions_file_n2_bp)
        top_functions_n2_mf, top_functions_to_info_n2_mf = read_functions_file(functions_file_n2_mf)
        #top_functions_n2 = set(top_functions_n2_bp) | set(top_functions_n2_mf)
        top_functions_n2 = set(top_functions_n2_bp)
        seeds_in_top_functions_n2 = seed_functions_n2 & top_functions_n2
        print('{} seed functions in the {} top functions of {}'.format(len(seeds_in_top_functions_n2), len(top_functions_n2), drug))

        functions_file_n2_seeds_bp = os.path.join(n2_dir, 'enrichment.GObp.fdr_bh.seeds.txt')
        functions_file_n2_seeds_mf = os.path.join(n2_dir, 'enrichment.GOmf.fdr_bh.seeds.txt')
        top_functions_seeds_n2_bp, top_functions_to_info_n2_seeds_bp = read_functions_file(functions_file_n2_seeds_bp)
        top_functions_seeds_n2_mf, top_functions_to_info_n2_seeds_mf = read_functions_file(functions_file_n2_seeds_mf)
        top_functions_seeds_n2 = set(top_functions_seeds_n2_bp) | set(top_functions_seeds_n2_mf)
        top_functions_without_seeds_n2 = top_functions_n2 - top_functions_seeds_n2

        # Calculate intersections
        intersection_top = top_nodes_n1 & top_nodes_n2
        intersection_seeds = seeds_n1 & seeds_n2
        intersection_without_seeds = intersection_top - intersection_seeds
        intersection_functions = set(top_functions_n1) & set(top_functions_n2)
        intersection_functions_seeds = set(top_functions_seeds_n1) & set(top_functions_seeds_n2)
        intersection_functions_without_seeds = intersection_functions - intersection_functions_seeds
        print('{}'.format(drug))
        print('{} common genes, {}, {}'.format(len(intersection_top), len(set(top_nodes_n1)-intersection_top), len(set(top_nodes_n2)-intersection_top)))
        print('{} common seeds '.format(len(intersection_seeds)))
        print('{} common functions, {}, {}'.format(len(intersection_functions), len(set(top_functions_n1)-intersection_functions), len(set(top_functions_n2)-intersection_functions)))
        print('{} common functions seeds'.format(len(intersection_functions_seeds)))
        print('{} common genes (without common seeds)'.format(len(intersection_without_seeds)))
        print('{} common functions (without common seed functions)'.format(len(intersection_functions_without_seeds)))

        # Calculate fisher test
        oddsratio, pvalue = calculate_fisher(len(intersection_top), len(top_nodes_n1), len(top_nodes_n2), len(ranking_n1))
        print('Intersection: odds ratio {:.3f}, p-value {:.1E}'.format(oddsratio, pvalue))
        oddsratio2, pvalue2 = calculate_fisher(len(intersection_seeds), len(seeds_n1), len(seeds_n2), len(ranking_n1))
        print('Intersection seeds: odds ratio {:.3f}, p-value {:.1E}'.format(oddsratio2, pvalue2))
        oddsratio3, pvalue3 = calculate_fisher(len(intersection_functions), len(top_functions_n1), len(top_functions_n2), num_functions_bp)
        print('Intersection functions: odds ratio {:.3f}, p-value {:.1E}'.format(oddsratio3, pvalue3))
        oddsratio4, pvalue4 = calculate_fisher(len(intersection_functions_seeds), len(top_functions_seeds_n1), len(top_functions_seeds_n2), num_functions_bp)
        print('Intersection functions seeds: odds ratio {:.3f}, p-value {:.1E}'.format(oddsratio4, pvalue4))
        oddsratio5, pvalue5 = calculate_fisher(len(intersection_without_seeds), len(top_nodes_n1 - seeds_in_top_nodes_n1), len(top_nodes_n2 - seeds_in_top_nodes_n2), len(ranking_n1))
        print('Intersection without seeds: odds ratio {:.3f}, p-value {:.1E}'.format(oddsratio5, pvalue5))
        oddsratio6, pvalue6 = calculate_fisher(len(intersection_functions_without_seeds), len(top_functions_n1 - top_functions_seeds_n1), len(top_functions_n2 - top_functions_seeds_n2), num_functions_bp)
        print('Intersection functions without seeds: odds ratio {:.3f}, p-value {:.1E}'.format(oddsratio6, pvalue6))

        # # Calculate hypergeometric
        # M = len(ranking_n1) # Sample size
        # n = len(top_nodes_n1) # Correct associations in sample
        # N = len(top_nodes_n2) # Associations of our prediction
        # k = len(intersection_top) # Number of correct associations in our prediction
        # #print('M: {}, n: {}, N: {}, k: {}'.format(M, n, N, k))
        # rv = hypergeom(M, n, N)
        # x = np.arange(0, n+1)
        # pmf = rv.pmf(x)
        # prb = hypergeom.pmf(k, M, n, N)
        # print('Intersection: hypergeom {:.1E}'.format(prb))

        # # Calculate hypergeometric without seeds
        # [k, M, n, N] = [len(intersection_without_seeds), len(ranking_n1), len(top_nodes_n1 - intersection_seeds), len(top_nodes_n2 - intersection_seeds)]
        # #print('M: {}, n: {}, N: {}, k: {}'.format(M, n, N, k))
        # rv = hypergeom(M, n, N)
        # x = np.arange(0, n+1)
        # pmf = rv.pmf(x)
        # prb = hypergeom.pmf(k, M, n, N)
        # print('Intersection without seeds: hypergeom {:.1E}'.format(prb))

        # # Calculate hypergeometric for functions
        # [k, M, n, N] = [len(intersection_functions), num_functions, len(top_functions_n1), len(top_functions_n2)]
        # #print('M: {}, n: {}, N: {}, k: {}'.format(M, n, N, k))
        # rv = hypergeom(M, n, N)
        # x = np.arange(0, n+1)
        # pmf = rv.pmf(x)
        # prb = hypergeom.pmf(k, M, n, N)
        # print('Intersection functions: hypergeom {:.1E}'.format(prb))

        # Insert the results in the dataframe
        row = ['{}'.format(len(intersection_top)), '{:.1E}'.format(pvalue), '{}'.format(len(intersection_seeds)), '{:.1E}'.format(pvalue2), '{}'.format(len(intersection_functions)), '{:.1E}'.format(pvalue3), '{}'.format(len(intersection_functions_seeds)), '{:.1E}'.format(pvalue4)]
        #row = [len(intersection_top), '{:.1E}'.format(pvalue), len(intersection_without_seeds), '{:.1E}'.format(pvalue5), len(intersection_seeds), '{:.1E}'.format(pvalue2), len(intersection_functions), '{:.1E}'.format(pvalue3), len(intersection_functions_without_seeds), '{:.1E}'.format(pvalue6), len(intersection_functions_seeds), '{:.1E}'.format(pvalue4)]
        df2 = pd.DataFrame([row], columns=columns, index=[drug])
        results = results.append(df2)

    results.to_csv(output_file, sep='\t')
    print(results)

    # Correct p-values
    results_corr = results.copy()
    g_pvalues = [float(x) for x in list(results['g-p'])]
    s_pvalues = [float(x) for x in list(results['s-p'])]
    f_pvalues = [float(x) for x in list(results['f-p'])]
    sf_pvalues = [float(x) for x in list(results['sf-p'])]

    results_corr['g-p'] = ['{:.1E}'.format(x) for x in multipletests(g_pvalues, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)[1]]
    results_corr['s-p'] = ['{:.1E}'.format(x) for x in multipletests(s_pvalues, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)[1]]
    results_corr['f-p'] = ['{:.1E}'.format(x) for x in multipletests(f_pvalues, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)[1]]
    results_corr['sf-p'] = ['{:.1E}'.format(x) for x in multipletests(sf_pvalues, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)[1]]

    results_corr.to_csv(output_file_corr, sep='\t')
    print(results_corr)

    return

#######################
#######################
# SECONDARY FUNCTIONS #
#######################
#######################

def calculate_fisher(n_common, n1, n2, n_node):
    contingency = [[n_common, n1 - n_common], [n2 - n_common, n_node - n1 - n2 + n_common]]
    oddsratio, pvalue = fisher_exact(contingency, alternative="greater")
    return oddsratio, pvalue


def read_seed_file(seed_file):
    seeds = set()
    with open(seed_file, 'r') as seed_fd:
        first_line = seed_fd.readline()
        for line in seed_fd:
            fields = line.strip().split('\t')
            biana_id, uniprot, gene, description, geneID, eq_entries = fields
            seeds.add(biana_id)
    return seeds


def read_scores_file(scores_file):
    ranking_nodes = []
    ranking_nodes_to_info = {}
    with open(scores_file, 'r') as scores_fd:
        first_line = scores_fd.readline()
        for line in scores_fd:
            fields = line.strip().split('\t')
            biana_id, uniprot, gene, seed, description, geneID, eq_entries, score, degree = fields
            ranking_nodes.append(biana_id)
            ranking_nodes_to_info[biana_id] = [biana_id, uniprot, gene, seed, description, geneID, eq_entries, score, degree]
    return ranking_nodes, ranking_nodes_to_info


def read_functions_file(functions_file):
    functions_top = []
    functions_top_to_info = {}
    with open(functions_file, 'r') as func_fd:
        for line in func_fd:
            if line[0] == '#':
                continue
            fields = line.strip().split('\t')
            #genes, total, odds, pval, pval_adj, go_id, go_name, go_type = fields
            go_id, go_name, genes, total, pval, pval_adj = fields
            if float(pval_adj) > 0.05:
                continue
            functions_top.append(go_id)
            functions_top_to_info[go_id] = [go_id, go_name, pval, pval_adj]
    return functions_top, functions_top_to_info


if  __name__ == "__main__":
    main()


