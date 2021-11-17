from scipy.stats import fisher_exact, hypergeom
from statsmodels.sandbox.stats.multicomp import multipletests
import numpy as np
import pandas as pd
import os, sys, re


def main():

    type_top = 'e' # p: percent / e: enrich
    top1 = 1
    top2 = 1

    pairs =[
        ['asthma', 'rheumatoid_arthritis'],
        ['asthma', 'breast_cancer_all'],
        ['rheumatoid_arthritis', 'breast_cancer_all'],
        ['asthma', 'breast_cancer'],
        ['rheumatoid_arthritis', 'breast_cancer']
    ]

    disease_to_top_enrich = {
        'asthma' : 181,
        'rheumatoid_arthritis' : 290,
        'breast_cancer_all' : 182,
        'breast_cancer' : 47
    }

    working_dir = '/Users/quim/Dropbox/UPF/PhD/Projects/GUILDify/case_studies_JMB/comorbidity'
    output_file = os.path.join(working_dir, 'results_comorbidity.tsv')
    output_file_corr = os.path.join(working_dir, 'results_comorbidity_corrected_pval.tsv')
    columns = ['g-n', 'g-p', 'gw-n', 'gw-p', 's-n', 's-p', 'f-n', 'f-p', 'fw-n', 'fw-p', 'sf-n', 'sf-p']
    results = pd.DataFrame(columns=columns)

    for pair in pairs:

        [disease1, disease2] = pair

        if type_top == 'e':
            top1=disease_to_top_enrich[disease1]
            top2=disease_to_top_enrich[disease2]
            top1_str = 'enrich_{}'.format(top1)
            top2_str = 'enrich_{}'.format(top2)
        else:
            top1_str = top1
            top2_str = top2


        #### DISEASE 1 ####

        n1_dir = os.path.join(working_dir, disease1)
        seed_file_n1 = os.path.join(n1_dir, 'seeds.txt')
        scores_file_n1 = os.path.join(n1_dir, 'guild_scores.txt')

        functions_seeds_file_n1_bp = os.path.join(n1_dir, 'enrichment.GObp.fdr_bh.seeds.txt')
        functions_seeds_file_n1_mf = os.path.join(n1_dir, 'enrichment.GOmf.fdr_bh.seeds.txt')
        functions_file_n1_bp = os.path.join(n1_dir, 'enrichment.GObp.fdr_bh.{}.txt'.format(top1_str))
        functions_file_n1_mf = os.path.join(n1_dir, 'enrichment.GOmf.fdr_bh.{}.txt'.format(top1_str))
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


        #### DISEASE 2 ####

        n2_dir = os.path.join(working_dir, '{}'.format(disease2))
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
        print('{} seeds in the {} top nodes of {}'.format(len(seeds_in_top_nodes_n2), len(top_nodes_n2), disease2))
        top_nodes_without_seeds_n2 = top_nodes_n2 - seeds_n2

        seed_functions_n2_bp, seed_functions_to_info_n2_bp = read_functions_file(functions_seeds_file_n2_bp)
        seed_functions_n2_mf, seed_functions_to_info_n2_mf = read_functions_file(functions_seeds_file_n2_mf)
        seed_functions_n2 = set(seed_functions_n2_bp) | set(seed_functions_n2_mf)
        top_functions_n2_bp, top_functions_to_info_n2_bp = read_functions_file(functions_file_n2_bp)
        top_functions_n2_mf, top_functions_to_info_n2_mf = read_functions_file(functions_file_n2_mf)
        #top_functions_n2 = set(top_functions_n2_bp) | set(top_functions_n2_mf)
        top_functions_n2 = set(top_functions_n2_bp)
        seeds_in_top_functions_n2 = seed_functions_n2 & top_functions_n2
        print('{} seed functions in the {} top functions of {}'.format(len(seeds_in_top_functions_n2), len(top_functions_n2), disease2))

        functions_file_n2_seeds_bp = os.path.join(n2_dir, 'enrichment.GObp.fdr_bh.seeds.txt')
        functions_file_n2_seeds_mf = os.path.join(n2_dir, 'enrichment.GOmf.fdr_bh.seeds.txt')
        top_functions_seeds_n2_bp, top_functions_to_info_n2_seeds_bp = read_functions_file(functions_file_n2_seeds_bp)
        top_functions_seeds_n2_mf, top_functions_to_info_n2_seeds_mf = read_functions_file(functions_file_n2_seeds_mf)
        top_functions_seeds_n2 = set(top_functions_seeds_n2_bp) | set(top_functions_seeds_n2_mf)
        top_functions_without_seeds_n2 = top_functions_n2 - top_functions_seeds_n2


        #### INTERSECTIONS ####

        intersection_top = top_nodes_n1 & top_nodes_n2
        intersection_seeds = seeds_n1 & seeds_n2
        intersection_without_seeds = intersection_top - intersection_seeds
        intersection_functions = set(top_functions_n1) & set(top_functions_n2)
        intersection_functions_seeds = set(top_functions_seeds_n1) & set(top_functions_seeds_n2)
        intersection_functions_without_seeds = intersection_functions - intersection_functions_seeds
        print('{}'.format(disease2))
        print('{} common genes, {}, {}'.format(len(intersection_top), len(set(top_nodes_n1)-intersection_top), len(set(top_nodes_n2)-intersection_top)))
        print('{} common seeds '.format(len(intersection_seeds)))
        print('{} common functions, {}, {}'.format(len(intersection_functions), len(set(top_functions_n1)-intersection_functions), len(set(top_functions_n2)-intersection_functions)))
        print('{} common functions seeds'.format(len(intersection_functions_seeds)))
        print('{} common genes (without common seeds)'.format(len(intersection_without_seeds)))
        print('{} common functions (without common seed functions)'.format(len(intersection_functions_without_seeds)))


        #### FISHER TEST ####

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

        # Insert the results in the dataframe
        row = ['{}'.format(len(intersection_top)), '{:.1E}'.format(pvalue), '{}'.format(len(intersection_without_seeds)), '{:.1E}'.format(pvalue5), '{}'.format(len(intersection_seeds)), '{:.1E}'.format(pvalue2), '{}'.format(len(intersection_functions)), '{:.1E}'.format(pvalue3), '{}'.format(len(intersection_functions_without_seeds)), '{:.1E}'.format(pvalue6), '{}'.format(len(intersection_functions_seeds)), '{:.1E}'.format(pvalue4)]
        df2 = pd.DataFrame([row], columns=columns, index=['{} - {}'.format(disease1, disease2)])
        results = results.append(df2)

    results.to_csv(output_file, sep='\t')
    print(results)

    # Correct p-values
    results_corr = results.copy()
    g_pvalues = [float(x) for x in list(results['g-p'])]
    gw_pvalues = [float(x) for x in list(results['gw-p'])]
    s_pvalues = [float(x) for x in list(results['s-p'])]
    f_pvalues = [float(x) for x in list(results['f-p'])]
    fw_pvalues = [float(x) for x in list(results['fw-p'])]
    sf_pvalues = [float(x) for x in list(results['sf-p'])]

    results_corr['g-p'] = ['{:.1E}'.format(x) for x in multipletests(g_pvalues, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)[1]]
    results_corr['gw-p'] = ['{:.1E}'.format(x) for x in multipletests(gw_pvalues, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)[1]]
    results_corr['s-p'] = ['{:.1E}'.format(x) for x in multipletests(s_pvalues, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)[1]]
    results_corr['f-p'] = ['{:.1E}'.format(x) for x in multipletests(f_pvalues, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)[1]]
    results_corr['fw-p'] = ['{:.1E}'.format(x) for x in multipletests(fw_pvalues, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)[1]]
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


