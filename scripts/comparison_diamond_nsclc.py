import sys, os, re
import math
import numpy as np
import scipy.stats as stats
from statsmodels.sandbox.stats.multicomp import multipletests

def main():
    """
    # non-small cell lung carcinoma DIAMOnD: 315682e4-b08d-462a-a5d3-a469515075d7 ==> non_small_cell_lung_carcinoma_diamond
    # erlotinib DIAMOnD:                     3ca9ddfb-79e5-4811-9de2-e927397b0cdc ==> erlotinib_diamond
    # afatinib DIAMOnD: 7083f6b5-b7e4-4f81-bc07-e1384afdee21
    # ceritinib DIAMOnD: 5580e02d-e22e-4adf-ae83-7030d782b3c3
    # crizotinib DIAMOnD: 429dbf00-2cce-4bdc-8be3-1b2c39794f5c
    # gefitinib DIAMOnD: af19313e-ddd3-488a-a58c-812e15f58443
    # palbociclib DIAMOnD: 2dc53bfb-3a69-467d-93c5-097a5e0357ed
    """
    sessions_dir='/var/www/html/guildify2/data/sessions'
    msigdb_all_file = '/home/quim/Databases/MSIgDB/c2.all.v6.2.entrez.gmt'
    msigdb_cp_file = '/home/quim/Databases/MSIgDB/c2.cp.v6.2.entrez.gmt'
    output_dir = '/home/quim/PHD/Projects/GUILDify/drug_repurposing/comparison_diamond_nsclc'
    create_directory(output_dir)
    correction='fdr_bh'
    type_top='percent'
    top_value=1

    sessions = ['non_small_cell_lung_carcinoma', 'afatinib', 'ceritinib', 'crizotinib', 'erlotinib', 'gefitinib', 'palbociclib']
    netscore_pvals=[]
    diamond_pvals=[]
    seeds_pvals=[]

    for session in sessions:

        session_ns = session
        session_dm = session+'_diamond'
        session_ns_dir = os.path.join(sessions_dir, session_ns)
        session_dm_dir = os.path.join(sessions_dir, session_dm)

        ####------------------------------------
        #### 1. Get pathways associated to seeds
        ####------------------------------------

        # get seeds
        seeds_file_ns = os.path.join(session_ns_dir, 'seeds.txt')
        seeds, seeds_to_info=read_seed_file(seeds_file_ns)

        # get ranking of nodes
        scores_file_ns = os.path.join(session_ns_dir, 'guild_scores.txt')
        scores_file_dm = os.path.join(session_dm_dir, 'guild_scores.txt')
        ranking_ns, node_to_info_ns = read_scores_file(scores_file_ns)
        ranking_dm, node_to_info_dm = read_scores_file(scores_file_dm)

        # get top-scoring nodes
        top_nodes_ns = get_top_scoring_nodes(ranking_ns, len(ranking_ns), top_value, type_top)
        top_nodes_dm = get_top_scoring_nodes(ranking_dm, len(ranking_ns), top_value, type_top)
        print('Top-nodes in {}: {}'.format(session_ns, len(top_nodes_ns)))
        print('Top-nodes in {}: {}'.format(session_dm, len(top_nodes_dm)))

        # get entrez top-scoring nodes
        top_entrez_ns=get_geneids_of_ranking_nodes(top_nodes_ns, node_to_info_ns)
        top_entrez_dm=get_geneids_of_ranking_nodes(top_nodes_dm, node_to_info_dm)
        print('Top-entrez in {}: {}'.format(session_ns,len(top_entrez_ns)))
        print('Top-entrez in {}: {}'.format(session_dm,len(top_entrez_dm)))

        # get entrez seeds
        entrez_seeds=get_geneids_of_seeds(seeds_to_info)

        # get all nodes in network as entrez genes
        entrez_genes_network=get_geneids_of_ranking_nodes(ranking_ns, node_to_info_ns)

        # parse MSIgDB and get only associations with nodes in network
        path_to_genes,gene_to_paths=parse_msigdb(msigdb_all_file, entrez_genes_network=entrez_genes_network)
        # print('Entrez genes in network: {}'.format(len(entrez_genes_network)))
        # print('Entrez genes in network with pathways: {}'.format(len(gene_to_paths)))
        # print('Pathways with genes in network: {}'.format(len(path_to_genes)))

        # get pathways associated with seeds
        pathways_seeds=get_pathways_associated_with_seeds(entrez_seeds, gene_to_paths)
        print('Pathways associated with {} seeds: {}'.format(session_ns, len(pathways_seeds)))

        ####---------------------------------------------
        #### 2. Get enriched pathways associated to seeds
        ####---------------------------------------------

        enriched_pathways_file = os.path.join(output_dir, 'enriched_pathways_seeds_{}_{}.txt'.format(session_ns, correction))
        if not fileExist(enriched_pathways_file):
            pathway_to_stats = calculate_enriched_pathways(pathways_seeds, path_to_genes, entrez_seeds, entrez_genes_network, correction=correction)
            write_enriched_pathways(pathway_to_stats, enriched_pathways_file)
        pathway_to_stats=read_enriched_pathways(enriched_pathways_file)
        print('Pathways significantly enriched associated with {} seeds: {}'.format(session_ns, len(pathway_to_stats)))

        ####------------------------------------------------------------
        #### 3. Calculate enrichment of top genes with enriched pathways
        ####------------------------------------------------------------

        # Get genes in enriched pathways
        genes_pathways_enriched = set()
        for pathway in pathway_to_stats:
            if pathway in path_to_genes:
                for gene in path_to_genes[pathway]:
                    genes_pathways_enriched.add(gene)

        # Enrichment
        oddsratio_ns, pvalue_ns = calculate_enrichment_top_genes(top_entrez_ns, genes_pathways_enriched, entrez_genes_network)
        oddsratio_dm, pvalue_dm = calculate_enrichment_top_genes(top_entrez_dm, genes_pathways_enriched, entrez_genes_network)
        oddsratio_sd, pvalue_sd = calculate_enrichment_top_genes(entrez_seeds, genes_pathways_enriched, entrez_genes_network)
        print('{} p-values. NetScore: {}; DIAMOnD: {}; Seeds: {}.'.format(session_ns, pvalue_ns, pvalue_dm, pvalue_sd))

        netscore_pvals.append(-math.log10(pvalue_ns))
        diamond_pvals.append(-math.log10(pvalue_dm))
        seeds_pvals.append(-math.log10(pvalue_sd))

    ####-----------------------
    #### 4. Plot the enrichment
    ####-----------------------

    from matplotlib import pyplot as plt
    import pylab

    plot_name = os.path.join(output_dir, 'enrichment_msigdb_nsclc.png')
    fig = pylab.figure(dpi=300)
    ax = pylab.axes()

    # Define the x locations of the groups
    ind = np.arange(len(netscore_pvals))  # the x locations for the groups
    width = 0.23  # the width of the bars (in this way, we let aprox 0.3 of white space)

    # Plot the groups
    rects1 = ax.bar(ind - width, netscore_pvals, width,
                    label='NetScore', color='SkyBlue')
    rects2 = ax.bar(ind, diamond_pvals, width,
                    label='DIAMOnD', color='Orange')
    rects3 = ax.bar(ind + width, seeds_pvals, width,
                    label='seed genes', color='IndianRed')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    labels = ['NSCLC']+sessions[1:]
    ax.set_ylabel('enrichment [- log p]')
    ax.set_xticks(ind)
    ax.set_xticklabels(labels)
    ax.legend()
    pylab.savefig(plot_name, format='png')


    return

#######################
#######################
# SECONDARY FUNCTIONS #
#######################
#######################

def fileExist(file):
    """
    Checks if a file exists AND is a file
    """
    return os.path.exists(file) and os.path.isfile(file)

def create_directory(directory):
    """
    Checks if a directory exists and if not, creates it
    """
    try:
        os.stat(directory)
    except:
        os.mkdir(directory)
    return

def parse_msigdb(input_file, entrez_genes_network=None):
    path_to_genes={}
    gene_to_paths={}
    with open(input_file, 'r') as input_fd:
        for line in input_fd:
            fields = line.strip().split('\t')
            path=fields[0]
            link=fields[1]
            genes=fields[2:]
            for gene in genes:
                if entrez_genes_network:
                    if gene in entrez_genes_network:
                        path_to_genes.setdefault(path, set()).add(gene)
                        gene_to_paths.setdefault(gene, set()).add(path)
                else:
                    path_to_genes.setdefault(path, set()).add(gene)
                    gene_to_paths.setdefault(gene, set()).add(path)
    return path_to_genes, gene_to_paths

def read_seed_file(seed_file):
    seeds = set()
    seed_to_info={}
    with open(seed_file, 'r') as seed_fd:
        first_line = seed_fd.readline()
        for line in seed_fd:
            fields = line.strip().split('\t')
            biana_id, uniprot, gene, description, geneID, eq_entries = fields
            seed_to_info[biana_id]=[uniprot, gene, description, geneID, eq_entries]
            seeds.add(biana_id)
    return seeds, seed_to_info

def get_geneids_of_seeds(seed_to_info):
    entrez_seeds = set()
    for seed in seed_to_info:
        entrez_string=seed_to_info[seed][3]
        for entrez in entrez_string.split('; '):
            entrez_seeds.add(entrez)
    return entrez_seeds

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

def get_top_scoring_nodes(ranking, num_nodes, top_value, type_top):
    if type_top == 'percent':
        length_top_nodes = int(round(top_value * float(num_nodes) / float(100)))
        top_nodes = set(ranking[0: length_top_nodes])
    elif type_top == 'enrichment':
        top_nodes = set(ranking[0: top_value])
    else:
        print('Unrecognised type of top: {}'.format(type_top))
        sys.exit(10)
    return top_nodes

def get_geneids_of_ranking_nodes(nodes, node_to_info):
    entrez_genes=set()
    for node in nodes:
        if node in node_to_info:
            entrez_string=node_to_info[node][5]
            for entrez in entrez_string.split('; '):
                entrez_genes.add(entrez)
    return entrez_genes

def get_pathways_associated_with_seeds(entrez_seeds, gene_to_paths):
    pathways_associated=set()
    for gene in entrez_seeds:
        if gene in gene_to_paths:
            for path in gene_to_paths[gene]:
                pathways_associated.add(path)
    return pathways_associated

def calculate_enriched_pathways(pathways_associated, path_to_genes, genes_seeds, genes_network, correction='fdr_bh'):

    #------------------------------------------------------------------------
    #                | genes seeds      | genes not seeds                   |
    #------------------------------------------------------------------------
    # genes path     | #common          | #path - #common                   |
    #------------------------------------------------------------------------
    # genes not path | #seeds - #common | #total - #seeds - #path - #common |
    #------------------------------------------------------------------------

    # Calculate Fisher's exact test

    initial_values=[]
    pvalues=[]

    for pathway in pathways_associated:

        genes_path=path_to_genes[pathway]
        genes_common=genes_path&genes_seeds

        n_common = len(genes_common)
        n_seeds = len(genes_seeds)
        n_path = len(genes_path)
        n_total = len(genes_network)
        #print(pathway, n_common, n_seeds, n_path, n_total)

        contingency_table = [
            [ n_common,         n_path-n_common                  ],
            [ n_seeds-n_common, n_total-n_seeds-n_path-n_common  ]
        ]
        oddsratio, pvalue = stats.fisher_exact(contingency_table)

        initial_values.append([pathway, n_common, n_path, n_seeds, oddsratio, pvalue])
        pvalues.append(pvalue)


    # Correct p-value
    pvals_corrected= multipletests(pvalues, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)[1]

    # Get final values
    pathway_to_stats={}
    for x in xrange(len(initial_values)):
        row=initial_values[x]
        [pathway, n_common, n_path, n_seeds, oddsratio, pvalue]=row
        pvalue_corrected=pvals_corrected[x]
        pathway_to_stats[pathway]=[n_common, n_path, n_seeds, oddsratio, pvalue, pvalue_corrected]

    return pathway_to_stats

def write_enriched_pathways(pathway_to_stats, output_file, cutoff=0.01):
    with open(output_file, 'w') as output_fd:
        for pathway, data in sorted(pathway_to_stats.items(), key=lambda x:x[1][5], reverse=False):
            [n_common, n_path, n_seeds, oddsratio, pvalue, pvalue_corrected]=pathway_to_stats[pathway]
            if pvalue_corrected < cutoff:
                output_fd.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(pathway, n_common, n_path, n_seeds, oddsratio, pvalue, pvalue_corrected))
    return

def read_enriched_pathways(input_file, cutoff=0.01):
    pathway_to_stats={}
    with open(input_file, 'r') as input_fd:
        for line in input_fd:
            pathway, n_common, n_path, n_seeds, oddsratio, pvalue, pvalue_corrected = line.strip().split('\t')
            if float(pvalue_corrected) < cutoff:
                pathway_to_stats[pathway]=[int(n_common), int(n_path), int(n_seeds), float(oddsratio), float(pvalue), float(pvalue_corrected)]
    return pathway_to_stats

def calculate_enrichment_top_genes(genes_top, genes_path, genes_network):

    #--------------------------------------------------------------------
    #               | genes path      | genes not path                  |
    #--------------------------------------------------------------------
    # genes top     | #common         | #top - #common                  |
    #--------------------------------------------------------------------
    # genes not top | #path - #common | #total - #path - #top - #common |
    #--------------------------------------------------------------------

    # Calculate Fisher's exact test
    genes_common=genes_top&genes_path

    n_common = len(genes_common)
    n_top = len(genes_top)
    n_path = len(genes_path)
    n_total = len(genes_network)
    print(n_common, n_top, n_path, n_total)

    contingency_table = [
        [ n_common,         n_top-n_common               ],
        [ n_path-n_common, n_total-n_path-n_top-n_common ]
    ]
    oddsratio, pvalue = stats.fisher_exact(contingency_table)

    return oddsratio, pvalue


if  __name__ == "__main__":
    main()