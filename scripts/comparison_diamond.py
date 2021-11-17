import sys, os, re
import math
import numpy as np
import scipy.stats as stats
from statsmodels.sandbox.stats.multicomp import multipletests

def main():
    """
    # asthma DIAMOnD:                   322189de-9f19-4b96-8374-fe594eff7769 ==> asthma_diamond
    # rheumatoid arthiritis DIAMOnD:    c3e61664-b0f7-49ee-af3e-c1cca6b65472 ==> rheumatoid_arthritis_diamond
    """
    sessions_dir='/var/www/html/guildify2/data/sessions'
    asthma_ns_dir = os.path.join(sessions_dir, 'asthma')
    asthma_dm_dir = os.path.join(sessions_dir, 'asthma_diamond')
    rarthritis_ns_dir = os.path.join(sessions_dir, 'rheumatoid_arthritis')
    rarthritis_dm_dir = os.path.join(sessions_dir, 'rheumatoid_arthritis_diamond')
    msigdb_all_file = '/home/quim/Databases/MSIgDB/c2.all.v6.2.entrez.gmt'
    msigdb_cp_file = '/home/quim/Databases/MSIgDB/c2.cp.v6.2.entrez.gmt'
    output_dir = '/home/quim/PHD/Projects/GUILDify/drug_repurposing/comparison_diamond'
    create_directory(output_dir)
    correction='fdr_bh'
    type_top='enrichment'
    top_value_asthma_ns=181
    top_value_asthma_dm=151
    top_value_rarthritis_ns=290
    top_value_rarthritis_dm=235


    ####------------------------------------
    #### 1. Get pathways associated to seeds
    ####------------------------------------

    # get seeds
    seeds_file_asthma_ns = os.path.join(asthma_ns_dir, 'seeds.txt')
    seeds_file_rarthritis_ns = os.path.join(rarthritis_ns_dir, 'seeds.txt')
    seeds_asthma, seeds_asthma_to_info=read_seed_file(seeds_file_asthma_ns)
    seeds_rarthritis, seeds_rarthritis_to_info=read_seed_file(seeds_file_rarthritis_ns)

    # get ranking of nodes
    scores_file_asthma_ns = os.path.join(asthma_ns_dir, 'guild_scores.txt')
    scores_file_asthma_dm = os.path.join(asthma_dm_dir, 'guild_scores.txt')
    scores_file_rarthritis_ns = os.path.join(rarthritis_ns_dir, 'guild_scores.txt')
    scores_file_rarthritis_dm = os.path.join(rarthritis_dm_dir, 'guild_scores.txt')
    ranking_asthma_ns, node_to_info_asthma_ns = read_scores_file(scores_file_asthma_ns)
    ranking_asthma_dm, node_to_info_asthma_dm = read_scores_file(scores_file_asthma_dm)
    ranking_rarthritis_ns, node_to_info_rarthritis_ns = read_scores_file(scores_file_rarthritis_ns)
    ranking_rarthritis_dm, node_to_info_rarthritis_dm = read_scores_file(scores_file_rarthritis_dm)

    # get top-scoring nodes
    top_nodes_asthma_ns = get_top_scoring_nodes(ranking_asthma_ns, top_value_asthma_ns, type_top)
    top_nodes_asthma_dm = get_top_scoring_nodes(ranking_asthma_dm, top_value_asthma_dm, type_top)
    top_nodes_rarthritis_ns = get_top_scoring_nodes(ranking_rarthritis_ns, top_value_rarthritis_ns, type_top)
    top_nodes_rarthritis_dm = get_top_scoring_nodes(ranking_rarthritis_dm, top_value_rarthritis_dm, type_top)
    top_nodes_intersection_ns = top_nodes_asthma_ns & top_nodes_rarthritis_ns
    top_nodes_intersection_dm = top_nodes_asthma_dm & top_nodes_rarthritis_dm
    print('Top-nodes in asthma (NetScore): {}'.format(len(top_nodes_asthma_ns)))
    print('Top-nodes in asthma (DIAMOnD): {}'.format(len(top_nodes_asthma_dm)))
    print('Top-nodes in r. arthritis (NetScore): {}'.format(len(top_nodes_rarthritis_ns)))
    print('Top-nodes in r. arthritis (DIAMOnD): {}'.format(len(top_nodes_rarthritis_dm)))
    print('Top-nodes in intersection (NetScore): {}'.format(len(top_nodes_intersection_ns)))
    print('Top-nodes in intersection (DIAMOnD): {}'.format(len(top_nodes_intersection_dm)))

    # get entrez top-scoring nodes
    top_entrez_asthma_ns=get_geneids_of_ranking_nodes(top_nodes_asthma_ns, node_to_info_asthma_ns)
    top_entrez_asthma_dm=get_geneids_of_ranking_nodes(top_nodes_asthma_dm, node_to_info_asthma_dm)
    top_entrez_rarthritis_ns=get_geneids_of_ranking_nodes(top_nodes_rarthritis_ns, node_to_info_rarthritis_ns)
    top_entrez_rarthritis_dm=get_geneids_of_ranking_nodes(top_nodes_rarthritis_dm, node_to_info_rarthritis_dm)
    top_entrez_intersection_ns = top_entrez_asthma_ns & top_entrez_rarthritis_ns
    top_entrez_intersection_dm = top_entrez_asthma_dm & top_entrez_rarthritis_dm
    print('Top-entrez in asthma (NetScore): {}'.format(len(top_entrez_asthma_ns)))
    print('Top-entrez in asthma (DIAMOnD): {}'.format(len(top_entrez_asthma_dm)))
    print('Top-entrez in r. arthritis (NetScore): {}'.format(len(top_entrez_rarthritis_ns)))
    print('Top-entrez in r. arthritis (DIAMOnD): {}'.format(len(top_entrez_rarthritis_dm)))
    print('Top-entrez in intersection (NetScore): {}'.format(len(top_entrez_intersection_ns)))
    print('Top-entrez in intersection (DIAMOnD): {}'.format(len(top_entrez_intersection_dm)))

    # get entrez seeds
    entrez_seeds_asthma=get_geneids_of_seeds(seeds_asthma_to_info)
    entrez_seeds_rarthritis=get_geneids_of_seeds(seeds_rarthritis_to_info)
    entrez_seeds_intersection=entrez_seeds_asthma&entrez_seeds_rarthritis

    # get all nodes in network as entrez genes
    entrez_genes_network=get_geneids_of_ranking_nodes(ranking_asthma_ns, node_to_info_asthma_ns)

    # parse MSIgDB and get only associations with nodes in network
    path_to_genes,gene_to_paths=parse_msigdb(msigdb_all_file, entrez_genes_network=entrez_genes_network)
    print('Entrez genes in network: {}'.format(len(entrez_genes_network)))
    print('Entrez genes in network with pathways: {}'.format(len(gene_to_paths)))
    print('Pathways with genes in network: {}'.format(len(path_to_genes)))

    # get pathways associated with seeds
    pathways_seeds_asthma=get_pathways_associated_with_seeds(entrez_seeds_asthma, gene_to_paths)
    pathways_seeds_rarthritis=get_pathways_associated_with_seeds(entrez_seeds_rarthritis, gene_to_paths)
    pathways_seeds_intersection=pathways_seeds_asthma&pathways_seeds_rarthritis
    print('Pathways associated with asthma seeds: {}'.format(len(pathways_seeds_asthma)))
    print('Pathways associated with rheumatoid arthritis seeds: {}'.format(len(pathways_seeds_rarthritis)))
    print('Pathways intersection: {}'.format(len(pathways_seeds_intersection)))

    ####---------------------------------------------
    #### 2. Get enriched pathways associated to seeds
    ####---------------------------------------------

    enriched_pathways_file_asthma = os.path.join(output_dir, 'enriched_pathways_seeds_asthma_{}.txt'.format(correction))
    enriched_pathways_file_rarthritis = os.path.join(output_dir, 'enriched_pathways_seeds_rheumatoid_arthritis_{}.txt'.format(correction))
    enriched_pathways_file_intersection = os.path.join(output_dir, 'enriched_pathways_seeds_intersection_{}.txt'.format(correction))
    if not fileExist(enriched_pathways_file_asthma):
        pathway_to_stats_asthma = calculate_enriched_pathways(pathways_seeds_asthma, path_to_genes, entrez_seeds_asthma, entrez_genes_network, correction=correction)
        write_enriched_pathways(pathway_to_stats_asthma, enriched_pathways_file_asthma)
    if not fileExist(enriched_pathways_file_rarthritis):
        pathway_to_stats_rarthritis = calculate_enriched_pathways(pathways_seeds_rarthritis, path_to_genes, entrez_seeds_rarthritis, entrez_genes_network, correction=correction)
        write_enriched_pathways(pathway_to_stats_rarthritis, enriched_pathways_file_rarthritis)
    if not fileExist(enriched_pathways_file_intersection):
        pathway_to_stats_intersection = calculate_enriched_pathways(pathways_seeds_intersection, path_to_genes, entrez_seeds_intersection, entrez_genes_network, correction=correction)
        write_enriched_pathways(pathway_to_stats_intersection, enriched_pathways_file_intersection)
    pathway_to_stats_asthma=read_enriched_pathways(enriched_pathways_file_asthma)
    pathway_to_stats_rarthritis=read_enriched_pathways(enriched_pathways_file_rarthritis)
    pathways_enriched_intersection=set(pathway_to_stats_asthma.keys())&set(pathway_to_stats_rarthritis.keys())
    pathways_to_stats_intersection=read_enriched_pathways(enriched_pathways_file_intersection)
    print('Pathways significantly enriched associated with asthma seeds: {}'.format(len(pathway_to_stats_asthma)))
    print('Pathways significantly enriched associated with r. arthritis seeds: {}'.format(len(pathway_to_stats_rarthritis)))
    print('Pathways significantly enriched associated with intersection seeds: {}'.format(len(pathways_to_stats_intersection)))
    print('Intersection of pathways significantly enriched: {}'.format(len(pathways_enriched_intersection)))

    ####------------------------------------------------------------
    #### 3. Calculate enrichment of top genes with enriched pathways
    ####------------------------------------------------------------

    # Get genes in enriched pathways
    genes_pathways_enriched_asthma = set()
    for pathway in pathway_to_stats_asthma:
        if pathway in path_to_genes:
            for gene in path_to_genes[pathway]:
                genes_pathways_enriched_asthma.add(gene)

    genes_pathways_enriched_rarthritis = set()
    for pathway in pathway_to_stats_rarthritis:
        if pathway in path_to_genes:
            for gene in path_to_genes[pathway]:
                genes_pathways_enriched_rarthritis.add(gene)

    genes_pathways_enriched_intersection = set()
    # for pathway in pathways_enriched_intersection:
    for pathway in pathways_to_stats_intersection:
        if pathway in path_to_genes:
            for gene in path_to_genes[pathway]:
                genes_pathways_enriched_intersection.add(gene)

    # Enrichment ASTHMA
    oddsratio_asthma_ns, pvalue_asthma_ns = calculate_enrichment_top_genes(top_entrez_asthma_ns, genes_pathways_enriched_asthma, entrez_genes_network)
    oddsratio_asthma_dm, pvalue_asthma_dm = calculate_enrichment_top_genes(top_entrez_asthma_dm, genes_pathways_enriched_asthma, entrez_genes_network)
    oddsratio_asthma_sd, pvalue_asthma_sd = calculate_enrichment_top_genes(entrez_seeds_asthma, genes_pathways_enriched_asthma, entrez_genes_network)
    print('ASTHMA p-values. NetScore: {}; DIAMOnD: {}; Seeds: {}.'.format(pvalue_asthma_ns, pvalue_asthma_dm, pvalue_asthma_sd))

    # Enrichment R.ARTHRITIS
    oddsratio_rarthritis_ns, pvalue_rarthritis_ns = calculate_enrichment_top_genes(top_entrez_rarthritis_ns, genes_pathways_enriched_rarthritis, entrez_genes_network)
    oddsratio_rarthritis_dm, pvalue_rarthritis_dm = calculate_enrichment_top_genes(top_entrez_rarthritis_dm, genes_pathways_enriched_rarthritis, entrez_genes_network)
    oddsratio_rarthritis_sd, pvalue_rarthritis_sd = calculate_enrichment_top_genes(entrez_seeds_rarthritis, genes_pathways_enriched_rarthritis, entrez_genes_network)
    print('R.ARTHRITIS p-values. NetScore: {}; DIAMOnD: {}; Seeds: {}.'.format(pvalue_rarthritis_ns, pvalue_rarthritis_dm, pvalue_rarthritis_sd))

    # Enrichment INTERSECTION
    oddsratio_intersection_ns, pvalue_intersection_ns = calculate_enrichment_top_genes(top_entrez_intersection_ns, genes_pathways_enriched_intersection, entrez_genes_network)
    oddsratio_intersection_dm, pvalue_intersection_dm = calculate_enrichment_top_genes(top_entrez_intersection_dm, genes_pathways_enriched_intersection, entrez_genes_network)
    oddsratio_intersection_sd, pvalue_intersection_sd = calculate_enrichment_top_genes(entrez_seeds_intersection, genes_pathways_enriched_intersection, entrez_genes_network)
    print('INTERSECTION p-values. NetScore: {}; DIAMOnD: {}; Seeds: {}.'.format(pvalue_intersection_ns, pvalue_intersection_dm, pvalue_intersection_sd))

    ####-----------------------
    #### 4. Plot the enrichment
    ####-----------------------

    from matplotlib import pyplot as plt
    import pylab

    plot_name = os.path.join(output_dir, 'enrichment_msigdb.png')
    fig = pylab.figure(dpi=300)
    ax = pylab.axes()

    # Define groups
    netscore_pvals = (-math.log10(pvalue_asthma_ns), -math.log10(pvalue_rarthritis_ns), -math.log10(pvalue_intersection_ns))
    diamond_pvals = (-math.log10(pvalue_asthma_dm), -math.log10(pvalue_rarthritis_dm), -math.log10(pvalue_intersection_dm))
    seeds_pvals = (-math.log10(pvalue_asthma_sd), -math.log10(pvalue_rarthritis_sd), -math.log10(pvalue_intersection_sd))

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
    ax.plot([ind[0] - width*1.5, ind[len(ind)-1] + width*1.5], [-math.log10(0.05), -math.log10(0.05)], '--', label='p-value=0.05', color='black')

    #ax.set_ylim(0,30)

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('enrichment [- log p]')
    ax.set_xticks(ind)
    ax.set_xticklabels(('Asthma', 'R.Arthritis', 'Intersection'))
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

def get_top_scoring_nodes(ranking, top_value, type_top):
    if type_top == 'percent':
        length_top_nodes = int(round(top_value * float(len(ranking)) / float(100)))
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