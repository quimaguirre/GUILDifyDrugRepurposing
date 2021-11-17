import sys, os, re
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.obo_parser import GODag
from goatools import associations
import math

def main():


    taxID = 9606
    top_genes_file = '/home/quim/PHD/Projects/GUILDify/enrichment_tests/alzheimer_omim/data/guild_scores.txt'
    all_genes_file = '/var/www/html/guildify2/data/BIANA/9606/all/node_info.txt'
    go_obo_file = '/var/www/html/guildify2/data/go/go-basic.obo'
    gene2go_filtered_file = '/var/www/html/guildify2/data/go/gene2go.{}.filtered'.format(taxID)
    output_file = '/home/quim/PHD/Projects/GUILDify/enrichment_tests/alzheimer_omim/goatools/output_enrichment.txt'
    temp_file = '/home/quim/PHD/Projects/GUILDify/enrichment_tests/alzheimer_omim/goatools/temp_enrichment.txt'
    threshold = 1
    evidence_codes_used = ['EXP', 'IDA', 'IMP', 'IGI', 'IEP', 'ISS', 'ISA', 'ISM', 'ISO']
    obodag = GODag(go_obo_file)
    enrichment_cut_off = 0.05

    #geneid2gos = associations.read_ncbi_gene2go(gene2go_file, taxids=[int(taxID)], evidence_set=evidence_codes_used)
    geneid2gos = associations.read_associations(gene2go_filtered_file, no_top=False)

    node_to_geneIDs, geneID_to_nodes = get_geneIDs_from_info_file(all_genes_file, seeds=None)
    user_entity_ranking, user_entity_to_geneids = parse_guild_scores_file(top_genes_file)
    all_nodes = set()
    for user_entity in user_entity_ranking:
        for geneID in user_entity_to_geneids[user_entity]:
            all_nodes.add(geneID)
    all_nodes = list(all_nodes)
    ntop=int(float(threshold)*len(user_entity_ranking)/100.0)
    top_nodes_biana = user_entity_ranking[0:ntop]
    top_nodes = set()
    for node in top_nodes_biana:
        if node in user_entity_to_geneids:
            for geneID in user_entity_to_geneids[node]:
                top_nodes.add(geneID)
    top_nodes = list(top_nodes)

    print(len(top_nodes))
    print(len(all_nodes))
    print(len(set(top_nodes)&set(all_nodes)))

    goeaobj = GOEnrichmentStudy(
            all_nodes, # List of background genes
            geneid2gos, # geneid/GO associations
            obodag, # Ontologies
            propagate_counts = False,
            alpha = 0.05, # default significance cut-off
            methods = ['fdr_bh']) # defult multipletest correction method #also we can use bonferroni, fdr_bh

    goea_results_all = goeaobj.run_study(top_nodes)
    goeaobj.wr_tsv(temp_file, goea_results_all)


    results = []
    with open(temp_file, 'r') as ft:
        first_line = ft.readline()
        for line in ft:
            fields = line.strip().split('\t')
            go = fields[0]
            type_func = fields[1]
            name = fields[3]
            ratio_in_study = fields[4]
            ratio_in_pop = fields[5]
            p_uncorrected = float(fields[6])
            depth = fields[7]
            study_count = fields[8]
            p_corrected = float(fields[9])

            # Calculating the Log of odds ratio
            n_in_study, all_in_study = ratio_in_study.split('/')
            n_in_pop, all_in_pop = ratio_in_pop.split('/')
            q = int(n_in_study) # num genes assoc with function in top
            m = int(n_in_pop) # num genes assoc with func in population
            k = int(all_in_study) # num genes in top
            t = int(all_in_pop) # num genes in population
            log_of_odds = calculate_log_of_odds_ratio(q, k, m, t)
            #print('q: {}, m: {}, k: {}, t: {}'.format(q, m, k, t))
            results.append([q, m, log_of_odds, p_uncorrected, p_corrected, go, name, type_func])

    # Writing the output file
    with open(output_file, 'w') as fo:
        fo.write('# num of genes\tnum of total genes\tLog of odds ratio\tP-value\tAdjusted p-value\tGO term ID\tGO term name\tGO type\n')

        for line in sorted(results, key=lambda x: x[4], reverse=False):
            q, m, log_of_odds, p_uncorrected, p_corrected, go, name, type_func = line
            if q <= 0: # Skip the functions that do not have at least 1 top gene associated (the log of odds ratio is -infinite!!)
                continue
            if log_of_odds < 0: # Skip the log of odds ratio that are negative (we will only use functions with log of odds ratio from 0 to inf)
                continue
            #if p_corrected > enrichment_cut_off:
            #    continue
            if type_func.upper() == 'BP' or type_func.upper() == 'MF':
                new_line = '{}\t{}\t{:.3f}\t{:.2E}\t{:.2E}\t{}\t{}\t{}\n'.format(q, m, log_of_odds, p_uncorrected, p_corrected, go, name, type_func)
                fo.write(new_line)

    # Removing the temporary file
    #command = "rm {}".format(temp_file)
    #os.system(command)



def parse_guild_scores_file(guild_scores_file):
    """
    Parse the output file (guild scores file) and get the values of interest
    """
    user_entity_ranking = []
    user_entity_to_geneids = {}
    with open(guild_scores_file, 'r') as output_file_fd:
        output_file_fd.readline() # skip header
        for line in output_file_fd:
            fields = line.strip().split("\t")
            if len(fields) == 8:
                [user_entity_id, entry_ids, gene_symbols, is_seed, descriptions, gene_ids, equivalent_entries, score] = fields
            elif len(fields) == 9:
                [user_entity_id, entry_ids, gene_symbols, is_seed, descriptions, gene_ids, equivalent_entries, score, degree] = fields
            user_entity_ranking.append(user_entity_id)
            user_entity_to_geneids.setdefault(user_entity_id, set())
            gene_ids = gene_ids.split(';')
            for geneid in gene_ids:
                user_entity_to_geneids[user_entity_id].add(geneid)
    return user_entity_ranking, user_entity_to_geneids


def get_geneIDs_from_info_file(info_file, seeds=None):
    """
    Parse the functional enrichment file created by GOAtools
    and provide every row of interest.
    """

    node_to_geneIDs = {}
    geneID_to_nodes = {}
    if seeds:
        seeds = set([str(seed) for seed in seeds])

    with open(info_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if fields[0] == 'BIANA ID':
                continue
            node = fields[0]
            geneIDs = fields[4].split(';')

            if seeds:
                if node in seeds:
                    node_to_geneIDs[node] = geneIDs
                    for geneID in geneIDs:
                        geneID_to_nodes.setdefault(geneID, [])
                        geneID_to_nodes[geneID].append(node)
            else:
                node_to_geneIDs[node] = geneIDs
                for geneID in geneIDs:
                    geneID_to_nodes.setdefault(geneID, [])
                    geneID_to_nodes[geneID].append(node)

    return node_to_geneIDs, geneID_to_nodes


def calculate_log_of_odds_ratio(q, k, m, t):
    """Calculates the log of the odds ratio"""
    # print("q: {}".format(q))
    # print("k: {}".format(k))
    # print("m: {}".format(m))
    # print("t: {}".format(t))
    odds_ratio = ( (float(q)/float(k)) / (float(m)/float(t)) )
    #odds_ratio = ( (float(q)/(float(k)-float(q))) / (float(m)/((float(t)-float(m)))) )
    if odds_ratio == 0:
        return -float('inf')
    else:
        return float(math.log(odds_ratio, 2))


if  __name__ == "__main__":
    main()

