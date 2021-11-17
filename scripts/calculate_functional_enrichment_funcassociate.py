import sys, os, re
import funcassociate.client as client

def main():

    taxID = 9606
    top_genes_file = '/home/quim/PHD/Projects/GUILDify/enrichment_tests/alzheimer_omim/data/guild_scores.txt'
    all_genes_file = '/var/www/html/guildify2/data/BIANA/9606/all/node_info.txt'
    go_obo_file = '/var/www/html/guildify2/data/go/go-basic.obo'
    gene2go_filtered_file = '/var/www/html/guildify2/data/go/gene2go.{}.filtered'.format(taxID)
    temp_file = '/home/quim/PHD/Projects/GUILDify/enrichment_tests/alzheimer_omim/funcassociate/output_enrichment.txt'
    threshold = 1
    evidence_codes_used = ['EXP', 'IDA', 'IMP', 'IGI', 'IEP', 'ISS', 'ISA', 'ISM', 'ISO']
    associations_file = '/home/quim/PHD/Projects/GUILDify/enrichment_tests/alzheimer_omim/data/funcassociate_go_associations.txt'

    if taxID == 9606:
        species = "Homo sapiens"

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

    # associations = []
    # attrib_dict = {}
    # gene_ids_in_associations = []
    # with open(associations_file, 'r') as inp_fd:
    #     for line in inp_fd:
    #         if line.startswith('#'):
    #             continue
    #         fields = line.strip().split('\t')
    #         if len(fields) == 3:
    #             go_id = fields[0]
    #             go_name = fields[1]
    #             gene_ids = fields[2].split(' ')
    #             gene_ids_in_associations=gene_ids_in_associations+gene_ids
    #             association = [go_id]+gene_ids
    #             associations.append(association)
    #             attrib_dict[go_id] = go_name
    #         elif len(fields) == 1:
    #             gene_ids = fields[0].split(' ')
    #             association = [''] + gene_ids
    #         else:
    #             print(fields)
    #             sys.exit(10)

    from goatools.obo_parser import GODag
    from goatools import associations
    obodag = GODag(go_obo_file)
    geneid2gos = associations.read_associations(gene2go_filtered_file, no_top=False)
    go2geneids = associations.get_b2aset(geneid2gos)
    associations = []
    attrib_dict = {}
    gene_ids_in_associations = []
    for go_id in go2geneids:
        gene_ids = list(go2geneids[go_id])
        association = [go_id]+gene_ids
        associations.append(association)
        gene_ids_in_associations=gene_ids_in_associations+gene_ids
        go_name = obodag[go_id].name
        attrib_dict[go_id] = go_name

    new_top_nodes = list(set(top_nodes)&set(gene_ids_in_associations))
    print(len(new_top_nodes))

    #functional_enrichment.check_functional_enrichment(subset_gene_ids=top_nodes, gene_weights=all_nodes, id_type='geneid', output_method=open(temp_file, 'w').write, species = species, mode = "unordered", request_info=False, tex_format=False, support=None, associations=associations)
    which = "over"
    mode = "unordered"
    cutoff = 0.05
    reps = 1000
    client_funcassociate = client.FuncassociateClient()
    response = client_funcassociate.functionate(query=new_top_nodes, associations=associations,
                                                  attrib_dict=attrib_dict, species=None, namespace=None,
                                                  genespace=None, mode=mode, which=which,
                                                  cutoff=cutoff, reps=reps)

    with open(temp_file, 'w') as out_fd:
        headers = [ "# of genes", "# of genes in the query", "# of total genes", "Log of odds ratio", "P-value", "Adjusted p-value", "GO term ID", "Go term name" ]
        out_fd.write('{}\n'.format('\t'.join(headers)))
        for row in response[which]:
            #[2, 5, 4, 3.54982620422173, 7.77737633501871e-07, 0.001, u'GO:0030826', u'regulation of cGMP biosynthetic process']
            n_genes, n_genes_query, n_genes_total, log_odds, p_value, p_adjusted, go_id, go_name = row
            out_fd.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(n_genes, n_genes_query, n_genes_total, log_odds, p_value, p_adjusted, go_id, go_name))

    print(response)

    return


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


if  __name__ == "__main__":
    main()

