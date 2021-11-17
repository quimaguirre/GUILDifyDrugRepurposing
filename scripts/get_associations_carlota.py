import sys, os

def main():

    selected_method = 'GOmf' # GObp / GOmf / Reactome
    #taxID = 9606 # human:9606, mouse:10090, rat: 10116, worm:6239, fly:7227, plant:3702, yeast:4932
    #taxIDs = [9606, 10090, 10116, 6239, 7227, 3702, 4932]
    taxIDs = [9606]
    evidence_codes_go = ['EXP', 'IDA', 'IMP', 'IGI', 'IEP', 'ISS', 'ISA', 'ISM', 'ISO']
    evidence_codes_reactome = ['IEA', 'TAS']
    reactome_file = '/Users/quim/Dropbox/UPF/PhD/Projects/GUILDify/case_studies_JMB/data/NCBI2Reactome.txt'
    data_dir = '/Users/quim/Dropbox/UPF/PhD/Projects/GUILDify/case_studies_JMB/data'
    go_dir = '/Users/quim/Documents/Projects/goatools'
    gene2go_file = os.path.join(go_dir, 'gene2go')
    taxID2species = {
        9606 : 'Homo sapiens',
        10090 : 'Mus musculus',
        10116 : 'Rattus norvegicus',
        6239 : 'Caenorhabditis elegans',
        7227 : 'Drosophila melanogaster',
        3702 : 'Arabidopsis thaliana',
        4932 : 'Saccharomyces cerevisiae'
    }

    for taxID in taxIDs: 

        network_file = os.path.join(data_dir, '{}/network_filtered.sif'.format(taxID))
        node_info_file = os.path.join(data_dir, '{}/node_info.txt'.format(taxID))
        node_to_geneIDs, geneID_to_nodes = get_geneIDs_from_info_file(node_info_file)
        nodes = get_nodes_from_network_file(network_file)
        geneIDs_in_network = set()
        for node in nodes:
            if node in node_to_geneIDs:
                for geneID in node_to_geneIDs[node]:
                    geneIDs_in_network.add(geneID)

        if selected_method == 'Reactome':
            # Evidence codes in Reactome: 'IEA', 'TAS'
            term2genes, term2name = parse_reactome(reactome_file, evidence_codes_reactome, taxID2species[taxID])
        elif selected_method == 'GObp':
            term2genes, term2name = parse_gene2go(go_dir, 'Process', evidence_codes_go, taxID)
        elif selected_method == 'GOmf':
            term2genes, term2name = parse_gene2go(go_dir, 'Function', evidence_codes_go, taxID)
        else:
            print('Incorrect method.')
            sys.exit(10)

        output_file = os.path.join(data_dir, '{}/{}2gene.{}.txt'.format(taxID, selected_method, taxID))
        with open(output_file, 'w') as out_fd:
            for termID in term2genes:
                term_name = term2name[termID]
                for geneID in term2genes[termID]:
                    if geneID in geneIDs_in_network:
                        out_fd.write('{}\t{}\t{}\n'.format(termID, term_name, geneID))

    return


def parse_reactome(reactome_file, evidence_codes, selected_species):
    reactome2genes = {}
    reactome2name = {}
    with open(reactome_file, 'r') as inp_fd:
        for line in inp_fd:
            fields = line.strip().split('\t')
            # 1 R-HSA-114608    https://reactome.org/PathwayBrowser/#/R-HSA-114608  Platelet degranulation  TAS Homo sapiens
            geneID, reactomeID, reactome_url, reactome_name, evidence_code, species_name = fields
            if evidence_code not in evidence_codes:
                continue
            if selected_species == species_name:
                reactome2genes.setdefault(reactomeID, set()).add(geneID)
                reactome2name[reactomeID] = reactome_name
    return reactome2genes, reactome2name


def parse_gene2go(go_dir, type_go, evidence_codes_used, taxID):

    if taxID != 4932:
        go2name = {}
        go2geneids = {}
        gene2go_file = os.path.join(go_dir, 'gene2go')
        #from goatools import associations
        #from goatools.associations import read_ncbi_gene2go
        #geneid2gos = read_ncbi_gene2go(gene2go_file, taxids=[int(taxID)], evidence_set=evidence_codes_used)
        #go2geneids = associations.get_b2aset(geneid2gos)
        with open(gene2go_file, 'r') as inp_fd:
            first_line = inp_fd.readline()
            for line in inp_fd:
                #tax_id GeneID  GO_ID   Evidence    Qualifier   GO_term PubMed  Category
                tax_id, geneID, goID, evidence, qualifier, go_name, pubmed, category  = line.strip().split('\t')
                if category != type_go or tax_id != str(taxID):
                    continue
                if evidence in evidence_codes_used:
                    go2geneids.setdefault(goID, set()).add(geneID)
                    go2name[goID] = go_name
    else:
        taxID_to_gaf_file = {
            9606 : 'goa_human.gaf',
            10090 : 'mgi.gaf',
            10116 : 'rgd.gaf',
            6239 : 'wb.gaf',
            7227 : 'fb.gaf',
            3702 : 'tair.gaf',
            4932 : 'sgd.gaf',
        } 
        taxID_to_geneinfo_file = {
            9606 : 'Homo_sapiens.gene_info',
            10090 : 'Mus_musculus.gene_info',
            10116 : 'Rattus_norvegicus.gene_info',
            6239 : 'Caenorhabditis_elegans.gene_info',
            7227 : 'Drosophila_melanogaster.gene_info',
            3702 : 'Arabidopsis_thaliana.gene_info',
            4932 : 'Saccharomyces_cerevisiae.gene_info',
        } 
        from goatools.obo_parser import GODag
        from goatools.associations import read_gaf
        go_obo_file = os.path.join(go_dir, 'go-basic.obo')
        obodag = GODag(go_obo_file)
        gaf_file = os.path.join(go_dir, 'goa_files/{}'.format(taxID_to_gaf_file[taxID]))
        geneinfo_file = os.path.join(go_dir, 'goa_files/{}'.format(taxID_to_geneinfo_file[taxID]))
        gene2gos = read_gaf(gaf_file, taxids=[int(taxID)], evidence_set=evidence_codes_used)
        gos = set()
        for gene in gene2gos:
            for go in gene2gos[gene]:
                gos.add(go)
        print('Genes: {}, GOs: {}'.format(len(gene2gos), len(gos)))
        #print(gene2gos)

        geneid2gos = {}
        go2name = {}
        with open(geneinfo_file, 'r') as inp_fd:
            first_line = inp_fd.readline()
            for line in inp_fd:
                fields = line.strip().split('\t')
                geneID = fields[1]
                gene = fields[5]
                #print(geneID, gene)
                if ':' in gene:
                    gene = gene.split(':')[1]
                if gene != '-' and geneID != '-':
                    if gene in gene2gos:
                        for go in gene2gos[gene]:
                            geneid2gos.setdefault(geneID, set()).add(go)
                            go_name = obodag[go].name
                            go2name[go] = go_name
                    else:
                        #print('Gene {} not found'.format(gene))
                        pass
        gos = set()
        for gene in geneid2gos:
            for go in geneid2gos[gene]:
                gos.add(go)
        print('GeneIDs: {}, GOs: {}'.format(len(geneid2gos), len(gos)))
        go2geneids = associations.get_b2aset(geneid2gos)
    return go2geneids, go2name


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


def get_nodes_from_network_file(network_file):
    nodes = set()
    with open(network_file, 'r') as net_fd:
        for line in net_fd:
            fields = line.strip().split('\t')
            if len(fields) == 3:
                node1 = fields[0]
                node2 = fields[2]
            else:
                node1 = fields[0]
                node2 = fields[1]
            nodes.add(node1)
            nodes.add(node2)
    return nodes


if  __name__ == "__main__":
    main()