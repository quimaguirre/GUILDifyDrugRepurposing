import os, sys, re
from goatools.associations import read_ncbi_gene2go, read_gaf

def main():
    """
    Read the NCBI gene2go file, parse it by taxID and evidence codes
    and print a simpler version of it, as follows:
    12  GO:0005615
    13  GO:0010898;GO:0019213;GO:0017171;GO:0004806;GO:0005789
    14  GO:0014909;GO:0010595;GO:0005829;GO:0008201;GO:0015630;GO:0045171;GO:0009986
    15  GO:0030187;GO:0007623;GO:0005829;GO:0071320;GO:0004059;GO:0006474;GO:0048471
    16  GO:0004813;GO:0002161;GO:0005829;GO:0006419;GO:0006400;GO:0002196
    """

    taxID = 9606 # human:9606, mouse:10090, rat: 10116, worm:6239, fly:7227, plant:3702, yeast:4932
    by_goa = False
    go_dir = '/var/www/html/guildify2/data/go'
    gene2go_file = os.path.join(go_dir, 'gene2go')
    gene2go_filtered_file = os.path.join(go_dir, 'gene2go.{}.filtered'.format(taxID))
    evidence_codes_used = ['EXP', 'IDA', 'IMP', 'IGI', 'IEP', 'ISS', 'ISA', 'ISM', 'ISO']

    if by_goa == False and taxID != 4932:

        geneid2gos = read_ncbi_gene2go(gene2go_file, taxids=[int(taxID)], evidence_set=evidence_codes_used)

        with open(gene2go_filtered_file, 'w') as out_fd:
            for geneid in geneid2gos:
                gos = ';'.join(geneid2gos[geneid])
                out_fd.write('{}\t{}\n'.format(geneid, gos))

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
                    else:
                        #print('Gene {} not found'.format(gene))
                        pass

        with open(gene2go_filtered_file, 'w') as out_fd:
            for geneid in geneid2gos:
                gos = ';'.join(geneid2gos[geneid])
                out_fd.write('{}\t{}\n'.format(geneid, gos))

        gos = set()
        for gene in geneid2gos:
            for go in geneid2gos[gene]:
                gos.add(go)
        print('GeneIDs: {}, GOs: {}'.format(len(geneid2gos), len(gos)))

    return

if  __name__ == "__main__":
    main()
