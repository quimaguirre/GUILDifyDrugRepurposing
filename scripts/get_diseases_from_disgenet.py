from biana.biana_commands import available_sessions, create_new_session
import sys, os, re
import mysql.connector
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pylab

def main():

    BIANA_USER = 'quim'
    BIANA_PASS = ''
    BIANA_HOST = 'localhost'
    BIANA_DATABASE = 'test_BIANA_MAY_2018'
    BIANA_UNIFICATION = 'geneID_seqtax_drugtarget'
    BIANA_ATTRIBUTES = ["genesymbol"]
    node_info_file = '/var/www/html/guildify2/data/BIANA/9606/all/node_info.txt'
    edge_scores_file = '/var/www/html/guildify2/data/BIANA/9606/all/edge_scores.txt'
    main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    repur_path = os.path.join(main_path, 'drug_repurposing')
    disease_names_file = os.path.join(repur_path, 'disease_names.txt')
    disease_associations_file = os.path.join(repur_path, 'disease_associations.txt')
    disease_associations_sources_file = os.path.join(repur_path, 'disease_associations_and_sources.txt')
    plot_distribution_disease_genes = os.path.join(repur_path, 'number_genes_per_disease.png')
    curated_sources = ['uniprot', 'psygenet', 'orphanet', 'hpo'] # We avoid ctd_human because it uses experimental phenotypes such as "liver cirrhosis, experimental"
    tax_id = 9606

    biana_cnx = mysql.connector.connect(user=BIANA_USER, 
                                        password=BIANA_PASS,
                                        host=BIANA_HOST,
                                        database=BIANA_DATABASE)


    # From DisGeNET, get the disease names that are associated with curated databases (except CTD_human which contains experimental terms).
    if not fileExist(disease_names_file):
        umls_to_name = parse_diseases_from_disgenet(biana_cnx, disease_names_file, curated_sources)
    else:
        umls_to_name = {}
        with open(disease_names_file, 'r') as inp_fd:
            for line in inp_fd:
                umls, name = line.strip().split('\t')
                umls_to_name[umls] = name
    print('Number of diseases: {}.'.format(len(umls_to_name)))


    # Create BIANA session
    biana_session_id = 'biana_session'
    success = create_new_session(sessionID=biana_session_id, dbname=BIANA_DATABASE, dbhost=BIANA_HOST, dbuser=BIANA_USER, dbpassword=BIANA_PASS, unification_protocol=BIANA_UNIFICATION)
    if success is None: # Mysql related error
        raise ValueError("biana")
    biana_session = available_sessions[biana_session_id]


    # Parse node info file
    user_entity_to_info = get_user_entities_info(node_info_file)
    print('Number of user entities in total: {}.'.format(len(user_entity_to_info)))


    # Parse the edge file
    user_entity_ids_in_network = set()
    with open(edge_scores_file) as edge_fd:
        for line in edge_fd:
            node1, _, node2 = line.strip().split()
            user_entity_ids_in_network.add(node1)
            user_entity_ids_in_network.add(node2)
    print('Number of user entities in the network: {}.'.format(len(user_entity_ids_in_network)))


    # Query GUILDify
    disease_to_genes_in_network = {}
    disease_to_genes = {}
    disease_to_gene_sources = {}
    #tool_diseases = {'C0036341':'schizophrenia', 'C0023890':'liver cirrhosis', 'C0014544':'epilepsy', 'C0023418':'leukemia', 'C0003873':'rheumatoid arthritis'}
    #for disease in tool_diseases:
    for disease in umls_to_name:
        associations, associations_in_network, association_sources = get_disease_gene_associations(disease, user_entity_to_info, user_entity_ids_in_network, biana_cnx, biana_session, BIANA_ATTRIBUTES, tax_id)
        disease_to_genes[disease] = associations
        disease_to_genes_in_network[disease] = associations_in_network
        disease_to_gene_sources[disease] = association_sources
        print('Keywords: {}. Num genes in network: {}. Num genes: {}.'.format(disease, len(disease_to_genes_in_network[disease]), len(disease_to_genes[disease])))
        print(disease_to_gene_sources[disease])


    # Write the disease and gene associations in a file sorted by DisGeNET associations
    with open(disease_associations_file, 'w') as out, open(disease_associations_sources_file, 'w') as out2:
        out.write('Disease\tDisease name\t# genes in disgenet and network\t# genes in network\t# genes in total\n')
        out2.write('Disease\tDisease name\t# genes in network\t# disgenet\t# omim\t# uniprot\t# go\n')
        for disease,sources in sorted(disease_to_gene_sources.items(), key=lambda x:x[1]['disgenet'], reverse=True):
            disease_name = '-'
            if disease in umls_to_name:
                disease_name = umls_to_name[disease]
            disgenet=0
            omim=0
            uniprot=0
            go=0
            if disease in disease_to_gene_sources:
                disgenet=disease_to_gene_sources[disease]['disgenet']
                omim=disease_to_gene_sources[disease]['omim']
                uniprot=disease_to_gene_sources[disease]['uniprot']
                go=disease_to_gene_sources[disease]['go']
            out.write('{}\t{}\t{}\t{}\t{}\n'.format(disease, disease_name, disgenet, len(disease_to_genes_in_network[disease]), len(disease_to_genes[disease])))
            out2.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(disease, disease_name, len(disease_to_genes_in_network[disease]), disgenet, omim, uniprot, go))


    # # Write the disease and gene associations in a file
    # with open(disease_associations_file, 'w') as out, open(disease_associations_sources_file, 'w') as out2:
    #     out.write('Disease\tDisease name\t# genes in network\t# genes in total\n')
    #     out2.write('Disease\tDisease name\t# genes in network\t# disgenet\t# omim\t# uniprot\t# go\n')
    #     for disease,genes_in_network in sorted(disease_to_genes_in_network.items(), key=lambda x:len(x[1]), reverse=True):
    #         disease_name = '-'
    #         if disease in umls_to_name:
    #             disease_name = umls_to_name[disease]
    #         disgenet=0
    #         omim=0
    #         uniprot=0
    #         go=0
    #         if disease in disease_to_gene_sources:
    #             disgenet=disease_to_gene_sources[disease]['disgenet']
    #             omim=disease_to_gene_sources[disease]['omim']
    #             uniprot=disease_to_gene_sources[disease]['uniprot']
    #             go=disease_to_gene_sources[disease]['go']
    #         out.write('{}\t{}\t{}\t{}\n'.format(disease, disease_name, len(genes_in_network), len(disease_to_genes[disease])))
    #         out2.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(disease, disease_name, len(genes_in_network), disgenet, omim, uniprot, go))

    genes_network = [len(x) for x in disease_to_genes_in_network.values()]
    n, bins, patches = plt.hist(np.array(genes_network), bins=50, weights=np.zeros_like(np.array(genes_network)) + 1. / np.array(genes_network).size, facecolor='r')
    plt.xlabel('Number of genes per disease')
    plt.ylabel('Relative frequency')
    plt.title('Distribution of the number of disease genes per disease')
    plt.savefig(plot_distribution_disease_genes, format='png', dpi=300)
    plt.clf()


    biana_cnx.close()

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

def parse_diseases_from_disgenet(biana_cnx, disease_names_file, sources):
    """
    Parse diseases from DisGeNET database and create a file with these diseases.
    """
    umls_to_name = {}
    cursor = biana_cnx.cursor()
    query1 = ("""SELECT G.value, GS.value, D.value, DN.value, DT.value, Sc.value, So.value 
                 FROM externalEntityRelation R, externalEntityRelationParticipant R1, externalEntityRelationParticipant R2, externalEntityGeneID G, externalEntityGeneSymbol GS, externalEntityUMLS_CUI D, externalEntityDisGeNET_disease_type DT, externalEntityName DN, externalEntityDisGeNET_score Sc, externalEntityDisGeNET_source So 
                 WHERE R.externalEntityRelationID = R1.externalEntityRelationID AND R.externalEntityRelationID = R2.externalEntityRelationID AND R1.externalEntityID != R2.externalEntityID AND R.type = 'gene_disease_association' AND R1.externalEntityID = G.externalEntityID AND G.externalEntityID = GS.externalEntityID AND R2.externalEntityID = D.externalEntityID AND D.externalEntityID = DN.externalEntityID AND D.externalEntityID = DT.externalEntityID AND R.externalEntityRelationID = Sc.externalEntityID AND R.externalEntityRelationID = So.externalEntityID""")
    cursor.execute(query1)
    for items in cursor:
        geneid, gene_symbol, umls, disease_name, disease_type, dg_score, dg_source = items
        if dg_source.lower() in sources and disease_type == 'disease':
            umls_to_name[umls] = disease_name
    cursor.close()

    # Write the disease names in a file
    with open(disease_names_file, 'w') as out_fd:
        for umls in sorted(umls_to_name):
            name = umls_to_name[umls]
            out_fd.write('{}\t{}\n'.format(umls, name))

    return umls_to_name

def get_disease_gene_associations(disease, user_entity_to_info, user_entity_ids_in_network, biana_cnx, biana_session, BIANA_ATTRIBUTES, tax_id):

    keywords = '"{}"'.format(disease)
    keywords = keywords.encode("utf-8")
    keywords, matched = tokenize_text(keywords)
    original_keywords = keywords

    # Use guildify_descriptions table to get associated gene_symbols
    keywords, gene_to_desc_id, desc_id_desc_and_src = get_associated_gene_symbols_and_descriptions_from_guildify_tables(biana_cnx, keywords, tax_id, matched)
    if not keywords:
        user_entity_ids_sorted_by_gene_symbol = []
    else:
        try:
            # Get user entitites corresponding to this gene symbols
            user_entity_ids = list(set(get_associated_user_entities(biana_session, keywords, BIANA_ATTRIBUTES, tax_id)))
            if not user_entity_ids:
                user_entity_ids_sorted_by_gene_symbol = []
            else:
                user_entity_to_values = {}
                for user_entity_id in user_entity_ids:
                    if user_entity_id in user_entity_to_info:
                        user_entity_to_values[user_entity_id] = user_entity_to_info[user_entity_id]
                if len(user_entity_to_values) == 0:
                    user_entity_ids_sorted_by_gene_symbol = []
                else:
                    # Sort user entities by gene symbol
                    user_entity_ids_sorted_by_gene_symbol = []
                    for user_entity_id, value in user_entity_to_values.iteritems():
                        entry_ids, gene_symbols, descriptions, gene_ids, other_ids = value 
                        user_entity_ids_sorted_by_gene_symbol.append((gene_symbols.split(";")[0], user_entity_id))
                    user_entity_ids_sorted_by_gene_symbol.sort()
                    user_entity_ids_sorted_by_gene_symbol = zip(*user_entity_ids_sorted_by_gene_symbol)[1]
        except:
            user_entity_ids_sorted_by_gene_symbol = []

    associations_in_network = set()
    associations = set()
    association_sources = {'disgenet':0, 'omim':0, 'uniprot':0, 'go':0}
    for user_entity_id in user_entity_ids_sorted_by_gene_symbol:
        entry_ids, gene_symbols, descriptions, gene_ids, other_ids = user_entity_to_values[user_entity_id]
        srcs = set()
        skip_flag = False
        for entry_id in entry_ids.split("; "):
            if entry_id == "-":
                skip_flag = True
            genes = gene_symbols.split('; ')
            gene = ''
            if len(genes)>1:
                # If more than one gene symbol, we check that one of them is among the keywords
                # If so, we use this gene to find the description
                for gene_symbol in genes:
                    if gene_symbol in keywords:
                        gene = gene_symbol
            else:
                gene = gene_symbols
            # We skip the gene if it is not among keywords (and therefore it won't have a description)
            if gene not in keywords:
                skip_flag = True
            if gene in gene_to_desc_id:
                desc_ids = gene_to_desc_id[gene]
                for desc_id in desc_ids:
                    for desc, src in desc_id_desc_and_src[desc_id]:
                        srcs.add(src.split(' ')[0]) # get the sources and if it is uniprot swissprot, put only uniprot
        if skip_flag:
            continue
        associations.add(user_entity_id)
        if user_entity_id in user_entity_ids_in_network:
            associations_in_network.add(user_entity_id)
            for source in srcs:
                association_sources[source]+=1
    return associations, associations_in_network, association_sources

def get_user_entities_info(file_name, user_entity_ids = None):
    """
    Get the information of the user entities from 'file_name'
    """
    user_entity_to_values = {}
    # Get information of these user entities
    if user_entity_ids is not None:
        user_entity_ids = set(user_entity_ids)
    header = True
    for line in open(file_name):
        if header:
            header = False
            continue
        #print(line)
        fields = line.strip().split('\t')
        #print(len(fields))
        user_entity_id, entry_ids, gene_symbols, descriptions, gene_ids, other_ids = line.strip().split("\t")

        # Remove isoforms
        entry_ids = remove_isoforms(entry_ids)

        if user_entity_ids is not None:
            if user_entity_id not in user_entity_ids:
                continue
        user_entity_to_values[user_entity_id] = [entry_ids, gene_symbols, descriptions, gene_ids, other_ids]
    return user_entity_to_values 

def remove_isoforms(entry_ids):
    """
    Remove the isoforms from a list of uniprot accessions
    """
    # Check how many non-isoforms are there
    non_isoforms = [ entry for entry in entry_ids.split('; ') if len(entry.split('-')) == 1 ]

    # Remove the isoforms
    old_entry_ids = entry_ids.split('; ')
    new_entry_ids = set()
    for entry in old_entry_ids:
        if len(non_isoforms) > 0:
            if len(entry.split('-')) == 1:
                new_entry_ids.add(entry)
        else:
            new_entry_ids.add(entry)
    entry_ids = '; '.join(list(new_entry_ids))
    return entry_ids

def tokenize_text(text):
    """
    Process text.
    """
    matched = True
    current = ""
    keywords = []
    for c in text:
        if c == '"':
            if matched:
                matched = False
                #current += "\""
            else:
                matched = True 
                #current += "\""
                keywords.append(current)
                current = ""
        elif c == " ":
            if matched:
                if current != "":
                    keywords.append(current)
                    current = ""
            else:
                current += c
        else:
            if c == "'":
                c = "\\'"
            current += c
    if current != "":
        if not matched:
            #current += "\""
            pass
        keywords.append(current)
    text = text.lstrip()
    if text[0] not in ("'", '"'):
        matched = False
    return keywords, matched

def get_associated_gene_symbols_and_descriptions_from_guildify_tables(biana_cnx, keywords, tax_id, matched=False):
    """
    Get the information of disease-gene associations from GUILDify tables of BIANA.
    """
    #print "keywords:", keywords
    cursor = biana_cnx.cursor() # Start cursor to MySQL

    if matched:
        query = "SELECT desc_id, description, source from (%s) T" % ("SELECT desc_id, description, source FROM guildify_descriptions WHERE MATCH(description) AGAINST('\"%s\"' IN BOOLEAN MODE)" % " ".join(keywords))
    else:
        query = "SELECT desc_id, description, source from (%s) T" % ("SELECT desc_id, description, source FROM guildify_descriptions WHERE MATCH(description) AGAINST(\"%s\" IN BOOLEAN MODE)" % " ".join(keywords))
    #print query

    cursor.execute(query)

    result = cursor.fetchall()

    if len(result) == 0: # No match for the keyword(s)
        return None, None, None
    desc_id_to_desc_and_src = {}
    for row in result:
        desc_id_to_desc_and_src.setdefault(row[0], []).append((row[1], row[2]))
    query = "SELECT gene_symbol, desc_id FROM guildify_genes WHERE tax_id = %s AND desc_id IN (%s)" % (tax_id, ",".join(map(str, desc_id_to_desc_and_src.keys()))) #zip(*result_desc)[0])))
    #print(query)
    cursor.execute (query)
    result = cursor.fetchall()
    gene_to_desc_id = {}
    for row in result:
        gene_to_desc_id.setdefault(row[0], []).append(row[1])
    genes = set(gene_to_desc_id.keys()) #set(zip(*result)[0])
    #print "genes: ", 
    #for gene in genes: print gene,
    #print ""

    # Add the keyword just in case that there is a gene symbol in the keyword
    for keyword in keywords:
        keyword = keyword.upper()
        if keyword not in gene_to_desc_id and keyword.lower() not in gene_to_desc_id:
            genes.add(keyword.upper())

    cursor.close()

    return genes, gene_to_desc_id, desc_id_to_desc_and_src


def get_associated_user_entities(biana_session, keywords, biana_attributes, tax_id):
    identifier_description_list = [ (attribute, keyword) for attribute in biana_attributes for keyword in keywords ]
    attribute_restriction_list = [] # [("taxid", tax_id)] # slowing down unnecessarily
    # Get group of biomolecules (user entity set) matching with the query
    user_entity_set_id = biana_session._get_next_uEs_id()
    user_entity_set = biana_session.create_new_user_entity_set( identifier_description_list=identifier_description_list, attribute_restriction_list=attribute_restriction_list, id_type="embedded", new_user_entity_set_id=user_entity_set_id, only_uniques=True )
    if user_entity_set is None:
        raise ValueError("biana")
    return map(lambda x: str(int(x)), user_entity_set.get_user_entity_ids()) 


if  __name__ == "__main__":
    main()