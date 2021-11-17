from biana.biana_commands import available_sessions, create_new_session
import sys, os, re
import mysql.connector
import matplotlib.pyplot as plt
import numpy as np
import pylab

def main():
    """
    python /home/quim/PHD/Projects/GUILDify/scripts/get_drugs_from_guildify.py
    """

    BIANA_USER = 'quim'
    BIANA_PASS = ''
    BIANA_HOST = 'localhost'
    BIANA_DATABASE = 'test_BIANA_MAY_2018'
    BIANA_UNIFICATION = 'geneID_seqtax_drugtarget'
    BIANA_ATTRIBUTES = ["genesymbol"]
    node_info_file = '/var/www/html/guildify2/data/BIANA/9606/all/node_info.txt'
    edge_scores_file = '/var/www/html/guildify2/data/BIANA/9606/all/edge_scores.txt'
    drug_info_file = '/var/www/html/guildify2/data/drug_info.txt'
    main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    repur_path = os.path.join(main_path, 'drug_repurposing')

    # Output files
    drugs_file = os.path.join(repur_path, 'drugs_list.txt')
    plot_distribution_targets = os.path.join(repur_path, 'number_targets_per_drug.png')
    target_cutoff = 1 # More than x targets
    drugs_target_cutoff_file = os.path.join(repur_path, 'drugs_list_more_than_{}_targets.txt'.format(target_cutoff))

    tax_id = 9606
    biana_cnx = mysql.connector.connect(user=BIANA_USER, 
                                        password=BIANA_PASS,
                                        host=BIANA_HOST,
                                        database=BIANA_DATABASE)


    # Parse node info file and get all gene symbols
    user_entity_to_info = get_user_entities_info(node_info_file)
    all_entry_ids = set()
    for user_entity_id in user_entity_to_info:
        [entry_ids, gene_symbols, descriptions, gene_ids, other_ids] = user_entity_to_info[user_entity_id]
        for entry_id in entry_ids.split('; '):
            all_entry_ids.add(entry_id.upper())
    print(len(all_entry_ids))

    # Parse the edge file and get the gene symbols
    all_entry_ids_in_network = set()
    user_entity_ids_in_network = set()
    with open(edge_scores_file) as edge_fd:
        for line in edge_fd:
            node1, _, node2 = line.strip().split()
            user_entity_ids_in_network.add(node1)
            user_entity_ids_in_network.add(node2)
            for node in [node1, node2]:
                if node in user_entity_to_info:
                    [entry_ids, gene_symbols, descriptions, gene_ids, other_ids] = user_entity_to_info[node]
                    for entry_id in entry_ids.split('; '):
                        all_entry_ids_in_network.add(entry_id.upper())
    print(len(all_entry_ids_in_network))

    # Parse drug info from GUILDify
    col_to_idx, drug_to_values, uniprot_to_drugs = get_drug_and_target_info(drug_info_file)
    drug_to_targets_in_network = {}
    drug_to_targets = {}
    drug_to_name = {}
    name_to_drug = {}
    for drug in drug_to_values:
        drug_name = drug_to_values[drug][0].lower()
        drug_to_name[drug] = drug_name
        name_to_drug[drug_name] = drug
        targets = drug_to_values[drug][5].split('; ')
        drug_to_targets[drug] = targets
        for target in targets:
            if target in all_entry_ids_in_network:
                drug_to_targets_in_network.setdefault(drug, set()).add(target)

    # Output all the drugs with targets in the network
    with open(drugs_file, 'w') as out_fd, open(drugs_target_cutoff_file, 'w') as out2_fd:
        for drugbank_id,values in sorted(drug_to_targets_in_network.items(), key=lambda x:len(x[1]), reverse=True):
            out_fd.write('{}\t{}\t{}\t{}\n'.format(drugbank_id, drug_to_name[drugbank_id], len(drug_to_targets_in_network[drugbank_id]), len(drug_to_targets[drugbank_id])))
            if len(drug_to_targets_in_network[drugbank_id]) > target_cutoff:
                out2_fd.write('{}\n'.format(drug_to_name[drugbank_id]))
            #print(drug_to_values[drugbank_id][5])

    # Output plot of targets per drug in the network
    targets = [len(x) for x in drug_to_targets_in_network.values()]
    n, bins, patches = plt.hist(np.array(targets), bins=50, weights=np.zeros_like(np.array(targets)) + 1. / np.array(targets).size, facecolor='green')
    plt.xlabel('Number of targets per drug')
    plt.ylabel('Relative frequency')
    plt.title('Distribution of the number of targets per drug')
    plt.savefig(plot_distribution_targets, format='png', dpi=300)
    plt.clf()

    # Parse Hetionet
    hetionet_input_file = os.path.join(main_path, '../DrugRepurposing/data/hetionet/edges.sif')
    hetionet_output_file = os.path.join(repur_path, 'hetionet.txt')
    hetionet_drugs_file = os.path.join(repur_path, 'hetionet_drugs_list.txt')
    hetionet_DOID_to_drugbank = parse_hetionet(hetionet_input_file, hetionet_output_file)
    hetionet_drug_to_targets = {}
    hetionet_drug_to_targets_in_network = {}
    for DOID in hetionet_DOID_to_drugbank:
        for drugbank in hetionet_DOID_to_drugbank[DOID]['all']:
            if drugbank in drug_to_targets_in_network:
                for target in drug_to_targets[drugbank]:
                    hetionet_drug_to_targets.setdefault(drugbank, set()).add(target)
                for target in drug_to_targets_in_network[drugbank]:
                    hetionet_drug_to_targets_in_network.setdefault(drugbank, set()).add(target)
    with open(hetionet_drugs_file, 'w') as out_fd:
        for drugbank_id,values in sorted(hetionet_drug_to_targets_in_network.items(), key=lambda x:len(x[1]), reverse=True):
            out_fd.write('{}\t{}\t{}\t{}\n'.format(drugbank_id, drug_to_name[drugbank_id], len(hetionet_drug_to_targets_in_network[drugbank_id]), len(hetionet_drug_to_targets[drugbank_id])))


    drug_associations_file = os.path.join(repur_path, 'drug_associations.txt')
    drug_associations_sources_file = os.path.join(repur_path, 'drug_associations_and_sources.txt')

    if not fileExist(drug_associations_file) or not fileExist(drug_associations_sources_file):

        # Create BIANA session
        biana_session_id = 'biana_session'
        success = create_new_session(sessionID=biana_session_id, dbname=BIANA_DATABASE, dbhost=BIANA_HOST, dbuser=BIANA_USER, dbpassword=BIANA_PASS, unification_protocol=BIANA_UNIFICATION)
        if success is None: # Mysql related error
            raise ValueError("biana")
        biana_session = available_sessions[biana_session_id]


        # Query GUILDify
        drug_to_genes_in_network = {}
        drug_to_genes_in_network_and_drugdbs = {}
        drug_to_genes = {}
        drug_to_gene_sources = {}
        #tool_drugs = ['calcium', 'copper', 'erlotinib', 'haloperidol']
        #for drug_name in tool_drugs:
        for drug_name in name_to_drug:
            associations, associations_in_network, associations_in_network_and_drugdbs, association_sources = get_disease_gene_associations(drug_name, user_entity_to_info, user_entity_ids_in_network, biana_cnx, biana_session, BIANA_ATTRIBUTES, tax_id)
            drug_to_genes[drug_name] = associations
            drug_to_genes_in_network[drug_name] = associations_in_network
            drug_to_genes_in_network_and_drugdbs[drug_name] = associations_in_network_and_drugdbs
            drug_to_gene_sources[drug_name] = association_sources
            print('Keywords: {}. Num genes in network: {}. Num genes: {}.'.format(drug_name, len(drug_to_genes_in_network[drug_name]), len(drug_to_gene_sources[drug_name])))
            print(drug_to_gene_sources[drug_name])


        # Write the disease and gene associations in a file sorted by DrugBank associations
        with open(drug_associations_file, 'w') as out:
            out.write('Drug name\tDrugBankID\t# genes in drugbank and network\t# genes in drug DBs and network\t# genes in network\t# genes in total\n')
            for drug_name,sources in sorted(drug_to_gene_sources.items(), key=lambda x:x[1]['drugbank'], reverse=True):
                drug = name_to_drug[drug_name]
                drugbank=0
                if drug_name in drug_to_gene_sources:
                    drugbank=drug_to_gene_sources[drug_name]['drugbank']
                out.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(drug_name, drug, drugbank, len(drug_to_genes_in_network_and_drugdbs[drug_name]), len(drug_to_genes_in_network[drug_name]), len(drug_to_genes[drug_name])))

        with open(drug_associations_sources_file, 'w') as out2:
            out2.write('Drug name\tDrugBankID\t# genes in drug DBs and network\t# drugbank\t# dgidb\t# chembl\t# drugcentral\n')
            for drug_name,associations in sorted(drug_to_genes_in_network_and_drugdbs.items(), key=lambda x:len(x[1]), reverse=True):
                drug = name_to_drug[drug_name]
                drugbank=0
                dgidb=0
                chembl=0
                drugcentral=0
                if drug_name in drug_to_gene_sources:
                    drugbank=drug_to_gene_sources[drug_name]['drugbank']
                    dgidb=drug_to_gene_sources[drug_name]['dgidb']
                    chembl=drug_to_gene_sources[drug_name]['chembl']
                    drugcentral=drug_to_gene_sources[drug_name]['drugcentral']
                out2.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(drug_name, drug, len(drug_to_genes_in_network_and_drugdbs[drug_name]), drugbank, dgidb, chembl, drugcentral))


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

def get_drug_and_target_info(drug_info_file):
    """
    Get the drug and target info from the drug_info.txt file.
    """
    uniprot_to_drugs = {}
    col_to_idx, drug_to_values = get_drug_info(drug_info_file)
    for drug, values in drug_to_values.iteritems():
        #print drug, values
        uniprots = values[col_to_idx["targets"]].split("; ")
        for uniprot in uniprots:
            uniprot_to_drugs.setdefault(uniprot, set()).add(drug)
    return col_to_idx, drug_to_values, uniprot_to_drugs

def get_drug_info(drug_info_file):
    drug_to_values = {}
    f = open(drug_info_file)
    header = f.readline().strip().split("\t")
    col_to_idx = dict((k, i) for i, k in enumerate(header[1:])) 
    for line in f:
        words = line.strip("\n").split("\t")
        drug_to_values[words[0]] = words[1:]
    return col_to_idx, drug_to_values

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
    associations_in_network_and_drugdbs = set()
    associations = set()
    association_sources = {'drugbank':0, 'dgidb':0, 'chembl':0, 'drugcentral':0}
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
                if source in association_sources:
                    association_sources[source]+=1
                    associations_in_network_and_drugdbs.add(user_entity_id)
    return associations, associations_in_network, associations_in_network_and_drugdbs, association_sources

def get_associated_user_entities(biana_session, keywords, biana_attributes, tax_id):
    identifier_description_list = [ (attribute, keyword) for attribute in biana_attributes for keyword in keywords ]
    attribute_restriction_list = [] # [("taxid", tax_id)] # slowing down unnecessarily
    # Get group of biomolecules (user entity set) matching with the query
    user_entity_set_id = biana_session._get_next_uEs_id()
    user_entity_set = biana_session.create_new_user_entity_set( identifier_description_list=identifier_description_list, attribute_restriction_list=attribute_restriction_list, id_type="embedded", new_user_entity_set_id=user_entity_set_id, only_uniques=True )
    if user_entity_set is None:
        raise ValueError("biana")
    return map(lambda x: str(int(x)), user_entity_set.get_user_entity_ids()) 

def parse_hetionet(input_file, output_file):
    """
    Parses Hetionet "edges.sif".
    Obtains a file containing the DOID, its corresponding disease name and Palliative/Treatment.
    """
    DOID_to_drugbank = {}
    if not fileExist(output_file):
        with open(input_file, 'r') as input_fd, open(output_file, 'w') as output_fd:
            #source metaedge    target
            first_line = input_fd.readline()
            output_fd.write('DOID\tDrugBankID\tPalliative/Treatment\n')
            for line in input_fd:
                source, metaedge, target = line.strip().split('\t')
                if metaedge == 'CtD': # Treatment
                    #Compound::DB00997  CtD Disease::DOID:363
                    drugbank = source.split('Compound::')[1].upper()
                    DOID = target.split('Disease::')[1].upper()
                    DOID_to_drugbank.setdefault(DOID, {})
                    DOID_to_drugbank[DOID].setdefault('treatment', set())
                    DOID_to_drugbank[DOID]['treatment'].add(drugbank)
                    DOID_to_drugbank[DOID].setdefault('all', set())
                    DOID_to_drugbank[DOID]['all'].add(drugbank)
                    output_fd.write('{}\t{}\t{}\n'.format(DOID, drugbank, 'treatment'))
                elif metaedge == 'CpD': # Palliative
                    #Compound::DB01175  CpD Disease::DOID:3312
                    drugbank = source.split('Compound::')[1].upper()
                    DOID = target.split('Disease::')[1].upper()
                    DOID_to_drugbank.setdefault(DOID, {})
                    DOID_to_drugbank[DOID].setdefault('palliative', set())
                    DOID_to_drugbank[DOID]['palliative'].add(drugbank)
                    DOID_to_drugbank[DOID].setdefault('all', set())
                    DOID_to_drugbank[DOID]['all'].add(drugbank)
                    output_fd.write('{}\t{}\t{}\n'.format(DOID, drugbank, 'palliative'))
    else:
        DOID_to_drugbank = parse_hetionet_output_file(output_file)
    return DOID_to_drugbank

def parse_hetionet_output_file(hetionet_output_file):
    """
    Parses the Hetionet output file and obtains a dictionary.
    """
    DOID_to_drugbank = {}
    with open(hetionet_output_file, 'r') as input_fd:
        first_line = input_fd.readline()
        for line in input_fd:
            DOID, drugbank, type_drug = line.strip().split('\t')
            DOID_to_drugbank.setdefault(DOID, {})
            DOID_to_drugbank[DOID].setdefault('type_drug', set())
            DOID_to_drugbank[DOID]['type_drug'].add(drugbank)
            DOID_to_drugbank[DOID].setdefault('all', set())
            DOID_to_drugbank[DOID]['all'].add(drugbank)
    return DOID_to_drugbank

if  __name__ == "__main__":
    main()