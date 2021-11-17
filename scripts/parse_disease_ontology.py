from biana.biana_commands import available_sessions, create_new_session
import mysql.connector
import numpy as np
import sys, os, re


def main():


    database_dir = '/home/quim/Databases/DiseaseOntology'
    results_dir = '/home/quim/PHD/Projects/GUILDify/analysis_genes'
    # database_dir = '/Users/quim/Documents/Databases/DiseaseOntology'
    # results_dir = '/Users/quim/Dropbox/UPF/PhD/Projects/GUILDify/analysis_genes'
    input_file = os.path.join(database_dir, 'doid.obo')
    diseases_file = os.path.join(results_dir, 'list_of_lowest_level_non_obsolete_disease_names.txt')
    tax_id = 9606
    BIANA_ATTRIBUTES = ["genesymbol"]

    GUILDIFY_VERSION = 1
    guildify_genes_file = os.path.join(results_dir, 'GUILDify{}_genes.txt'.format(str(GUILDIFY_VERSION)))
    guildify1_genes_file = os.path.join(results_dir, 'GUILDify1_genes.txt')
    guildify2_genes_file = os.path.join(results_dir, 'GUILDify2_genes.txt')
    if GUILDIFY_VERSION == 1:
        BIANA_USER = 'jgarcia'
        BIANA_PASS = ''
        BIANA_HOST = 'sbi.upf.edu'
        BIANA_DATABASE = 'test_GUILDIFY_2013'
        BIANA_UNIFICATION = 'uniprot_seq_geneID'
        node_info_file = '/home/quim/PHD/Projects/GUILDify/analysis_genes/node_info_GUILDIFY1.txt'
    else:
        BIANA_USER = 'quim'
        BIANA_PASS = ''
        BIANA_HOST = 'localhost'
        BIANA_DATABASE = 'test_BIANA_MAY_2018'
        BIANA_UNIFICATION = 'geneID_seqtax_drugtarget'
        node_info_file = '/var/www/html/GUILDify2/data/BIANA/9606/all/node_info.txt'


    # Create BIANA connection
    biana_cnx = mysql.connector.connect(user=BIANA_USER, 
                                        password=BIANA_PASS,
                                        host=BIANA_HOST,
                                        database=BIANA_DATABASE)


    # Create BIANA session
    biana_session_id = 'biana_session'
    success = create_new_session(sessionID=biana_session_id, dbname=BIANA_DATABASE, dbhost=BIANA_HOST, dbuser=BIANA_USER, dbpassword=BIANA_PASS, unification_protocol=BIANA_UNIFICATION)
    if success is None: # Mysql related error
        raise ValueError("biana")
    biana_session = available_sessions[biana_session_id]

    # Parse node info file
    user_entity_to_info = get_user_entities_info(node_info_file)

    if not fileExist(diseases_file):

        # Parse Disease Ontology
        term_id_to_info, is_a_dict = parse_disease_ontology(input_file)

        # Get non-obsolete lowest level disease names
        all_term_ids = set(term_id_to_info.keys())
        higher_level_term_ids = set(is_a_dict.values())
        print('Number of term IDs: {}'.format(len(all_term_ids)))
        non_obsolete_term_ids = set([term for term in all_term_ids if term_id_to_info[term]['term_is_obsolete'] == False])
        print('Number of non-obsolete term IDs: {}'.format(len(non_obsolete_term_ids)))
        lowest_level_term_ids = set([term for term in all_term_ids if term not in higher_level_term_ids])
        print('Number of term IDs in the lowest level: {}'.format(len(lowest_level_term_ids)))
        lowest_level_non_obsolete_term_ids = lowest_level_term_ids & non_obsolete_term_ids
        print('Number of non-obsolete term IDs in the lowest level: {}'.format(len(lowest_level_non_obsolete_term_ids)))
        #print(lowest_level_non_obsolete_term_ids)
        lowest_level_non_obsolete_names = sorted([term_id_to_info[term]['term_name'].lower() for term in lowest_level_non_obsolete_term_ids])
        #print(lowest_level_non_obsolete_names)

        # Write output file with names
        with open(diseases_file, 'w') as diseases_fd:
            for name in lowest_level_non_obsolete_names:
                diseases_fd.write('{}\n'.format(name))

    else:

        # Get non-obsolete lowest level disease names from file
        lowest_level_non_obsolete_names = []
        with open(diseases_file, 'r') as diseases_fd:
            for line in diseases_fd:
                lowest_level_non_obsolete_names.append(line.strip())


    # Query GUILDify
    if not fileExist(guildify_genes_file):
        with open(guildify_genes_file, 'w') as guildify_fd:
            guildify_fd.write('Keywords\tNum genes\n')
            for disease in lowest_level_non_obsolete_names:
                keywords = '"{}"'.format(disease)
                keywords = keywords.encode("utf-8")
                keywords, matched = tokenize_text(keywords)
                original_keywords = keywords
                keywords, gene_to_desc_id, desc_id_desc_and_src = get_associated_gene_symbols_and_descriptions_from_guildify_tables(biana_cnx, keywords, tax_id, matched)
                if not keywords:
                    genes = []
                else:
                    try:
                        user_entity_ids = list(set(get_associated_user_entities(biana_session, keywords, BIANA_ATTRIBUTES, tax_id)))
                        if not user_entity_ids:
                            genes = []
                        else:
                            user_entity_to_values = {}
                            for user_entity_id in user_entity_ids:
                                if user_entity_id in user_entity_to_info:
                                    user_entity_to_values[user_entity_id] = user_entity_to_info[user_entity_id]
                            if len(user_entity_to_values) == 0:
                                genes = []
                            else:
                                genes = sorted(user_entity_to_values.keys())
                                #print(user_entity_to_values)
                    except:
                        genes = []
                guildify_fd.write('{}\t{}\n'.format(original_keywords[0], len(genes)))
                print('Keywords: {}. Num genes: {}'.format(original_keywords, len(genes)))


    trial = 'oculocutaneous albinism type iv'
    keywords = '"{}"'.format(trial)
    keywords = keywords.encode("utf-8")
    keywords, matched = tokenize_text(keywords)
    original_keywords = keywords
    keywords, gene_to_desc_id, desc_id_desc_and_src = get_associated_gene_symbols_and_descriptions_from_guildify_tables(biana_cnx, keywords, tax_id, matched)
    user_entity_ids = list(set(get_associated_user_entities(biana_session, keywords, BIANA_ATTRIBUTES, tax_id)))
    user_entity_to_values = get_user_entities_info(node_info_file, user_entity_ids=user_entity_ids)
    print(original_keywords)
    print(keywords)
    print(gene_to_desc_id)
    #print(desc_id_desc_and_src)
    print(user_entity_to_values)
    for ue in user_entity_to_values:
        print(ue, user_entity_to_values[ue])
    print(len(user_entity_to_values.keys()))

    # Close BIANA connexion
    biana_cnx.close()


    # Gather information from queries in GUILDify v1 and v2
    if fileExist(guildify1_genes_file) and fileExist(guildify2_genes_file):

        guildify1_dict = read_guildify_genes_file(guildify1_genes_file)
        guildify2_dict = read_guildify_genes_file(guildify2_genes_file)

        diseases = set(guildify1_dict.keys()) | set(guildify2_dict.keys())
        
        # Compare two GUILDify
        comparison_numbers_file = os.path.join(results_dir, 'comparison_numbers.txt')
        with open(comparison_numbers_file, 'w') as comparison_fd:
            for disease in diseases:
                #comparison_fd.write('GUILDify V1: {}\tGUILDify V2: {}\tDisease: {}\n'.format(guildify1_dict[disease], guildify2_dict[disease], disease))
                comparison_fd.write('{}\t{}\t{}\n'.format(disease, guildify1_dict[disease], guildify2_dict[disease]))

        # Generate some numbers
        diseases_with_genes = set()
        guildify1_diseases = set()
        guildify2_diseases = set()
        guildify1_2_diseases = set()
        guildify1_genes = 0
        guildify2_genes = 0
        more_in_2 = 0
        more_in_1 = 0
        equal = 0
        guildify1_numbers = []
        guildify2_numbers = []
        for disease in diseases:
            if guildify1_dict[disease] > 0:
                guildify1_diseases.add(disease)
                guildify1_numbers.append(guildify1_dict[disease])
            if guildify2_dict[disease] > 0:
                guildify2_diseases.add(disease)
                guildify2_numbers.append(guildify2_dict[disease])
            if guildify1_dict[disease] > 0 or guildify2_dict[disease] > 0:
                diseases_with_genes.add(disease)
                guildify1_genes += guildify1_dict[disease]
                guildify2_genes += guildify2_dict[disease]
                if guildify2_dict[disease] > guildify1_dict[disease]:
                    more_in_2 += 1
                elif guildify2_dict[disease] < guildify1_dict[disease]:
                    more_in_1 += 1
                else:
                    equal += 1
            if guildify1_dict[disease] > 0 and guildify2_dict[disease] > 0:
                guildify1_2_diseases.add(disease)
        print('Total of diseases: {}'.format(len(diseases)))
        print('Total of diseases with genes: {}'.format(len(diseases_with_genes)))
        print('Total of diseases with genes in GUILDify V1: {}'.format(len(guildify1_diseases)))
        print('Total of diseases with genes in GUILDify V2: {}'.format(len(guildify2_diseases)))
        print('Total of diseases with genes in both GUILDifies: {}'.format(len(guildify1_2_diseases)))
        print('GUILDify V1 - num genes: {}'.format(guildify1_genes))
        print('GUILDify V2 - num genes: {}'.format(guildify2_genes))
        print('GUILDify V1 - median: {} - mean: {}'.format(np.median(guildify1_numbers), np.mean(guildify1_numbers)))
        print('GUILDify V2 - median: {} - mean: {}'.format(np.median(guildify2_numbers), np.mean(guildify2_numbers)))
        print('More in 1: {} / More in 2: {} / Equal: {}'.format(more_in_1, more_in_2, equal))

    return


def parse_disease_ontology(input_file):

    # Start variables
    term_id_to_info = {}
    is_a_dict = {}
    term_id = None
    term_name = None
    term_def = None
    term_namespace = None 
    term_is_a = []
    term_part_of = []
    term_exact_synonyms = []
    term_related_synonyms = []
    term_broad_synonyms = []
    term_narrow_synonyms = []
    term_alt_id = []
    term_is_obsolete = False

    with open(input_file, 'r') as input_file_fd:
        for line in input_file_fd:

            # Quim Aguirre: I have included to recognise [Typedef], so that the [Term] entries are recorded well when they are finished and there is a [Typedef] afterwards
            if re.search("\[Term\]",line) or re.search("\[Typedef\]",line):

                # New term
                if term_id is not None:

                    # insert previous
                    #print(term_id)
                    #print(term_name)
                    term_id_to_info.setdefault(term_id, {})
                    term_id_to_info[term_id]['term_name'] = term_name
                    term_id_to_info[term_id]['term_def'] = term_def
                    term_id_to_info[term_id]['term_namespace'] = term_namespace
                    term_id_to_info[term_id]['term_is_a'] = term_is_a
                    term_id_to_info[term_id]['term_part_of'] = term_part_of
                    term_id_to_info[term_id]['term_exact_synonyms'] = term_exact_synonyms
                    term_id_to_info[term_id]['term_related_synonyms'] = term_related_synonyms
                    term_id_to_info[term_id]['term_broad_synonyms'] = term_broad_synonyms
                    term_id_to_info[term_id]['term_narrow_synonyms'] = term_narrow_synonyms
                    term_id_to_info[term_id]['term_alt_id'] = term_alt_id
                    term_id_to_info[term_id]['term_is_obsolete'] = term_is_obsolete

                    # Annotate the is_a in the dictionary
                    for is_a in term_is_a:
                        is_a_dict[term_id] = is_a

                # Restart variables
                term_id = None
                term_name = None
                term_def = None
                term_namespace = None 
                term_is_a = []
                term_part_of = []
                term_exact_synonyms = []
                term_related_synonyms = []
                term_broad_synonyms = []
                term_narrow_synonyms = []
                term_alt_id = []
                term_is_obsolete = False

                if re.search("\[Typedef\]",line):
                    typedef = True
                else:
                    typedef = False


            elif re.search("^id\:",line):
                if typedef == True: # If typedef tag is true, we do not want to record anything
                    continue

                temp = re.search("(DOID\:\d+)",line)
                
                if temp:
                    term_id = temp.group(1)

            elif re.search("^name\:",line):
                if typedef == True:
                    continue

                temp = re.search("name:\s+(.+)",line)
                term_name = temp.group(1)

            elif re.search("^namespace\:",line):
                if typedef == True:
                    continue

                temp = re.search("namespace:\s+(.+)",line)
                term_namespace = temp.group(1)

            elif re.search("^def\:",line):
                if typedef == True:
                    continue

                temp = re.search("\"(.+)\"",line)
                term_def = temp.group(1)

            elif re.search("synonym\:",line):
                if typedef == True:
                    continue

                temp = re.search("\"(.+)\"\s+(\w+)",line)
                if temp.group(2) == "EXACT":
                    term_exact_synonyms.append(temp.group(1))
                elif temp.group(2) == "RELATED":
                    term_related_synonyms.append(temp.group(1))
                # Quim Aguirre: I have added the broad and narrow synonyms
                elif temp.group(2) == "BROAD":
                    term_broad_synonyms.append(temp.group(1))
                elif temp.group(2) == "NARROW":
                    term_narrow_synonyms.append(temp.group(1))

            elif re.search("^alt_id\:",line):
            # Quim Aguirre: Recognison of the "alt_id" tags
            # Example --> alt_id: DOID:267
                if typedef == True:
                    continue

                temp = re.search("(DOID\:\d+)",line)
                
                if temp:
                    term_alt_id.append(temp.group(1))

            elif re.search("is_a\:",line):
                if typedef == True:
                    continue

                temp = re.search("(DOID\:\d+)",line)
                if temp is not None:
                    #print "??:", line # malformation --> is_a: regulates ! regulates
                    term_is_a.append(temp.group(1))

            elif re.search("relationship\:",line):
                if typedef == True:
                    continue

                if( re.search("part_of",line) ):
                    temp = re.search("part_of\s+(DOID\:\d+)",line)
                    if temp is not None:
                        term_part_of.append(temp.group(1))

            elif re.search("is_obsolete\:",line):
                if typedef == True:
                    continue

                temp = re.search("is_obsolete:\s+(.+)",line)
                is_obsolete = temp.group(1)
                if is_obsolete == 'true':
                    term_is_obsolete = True

    # Insert last term
    if term_id is not None:

        #print(term_id)
        #print(term_name)
        term_id_to_info.setdefault(term_id, {})
        term_id_to_info[term_id]['term_name'] = term_name
        term_id_to_info[term_id]['term_def'] = term_def
        term_id_to_info[term_id]['term_namespace'] = term_namespace
        term_id_to_info[term_id]['term_is_a'] = term_is_a
        term_id_to_info[term_id]['term_part_of'] = term_part_of
        term_id_to_info[term_id]['term_exact_synonyms'] = term_exact_synonyms
        term_id_to_info[term_id]['term_related_synonyms'] = term_related_synonyms
        term_id_to_info[term_id]['term_broad_synonyms'] = term_broad_synonyms
        term_id_to_info[term_id]['term_narrow_synonyms'] = term_narrow_synonyms
        term_id_to_info[term_id]['term_alt_id'] = term_alt_id
        term_id_to_info[term_id]['term_is_obsolete'] = term_is_obsolete

        # Annotate the is_a in the dictionary
        for is_a in term_is_a:
            is_a_dict[term_id] = is_a


    return term_id_to_info, is_a_dict


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


def read_guildify_genes_file(input_file):
    """
    Read the file containing the number of genes per disease from GUILDify.
    """
    disease_to_num_genes = {}
    with open(input_file, 'r') as input_file_fd:
        first_line = input_file_fd.readline()
        for line in input_file_fd:
            keywords, num_genes = line.strip().split('\t')
            disease_to_num_genes[keywords] = int(num_genes)
    return disease_to_num_genes


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


if  __name__ == "__main__":
    main()
