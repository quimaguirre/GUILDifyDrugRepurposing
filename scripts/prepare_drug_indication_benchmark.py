import sys, os, re
import mysql.connector
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pylab

def main():

    main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    repur_path = os.path.join(main_path, 'drug_repurposing')
    disgenet_mappings_file = '/home/quim/Databases/disgenet/disease_mappings.tsv'
    hetionet_input_file = os.path.join(main_path, '../DrugRepurposing/data/hetionet/edges.sif')
    proximity_input_file = os.path.join(main_path, '../DrugRepurposing/data/medi/proximity.dat')
    doid_input_file = os.path.join(main_path, '../DrugRepurposing/data/DOID/DOID.csv')
    mesh_input_file = '/home/quim/Databases/MESH/MESH.csv'
    drug_field_to_num = {'total':2,'drugbank':3,'dgidb':4,'chembl':5,'drugcentral':6}
    disease_field_to_num = {'total':2,'disgenet':3,'omim':4,'uniprot':5,'go':6}
    drug_cutoff=10
    drug_cutoff_field='total'
    disease_cutoff=10
    disease_cutoff_field='disgenet'

    # Parsing of drug queries to GUILDify:
    # This file contains the queries of drug names to GUILDify and its results.
    #drug_associations_file = os.path.join(repur_path, 'drug_associations.txt')
    drug_associations_file = os.path.join(repur_path, 'drug_associations_and_sources.txt')
    #drugname_to_assoc, drugname_to_drugbank = parse_disease_associations(drug_associations_file, cutoff=drug_cutoff)
    drugname_to_assoc, drugname_to_drugbank, drugbank_to_drugname = parse_associations_file(drug_associations_file, cutoff=drug_cutoff, cutoff_field_num=drug_field_to_num[drug_cutoff_field])
    drugs = set(drugname_to_assoc.keys())
    print('GUILDify: {} drugs selected. Cut-off: {}'.format(len(drugs), drug_cutoff))

    # Parsing of disease queries to GUILDify:
    # This file contains the queries of disease CUIs from DisGeNET to GUILDify
    # and its results.
    #disease_associations_file = os.path.join(repur_path, 'disease_associations.txt')
    disease_associations_file = os.path.join(repur_path, 'disease_associations_and_sources.txt')
    #cui_to_assoc, cui_to_name = parse_disease_associations(disease_associations_file, cutoff=cutoff)
    cui_to_assoc, cui_to_name, name_to_cui = parse_associations_file(disease_associations_file, cutoff=disease_cutoff, cutoff_field_num=disease_field_to_num[disease_cutoff_field])
    disgenet_CUIs = set(cui_to_assoc.keys())
    print('GUILDify: {} CUIs selected. Cut-off: {}'.format(len(disgenet_CUIs), disease_cutoff))

    # Parsing of DisGeNET mappings:
    # This file contains the mappings of the CUIs (UMLS identifiers) in DisGeNET
    # to identifiers of other databases (MESH, DOID)
    CUI_to_MESH, MESH_to_CUI, CUI_to_DOID, DOID_to_CUI = parse_disease_mappings_disgenet(disgenet_mappings_file)

    # Parsing of Hetionet: 
    # Hetionet is a drug-disease indication DB with drugs in DrugBankID format 
    # and diseases in Disease Ontology ID (DOID) format.
    hetionet_output_file = os.path.join(repur_path, 'hetionet.txt')
    hetionet_DOID_to_drugbank = parse_hetionet(hetionet_input_file, hetionet_output_file)

    # Translation of Hetionet entities:
    # To find the entities of Hetionet in a GUILDify, we need to translate them.
    # The diseases have to be translated into CUI (using DOID_to_CUI).
    # The drugs are already in DrugBank, but we get the drug names from our 
    # mapping in GUILDify.

    hetionet_CUIs = set()
    for DOID in hetionet_DOID_to_drugbank:
        if DOID in DOID_to_CUI:
            for CUI in DOID_to_CUI[DOID]:
                hetionet_CUIs.add(CUI)
    print('Hetionet: {} DOID; {} CUIs.'.format(len(hetionet_DOID_to_drugbank), len(hetionet_CUIs)))

    hetionet_CUIs_selected = set()
    hetionet_guildify_diseases_file = os.path.join(repur_path, 'hetionet_guildify_diseases.txt')
    with open(hetionet_guildify_diseases_file, 'w') as output_fd:
        output_fd.write('#CUI\tDisease name\t# associations\n')
        #for cui, num_assoc in sorted(cui_to_assoc.items(), key=lambda x:x[1], reverse=True):
        for cui, data in sorted(cui_to_assoc.items(), key=lambda x:x[1][disease_field_to_num[disease_cutoff_field]], reverse=True):
            if cui in hetionet_CUIs:
                num_assoc = data[disease_field_to_num[disease_cutoff_field]]
                output_fd.write('{}\t{}\t{}\n'.format(cui, cui_to_name[cui], num_assoc))
                hetionet_CUIs_selected.add(cui)
    print('Hetionet: {} CUIs selected.'.format(len(hetionet_CUIs_selected)))

    mapping_file = os.path.join(repur_path, 'job_to_disease_mapping_{}.txt'.format(disease_cutoff))
    new_mapping_file = os.path.join(repur_path, 'job_to_disease_hetionet_{}.txt'.format(disease_cutoff))
    if fileExist(mapping_file) and not fileExist(new_mapping_file):
        id_to_job = parse_job_to_identifier_mapping(mapping_file)
        with open(new_mapping_file, 'w') as output_fd:
            for identifier in hetionet_CUIs_selected:
                if identifier in id_to_job:
                    job_id=id_to_job[identifier]
                    output_fd.write('{}\t{}\n'.format(job_id, identifier))

    hetionet_drugs = set()
    for DOID in hetionet_DOID_to_drugbank:
        for drugbank in hetionet_DOID_to_drugbank[DOID]['all']:
            hetionet_drugs.add(drugbank.upper())
    print('Hetionet: {} DrugBankIDs'.format(len(hetionet_drugs)))

    hetionet_drugbanks_selected = set()
    hetionet_drugs_selected = set()
    hetionet_guildify_drugs_file = os.path.join(repur_path, 'hetionet_guildify_drugs.txt')
    with open(hetionet_guildify_drugs_file, 'w') as output_fd:
        output_fd.write('#Drug name\tDrugbankID\t# associations\n')
        for drug_name, data in sorted(drugname_to_assoc.items(), key=lambda x:x[1][drug_field_to_num[drug_cutoff_field]], reverse=True):
            drugbank = drugname_to_drugbank[drug_name].upper()
            if drugbank in hetionet_drugs:
                num_assoc = data[drug_field_to_num[drug_cutoff_field]]
                output_fd.write('{}\t{}\t{}\n'.format(drug_name.lower(), drugbank, num_assoc))
                hetionet_drugbanks_selected.add(drugbank)
                hetionet_drugs_selected.add(drug_name)
    print('Hetionet: {} drugs selected.'.format(len(hetionet_drugbanks_selected)))

    mapping_file = os.path.join(repur_path, 'job_to_drug_mapping_{}.txt'.format(drug_cutoff))
    new_mapping_file = os.path.join(repur_path, 'job_to_drug_hetionet_{}.txt'.format(drug_cutoff))
    if fileExist(mapping_file) and not fileExist(new_mapping_file):
        id_to_job = parse_job_to_identifier_mapping(mapping_file)
        with open(new_mapping_file, 'w') as output_fd:
            for identifier in hetionet_drugs_selected:
                if identifier in id_to_job:
                    job_id=id_to_job[identifier]
                    output_fd.write('{}\t{}\n'.format(job_id, identifier))

    hetionet_indications_selected = []
    hetionet_guildify_indications_file = os.path.join(repur_path, 'hetionet_guildify_indications_diseases{}_drugs{}.txt'.format(disease_cutoff, drug_cutoff))
    with open(hetionet_guildify_indications_file, 'w') as output_fd:
        output_fd.write('#CUI\tDisease name\tDrug name\tDrugbankID\n')
        for DOID in hetionet_DOID_to_drugbank:
            if DOID in DOID_to_CUI:
                for CUI in DOID_to_CUI[DOID]:
                    if CUI in hetionet_CUIs_selected:
                        for drugbank in hetionet_DOID_to_drugbank[DOID]['all']:
                            if drugbank in hetionet_drugbanks_selected:
                                disease_name = cui_to_name[CUI]
                                drug_name = drugbank_to_drugname[drugbank]
                                indication = [CUI, drug_name]
                                hetionet_indications_selected.append(indication)
                                output_fd.write('{}\t{}\t{}\t{}\n'.format(CUI, disease_name.lower(), drug_name.lower(), drugbank.upper()))
    print('Hetionet: {} indications selected.'.format(len(hetionet_indications_selected)))

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

def parse_disease_to_name(disease_to_name_file):
    """
    Obtain a dictionary mapping the CUIs to the disease names.
    """
    cui_to_name = {}
    with open(disease_to_name_file, 'r') as assoc_fd:
        for line in assoc_fd:
            fields = line.strip().split('\t')
            cui = fields[0]
            disease_name = fields[1]
            cui_to_name[cui] = disease_name
    return cui_to_name

def parse_disease_associations(disease_associations_file, cutoff=20):
    """
    Get the CUIs of the diseases. cutoff is the minimum number of 
    disease-gene associations considered.
    """
    cui_to_name = {}
    cui_to_assoc = {}
    with open(disease_associations_file, 'r') as assoc_fd:
        first_line = assoc_fd.readline()
        for line in assoc_fd:
            fields = line.strip().split('\t')
            cui = fields[0].upper()
            disease_name = fields[1]
            num_assoc = int(fields[2])
            if num_assoc < cutoff:
                break
            cui_to_name[cui] = disease_name.lower()
            cui_to_assoc[cui] = num_assoc
    return cui_to_assoc, cui_to_name

def parse_associations_file(associations_file, cutoff=1, cutoff_field_num=2):
    """
    Parse an associations file.
    """
    id_to_assoc = {}
    id_to_id2 = {}
    id2_to_id = {}
    with open(associations_file, 'r') as assoc_fd:
        first_line = assoc_fd.readline()
        for line in assoc_fd:
            fields = line.strip().split('\t')
            id1 = fields[0]
            id2 = fields[1]
            data = [id1,id2]+[int(x) for x in fields[2:]]
            if data[cutoff_field_num] >= cutoff:
                id_to_assoc[id1] = data
                id_to_id2[id1] = id2
                id2_to_id[id2] = id1
    return id_to_assoc, id_to_id2, id2_to_id

def parse_job_to_identifier_mapping(mapping_file):
    """
    Get the mapping between the GUILDify job IDs and the identifiers of the diseases/drugs.
    """
    id_to_job = {}
    with open(mapping_file, 'r') as assoc_fd:
        for line in assoc_fd:
            fields = line.strip().split('\t')
            job_id = fields[0]
            identifier = fields[1]
            id_to_job[identifier] = job_id
    return id_to_job

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

def parse_proximity_data(input_file, output_file):
    """
    Parses Proximity data "proximity.dat".
    Obtains a file containing the DOID, its corresponding disease name and Palliative/Treatment.
    """
    disease_to_drugbank = {}
    if not fileExist(output_file):
        with open(input_file, 'r') as input_fd, open(output_file, 'w') as output_fd:
            #group group.name disease flag n.target n.disease n.overlap j.overlap d.target d.disease symptomatic ra re z.pathway z.side d z pval pval.adj
            first_line = input_fd.readline()
            output_fd.write('DiseaseName\tDrugBankID\tZScore\tAdjPvalue\n')
            for line in input_fd:
                fields = line.strip().split(' ')
                drugbank = fields[0].upper()
                drug_name = fields[1].lower()
                disease_name = fields[2].lower()
                flag = fields[3].lower()
                zscore = float(fields[16])
                adj_pval = float(fields[18])
                if flag == 'true':
                    output_fd.write('{}\t{}\t{:.3f}\t{:.1E}\n'.format(disease_name, drugbank, zscore, adj_pval))
                    disease_to_drugbank.setdefault(disease_name, set()).add(drugbank)
    else:
        with open(output_file, 'r') as input_fd:
            first_line = input_fd.readline()
            for line in input_fd:
                disease_name, drugbank, zscore, adj_pval = line.strip().split('\t')
                disease_to_drugbank.setdefault(disease_name, set()).add(drugbank)
    return disease_to_drugbank

def parse_DOID(input_file, output_file):
    """
    Parses DOID and obtains its mapping to UMLS CUI
    """
    DOID_to_UMLS = {}
    DOID_to_MESH = {}
    if not fileExist(output_file):
        df_doid = pd.read_csv(input_file, sep=',')
        for index, row in df_doid.iterrows():

            # Obtain DOID id from 'http://purl.obolibrary.org/obo/DOID_8986'
            DOID = str(row['Class ID'])
            DOID_split = DOID.split('/')
            DOID = DOID_split[len(DOID_split)-1]
            DOID = DOID.split('DOID_')[1]
            DOID = 'DOID:'+DOID
            DOID = DOID.upper()

            # Obtain cross-references
            #EFO:0000614|ICD10CM:G47.419|MESH:D009290|ICD10CM:G47.41|OMIM:612851|SNOMEDCT_US_2016_03_01:267702006|SNOMEDCT_US_2016_03_01:60380001|SNOMEDCT_US_2016_03_01:155059003|OMIM:614223|OMIM:605841|ICD9CM:347.0|OMIM:161400|UMLS_CUI:C0027404|OMIM:609039|OMIM:614250|OMIM:612417|NCI:C84489|ORDO:2073
            crossreferences = str(row['database_cross_reference'])
            if crossreferences != 'nan':
                cross_split = crossreferences.split('|')
                for cross in cross_split:
                    cr_split = cross.split(':')
                    if len(cr_split) == 2:
                        type_cross, content = cr_split
                        identifier = content.upper()
                        if type_cross == 'UMLS_CUI':
                            DOID_to_UMLS.setdefault(DOID, set()).add(identifier)
                        if type_cross == 'MESH':
                            DOID_to_MESH.setdefault(DOID, set()).add(identifier)
    return DOID_to_UMLS, DOID_to_MESH

def parse_MESH(input_file, output_file):
    """
    Parses MESH and obtains its mapping to UMLS CUI and disease names
    """
    MESH_to_name = {}
    name_to_MESH = {}
    MESH_to_UMLS = {}
    if not fileExist(output_file):
        with open(output_file, 'w') as output_fd:
            output_fd.write('MESH\tMESH name\tCUIs\n')
            df_mesh = pd.read_csv(input_file, sep=',')
            for index, row in df_mesh.iterrows():

                # Obtain the identifier from http://purl.bioontology.org/ontology/MESH/C585345
                MESH = str(row['Class ID'])
                MESH_split = MESH.split('/')
                MESH = MESH_split[len(MESH_split)-1]
                MESH = MESH.upper()

                # Obtain the name
                MESH_name = str(row['Preferred Label']).lower()
                MESH_name = '.'.join('.'.join('.'.join('.'.join(MESH_name.split(', ')).split(' ')).split(',')).split('-'))
                MESH_to_name[MESH] = MESH_name
                name_to_MESH[MESH_name] = MESH

                # Obtain CUI
                CUIs = str(row['CUI'])
                if CUIs != 'nan':
                    CUIs = CUIs.split('|')
                    output_fd.write('{}\t{}\t{}\n'.format(MESH, MESH_name, ';'.join(CUIs)))
                    for CUI in CUIs:
                        MESH_to_UMLS.setdefault(MESH, set()).add(CUI)
    else:
        with open(output_file, 'r') as input_fd:
            first_line = input_fd.readline()
            for line in input_fd:
                MESH, MESH_name, CUIs = line.strip().split('\t')
                CUIs = CUIs.split(';')
                MESH_to_name[MESH] = MESH_name
                name_to_MESH[MESH_name] = MESH
                for CUI in CUIs:
                    MESH_to_UMLS.setdefault(MESH, set()).add(CUI)
    return MESH_to_name, name_to_MESH, MESH_to_UMLS

def parse_disease_mappings_disgenet(input_file):
    """
    Parse the disease_mappings.tsv file from DisGeNET.
    """
    CUI_to_MESH = {}
    MESH_to_CUI = {}
    CUI_to_DOID = {}
    DOID_to_CUI = {}
    with open(input_file, 'r') as input_fd:
        first_line = input_fd.readline()
        #diseaseId  name    vocabulary  code    vocabularyName
        for line in input_fd:
            CUI, name, vocabulary, vocab_code, vocab_name = line.strip().split('\t')
            if vocabulary.upper() == 'MSH':
                CUI_to_MESH.setdefault(CUI.upper(), set()).add(vocab_code.upper())
                MESH_to_CUI.setdefault(vocab_code.upper(), set()).add(CUI.upper())
            elif vocabulary.upper() == 'DO':
                vocab_code = 'DOID:'+str(int(vocab_code))
                CUI_to_DOID.setdefault(CUI.upper(), set()).add(vocab_code.upper())
                DOID_to_CUI.setdefault(vocab_code.upper(), set()).add(CUI.upper())
    return CUI_to_MESH, MESH_to_CUI, CUI_to_DOID, DOID_to_CUI


if  __name__ == "__main__":
    main()