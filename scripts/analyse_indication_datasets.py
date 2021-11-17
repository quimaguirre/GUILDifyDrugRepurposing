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
    cut_off = 20

    # Parse DisGeNET mappings file
    CUI_to_MESH, MESH_to_CUI, CUI_to_DOID, DOID_to_CUI = parse_disease_mappings_disgenet(disgenet_mappings_file)
    disgenet_all_MESH = set()
    for CUI in CUI_to_MESH:
        for MESH in CUI_to_MESH[CUI]:
            disgenet_all_MESH.add(MESH)
    disgenet_all_DOID = set()
    for CUI in CUI_to_DOID:
        for DOID in CUI_to_DOID[CUI]:
            disgenet_all_DOID.add(DOID)
    print('DisGeNET: {} MESH. {} DOID.'.format(len(disgenet_all_MESH), len(disgenet_all_DOID)))

    cui_to_assoc, cui_to_name = parse_disease_associations(repur_path, cut_off=cut_off)
    disgenet_UMLS = set(cui_to_assoc.keys())
    print('DisGeNET: {} CUIs selected. Cut-off: {}'.format(len(disgenet_UMLS), cut_off))

    # Parse Hetionet
    hetionet_output_file = os.path.join(repur_path, 'hetionet.txt')
    hetionet_DOID_to_drugbank = parse_hetionet(hetionet_input_file, hetionet_output_file)

    # Parse Proximity data
    proximity_output_file = os.path.join(repur_path, 'proximity.txt')
    proximity_disease_to_drugbank = parse_proximity_data(proximity_input_file, proximity_output_file)

    # Parse DOID
    doid_output_file = os.path.join(repur_path, 'DOID.txt')
    DOID_to_UMLS, DOID_to_MESH = parse_DOID(doid_input_file, doid_output_file)

    # Parse MESH
    mesh_output_file = os.path.join(repur_path, 'MESH.txt')
    MESH_to_name, name_to_MESH, MESH_to_UMLS = parse_MESH(mesh_input_file, mesh_output_file)

    # Transform Proximity to CUI
    print('Proximity: {} MESH names.'.format(len(proximity_disease_to_drugbank)))
    proximity_MESH = set()
    for disease in proximity_disease_to_drugbank:
        if disease in name_to_MESH:
            MESH = name_to_MESH[disease]
            proximity_MESH.add(MESH)
        else:
            print('Wrong: {}'.format(disease))
    print('Proximity: {} MESH IDs.'.format(len(proximity_MESH)))
    proximity_UMLS_per_MESH = set()
    for MESH in proximity_MESH:
        if MESH in MESH_to_UMLS:
            proximity_UMLS_per_MESH.add(';'.join(MESH_to_UMLS[MESH]))
    print('Proximity: {} CUIs.'.format(len(proximity_UMLS_per_MESH)))
    proximity_UMLS = set()
    for CUI in proximity_UMLS_per_MESH:
        proximity_UMLS.add(CUI)
    print('Proximity: {} total CUIs.'.format(len(proximity_UMLS)))

    # Transform Hetionet to CUI
    print('Hetionet: {} DOID names.'.format(len(hetionet_DOID_to_drugbank)))
    hetionet_UMLS_per_DOID = set()
    for DOID in hetionet_DOID_to_drugbank:
        if DOID in DOID_to_UMLS:
            hetionet_UMLS_per_DOID.add(';'.join(DOID_to_UMLS[DOID]))
        else:
            print('Wrong DOID: {}.'.format(DOID))
    print('Hetionet: {} CUIs.'.format(len(hetionet_UMLS_per_DOID)))
    hetionet_UMLS = set()
    for CUI in hetionet_UMLS_per_DOID:
        hetionet_UMLS.add(CUI)
    print('Hetionet: {} total CUIs.'.format(len(hetionet_UMLS)))

    # Overlap between DisGeNET and Proximity
    overlap_prox = disgenet_UMLS & proximity_UMLS
    print('Overlap DisGeNET & Proximity: {} CUIs.'.format(len(overlap_prox)))

    # Overlap between DisGeNET and Hetionet
    overlap_het = disgenet_UMLS & hetionet_UMLS
    print('Overlap DisGeNET & Hetionet: {} CUIs.'.format(len(overlap_het)))

    # Overlap between DisGeNET mappings and Proximity
    overlap_prox = disgenet_all_MESH & set(proximity_MESH)
    print('Overlap DisGeNET mappings & Proximity: {} MESH.'.format(len(overlap_prox)))

    # Overlap between DisGeNET mappings and Hetionet
    overlap_het = disgenet_all_DOID & set(hetionet_DOID_to_drugbank.keys())
    print('Overlap DisGeNET mappings & Hetionet: {} DOID.'.format(len(overlap_het)))

    disgenet_MESH = set()
    for CUI in cui_to_assoc:
        if CUI in CUI_to_MESH:
            for MESH in CUI_to_MESH[CUI]:
                disgenet_MESH.add(MESH)

    disgenet_DOID = set()
    for CUI in cui_to_assoc:
        if CUI in CUI_to_DOID:
            for DOID in CUI_to_DOID[CUI]:
                disgenet_DOID.add(DOID)

    overlap_prox = disgenet_MESH & set(proximity_MESH)
    print('Overlap DisGeNET & Proximity: {} MESH.'.format(len(overlap_prox)))
    overlap_het = disgenet_DOID & set(hetionet_DOID_to_drugbank.keys())
    print('Overlap DisGeNET & Hetionet: {} DOID.'.format(len(overlap_het)))

    proximity_UMLS_per_MESH_from_mappings = set()
    proximity_UMLS_from_mappings = set()
    for MESH in proximity_MESH:
        if MESH in MESH_to_CUI:
            proximity_UMLS_per_MESH_from_mappings.add(';'.join(MESH_to_CUI[MESH]))
            for CUI in MESH_to_CUI[MESH]:
                proximity_UMLS_from_mappings.add(CUI)
    print('Proximity: {} CUIs from mappings.'.format(len(proximity_UMLS_per_MESH_from_mappings)))
    print('Proximity: {} total CUIs from mappings.'.format(len(proximity_UMLS_from_mappings)))

    hetionet_UMLS_per_DOID_from_mappings = set()
    hetionet_UMLS_from_mappings = set()
    for DOID in hetionet_DOID_to_drugbank:
        if DOID in DOID_to_CUI:
            hetionet_UMLS_per_DOID_from_mappings.add(';'.join(DOID_to_CUI[DOID]))
            for CUI in DOID_to_CUI[DOID]:
                hetionet_UMLS_from_mappings.add(CUI)
    print('Hetionet: {} CUIs from mappings.'.format(len(hetionet_UMLS_per_DOID_from_mappings)))
    print('Hetionet: {} total CUIs from mappings.'.format(len(hetionet_UMLS_from_mappings)))

    overlap_prox = disgenet_UMLS & proximity_UMLS_from_mappings
    print('Overlap DisGeNET & Proximity from mappings: {} CUIs.'.format(len(overlap_prox)))
    overlap_het = disgenet_UMLS & hetionet_UMLS_from_mappings
    print('Overlap DisGeNET & Hetionet from mappings: {} CUIs.'.format(len(overlap_het)))

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

def parse_disease_to_name(drug_repurposing_dir):
    """
    Obtain a dictionary mapping the CUIs to the disease names.
    """
    disease_to_name_file = os.path.join(drug_repurposing_dir, 'disease_names.txt')
    cui_to_name = {}
    with open(disease_to_name_file, 'r') as assoc_fd:
        for line in assoc_fd:
            fields = line.strip().split('\t')
            cui = fields[0]
            disease_name = fields[1]
            cui_to_name[cui] = disease_name
    return cui_to_name

def parse_disease_associations(drug_repurposing_dir, cut_off=20):
    """
    Get the CUIs of the diseases. cut_off is the minimum number of 
    disease-gene associations considered.
    """
    disease_associations_file = os.path.join(drug_repurposing_dir, 'disease_associations.txt')
    cui_to_name = {}
    cui_to_assoc = {}
    with open(disease_associations_file, 'r') as assoc_fd:
        first_line = assoc_fd.readline()
        for line in assoc_fd:
            fields = line.strip().split('\t')
            cui = fields[0].upper()
            disease_name = fields[1]
            num_assoc = int(fields[2])
            if num_assoc < cut_off:
                break
            cui_to_name[cui] = disease_name.lower()
            cui_to_assoc[cui] = num_assoc
    return cui_to_assoc, cui_to_name

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