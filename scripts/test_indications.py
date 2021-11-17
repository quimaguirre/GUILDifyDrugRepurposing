import sys, os, re
import mysql.connector
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pylab

def main():

    repur_path = '/var/www/html/sbi/guildify2/data/drug_repurposing'
    sessions_dir = '/var/www/html/sbi/guildify2/data/sessions'
    #jobs_dir = os.path.join(repur_path, 'jobs_diseases')
    jobs_dir = os.path.join(repur_path, 'jobs_drugs')
    drug_field_to_num = {'total':2,'drugbank':3,'dgidb':4,'chembl':5,'drugcentral':6}
    disease_field_to_num = {'total':2,'disgenet':3,'omim':4,'uniprot':5,'go':6}
    drug_cutoff=10
    drug_cutoff_field='total'
    disease_cutoff=20
    disease_cutoff_field='disgenet'

    disease_associations_file = os.path.join(repur_path, 'disease_associations_and_sources.txt')
    cui_to_assoc, cui_to_name, name_to_cui = parse_associations_file(disease_associations_file, cutoff=disease_cutoff, cutoff_field_num=disease_field_to_num[disease_cutoff_field])
    drug_associations_file = os.path.join(repur_path, 'hetionet_guildify_drugs.txt')
    drugname_to_assoc, drugname_to_drugbank, drugbank_to_drugname = parse_associations_file(drug_associations_file, cutoff=drug_cutoff, cutoff_field_num=drug_field_to_num[drug_cutoff_field])

    mapping_file = os.path.join(repur_path, 'job_to_disease_mapping.txt')
    cui_to_job = parse_job_to_identifier_mapping(mapping_file)
    mapping_file = os.path.join(repur_path, 'job_to_drug_all_targets.txt')
    drugname_to_job = parse_job_to_identifier_mapping(mapping_file)

    CUIs = cui_to_job.keys()

    for identifier in id_to_assoc:
        identifier=identifier.split(' ')[0]
        if identifier not in id_to_job:
            print('Identifier missing: {}'.format(identifier))
            sys.exit(10)
        job = id_to_job[identifier]
        input_dir = os.path.join(sessions_dir, job)
        output_dir = os.path.join(jobs_dir, job)
        if not os.path.exists(input_dir):
            print('{} and {} does not exists'.format(identifier, job))
        if not os.path.exists(output_dir):
            command = 'cp -rf {} {}'.format(input_dir, output_dir)
            os.system(command)

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

def parse_drug_repurposing_table(genes_table, functions_table, identifiers):
    """
    Create a table containing all the diseases and its genes and functions
    """
    if not os.path.exists(genes_table) or not os.path.exists(functions_table): 
        with open(genes_table, 'w') as out1, open(functions_table, 'w') as out2:
            out1.write('JobID\tCUI\tDisease name\tTop1\tTop2\tTopEnrich\n')
            out2.write('JobID\tCUI\tDisease name\tTop1-GObp-fdr_bh\tTop1-GObp-bonferroni\tTop1-GOmf-fdr_bh\tTop1-GOmf-bonferroni\tTop1-Reactome-fdr_bh\tTop1-Reactome-bonferroni\tTop2-GObp-fdr_bh\tTop2-GObp-bonferroni\tTop2-GOmf-fdr_bh\tTop2-GOmf-bonferroni\tTop2-Reactome-fdr_bh\tTop2-Reactome-bonferroni\tTopEnrich-GObp-fdr_bh\tTopEnrich-GObp-bonferroni\tTopEnrich-GOmf-fdr_bh\tTopEnrich-GOmf-bonferroni\tTopEnrich-Reactome-fdr_bh\tTopEnrich-Reactome-bonferroni\n')
            for identifier in identifiers:
                if identifier in cui_to_name and cui in cui_to_job:
                    disease_name = cui_to_name[cui]
                    job_id = cui_to_job[cui]
                    check_retrieve = GUILDifier(job_id, species = None, settings = self.settings, retrieve=True)
                    if not check_retrieve.retrieved:
                        print('Error in JOB {}, CUI {}: job not finished '.format(job_id, cui))
                        continue
                    disease_result = GUILDifier(job_id, species = None, settings = self.settings, tissue = None)
                    ready_flag = disease_result.create_output_file()
                    if not ready_flag:
                        print('Error in JOB {}, CUI {}: job not finished '.format(job_id, cui))
                        continue
                    parser = GUILDScoreParser(disease_result.output_file, disease_result.edge_scores_file)
                    # Calculate the biological validation
                    position_values, enrichment_values, seed_values, validated_nodes, non_validated_nodes, rank_range, enrichment_cutoff = disease_result.calculate_biological_validation_using_Carlota(disease_result.species)
                    validation_div = disease_result.plot_biological_validation_with_plotly(position_values, enrichment_values, seed_values, validated_nodes, non_validated_nodes, rank_range)
                    # Calculate the top genes and functional enrichment
                    selected_methods = disease_result.selected_method_to_file.keys()
                    criteria_list=['fdr_bh', 'bonferroni']
                    top_enrichment = 'enrich_'+str(enrichment_cutoff)
                    top_percentages = [1, 2, top_enrichment]
                    disease_result.calculate_enrichment_of_top_nodes(parser, disease_result.species, top_percentages)
                    #seed_enrichment_file = os.path.join(disease_result.session_dir, 'enrichment.{}.{}.seeds.txt'.format(selected_method, criteria))
                    #seed_go_rows_bp = disease_result.parse_enrichment_file_Carlota(seed_enrichment_file)
                    top_to_genes = {}
                    top_to_functions = {}
                    for top in top_percentages:
                        user_entity_id_to_values = parser.get_top_scoring_ueids_and_values(top, disease_result.diamond)
                        top_genes = ';'.join(sorted(user_entity_id_to_values.keys()))
                        print(job_id, top, top_genes)
                        top_to_genes[top] = top_genes
                        for selected_method in selected_methods:
                            for criteria in criteria_list:
                                enrichment_file = os.path.join(disease_result.session_dir, 'enrichment.{}.{}.{}.txt'.format(selected_method, criteria, str(top)))
                                top_go_rows = disease_result.parse_enrichment_file_Carlota(enrichment_file)
                                go_ids = ';'.join(sorted([row[0] for row in top_go_rows]))
                                top_to_functions.setdefault(top, {})
                                top_to_functions[top].setdefault(selected_method, {})
                                top_to_functions[top][selected_method][criteria] = go_ids
                    out1.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(job_id, cui, disease_name, top_to_genes[1], top_to_genes[2], top_to_genes[top_enrichment]))
                    out2.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(job_id, cui, disease_name, top_to_functions[1]['GObp']['fdr_bh'], top_to_functions[1]['GObp']['bonferroni'], top_to_functions[1]['GOmf']['fdr_bh'], top_to_functions[1]['GOmf']['bonferroni'], top_to_functions[1]['Reactome']['fdr_bh'], top_to_functions[1]['Reactome']['bonferroni'],
 


if  __name__ == "__main__":
    main()