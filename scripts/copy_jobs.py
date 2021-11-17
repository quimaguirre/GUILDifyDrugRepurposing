import sys, os, re

def main():

    repur_path = '/var/www/html/sbi/guildify2/data/drug_repurposing'
    sessions_dir = '/var/www/html/sbi/guildify2/data/sessions'
    #jobs_dir = os.path.join(repur_path, 'jobs_diseases')
    jobs_dir = os.path.join(repur_path, 'jobs_drugs')
    drug_field_to_num = {'total':2,'drugbank':3,'dgidb':4,'chembl':5,'drugcentral':6}
    disease_field_to_num = {'total':2,'disgenet':3,'omim':4,'uniprot':5,'go':6}
    drug_cutoff=10
    drug_cutoff_field='total'
    disease_cutoff=10
    disease_cutoff_field='disgenet'

    #disease_associations_file = os.path.join(repur_path, 'disease_associations_and_sources.txt')
    #id_to_assoc, id_to_id2, id2_to_id = parse_associations_file(disease_associations_file, cutoff=disease_cutoff, cutoff_field_num=disease_field_to_num[disease_cutoff_field])
    drug_associations_file = os.path.join(repur_path, 'drug_associations_and_sources.txt')
    id_to_assoc, id_to_id2, id2_to_id = parse_associations_file(drug_associations_file, cutoff=drug_cutoff, cutoff_field_num=drug_field_to_num[drug_cutoff_field])

    #mapping_file = os.path.join(repur_path, 'job_to_disease_mapping_10.txt')
    mapping_file = os.path.join(repur_path, 'job_to_drug_mapping_10.txt')
    id_to_job = parse_job_to_identifier_mapping(mapping_file)

    for identifier in id_to_assoc:
        if identifier not in id_to_job:
            print('Identifier missing: {}'.format(identifier))
            continue
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

def parse_disease_associations(disease_associations_file, cut_off=20):
    """
    Get the CUIs of the diseases. cut_off is the minimum number of 
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
            if num_assoc < cut_off:
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

if  __name__ == "__main__":
    main()

