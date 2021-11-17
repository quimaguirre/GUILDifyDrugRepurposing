import sys, os, re
import argparse
import time

def main():
    """
    Get job IDs and create links with the names of the corresponding keywords. Save the results into another folder. 
    Example:
    python /home/quim/PHD/Projects/DIANA/GUILDIFY_test/scripts/copy_guildify_results_from_aleph.py
    """

    # Directory with stored results from GUILDify (in Aleph)
    guildify_sessions = '/var/www/html/sbi/guildify2/data/sessions'
    job_list = '/var/www/html/sbi/guildify2/data/drug_repurposing/job_to_drug_mapping_5.txt'
    output_path = '/var/www/html/sbi/guildify2/data/drug_repurposing_additional_results'

    # Read keywords and job ids from job list and fetch the files of interest
    create_directory(output_path)
    if fileExist(job_list):
        with open(job_list, 'r') as input_fd:
            for line in input_fd:
                fields = line.strip().split('\t')
                job_id = fields[0]
                keyword = fields[1]
                if job_id != keyword:
                    job_dir = os.path.join(guildify_sessions, job_id)
                    link_dir = os.path.join(guildify_sessions, keyword)
                    # Create link with keyword name
                    if not os.path.exists(link_dir):
                        command = 'ln -s {} "{}"'.format(job_dir, link_dir)
                        print(command)
                        os.system(command)
                    # Copy results to a new folder
                    new_dir = os.path.join(output_path, keyword)
                    if not os.path.exists(new_dir):
                        os.system('cp -r {} "{}"'.format(job_dir, new_dir))
    else:
        print('Job list {} does not exist!'.format(job_list))
        sys.exit(10)

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

if  __name__ == "__main__":
    main()