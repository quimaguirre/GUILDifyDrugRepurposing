import sys, os, re
import argparse
import time

def main():
    """
    Retrieve jobs from GUILDify using the R package guildifyR
    Example:
    python /home/quim/PHD/Projects/DIANA/GUILDIFY_test/scripts/retrieve_jobs_from_guildify.py -j /home/quim/PHD/Projects/DIANA/GUILDIFY_test/outputs/Rscripts_drugs_submit/job_ids.txt -r /home/quim/PHD/Projects/DIANA/GUILDIFY_test/outputs/Rscripts_drugs_retrieve
    """
    options = parse_user_arguments()
    send_jobs_to_guildify(options)
    return

def parse_user_arguments(*args, **kwds):
    """
    Parses the arguments of the program
    """

    parser = argparse.ArgumentParser(
        description = "Send jobs to GUILDify using the guildifyR package",
        epilog      = "@oliva's lab 2020")
    parser.add_argument('-j','--job_list',dest='job_list',action = 'store',
                        help = """ File containing the mapping of keywords to job IDs: keyword<tab>job_id. """)
    parser.add_argument('-r','--r_path',dest='r_path',action = 'store',
                        help = """ Folder to store the RScripts that will be created. """)

    options=parser.parse_args()

    return options

#################
#################
# MAIN FUNCTION #
#################
#################

def send_jobs_to_guildify(options):
    """
    Send jobs to GUILDify using guildifyR package.
    """
    # Directory with stored results from GUILDify
    guildify_sessions = '/var/www/html/guildify2/data/sessions'

    # Read keywords and job ids from input list
    keyword_to_job_id = {}
    job_id_to_keyword = {}
    if fileExist(options.job_list):
        with open(options.job_list, 'r') as input_fd:
            for line in input_fd:
                fields = line.strip().split('\t')
                job_id = fields[0]
                keyword = fields[1]
                keyword_to_job_id[keyword] = job_id
                job_id_to_keyword[job_id] = keyword
    else:
        print('Job list {} does not exist!'.format(options.job_list))
        sys.exit(10)

    # Create R directory
    if not options.r_path:
        print('Introduce an R directory to store the scripts.')
        sys.exit(10)
    create_directory(options.r_path)
    
    # Create R script to fetch results
    rscript = os.path.join(options.r_path, 'retrieve_jobs.R')
    log_file = os.path.join(options.r_path, 'retrieve_jobs.log')
    with open(rscript, 'w') as rscript_fd:
        rscript_fd.write('library(guildifyR)\n')
        for job_id in job_id_to_keyword:
            rscript_fd.write('retrieve.job("{}", n.top = NULL, fetch.files = F, output.dir = NULL)\n'.format(job_id))
    os.system('Rscript {} &> {}'.format(rscript, log_file))


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

def create_rscript(rscript, keyword, taxid, tissue, network, databases, scoring, repetitions=3, iterations=2):
    """
    Create an R script that sends a GUILDify job
    """
    databases_mod = ['"{}"'.format(database) for database in databases]
    with open(rscript, 'w') as rscript_fd:

        rscript_fd.write('library(guildifyR)\n')
        rscript_fd.write('library(tidyr)\n')
        rscript_fd.write('drug_databases = c({})\n'.format(', '.join(databases_mod)))
        rscript_fd.write('result.table = query("{}", species="{}", tissue="{}", network.source = "{}")\n'.format(keyword, taxid, tissue, network))
        rscript_fd.write('sep.result.table = separate_rows(result.table, source, sep=", ")\n')
        rscript_fd.write('filt.result.table = sep.result.table[ which( sep.result.table$source %in% drug_databases & sep.result.table$in.network==1 ), ]\n')

        if scoring in ['netscore', 'netcombo']:
            rscript_fd.write('submit.job(filt.result.table, species="{}", tissue="{}", network.source="{}", scoring.options = list({}=T, repetitionSelector={}, iterationSelector={}))\n'.format(taxid, tissue, network, scoring, repetitions, iterations))
        elif scoring == netzcore:
            rscript_fd.write('submit.job(filt.result.table, species="{}", tissue="{}", network.source="{}", scoring.options = list(netzcore=T, iterationSelector={}))\n'.format(taxid, tissue, network, iterations))
        else:
            rscript_fd.write('submit.job(filt.result.table, species="{}", tissue="{}", network.source="{}", scoring.options = list({}=T))\n'.format(taxid, tissue, network, scoring))
        return

def extract_job_id_from_log_file(log_file):
    """
    Obtain the job id from the log file of the guildifyR submission
    """
    job_id = None
    with open(log_file, 'r') as log_fd:
        for line in log_fd:
            if 'Job id: ' in line:
                # [1] "Job id: 4682695b-d09c-4769-822e-c4549d2dbd2e"
                # ['[1] "', '4682695b-d09c-4769-822e-c4549d2dbd2e"']
                fields = line.strip().split('Job id: ') 
                job_id = fields[1].rstrip('"')
    if job_id:
        return job_id
    else:
        print('No job ID in log file: {}'.format(log_file))
        sys.exit(10)

def get_keyword_filename(keyword):
    """
    If a keyword is composed by multiple words, join them.
    """
    keyword_parts = keyword.split(' ')
    if len(keyword_parts) > 1:
        keyword_filename = '_'.join(keyword_parts)
        return keyword_filename
    else:
        return keyword

if  __name__ == "__main__":
    main()