import sys, os, re
import argparse
import time

def main():
    """
    Sends jobs to GUILDify using the R package guildifyR
    Example:
    python /home/quim/PHD/Projects/DIANA/GUILDIFY_test/scripts/send_jobs_to_guildify.py -i /home/quim/PHD/Projects/DIANA/GUILDIFY_test/data/Cheng_NatCom19_Drugs_More1Targets.txt -r /home/quim/PHD/Projects/DIANA/GUILDIFY_test/outputs/Rscripts_drugs_submit -t 9606 -n biana -d drugbank,dgidb,drugcentral,chembl -s netscore -re 3 -it 2 -ch
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
    parser.add_argument('-i','--input_list',dest='input_list',action = 'store',
                        help = """ File containing the list of keywords to be sent to GUILDify. """)
    parser.add_argument('-r','--r_path',dest='r_path',action = 'store',
                        help = """ Folder to store the RScripts that will be created. """)
    parser.add_argument('-t','--taxid',dest='taxid',action = 'store', default='9606',
                        help = """ Taxonomy ID of the species of the network: 9606 (human) / 10090 (mouse) / 3702 (plant) / 6239 (worm) /
                                   7227 (fly) / 10116 (rat) / 4932 (yeast) """)
    parser.add_argument('-n','--network',dest='network',action = 'store', default='biana',
                        help = """ Name of the network: biana / consensuspathdb / hippie / i2d / inbio_map / string """)
    parser.add_argument('-d','--databases',dest='databases',action = 'store', default='disgenet,omim,uniprot,go,drugbank,dgidb,drugcentral,chembl',
                        help = """ Name of the databases (separated by , if more than one):
                                   For drugs: drugbank,dgidb,drugcentral,chembl
                                   For phenotypes: disgenet,omim,uniprot,go""")
    parser.add_argument('-s','--scoring',dest='scoring',action = 'store',
                        help = """ Scoring function name: netscore / netzcore / netshort / netcombo / diamond. """)
    parser.add_argument('-re','--repetitions',dest='repetitions',action = 'store', default='3',
                        help = """ Number of repetitions (for netscore and netcombo) """)
    parser.add_argument('-it','--iterations',dest='iterations',action = 'store', default='2',
                        help = """ Number of iterations (for netscore, netzcore and netcombo) """)
    parser.add_argument('-ch','--check_if_exists',dest='check_if_exists',action = 'store_true',
                        help = """ Check if the keyword exists as a job ID (only available if executed locally) """)

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
    # Count number of jobs submitted
    # If exceeds X number, sleep for X minutes
    counter = 0
    job_limit = 25
    sleep_time = 600 # 10 minutes

    # Tissue will be all for now
    tissue = 'all'

    # Directory with stored results from GUILDify
    guildify_sessions = '/var/www/html/guildify2/data/sessions'

    # Read keywords from input list
    keywords = []
    if fileExist(options.input_list):
        with open(options.input_list, 'r') as input_fd:
            for line in input_fd:
                keyword = line.strip()
                keywords.append(keyword)
    else:
        print('Input list {} does not exist!'.format(options.input_list))
        sys.exit(10)

    # Check databases
    databases_available = ['disgenet','omim','uniprot','go','drugbank','dgidb','drugcentral','chembl']
    databases = options.databases.split(',')
    if len(databases) > 0:
        for database in databases:
            if database not in databases_available:
                print('Incorrect database: {}. Introduce databases among the following list: {}'.format(database, ', '.join(databases_available)))
                sys.exit(10)
    else:
        print('Incorrect databases: {}. Introduce databases among the following list: {}'.format(options.databases, ', '.join(databases_available)))
        sys.exit(10)

    # Check network
    networks_available = {'biana':'BIANA', 'inbio_map':'InBio_Map', 'i2d':'I2D', 'consensuspathdb':'ConsensusPathDB', 'hippie':'HIPPIE', 'string':'STRING'}
    if options.network.lower() in networks_available:
        network = networks_available[options.network.lower()]
    else:
        print('Network {} is not available. Introduce any of the following networks: {}'.format(options.network.lower(), ', '.join(networks_available.keys())))
        sys.exit(10)

    # Check scoring
    scoring_functions = ['netcombo', 'netscore', 'netzcore', 'netshort', 'diamond']
    if options.scoring.lower() in scoring_functions:
        scoring = options.scoring.lower()
        if scoring == 'netcombo':
            sleep_time = 1200
    else:
        print('Scoring function {} does not exist. Introduce any of the following scoring functions: {}'.format(options.scoring.lower(), ', '.join(scoring_functions)))
        sys.exit(10)

    # Create R directory
    if not options.r_path:
        print('Introduce an R directory to store the scripts.')
        sys.exit(10)
    create_directory(options.r_path)

    # Create RScripts and send them (skip if the log file exists!)
    # Save the job IDs in a file
    job_ids_file = os.path.join(options.r_path, 'job_ids.txt')
    with open (job_ids_file, 'w') as job_fd:
        for keyword in keywords:
            keyword_filname = get_keyword_filename(keyword)
            rscript = os.path.join(options.r_path, '{}.R'.format(keyword_filname))
            log_file = os.path.join(options.r_path, '{}.log'.format(keyword_filname))
            if fileExist(log_file):
                job_id = extract_job_id_from_log_file(log_file)
                if job_id:
                    job_fd.write('{}\t{}\n'.format(job_id, keyword))
                    print('JOB {} EXISTS with job ID {}'.format(keyword, job_id))
                continue
            else:
                if options.check_if_exists:
                    # Check if a job with the keyword name already exists and skip it if so
                    job_dir = os.path.join(guildify_sessions, keyword)
                    if os.path.exists(job_dir):
                        print('JOB {} EXISTS!'.format(keyword))
                        job_fd.write('{}\t{}\n'.format(keyword, keyword))
                        continue

                create_rscript(rscript, keyword, options.taxid, tissue, network, databases, scoring, repetitions=options.repetitions, iterations=options.iterations)
                print('Submitting Rscript for keyword: {}'.format(keyword))
                os.system('Rscript {} &> {}'.format(rscript, log_file))
                time.sleep(20) # Wait for the job to be submitted
                job_id = extract_job_id_from_log_file(log_file)
                if job_id:
                    job_fd.write('{}\t{}\n'.format(job_id, keyword))
                    print('Keyword "{}" sumbitted with job ID "{}"'.format(keyword, job_id))
                else:
                    continue
                # If counter exceeds job limit, sleep X time
                counter += 1
                if counter >= job_limit:
                    counter = 0
                    time.sleep(sleep_time)

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
        elif scoring == "netzcore":
            rscript_fd.write('submit.job(filt.result.table, species="{}", tissue="{}", network.source="{}", scoring.options = list(netzcore=T, iterationSelector={}))\n'.format(taxid, tissue, network, iterations))
        else:
            rscript_fd.write('submit.job(filt.result.table, species="{}", tissue="{}", network.source="{}", scoring.options = list({}=T))\n'.format(taxid, tissue, network, scoring))
        return

def extract_job_id_from_log_file(log_file):
    """
    Obtain the job id from the log file of the guildifyR submission
    """
    job_id = None
    if fileExist(log_file):
        with open(log_file, 'r') as log_fd:
            for line in log_fd:
                if 'Job id: ' in line:
                    # [1] "Job id: 4682695b-d09c-4769-822e-c4549d2dbd2e"
                    # ['[1] "', '4682695b-d09c-4769-822e-c4549d2dbd2e"']
                    fields = line.strip().split('Job id: ') 
                    job_id = fields[1].rstrip('"')
        if job_id:
            return job_id
    print('NO JOB ID IN LOG FILE: {}'.format(log_file))
    return None

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