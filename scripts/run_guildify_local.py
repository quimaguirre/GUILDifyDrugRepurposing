import sys, os, re
import argparse
import time



def main():

    options = parse_user_arguments()
    run_guildify_local(options)

def parse_user_arguments(*args, **kwds):
    """
    Parses the arguments of the program.
    Examples:
    python /home/quim/GUILDifyTools/scripts/run_guildify_local.py -i "TP53; BRCA1; BRCA2" -t 9606 -f /var/www/html/sbi/guildp/data/BIANA_phy/9606/all/node_info.txt -o /home/quim/data/guildify_local_outputs/example
    python /home/quim/GUILDifyTools/scripts/run_guildify_local.py -i BRCA1 -t 9606 -f /var/www/html/sbi/guildp/data/BIANA_phy/9606/all/node_info.txt -o /home/quim/data/guildify_local_outputs/example_BRCA1
    """

    parser = argparse.ArgumentParser(
        description = "Generate the profiles of GUILDify",
        epilog      = "@oliva's lab 2023")
    parser.add_argument('-i','--input',dest='input',action = 'store',
                        help = """Input seeds (in gene symbol notation). If there is more than one seed, separate them using semicolons (e.g. "TP53; BRCA1; BRCA2")""")
    parser.add_argument('-n','--network_file',dest='network_file',action = 'store', default='biana',
                        help = """ Network file. """)
    parser.add_argument('-f','--node_info_file',dest='node_info_file',action = 'store', default='biana',
                        help = """ Node information file """)
    parser.add_argument('-d','--databases',dest='databases',action = 'store', default='disgenet,omim,uniprot,go,drugbank,dgidb,drugcentral,chembl',
                        help = """ Name of the databases (separated by , if more than one):
                                   For drugs: drugbank,dgidb,drugcentral,chembl
                                   For phenotypes: disgenet,omim,uniprot,go""")
    parser.add_argument('-t','--taxid',dest='taxid',action = 'store', default='9606',
                        help = """ Taxonomy ID of the species of the network: 9606 (human) / 10090 (mouse) / 3702 (plant) / 6239 (worm) /
                                   7227 (fly) / 10116 (rat) / 4932 (yeast) """)
    parser.add_argument('-ti','--tissue',dest='tissue',action = 'store', default='all',
                        help = """ Type of tissue. If all of them considered, enter 'all' """)
    parser.add_argument('-ns','--network_source',dest='network_source',action = 'store', default='BIANA_phy',
                        help = """ Name of the type of network considered (e.g. 'BIANA_phy' or 'BIANA_func') """)
    parser.add_argument('-s','--scoring',dest='scoring',action = 'store',
                        help = """ Scoring function name: netscore / netzcore / netshort / netcombo / diamond. """)
    parser.add_argument('-re','--repetitions',dest='repetitions',action = 'store', default='3',
                        help = """ Number of repetitions (for netscore and netcombo) """)
    parser.add_argument('-it','--iterations',dest='iterations',action = 'store', default='2',
                        help = """ Number of iterations (for netscore, netzcore and netcombo) """)
    parser.add_argument('-o','--output_dir',dest='output_dir',action = 'store',
                        help = """ Output directory. """)

    options=parser.parse_args()

    return options


#################
#################
# MAIN FUNCTION #
#################
#################

def run_guildify_local(options):
    """
    Generates GUILDify profiles
    """
    # Example of outputs
    #ls /var/www/html/sbi/guildify2/data/sessions/8075b9ad-3cec-4138-83c5-8a03d3a2fab8/

    #-------------------#
    #  PROCESS INPUTS   #
    #-------------------#

    # Create a directory for GUILD results
    output_dir = options.output_dir
    create_directory(output_dir)

    # Get associated gene_symbols
    words = options.input.split(";")
    if len(words) > 1:
        # Gene symbol list if query contains ";" 
        keywords = [ word.strip() for word in words if word != "" ]
    else:
        keywords = words
    print(keywords)

    # Read node info file
    user_entity_to_values = get_user_entities_info(options.node_info_file)
    gene_symbol_to_user_entity = get_gene_symbol_to_user_entity_dict(user_entity_to_values)

    # Create seed info file
    seed_info_file = os.path.join(output_dir, "seeds.txt")
    create_seed_info_file(keywords, user_entity_to_values, gene_symbol_to_user_entity, seed_info_file)
    sys.exit(0)

    # Save the species
    species = str(options.taxid).lower()
    species_file = os.path.join(output_dir, "species.txt")
    open(species_file, "w").write("{}\n".format(species))

    # Save the tissue
    tissue = str(options.tissue).lower()
    tissue_file = os.path.join(output_dir, "tissue.txt")
    open(tissue_file, "w").write("{}\n".format(tissue))

    # Save the network source
    network_source = str(options.network_source)
    network_source_file = os.path.join(output_dir, "network_source.txt")
    open(network_source_file, "w").write("{}\n".format(network_source))


    #-------------------------------------#
    #  SCORE NODES IN THE NETWORK (GUILD) #
    #-------------------------------------#

    # Run GUILD
    pvalue_file = os.path.join(guild_output_dir, 'output_scores.sif.netcombo.pval')
    if not fileExist(pvalue_file):

        guild_command = 'python {} {} {} {} {} {} {}'.format( os.path.join(toolbox_dir, 'run_guild.py'), drug_dir, network_targets_file, options.sif, guild_output_dir, random_networks_dir, config.get('Paths', 'guild_path') )
        os.system(guild_command)
        print('  DIANA INFO:\tGUILD has finished.\n')

    else:
        print('  DIANA INFO:\tThe scoring of the network with GUILD for {} was already done and it has been skipped.\n'.format(options.drug_name))

    # Creating an instance of the file generated by GUILD
    guild_profile_instance = network_analysis.GUILDProfile(pvalue_file, network_instance.type_id, 100)

    # Translate the NODE profile to protein_type_id if the type of id is 'biana'
    if guild_profile_instance.type_id == 'biana' and translation_file:
        output_file = os.path.join(guild_output_dir, 'node_profile_top_100_{}.txt'.format(options.proteins_type_id))
        guild_profile_geneid = guild_profile_instance.translate_pvalue_file(translation_file, options.proteins_type_id, output_file, verbose=False)



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
    Checks if a directory exists and if not, creates it.
    """
    try:
        os.stat(directory)
    except:
        os.mkdir(directory)
    return

def get_user_entities_info(file_name, user_entity_ids = None):
    """
    Get the information of the user entities from 'file_name'.
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
        user_entity_to_values[user_entity_id] = [entry_ids, gene_symbols.strip(), descriptions, gene_ids, other_ids]
    return user_entity_to_values 

def remove_isoforms(entry_ids):
    """
    Remove the isoforms from a list of uniprot accessions.
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

def get_gene_symbol_to_user_entity_dict(user_entity_to_values):
    """
    Get a dictionary that maps gene symbols to BIANA user entities.
    """
    gene_symbol_to_user_entity = {}
    for user_entity_id in user_entity_to_values:
        [entry_ids, gene_symbols, descriptions, gene_ids, other_ids] = user_entity_to_values[user_entity_id]
        for gene_symbol in gene_symbols.split('; '):
            gene_symbol_to_user_entity.setdefault(gene_symbol, set()).add(user_entity_id)
    return gene_symbol_to_user_entity

def create_seed_info_file(keywords, user_entity_to_values, gene_symbol_to_user_entity, output_file):
    """
    Create file containing the seeds and the information about the seeds.
    """
    # Get the user entities of the seeds
    user_entity_ids = set()
    for keyword in keywords:
        print(keyword)
        if keyword in gene_symbol_to_user_entity:
            print(keyword)
            user_entities_keyword = gene_symbol_to_user_entity[keyword]
            print(user_entities_keyword)
            for user_entity in user_entities_keyword:
                user_entity_ids.add(user_entity)
    # Write the output file containing information about each seed
    with open(output_file, "w") as out_fd:
        out_fd.write("BIANA ID\tUniProt ID\tGene Symbol\tDescription\tGene ID\tEquivalent Entries\n") # Write header
        for user_entity_id in user_entity_ids:
            entry_ids, gene_symbols, descriptions, gene_ids, equivalent_entries = user_entity_to_values[user_entity_id]
            if entry_ids == "-":
                continue
            entry_ids = entry_ids.strip("; ")
            inner_values = [user_entity_id, entry_ids, gene_symbols, descriptions, gene_ids, equivalent_entries]
            out_fd.write("{}\n".format("\t".join(inner_values)))
    return

if  __name__ == "__main__":
    main()