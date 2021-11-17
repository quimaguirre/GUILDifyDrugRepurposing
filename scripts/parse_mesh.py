import sys, os, re
import pandas as pd

def main():

    BIANA_USER = 'quim'
    BIANA_PASS = ''
    BIANA_HOST = 'localhost'
    BIANA_DATABASE = 'test_BIANA_MAY_2018'
    BIANA_UNIFICATION = 'geneID_seqtax_drugtarget'
    mesh_file = '/home/quim/Databases/MESH/MESH.csv'
    main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

    df = pd.read_csv(mesh_file, sep=',')
    cuis_from_df = list(df['CUI'])
    cuis = set()
    for cui in cuis_from_df:
        if str(cui) != 'nan':
            for sub_cui in cui.split('|'):
                cuis.add(sub_cui)
    print(len(cuis))


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