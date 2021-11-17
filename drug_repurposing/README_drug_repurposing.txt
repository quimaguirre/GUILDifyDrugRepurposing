######################
#  DRUG REPURPOSING  #
######################

#---------------------#
# 1. LIST OF DISEASES #
#---------------------#

1. Execute get_diseases_from_disgenet.py
    
    Location: /home/quim/PHD/Projects/GUILDify/scripts/get_diseases_from_disgenet.py

    1) It queries BIANA and retrieves all "disease" terms from DisGeNET that 
       have associations from the curated databases: 'uniprot', 'psygenet', 
       'orphanet', 'hpo'. We avoid 'ctd_human' because it contains experimental 
       terms such as "liver cirrhosis, experimental" which are only used in the 
       laboratory. The output file is called 'disease_names.txt'.

    2) It uses the CUIs of the diseases retrieved in the previous step to query
       GUILDify (from the command line) and obtain the number of disease-gene
       associations from DisGeNET. The output file is called 'disease_associations.txt'.

2. We define a cut-off of associations so that we only use the disease terms that
   have more associations and are more reliable.
   We used the cut-off 20, and manually created a new file where we have stored
   390 diseases: input_diseases_20.txt
   For a cut-off of 10, there are 757 diseaes: input_diseases_10.txt



#----------------------------#
# 2. SEND JOBS AUTOMATICALLY #
#----------------------------#

1. To query GUILDify automatically, we use the script: create.rscript_KW.pl.
   This script will automatically generate RScripts which will send jobs to GUILDify
   every certain time. We only have to modify the script to specify the folder where
   we will store those Rscripts. We also have to specify the number of jobs sent 
   every period and the time between periods. And then we execute it using as 
   input the file generated before "input_diseases_20.txt":

   nohup perl /home/quim/PHD/Projects/GUILDify/scripts/create.rscript_KW.pl /home/quim/PHD/Projects/GUILDify/drug_repurposing/input_diseases_20.txt >> jobs.ids.log &

   nohup perl /home/quim/PHD/Projects/GUILDify/scripts/create.rscript_KW.pl /home/quim/PHD/Projects/GUILDify/drug_repurposing/input_diseases_10.txt >> jobs.ids.log &

   my $limit   = 15; 
   my $sleep   = 800; 
   my $Rscriptspath = "/home/quim/PHD/Projects/GUILDify/drug_repurposing/RScripts";

   WARNING: Sometimes, the program does not run correctly (Execution halted) or
   it stops and does not send the job. In does cases, remove the R script of 
   the unsent CUI and send it again.


2. When all the jobs have been sent, we have to parse the file "jobs.ids.log" which
   will contain the job IDs generated, and create a file mapping each disease with
   its Job ID. We use the following command:

   perl -e 'while(<>){if ($_ =~ /BIANA/){@t=split(/: /); $id=$t[$#t];} if ($_ =~ /^\[1\]\s\"Job id\:/){@t=split; $t[$#t]=~s/\"$//; print $t[$#t]."\t". $id;}}' < jobs.ids.log > /home/quim/PHD/Projects/GUILDify/drug_repurposing/job_to_disease_mapping.txt

   perl -e 'while(<>){if ($_ =~ /BIANA/){@t=split(/: /); $id=$t[$#t];} if ($_ =~ /^\[1\]\s\"Job id\:/){@t=split; $t[$#t]=~s/\"$//; print $t[$#t]."\t". $id;}}' < jobs.ids.log > /home/quim/PHD/Projects/GUILDify/drug_repurposing/job_to_disease_mapping_10.txt

   The mapping file will be the file "job_to_disease_mapping.txt" (or "job_to_disease_mapping_10.txt")


3. When all the jobs have finished, we have to retrieve them. This is done 
   executing the following script:

   perl /home/quim/PHD/Projects/GUILDify/scripts/retrieve_expansions_KW.pl /home/quim/PHD/Projects/GUILDify/drug_repurposing/job_to_disease_mapping.txt

   my $Rscriptdir = "/home/quim/PHD/Projects/GUILDify/drug_repurposing/RScripts_retrieve";


4. We have to copy some files in GUILDify drug_repurposing data folder of the web server:

   Folder: /var/www/html/guildify2/data/drug_repurposing

   Files:
   - "disease_names.txt"
   - "job_to_disease_mapping.txt"

   And execute GUILDify drug repurposing application. The first time it will 
   parse all the jobs and obtain a mapping of the genes and functions of all 
   the diseases.


#------------------------------#
# 3. SAME PROCEDURE WITH DRUGS #
#------------------------------#

/home/quim/PHD/Projects/GUILDify/scripts/prepare_drug_indication_benchmark.py ==> hetionet_guildify_drugs.txt

nohup perl /home/quim/PHD/Projects/GUILDify/scripts/create.rscript_drug.pl /home/quim/PHD/Projects/GUILDify/drug_repurposing/hetionet_guildify_drugs.txt >> /home/quim/PHD/Projects/GUILDify/drug_repurposing/RScripts_drugs/jobs.ids.log &

perl -e 'while(<>){if ($_ =~ /BIANA/){@t=split(/: /); $id=$t[$#t];} if ($_ =~ /^\[1\]\s\"Job id\:/){@t=split; $t[$#t]=~s/\"$//; print $t[$#t]."\t". $id;}}' < jobs.ids.log > /home/quim/PHD/Projects/GUILDify/drug_repurposing/job_to_drug_hetionet.txt

# With all drugs, not only hetionet

/home/quim/PHD/Projects/GUILDify/scripts/get_drugs_from_guildify.py ==> drug_associations_and_sources.txt ==> input_drugs_10.txt

nohup perl /home/quim/PHD/Projects/GUILDify/scripts/create.rscript_drug.pl /home/quim/PHD/Projects/GUILDify/drug_repurposing/input_drugs_10.txt >> /home/quim/PHD/Projects/GUILDify/drug_repurposing/RScripts_drugs/jobs.ids.log &

perl -e 'while(<>){if ($_ =~ /BIANA/){@t=split(/: /); $id=$t[$#t];} if ($_ =~ /^\[1\]\s\"Job id\:/){@t=split; $t[$#t]=~s/\"$//; print $t[$#t]."\t". $id;}}' < jobs.ids.log > /home/quim/PHD/Projects/GUILDify/drug_repurposing/job_to_drug_mapping_10.txt


