# Load variables
version <- 'v7'
type_analysis <- 'tpm'

# Define file names
if (version == 'v7') {
  setwd('/Users/quim/Documents/Databases/GTEx/v7')
  counts_file = '/Users/quim/Documents/Databases/GTEx/v7/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct'
  tpm_file = '/Users/quim/Documents/Databases/GTEx/v7/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct'
  subjects_file = '/Users/quim/Documents/Databases/GTEx/v7/GTEx_v7_Annotations_SubjectPhenotypesDS.txt'
  samples_file = '/Users/quim/Documents/Databases/GTEx/v7/GTEx_v7_Annotations_SampleAttributesDS.txt'
  counts_filtered_file = '/Users/quim/Documents/Databases/GTEx/v7/gene_counts_filtered_by_tissue.gct'
  counts_median_file = '/Users/quim/Documents/Databases/GTEx/v7/gene_counts_filtered_by_tissue_median.gct'
  tpm_filtered_file = '/Users/quim/Documents/Databases/GTEx/v7/gene_tpm_filtered_by_tissue.gct'
  tpm_median_file = '/Users/quim/Documents/Databases/GTEx/v7/gene_tpm_filtered_by_tissue_median.gct'
  tpm_median_unified_file = '/Users/quim/Documents/Databases/GTEx/v7/gene_tpm_filtered_by_tissue_median_unified.gct'
} else if (version == 'v6p') {
  setwd('/Users/quim/Documents/Databases/GTEx/v6p')
  counts_file = '/Users/quim/Documents/Databases/GTEx/v6p/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct'
  tpm_file = '/Users/quim/Documents/Databases/GTEx/v6p/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct'
  subjects_file = '/Users/quim/Documents/Databases/GTEx/v6p/GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt'
  samples_file = '/Users/quim/Documents/Databases/GTEx/v6p/GTEx_Data_V6_Annotations_SampleAttributesDS.txt'
  counts_filtered_file = '/Users/quim/Documents/Databases/GTEx/v6p/gene_counts_filtered_by_tissue.gct'
  counts_median_file = '/Users/quim/Documents/Databases/GTEx/v6p/gene_counts_filtered_by_tissue_median.gct'
  tpm_filtered_file = '/Users/quim/Documents/Databases/GTEx/v6p/gene_tpm_filtered_by_tissue.gct'
  tpm_median_file = '/Users/quim/Documents/Databases/GTEx/v6p/gene_tpm_filtered_by_tissue_median.gct'
  tpm_median_unified_file = '/Users/quim/Documents/Databases/GTEx/v6p/gene_tpm_filtered_by_tissue_median_unified.gct'
} else if (version == 'v6') {
  setwd('/Users/quim/Documents/Databases/GTEx/v6')
  counts_file = '/Users/quim/Documents/Databases/GTEx/v6/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct'
  tpm_file = '/Users/quim/Documents/Databases/GTEx/v6/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct'
  subjects_file = '/Users/quim/Documents/Databases/GTEx/v6/GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt'
  samples_file = '/Users/quim/Documents/Databases/GTEx/v6/GTEx_Data_V6_Annotations_SampleAttributesDS.txt'
  counts_filtered_file = '/Users/quim/Documents/Databases/GTEx/v6/gene_counts_filtered_by_tissue.gct'
  counts_median_file = '/Users/quim/Documents/Databases/GTEx/v6/gene_counts_filtered_by_tissue_median.gct'
  tpm_filtered_file = '/Users/quim/Documents/Databases/GTEx/v6/gene_tpm_filtered_by_tissue.gct'
  tpm_median_file = '/Users/quim/Documents/Databases/GTEx/v6/gene_tpm_filtered_by_tissue_median.gct'
  tpm_median_unified_file = '/Users/quim/Documents/Databases/GTEx/v6/gene_tpm_filtered_by_tissue_median_unified.gct'
}

# Read subjects dataframe and filter subjects by traumatic injury
subjects_df = read.csv(subjects_file, header = TRUE, sep = "\t")
#subjects_filtered_df <- subjects_df[(subjects_df$DTHHRDY==1 | subjects_df$DTHHRDY==2),]
subjects_filtered_df <- subjects_df[subjects_df$DTHHRDY==1,] # Filter by DTHHRDY=1, traumatic injury
subjects_filtered <- unique(subjects_filtered_df["SUBJID"])
subjects_filtered <- subjects_filtered[!is.na(subjects_filtered)]

# Read samples dataframe
samples_df = read.csv(samples_file, header = TRUE, sep = "\t")

# Add column with subject ID 
# sample ID --> GTEX-14753-1626-SM-5NQ9L
# subject ID --> GTEX-14753
# info --> https://sites.google.com/broadinstitute.org/gtex-faqs/home
donor_ids <- list()
i = 1
for(sample_id in samples_df$SAMPID){
  split_list <- unlist(strsplit(sample_id, "-", fixed = TRUE))
  donor_id <- paste(split_list[1], split_list[2], sep='-')
  donor_ids[i] <- donor_id
  i<- i+1
}
samples_df["SUBJID"] <- unlist(donor_ids)
#samples_df[c("SAMPID", "SUBJID")] # Print columns

# Merge samples dataframe with subjects dataframe
samples_df <- merge(x = samples_df, y = subjects_df, by = "SUBJID", all = TRUE)
#samples_df[c("SAMPID", "SUBJID", "DTHHRDY")] # Print columns

# Filter samples by traumatic injury (DTHHRDY=1)
samples_filtered_df <- samples_df[samples_df$DTHHRDY==1,]
samples_filtered <- unique(samples_filtered_df["SAMPID"])
samples_filtered <- samples_filtered[!is.na(samples_filtered)]

# Get tissues and filter tissues with more than five samples 
library(plyr)
tissue_count <- count(samples_filtered_df, "SMTSD")
tissues <- tissue_count$SMTSD[tissue_count$freq >= 5] # Filter by tissues with more than 5 samples
tissues <- tissues[!is.na(tissues)]
print(length(tissues))
samples_filtered_tissues_df <- samples_filtered_df[samples_filtered_df$SMTSD %in% tissues,]
samples_filtered_tissues <- unique(samples_filtered_tissues_df["SAMPID"])
samples_filtered_tissues <- samples_filtered_tissues[!is.na(samples_filtered_tissues)]
print(length(samples_filtered_tissues))
#samples_filtered_tissues_df[c("SAMPID", "SUBJID", "DTHHRDY","SMTSD")] # Print columns

# Define the expression file (counts / tpm)
if (type_analysis == 'counts'){
  expression_file <- counts_file
  expression_filtered_file <- counts_filtered_file
  expression_median_file <- counts_median_file
  expression_median_unified_file <- counts_median_unified_file
} else if (type_analysis == 'tpm'){
  expression_file <- tpm_file
  expression_filtered_file <- tpm_filtered_file
  expression_median_file <- tpm_median_file
  expression_median_unified_file <- tpm_median_unified_file
}

# Create the expression file filtered by the samples of interest
if (!file.exists(expression_filtered_file)){
  # Read expression dataframe
  expression_df = read.csv(expression_file, header = TRUE, sep = "\t", skip=2)
  column_names = colnames(expression_df)
  new_columns <- list()
  for(i in 1:length(column_names)){
    column_name = column_names[i]
    new_name = gsub("\\.", "-", column_name)
    new_columns[i] <- new_name
  }
  colnames(expression_df) <- new_columns
  
  # Get number of filtered samples in the expression file 
  print(length(new_columns[new_columns %in% samples_filtered_tissues]))
  
  # Filter expression data by the filtered samples
  expression_filtered_columns <- c(c("Name", "Description"), samples_filtered_tissues)
  expression_filtered_df <- expression_df[colnames(expression_df) %in% expression_filtered_columns]
  
  # Write expression file filtered
  write.table(expression_filtered_df, file = expression_filtered_file, sep = "\t", row.names = FALSE)
}

# Read expression dataframe filtered
expression_filtered_df = read.csv(expression_filtered_file, header = TRUE, sep = "\t")

# Put SAMPID names correctly (with '-')
new_columns <- list()
for(i in 1:length(colnames(expression_filtered_df))){
  column_name = colnames(expression_filtered_df)[i]
  new_name = gsub("\\.", "-", column_name)
  new_columns[i] <- new_name
}
colnames(expression_filtered_df) <- new_columns
print(colnames(expression_filtered_df))



# Create new dataframe with the median of the samples expression for each tissue
if (!file.exists(expression_median_file)){
  expression_median_df <- expression_filtered_df[,c("Name", "Description")]

  for (tissue in tissues){
    print(tissue)
    samples <- samples_filtered_tissues_df[samples_filtered_tissues_df$SMTSD==tissue,]$SAMPID
    #print(samples)
    #print(expression_filtered_df[samples,])
    median_values <- apply(expression_filtered_df[colnames(expression_filtered_df) %in% samples], 1, median)
    print(median_values)
    expression_median_df[,tissue] <- median_values
    #for(i in 1:length(row.names(expression_filtered_df))){
    #  row <- expression_filtered_df[i,]
    #  expression <- row[colnames(row) %in% samples]
    #  print(unlist(expression))
    #  print(median(unlist(expression)))
    #}
  }
  
  # Write the dataframe
  write.table(expression_median_df, file = expression_median_file, sep = "\t", row.names = FALSE)
}

expression_median_df = read.csv(expression_median_file, header = TRUE, sep = "\t")

# Unify subtissues in a main tissue
print(colnames(expression_median_df))
expression_median_unified_df <- expression_median_df[,c("Name", "Description")]

adipose_tissues <- c("Adipose...Subcutaneous", "Adipose...Visceral..Omentum.")
adrenal_tissues <- c("Adrenal.Gland")
artery_tissues <- c("Artery...Aorta", "Artery...Coronary", "Artery...Tibial")
brain_tissues <- c("Brain...Amygdala", "Brain...Anterior.cingulate.cortex..BA24.", "Brain...Caudate..basal.ganglia.", "Brain...Cerebellar.Hemisphere", "Brain...Cerebellum", "Brain...Cortex", "Brain...Frontal.Cortex..BA9.", "Brain...Hippocampus", "Brain...Hypothalamus", "Brain...Nucleus.accumbens..basal.ganglia.", "Brain...Putamen..basal.ganglia.", "Brain...Spinal.cord..cervical.c.1.")
breast_tissues <- c("Breast...Mammary.Tissue")
cells_tissues <- c("Cells...Transformed.fibroblasts")
colon_tissues <- c("Colon...Sigmoid")
esophagus_tissues <- c("Esophagus...Gastroesophageal.Junction", "Esophagus...Mucosa", "Esophagus...Muscularis")
heart_tissues <- c("Heart...Atrial.Appendage", "Heart...Left.Ventricle")
liver_tissues <- c("Liver")
lung_tissues <- c("Lung")
muscle_tissues <- c("Muscle...Skeletal")
nerve_tissues <- c("Nerve...Tibial")
ovary_tissues <- c("Ovary")
pituitary_tissues <- c("Pituitary")
prostate_tissues <- c("Prostate")
skin_tissues <- c("Skin...Not.Sun.Exposed..Suprapubic.", "Skin...Sun.Exposed..Lower.leg.")
testis_tissues <- c("Testis")
thyroid_tissues <- c("Thyroid")
uterus_tissues <- c("Uterus")
vagina_tissues <- c("Vagina")
blood_tissues <- c("Whole.Blood")

expression_median_unified_df[,'Adipose'] <- apply(expression_median_df[colnames(expression_median_df) %in% adipose_tissues], 1, max)
expression_median_unified_df[,'Adrenal Gland'] <- apply(expression_median_df[colnames(expression_median_df) %in% adrenal_tissues], 1, max)
expression_median_unified_df[,'Artery'] <- apply(expression_median_df[colnames(expression_median_df) %in% artery_tissues], 1, max)
expression_median_unified_df[,'Brain'] <- apply(expression_median_df[colnames(expression_median_df) %in% brain_tissues], 1, max)
expression_median_unified_df[,'Breast - Mammary Tissue'] <- apply(expression_median_df[colnames(expression_median_df) %in% breast_tissues], 1, max)
expression_median_unified_df[,'Cells - Transformed fibroblasts'] <- apply(expression_median_df[colnames(expression_median_df) %in% cells_tissues], 1, max)
expression_median_unified_df[,'Colon - Sigmoid'] <- apply(expression_median_df[colnames(expression_median_df) %in% colon_tissues], 1, max)
expression_median_unified_df[,'Esophagus'] <- apply(expression_median_df[colnames(expression_median_df) %in% esophagus_tissues], 1, max)
expression_median_unified_df[,'Heart'] <- apply(expression_median_df[colnames(expression_median_df) %in% heart_tissues], 1, max)
expression_median_unified_df[,'Liver'] <- apply(expression_median_df[colnames(expression_median_df) %in% liver_tissues], 1, max)
expression_median_unified_df[,'Lung'] <- apply(expression_median_df[colnames(expression_median_df) %in% lung_tissues], 1, max)
expression_median_unified_df[,'Muscle - Skeletal'] <- apply(expression_median_df[colnames(expression_median_df) %in% muscle_tissues], 1, max)
expression_median_unified_df[,'Nerve - Tibial'] <- apply(expression_median_df[colnames(expression_median_df) %in% nerve_tissues], 1, max)
expression_median_unified_df[,'Ovary'] <- apply(expression_median_df[colnames(expression_median_df) %in% ovary_tissues], 1, max)
expression_median_unified_df[,'Pituitary'] <- apply(expression_median_df[colnames(expression_median_df) %in% pituitary_tissues], 1, max)
expression_median_unified_df[,'Prostate'] <- apply(expression_median_df[colnames(expression_median_df) %in% prostate_tissues], 1, max)
expression_median_unified_df[,'Skin'] <- apply(expression_median_df[colnames(expression_median_df) %in% skin_tissues], 1, max)
expression_median_unified_df[,'Testis'] <- apply(expression_median_df[colnames(expression_median_df) %in% testis_tissues], 1, max)
expression_median_unified_df[,'Thyroid'] <- apply(expression_median_df[colnames(expression_median_df) %in% thyroid_tissues], 1, max)
expression_median_unified_df[,'Uterus'] <- apply(expression_median_df[colnames(expression_median_df) %in% uterus_tissues], 1, max)
expression_median_unified_df[,'Vagina'] <- apply(expression_median_df[colnames(expression_median_df) %in% vagina_tissues], 1, max)
expression_median_unified_df[,'Whole Blood'] <- apply(expression_median_df[colnames(expression_median_df) %in% blood_tissues], 1, max)

write.table(expression_median_unified_df, file = expression_median_unified_file, sep = "\t", row.names = FALSE)

