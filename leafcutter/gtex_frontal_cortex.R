library(data.table)
library(stringr)
# can we use leafcutter with the GTEX samples to boost our controls?


# example - frontal cortex samples

gtex_fc <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/GTEX/Brain/frontal_cortex.Rdata"
gtex_samples="/SAN/vyplab/HuRNASeq/opthalmology_work/GTEx_Data_V6_Annotations_SampleAttributesDS.txt"
gtex_subjects="/SAN/vyplab/HuRNASeq/opthalmology_work/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct"
gtex_covariates="/SAN/vyplab/HuRNASeq/opthalmology_work/GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt"

gtex_junctions <- list.files("/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/GTEX/Brain/gtex_junc_files", full.names =T)
gtex_sra <- gsub( ".junc", "", str_split_fixed( gtex_junctions, "_", 12 )[,12] )


#load(gtex_fc)
gtex_samples <- fread(gtex_samples)

brain_samples <- gtex_samples[ gtex_samples$SMTS == "Brain", ]

fc_samples <- subset(gtex_samples, SAMPID %in% names(expression) ) 

# SMRDLGTH is the maximum read lenght. Ideally this should match my ALS brain data
table(fc_samples$SMRDLGTH)
# all samples are 76bp paired end

# DTHHRDY is how they died
# how to exclude neurodegenerative disease samples?


# leafcutter junctions are named by the SRA id. Pull out from the SRA run table for all 9000 samples
sra_run_table <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/GTEX/gtex_samples_sra_bigger.txt"

sra_run_table <- as.data.frame(fread(sra_run_table))
#sra_run_table <- 

sra_brains <- subset( sra_run_table, Sample_Name_s %in% brain_samples$SAMPID )$Run_s

sra.brain.list <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/GTEX/Brain_SRA_samples"

writeLines( sra_brains, sra.brain.list)

present <- sapply( sra_fc, FUN = function(x) gtex_junctions[ grepl(x, gtex_junctions) ] )