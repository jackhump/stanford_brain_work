# using all the GTEX data, extract the TARDBP expression levels for each brain sample. Correlate TARDBP expression with the expression of every other gene and create QQ-plots of the Spearman and Pearson correlation tests.


library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(qqman)
options(echo=T)
iFolder <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/GTEX"
oFolder <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/GTEX/Brain"
plotFolder <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/GTEX/plots"
if(!dir.exists(plotFolder)){dir.create(plotFolder) }
if(!dir.exists(oFolder)){dir.create(oFolder) }
# input files
annotation_file="/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Human_hg38/biomart_annotations_human.tab"
gtex_samples="/SAN/vyplab/HuRNASeq/opthalmology_work/GTEx_Data_V6_Annotations_SampleAttributesDS.txt"
gtex_subjects="/SAN/vyplab/HuRNASeq/opthalmology_work/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct"
gtex_covariates="/SAN/vyplab/HuRNASeq/opthalmology_work/GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt"

gtex_samples <- fread(gtex_samples, header=T )
gtex_cov <- fread(gtex_covariates, header=T )

all_brain <- subset(gtex_samples, SMTS == "Brain")

if(length(list.files(oFolder)) ==0 ){
  # load gtex - this takes ~ 2.5 minutes
  gtex <- as.data.frame( fread(gtex_subjects, header=T ) )

  # save just the gene IDs and names
  gtex_sample_key <- gtex[, 1:2]
  write.table( gtex_sample_key, paste0(oFolder, "/", "all_genes_sample_key.tab"), sep = "\t", col.names = T, quote = F )


  # extract each brain region and save expression data separately
  for( region in names(table(all_brain$SMTSD)) ){

    samples <- subset(gsample, gsample$SMTSD == region)$SAMPID

    expression <- gtex[, which(names(gtex) %in% c("Name", "Description",samples) ) ]

    outCode <- str_split_fixed(region, " - ", 2)[,2]
    outCode <- gsub( " \\(.*\\)", "", outCode )
    outCode <- gsub( " ", "_", outCode )
    outCode <- tolower( outCode )
   
    outFile <- paste0( oFolder, "/", outCode,".Rdata" )

    save( expression, file = outFile )

  }
}


# for each dataset, correlate TARDBP expression with every other gene
# to get a numeric vector, ignore the first two rows

allRegions <- list.files(oFolder, pattern = "Rdata",full.names = TRUE)

for( region in allRegions ){

  load(region)

  # file is "expression"
  row.names(expression) <- expression$Name
  expression <- select( expression, -Name, -Description)

  # TARDBP is ENSG00000120948.11 

 #remove genes with zero expression
  meanExpression <- apply( expression, MAR = 1, FUN = mean)
  expression <- expression[ meanExpression > 1,]

  tdpExpression <- as.numeric( subset(expression, row.names(expression) == "ENSG00000120948.11") )

  # remove TARDBP from the list of genes
  expression <- subset( expression, row.names(expression) != "ENSG00000120948.11")

  correlation <- apply( expression, MAR = 1, FUN = function(x) cor.test(x, tdpExpression, method = "pearson")$p.value )

  # create a null distribution where the correlation test where the TDP expression values have been shuffled

  nullExpression <- sample(x = tdpExpression, replace = FALSE )

  nullCorrelation <- apply( expression, MAR = 1, FUN = function(x) cor.test(x, nullExpression, method = "pearson")$p.value )

  # quick and dirty with qqman

  outCode <- basename(region)
  outCode <- gsub(".Rdata", "", outCode)

  pdf( paste0( plotFolder, "/", outCode,"_correlation_qq.pdf"))
  # side by side
  par(mfrow=c(1,2))
  qq( correlation, main = paste0( outCode, "\n TARDBP Pearson correlation"))
  qq( nullCorrelation, main = paste0( outCode, "\n null correlation "))

  dev.off()
}

for( region in allRegions ){
  # find all the covariates for each sample
  # age and sex require the subject ID
  subjectID <- paste( str_split_fixed( names(expression), "-", 3 )[,1], str_split_fixed( names(expression), "-", 3 )[,2], sep = "-" )

  sampleAge <- gtex_cov$AGE[ match( subjectID, gtex_cov$SUBJID )]
  sampleSex <- gtex_cov$GENDER[ match( subjectID, gtex_cov$SUBJID )]
  # all the others can be done with the full sample name
  sampleBatch <- gtex_samples$SMNABTCHD[ match( names(expression), gtex_samples$SAMPID )]
  sampleRIN <- gtex_samples$SMRIN[ match( names(expression), gtex_samples$SAMPID )]
  sampleLibSize <- gtex_samples$SMESTLBS[ match( names(expression), gtex_samples$SAMPID )]

  # regress gene expression against TDP expression for each genes.
  # apply row-wise and return list of lm objects

  lmFull <- apply( expression, MAR = 1, FUN = function(x){
  lm( as.numeric(x) ~ tdpExpression + sampleAge + sampleSex + sampleBatch + sampleRIN + sampleLibSize)
  })
  
  lmReduced <- apply( expression, MAR = 1, FUN = function(x){ 
  lm( as.numeric(x) ~ sampleAge + sampleSex + sampleBatch + sampleRIN + sampleLibSize)
  })

  lmNull <- apply( expression, MAR = 1, FUN = function(x){ 
  lm( as.numeric(x) ~ nullExpression + sampleAge + sampleSex + sampleBatch + sampleRIN + sampleLibSize)
  })

  # mapply is multiple list apply
  anovaPFull <- mapply( function(X,Y) { unlist( anova( X, Y ) )[ "Pr(>F)2" ]  }, SIMPLIFY = TRUE, X = lmFull, Y = lmReduced)
  anovaPNull <- mapply( function(X,Y) { unlist( anova( X, Y ) )[ "Pr(>F)2" ]  }, SIMPLIFY = TRUE, X = lmNull, Y = lmReduced)

  pdf( paste0( plotFolder, "/", outCode,"_confounders_qq.pdf"))
  # side by side
  par(mfrow=c(1,2))
  qq( anovaPFull, main = paste0( outCode, "\n TARDBP confounders"))
  qq( anovaPNull, main = paste0( outCode, "\n null confounders "))

  dev.off()


}

