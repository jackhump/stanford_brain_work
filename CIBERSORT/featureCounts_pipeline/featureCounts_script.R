# Take a bam file and assign reads to genes in a GTF

library(Rsubread)
#library(stringr)
library(optparse)
options(echo=T)
########################## read arguments

option_list <- list(
	make_option(c('--bamFile'), help = '', default = '/SAN/vyplab/HuRNASeq/misc/mouse_chr8_onegene_test.bam'),
    make_option(c('--paired'), help = '', default = 'yes'),
    make_option(c('--countStrand'), help='', default='no'),
    make_option(c('--GTF '), help='', default = "/SAN/vyplab/HuRNASeq/misc/mouse_chr8_test.gtf"),
    make_option(c('--outFile'), help='', default="/SAN/vyplab/HuRNASeq/misc/mouse_chr8_test_counts.tab")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

bamFile <- opt$bamFile
paired <- opt$paired
countStrand <- opt$countStrand
GTF <- opt$GTF
outFile <- opt$outFile

# check if GTF exists
if( !file.exists(GTF) ){
  stop(paste0(GTF," does not exist!"))
}


# get strandedness correct
strandedness <- 0
if(countStrand == "yes"){
	strandedness <- 1
}
if(countStrand == "reverse"){
	strandedness <- 2
}

isPaired <- FALSE
if(paired == "yes"){
	isPaired <- TRUE
}

counts <- featureCounts(bamFile,
	reportReads=FALSE,
	annot.ext = GTF,
	isGTFAnnotationFile=TRUE,
	GTF.featureType = "exon",
  GTF.attrType = "gene_id",
	strandSpecific = strandedness,
	countMultiMappingReads=FALSE,
	useMetaFeatures=TRUE,
	allowMultiOverlap=FALSE,
	fraction=FALSE,
	ignoreDup=TRUE,
	isPairedEnd=isPaired
)

write.table(counts$counts, file = outFile, sep = "\t", col.names = FALSE, row.names = TRUE, quote=F)



