#!/usr/bin/env Rscript

library(DESeq2) 
library(optparse)
library(data.table)
library(stringr)
options(echo=T)
########################## read arguments

GTF.file <- "/SAN/vyplab/HuRNASeq/GENCODE/gencode.v25.annotation.gtf"
featureCounts.list <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/CIBERSORT/Prudencio_CBL/Prudencio_CBL_counts_list.txt"
support.frame <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/CIBERSORT/Prudencio_CBL/Prudencio_CBL_support.tab"
outFolder <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/CIBERSORT/Prudencio_CBL/"
code <- "Prudencio_CBL"


option_list <- list(
    make_option(c('--support.frame'), help='', default=''),
    make_option(c('--code'), help='', default = ""),
    make_option(c('--GTF'), help='', default=""),
    make_option(c('--outFolder'), help='', default=""),
    make_option(c('--featureCounts.list'))
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

support.frame <- opt$support.frame
code <- opt$code
GTF.file <- opt$GTF
outFolder <- opt$outFolder
featureCounts.list <- opt$featureCounts.list

# do checks!

for( file in c(support.frame, featureCounts.list, outFolder, GTF.file)){
  if( ! file.exists(file) ){
    stop( paste0( file, "doesn't exist!"))
  }
}

# create vector of count file paths
count.files <- readLines( featureCounts.list)
if( any( !file.exists(count.files))){
  print(subset(count.files, ! file.exists(count.files)))
  stop( "some count files are missing!")
}



######## outputs are RPKM file, size factors and raw counts
FPKM.file <- paste0(outFolder,"/", code, "_all_FPKMs.txt")
FPKM.HUGO.file <- paste0( outFolder, "/", code, "_all_FPKMs_HUGO.txt")
FPKM.ensembl.file <- paste0( outFolder, "/", code, "_all_FPKMs_ensembl.txt" )
raw.counts <- paste0(outFolder,"/", code, "_raw_counts.txt")
size.factors <- paste0(outFolder,"/",code,"_size_factors.txt")


# bring in the GTF. This will provide the gene name and the gene length too!

GTF <- fread(GTF.file)
GTF <- GTF[ GTF$V3 == "gene", ]
GTF$ensemblID <- sub( "\"", "" ,
                  sub( "[a-z].* \"", "", 
                      str_split_fixed( GTF$V9, "; ", 6)[,1] 
                      )
                  )
# find the "gene_name" field in column 9 as it is different between mouse and human, stupidly

col9 <- str_split_fixed( as.character( head(GTF$V9,100) ), "; ", n = 30 ) 

gene_position <- median( apply( col9, MAR = 1, FUN = function(x){
  which(grepl("gene_name", x) )
  }) )
# should return 4 for human and 3 for mouse


GTF$gene <- sub( "\"", "" ,
                  sub( "[a-z].* \"", "", 
                      str_split_fixed( GTF$V9, "; ", 6)[, gene_position] 
                      )
                  )
GTF$V9 <- NULL




###check input files and data frame
message('Now reading ', support.frame)
support <- read.table(support.frame, header = FALSE, stringsAsFactors = FALSE)



sample.ids <- support$V2




sampleTable <- data.frame( sample = sample.ids , count.file = count.files, condition = c( rep(1, floor(length(count.files) / 2) ), rep(2, ceiling(length(count.files) / 2) ) ) )


dds <- DESeqDataSetFromHTSeqCount( sampleTable = sampleTable, design = ~ condition, directory = "")
# get size factors
dds <- estimateSizeFactors(dds)
sizeFactors <- sizeFactors(dds)

GTF$length <- GTF$V5 - GTF$V4

genes.counts <- counts(dds)

gene.lengths <- GTF$length[ match( rownames( genes.counts), GTF$ensemblID) ]




rpkms <- genes.counts

average.depth <- sum(as.numeric(rpkms))/(10^6*ncol(rpkms))

########## Now compute the length of each feature
compute.lengths <- data.frame( id = dimnames(rpkms)[[1]], length = NA)

compute.lengths$length <- GTF$length[ match( compute.lengths$id, table = GTF$ensemblID) ]



###### Now normalize for everything
for (i in 1:ncol(rpkms)) {

  rpkms[, i ] <- rpkms[, i ] /( (compute.lengths$length/1000)* average.depth * sizeFactors [ i ])
}
my.median <- apply(rpkms, MAR = 1, FUN = median)
rpkms <- as.data.frame(rpkms [ order( my.median, decreasing = TRUE), ]) ##reorder

samples.names <- names(rpkms)

####### Now fix the gene names
rpkms$external_gene_id <- as.character(GTF$gene [ match(row.names(rpkms), GTF$ensemblID) ])



########## finalize the computation
rpkms$ensemblID <- row.names( rpkms ) 
rpkms <- rpkms[, c('external_gene_id', 'ensemblID', samples.names) ]

rpkms_HUGO <- rpkms[, c('external_gene_id', samples.names )]
rpkms_ensembl <- rpkms[ , c('ensemblID', samples.names )]

########## write out all the files
write.table( x = rpkms, file = FPKM.file, row.names = FALSE, quote = FALSE, sep = "\t")
write.table( x = rpkms_HUGO, file = FPKM.HUGO.file, row.names = FALSE, quote = FALSE, sep = "\t" )
write.table( x = rpkms_ensembl, file = FPKM.ensembl.file, row.names = FALSE, quote = FALSE, sep = "\t" )


write.table( x = data.frame( sample =  names(sizeFactors), sizeFactors = sizeFactors) ,
           file = size.factors, row.names = FALSE, quote = TRUE, sep = "\t" )
write.table( x = genes.counts, file = raw.counts, row.names = TRUE, quote = TRUE)

message("Done printing FPKM values")



