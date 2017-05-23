# wrangle leafcutter output into something interpretable
library(data.table)
library(leafcutter)
library(stringr)
library(dplyr)
library(optparse)
library(readr)

FDR_limit <- 0.05

options(echo=TRUE)


outFolder <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/leafcutter/FTD_C9/Individual_samples"   
species <- "human"                                                
groups_file <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/leafcutter/FTD_C9/Individual_samples/A075_ds_support.tab"
counts_file <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/leafcutter/FTD_C9/Individual_samples/A075_perind_numers.counts.gz"
annotation_code <- "/SAN/vyplab/HuRNASeq/leafcutter/leafcutter/data/gencode_hg38"
code <- "A075"



## testing
#setwd("/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/leafcutter/F210I_norm")
#species <- "mouse"
#outFolder <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/leafcutter/F210I_norm"
#code <- "F210I_norm"
#leafcutter_dir <- "/SAN/vyplab/HuRNASeq/leafcutter/"
#annotation_code <- "/SAN/vyplab/HuRNASeq/leafcutter/leafcutter/data/gencode_mm10"
#counts_file <- paste0(outFolder, "/", code, "_perind_numers.counts.gz")
#groups_file <- paste0( outFolder, "/", code, "_ds_support.tab")
#
## more testing
#outFolder <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/leafcutter/FTD_C9"                
#species <- "human"               
#groups_file <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/leafcutter/FTD_C9/FTD_C9_ds_support.tab"                
#counts_file <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/leafcutter/FTD_C9/FTD_C9_perind_numers.counts.gz"               
#annotation_code <- "/SAN/vyplab/HuRNASeq/leafcutter/leafcutter/data/gencode_hg38"
#
#
opt <- parse_args(OptionParser(option_list=list(
  make_option( c("-o","--outFolder") ),
  make_option( "--species", default=NULL, help="which species to use. Either mm10 or hg38 are supported, for now." ),
  make_option( "--groups_file", default=NULL, help="The support file used in the differential splicing analysis. Columns should be file name and condition"),
  make_option( "--counts_file", default = NULL, help = "the perind.counts.gz file created by the cluster discovery step"),
  make_option( "--annotation_code", default=NULL, help = "a path to the annotation files of exons, introns and splice sites"),
  make_option( "--code", default=NULL, help = "the same dataset-specific code used throughout the pipeline"))
	)
)
outFolder <- opt$outFolder
species <- opt$species
groups_file <- opt$groups_file
counts_file <- opt$counts_file
annotation_code <- opt$annotation_code
code <- opt$code

print(outFolder)

print(annotation_code)

for( file in c("_all_introns.bed","_threeprime.bed","_fiveprime.bed", "_all_exons.txt.gz" )){
  if( ! file.exists( 
    paste0( annotation_code, file) 
    ) 
  ){
    stop( paste0( annotation_code, file , " does not exist!"))
  }
}

# options

# outFolder 
# species
# ds_support
# exon_file


resultsFolder <- paste0(outFolder,"/results")

if( !dir.exists(resultsFolder)){
  dir.create(resultsFolder)
}

effect.sizes.file <- paste0(outFolder,"/",code,"_ds_effect_sizes.txt") 
results.file <- paste0(outFolder, "/",code,"_ds_cluster_significance.txt")

# check if the required results files exist
for( file in c(effect.sizes.file, results.file)){
  if( !file.exists(file)){
    stop( paste0(file, " does not exist"))
  }
}

effectSizes <- fread(effect.sizes.file, stringsAsFactors = FALSE )
effectSizesSplit <-  as.data.frame(str_split_fixed(effectSizes$intron, ":", 4), stringsAsFactors = FALSE )
names(effectSizesSplit) <- c("chr","start","end","clusterID")

print(head(effectSizes))

print(head(effectSizesSplit))



effectSizes <- cbind( effectSizes, effectSizesSplit)
effectSizes$cluster <- paste(effectSizesSplit$chr, effectSizesSplit$clusterID, sep = ":")

results <- fread(results.file, stringsAsFactors = F)

results$FDR <- p.adjust( results$p, method = "fdr")

all <- merge(x = results, y = effectSizes, by = "cluster")
all <- all[ order(all$FDR),]

all <- subset( all, FDR <= FDR_limit )
all$start <- as.numeric(all$start)
all$end <- as.numeric(all$end)


# for each splice site write out a bed file  
all.fiveprime <- data.frame( chr = all$chr,
                             start = all$start,
                             end = as.numeric( as.character(all$start) ) + 1,
                             clusterID = all$clusterID)
all.threeprime <- data.frame( chr = all$chr,
                             start = all$end,
                             end = as.numeric( as.character(all$end) ) + 1,
                             clusterID = all$clusterID)
all.fiveprime.file <- paste0(resultsFolder, "/all.fiveprime.bed")
all.threeprime.file <- paste0(resultsFolder, "/all.threeprime.bed")

write.table( all.threeprime, all.threeprime.file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t" )
write.table( all.fiveprime, all.fiveprime.file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t" )

# intersect with bedtools to find the annotations of each splice site
threeprime.cmd <- paste0( "bedtools intersect -a ", all.threeprime.file, " -b ", annotation_code,"_threeprime.bed", " -wa -wb -loj -f 1" )

threeprime_intersect <- fread(threeprime.cmd)

fiveprime.cmd <- paste0( "bedtools intersect -a ", all.fiveprime.file, " -b ", annotation_code,"_fiveprime.bed", " -wa -wb -loj -f 1" )

fiveprime_intersect <- fread(fiveprime.cmd)

rm.cmd <- paste("rm ", all.fiveprime.file, all.threeprime.file) 
system(rm.cmd)

# now I have two lists of splice site annotation
# using these, can I annotate each intron in each cluster?
# for testing
#cluster <- all[ all$clusterID == "clu_4879" , ]


verdict.list <- list()
coord.list <- list()
gene.list <- list()
ensemblID.list <- list()
transcripts.list <- list()
constitutive.list <- list()

clusters <- unique( all$clusterID ) 
for( clu in clusters ){
  # for each intron in the cluster, check for coverage of both
  # output a vector of string descriptions 
  cluster <- all[ all$clusterID == clu , ]
  
      # for each intron in the cluster:
  #   create vector of overlapping splice sites, indexed by the row of the intersect
  # five prime splice sites
  fprime <- apply( cluster, MAR = 1, FUN = function(x) {
    chr <- which( names(cluster) == "chr" )
    start <- which( names(cluster) == "start" )
    fiveprime_intersect[   
      fiveprime_intersect$V1 == x[chr] & 
      fiveprime_intersect$V2 == as.numeric( x[start] ),]
  } )
  # three prime splice sites
  tprime <- apply( cluster, MAR = 1, FUN = function(x) {
    chr <- which( names(cluster) == "chr" )
    end <- which( names(cluster) == "end" )
    threeprime_intersect[   
      threeprime_intersect$V1 == x[chr] & 
      threeprime_intersect$V2 == as.numeric( x[end] ),]
  } )

  # find gene and ensemblID by the most represented gene among all the splice sites
  cluster_genes <- names(sort(table(do.call( what = rbind, tprime )$V8), decreasing = TRUE ))

  cluster_gene <- cluster_genes[ cluster_genes != "." ][1]
  # if no cluster gene found then leave as "."
  if( length(cluster_gene) == 0){
    cluster_gene == "."
  }


  cluster_ensemblIDs <- names(sort(table(do.call( what = rbind, tprime )$V9), decreasing = TRUE ))
  cluster_ensemblID <- cluster_ensemblIDs[ cluster_ensemblIDs != "." ][1]
  if( length( cluster_ensemblID ) == 0 ){
    cluster_ensemblID == "."
  }


  verdict <- c()
  coord <- c()
  gene <- c()
  ensemblID <- c()
  transcripts <- list() 
  
  for( intron in 1:nrow(cluster) ){
    coord[intron] <- paste(cluster[intron]$chr,cluster[intron]$start, cluster[intron]$end )
    # record all transcripts that use the splice sites for each intron

    # tgene <- names(sort(table( tprime[[intron]]$V8 ), decreasing = TRUE)[1])
    # fgene <- names(sort(table( fprime[[intron]]$V8 ), decreasing = TRUE)[1])

    # tensemblID <- names(sort(table( tprime[[intron]]$V8 ), decreasing = TRUE)[1])
    # fensemblID <- names(sort(table( fprime[[intron]]$V8 ), decreasing = TRUE)[1])


    # gene[intron] <- ifelse( tgene == ".",  no = fgene, yes = tgene )
    # ensemblID[intron]<- ifelse( tensemblID == ".", no = fensemblID, yes = tensemblID )

    gene[intron] <- cluster_gene
    ensemblID[intron] <- cluster_ensemblID

    # for each intron create vector of all transcripts that contain both splice sites
    transcripts[[intron]] <- unique( intersect( tprime[[intron]]$V10, fprime[[intron]]$V10 ) ) 


    verdict[intron] <- "error"
    if(
    all( tprime[[intron]]$V5 == ".") & all( fprime[[intron]]$V5 == "." )
    ){ verdict[intron] <- "cryptic_unanchored"
    }
    if(
    all( tprime[[intron]]$V5 == ".") & all( fprime[[intron]]$V5 != "." )
    ){ verdict[intron] <- "cryptic_threeprime"
    }
    if(
    all( tprime[[intron]]$V5 != ".") & all( fprime[[intron]]$V5 == "." )
    ){ verdict[intron] <- "cryptic_fiveprime"
    }
    if(
      all( tprime[[intron]]$V5 != "." ) & all( fprime[[intron]]$V5 != "." )
    ){ 
      # test if the splice sites are paired in a known intron
      tp <- paste( tprime[[intron]]$V9, tprime[[intron]]$V10 )
      fp <- paste( fprime[[intron]]$V9, fprime[[intron]]$V10 )
      if( length( intersect(fp,tp) ) > 0 ){
        verdict[intron] <- "annotated"
      }else{
        verdict[intron] <- "skiptic"
      }
    }
    verdict.list[[clu]] <- verdict
    coord.list[[clu]] <- coord
    gene.list[[clu]] <- gene
    ensemblID.list[[clu]] <- ensemblID
    #transcripts.list[[clu]] <- transcripts

    # once all the transcripts for all the introns are found, go back and work out how many constitutive each junction is. Does the junction appear in every transcript? 

    if( intron == nrow(cluster)){ # only on final intron
      all_transcripts <- unique( unlist( transcripts ) )
      # remove "." - non-existent transcripts
      all_transcripts <- all_transcripts[ all_transcripts != "." ]
 
      constitutive <- lapply( transcripts, FUN = function(x) {
        # for each intron how many transcripts is it seen in?
        x <- x[ x != "." ]
        length(x) / length( all_transcripts)

        })

      constitutive.list[[clu]] <- constitutive

      # collapse all transcripts for each intron into a single string
      transcripts.list[[clu]] <- lapply(transcripts, FUN = function(x) paste( x, collapse = "+" ) )

    }

  }
}

# verdict.list and coord.list are not equal lengths. why?

# need to match ideally
all$verdict <- unlist(verdict.list)[ match( paste( all$chr, all$start, all$end ), unlist(coord.list)) ]

all$gene <- unlist(gene.list)[ match( paste( all$chr, all$start, all$end ), unlist(coord.list)) ]

all$ensemblID <- unlist(ensemblID.list)[ match( paste( all$chr, all$start, all$end ), unlist(coord.list)) ]

all$transcripts <- unlist( transcripts.list )[ match( paste( all$chr, all$start, all$end ), unlist(coord.list)) ]

all$constitutive.score <-  unlist( constitutive.list )[ match( paste( all$chr, all$start, all$end ), unlist(coord.list)) ]

# replace NA values with "."
all$gene[ is.na( all$gene) ] <- "."
all$ensemblID[ is.na( all$ensemblID) ] <- "."
# replace missing transcripts with "."
all[ all$transcripts == "", ]$transcripts <- "."
all$constitutive.score <- signif(all$constitutive.score, digits = 2)


# prepare results
results$clusterID <- str_split_fixed(results$cluster, ":", 2)[,2]
results$N <- results$df + 1
sig <- subset(results, FDR < FDR_limit)
sig$clusterID <- str_split_fixed(sig$cluster, ":", 2)[,2]

# testing
#clu <- "clu_4879"

sig.annotated <- lapply(sig$clusterID, FUN = function(clu){
  cluster <- all[ all$clusterID == clu, ]
  chr <- unique( cluster$chr )[1] # this should always be one number
  start <- min( cluster$start )
  end <- max( cluster$end )
  # get most common gene name that is not "."
  gene <- names( sort( table( unique(cluster$gene) ), decreasing = TRUE ) )[1]
  ensemblID <- names( sort( table( unique(cluster$ensemblID) ), decreasing = TRUE ) )[1]
  annotation <- "annotated"
  if( any(grepl( "cryptic", cluster$verdict)) | any( grepl("skiptic", cluster$verdict)) ){
    annotation <- "cryptic"
  }  
  return( data.frame( clusterID = clu, chr = chr, start = start, end = end, gene = gene, ensemblID = ensemblID, annotation = annotation ) )
  })

sig.annotated <- do.call( what = rbind, args = sig.annotated)

sig.annotated$FDR  <- results$FDR[ match( sig.annotated$clusterID, results$clusterID)]
sig.annotated$FDR <- signif( sig.annotated$FDR, digits = 3)
sig.annotated$N  <- results$N[ match( sig.annotated$clusterID, results$clusterID)]
# fudge it by renaming
all.clusters <- sig.annotated
all.introns <- all

# write out the sig.annotated and all.clusters
cluster_results <- paste0(resultsFolder,"/per_cluster_results.tab")
intron_results <- paste0(resultsFolder, "/per_intron_results.tab")

write.table( sig.annotated, cluster_results, quote = FALSE, row.names = FALSE, sep = "\t" )
write.table( all, intron_results, quote = FALSE, row.names = FALSE, sep = "\t" )

# make png graphs of all clusters for visualising

#counts_file <- paste0( outFolder, "/", code, "_perind_numers.counts.gz" )


cat("Loading counts from",counts_file,"\n")
if (!file.exists(counts_file)) stop("File ",counts_file," does not exist")
counts=read.table(counts_file)

cat("Loading metadata from",groups_file,"\n")
if (!file.exists(groups_file)) stop("File ",groups_file," does not exist")
meta=read.table(groups_file, header=F, stringsAsFactors = F)
colnames(meta)=c("sample","group")

exon_file <- paste0(annotation_code, "_all_exons.txt.gz")

exon_table=if (!is.null( exon_file )) {
  cat("Loading exons from",exon_file,"\n")
  if (!file.exists(exon_file)) stop("File ",exon_file," does not exist")
  read_table(exon_file)
} else {
  cat("No exon_file provided.\n")
  NULL
}


counts=counts[,meta$sample]

# PCA of the counts matrix


meta$group=as.factor(meta$group)
group_names=levels(meta$group)

introns=leafcutter:::get_intron_meta(rownames(counts))
cluster_ids=paste(introns$chr,introns$clu,sep = ":")

plotFolder <- paste0(resultsFolder, "/plots")
if( ! dir.exists(plotFolder)){ dir.create(plotFolder)}


save( all.introns, all.clusters, counts, meta, exon_table, 
  file = paste0( resultsFolder, "/results.Rdata") )



quit()

# for all signiicant clusters, make a png of the visualisation
for( i in 1:nrow(sig.annotated) ){

  cluster_to_plot <- paste( sig.annotated$chr[i], sig.annotated$clusterID[i], sep = ":")
  cluster_gene <- sig.annotated$gene[i]
  if( cluster_gene == "."){ cluster_gene <- "Null"}

  out_png <- paste0( plotFolder, "/",sig.annotated$clusterID[i], ".png")

  png(out_png,width = 1024, height = 1024, units = "px", type = "cairo-png", pointsize = 40)
  y=t(counts[ cluster_ids==cluster_to_plot, ])
  make_differential_splicing_plot(y, meta$group, exons_table=exon_table)
  dev.off()
}


