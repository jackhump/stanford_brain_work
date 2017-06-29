#### Jack Humphrey 2017
#### wrangle leafcutter differential splicing output into something interpretable

library(data.table)
library(leafcutter)
library(stringr)
library(dplyr)
library(optparse)
library(readr)

FDR_limit <- 0.05

mode <- "differential_splicing"

#options(echo=TRUE)


### for testing
# outFolder <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/leafcutter/ENCODE/TARDBP_K562"   
# species <- "human"                                                
# groups_file <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/leafcutter/ENCODE_TARDBP_K562/TARDBP_K562_ds_support.tab"
# counts_file <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/leafcutter/ENCODE/TARDBP_K562/TARDBP_K562_perind_numers.counts.gz"
# annotation_code <- "/SAN/vyplab/HuRNASeq/leafcutter/leafcutter/data/gencode_hg38"
# code <- "TARDBP_K562"



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


opt <- parse_args(
  OptionParser(
    option_list=list(
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

# script depends on a set of splice junction annotation files
# check that the annotation files have been created

for( file in c("_all_introns.bed","_threeprime.bed","_fiveprime.bed", "_all_exons.txt.gz" )){
  if( ! file.exists( paste0( annotation_code, file ))){
    stop( paste0( annotation_code, file , " does not exist!"))
  }
}

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

print( "BedTools intersect junctions with list of known splice sites")

# intersect with bedtools to find the annotations of each splice site
threeprime.cmd <- paste0( "bedtools intersect -a ", all.threeprime.file, " -b ", annotation_code,"_threeprime.bed", " -wa -wb -loj -f 1" )

threeprime_intersect <- fread(threeprime.cmd)

fiveprime.cmd <- paste0( "bedtools intersect -a ", all.fiveprime.file, " -b ", annotation_code,"_fiveprime.bed", " -wa -wb -loj -f 1" )

fiveprime_intersect <- fread(fiveprime.cmd)

# remove temporary files
rm.cmd <- paste("rm ", all.fiveprime.file, all.threeprime.file) 
system(rm.cmd)

# now I have two lists of splice site annotation
# for testing
#cluster <- all[ all$clusterID == "clu_4879" , ]

print("Annotating junctions")

verdict.list <- list()
coord.list <- list()
gene.list <- list()
ensemblID.list <- list()
transcripts.list <- list()
constitutive.list <- list()
classification.list <- list()


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
  # do the same for EnsemblID
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

  # predicting the event type from the shape of the junctions
  # easy start - cassette exons
  #print(clu)

  if( nrow(cluster) != 3){ 
    classification.list[[clu]] <- "." 
    next
  }else{
    classification.list[[clu]] <- "."

    tab <- select(cluster, start, end)
    
    # the junctions are sorted by start and end coordinates

    # check for the presence of a junction that spans the entire length of the cluster
    if( !any(  which( tab$start == min(tab$start) ) %in% which( tab$end == max(tab$end) )  ) ){
      classification.list[[clu]] <- "."
      next
    }

    # therefore for a cassette exon arrangement the longest junction always comes second 
    if( which( tab$start ==  min(tab$start) & tab$end == max(tab$end ) ) != 2 ){
     classification.list[[clu]] <- "." 
     next 
    }

    # now we know that junction 2 is the parent, junction 1 is the left most child and junction 3 is the right most
    # check that the end of junction 1 comes before the start of junction 3

    if( tab[1,"end"] > tab[3,"start"] ){
      classification.list[[clu]] <- "."
      next
    }

    # double check the starts and ends
    if( tab[1, "start"] != tab[2,"start"] | tab[3,"end"] != tab[2,"end"] ){
      classification.list[[clu]] <- "."
      next
    }

    # work out direction of change
    if( cluster[1, "deltapsi"] > 0 & cluster[3, "deltapsi"] > 0 & cluster[2,"deltapsi"] < 0){
      classification.list[[clu]] <- "cassette exon - increased"
    }
    if( cluster[1, "deltapsi"] < 0 & cluster[3, "deltapsi"] < 0 & cluster[2,"deltapsi"] > 0){
      classification.list[[clu]] <- "cassette exon - decreased"
    }

    # work out annotation status
    if( all( verdict.list[[clu]] == "annotated") ){
      classification.list[[clu]] <- paste0( classification.list[[clu]], " - annotated")
    }

    if( verdict.list[[clu]][2] == "annotated" & verdict.list[[clu]][1] != "annotated" & verdict.list[[clu]][3] != "annotated"  ){
      classification.list[[clu]] <- paste0( classification.list[[clu]], " - cryptic")
    }

    if( verdict.list[[clu]][2] != "annotated" & verdict.list[[clu]][1] == "annotated" & verdict.list[[clu]][3] == "annotated"  ){
      classification.list[[clu]] <- paste0( classification.list[[clu]], " - skiptic")
    }


  # print(n_parent)
  # print(children_list)

  }
  
  # print(clu)
  # # print(n_parent)
  # # print(children_list)
  # print(classification.list[[clu]])

}

print("Preparing results")

# match all the lists together
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
  return( 
    data.frame( 
      clusterID = clu, 
      chr = chr, 
      start = start, 
      end = end, 
      gene = gene, 
      ensemblID = ensemblID, 
      annotation = annotation ) )
  })
sig.annotated <- do.call( what = rbind, args = sig.annotated)

sig.annotated$FDR  <- results$FDR[ match( sig.annotated$clusterID, results$clusterID)]
sig.annotated$FDR <- signif( sig.annotated$FDR, digits = 3)
sig.annotated$N  <- results$N[ match( sig.annotated$clusterID, results$clusterID)]

# add classification 
sig.annotated$verdict <- unlist(classification.list)[ match(sig.annotated$clusterID, names(classification.list))]

# fudge it by renaming
all.clusters <- sig.annotated
all.introns <- all

# write out the sig.annotated and all.clusters
cluster_results <- paste0(resultsFolder,"/per_cluster_results.tab")
intron_results <- paste0(resultsFolder, "/per_intron_results.tab")

write.table( sig.annotated, cluster_results, quote = FALSE, row.names = FALSE, sep = "\t" )
write.table( all, intron_results, quote = FALSE, row.names = FALSE, sep = "\t" )

counts <- counts[,meta$sample]

# PCA of the counts matrix

meta$group <- as.factor(meta$group)
group_names <- levels(meta$group)

# introns=leafcutter:::get_intron_meta(rownames(counts))
# cluster_ids=paste(introns$chr,introns$clu,sep = ":")

make_pca <- function(counts,meta){
  dev <- apply( counts, MAR = 1, FUN = sd )
  # remove rows with 0 variance
  counts <- counts[ dev != 0, ]
  pca <- prcomp( t(counts), scale = TRUE )
  importance <- signif( summary(pca)$importance[2,], digits = 2) * 100
  pca <- as.data.frame(pca$x)
  #names(pca) <- paste0( names(pca), " (", signif( importance, digits =  2) * 100, "%)" )
  pca$groups <- meta$group[ match( rownames(pca), meta$sample )]
  return(list( pca, importance) )
}

# sort out clusters table
# use on all.clustersv
fix_clusters <- function(clusters){
  clusters$FDR <- signif( clusters$FDR, digits = 3)
  clusters$coord <- paste0( clusters$chr, ":", clusters$start, "-", clusters$end)
  clusters <- clusters[ order(clusters$FDR, decreasing = FALSE),]
  # removed ensemblID - this could be an option?
  #clusters <- select( clusters, clusterID, N, coord, gene, annotation, FDR, verdict)
  clusters <- select( 
    clusters, 
    clusterID, 
    N, 
    coord, 
    gene, 
    annotation, 
    FDR)
  clusters$gene <- paste0("<i>",clusters$gene,"</i>")
return(clusters)
}

# use on all.introns
fix_introns <- function(introns){
  introns <- select(introns, 
    clusterID, 
    gene, 
    ensemblID, 
    chr, 
    start, 
    end, 
    verdict, 
    deltapsi, 
    #constitutive.score, 
    transcripts)
  #introns$constitutive.score <- signif(introns$constitutive.score, digits = 3)
  introns$deltapsi<- round(introns$deltapsi, digits = 3)
return(introns)
}


cluster_summary <- function(clusters){
  summary <- data.frame( 
              Results = c(
                paste0("Number of differentially spliced clusters at FDR = 0.05 ") , 
                        "Fully annotated",
                        "Contain unannotated junctions"),
                n = c( nrow(clusters),
                       nrow( clusters[ clusters$annotation == "annotated", ]),
                       nrow( clusters[ clusters$annotation == "cryptic", ]) 
                       ) 
                )
  return(summary)
}

intron_summary <- function(all.introns){
    summary <- data.frame( 
                Results = c( "Number of fully annotated junctions",
                             "Number of junctions with cryptic 5' splice site", 
                              "Number of junctions with cryptic 3' splice site",  
                              "Number of junctions with two cryptic splice sites",
                              "Number of novel junctions that connect two annotated splice sites"),

                  n = c( nrow(all.introns[ all.introns$verdict == "annotated",]),
                         nrow(all.introns[ all.introns$verdict == "cryptic_fiveprime",]),
                         nrow(all.introns[ all.introns$verdict == "cryptic_threeprime",]),
                         nrow(all.introns[ all.introns$verdict == "cryptic_unanchored",]),
                         nrow(all.introns[ all.introns$verdict == "skiptic",])  
                  )
                )
    return( summary )
}

# create all the objects for visualisation
pca <- make_pca(counts, meta = meta)
clusters <- fix_clusters(all.clusters)
introns <- fix_introns(all.introns)
intron_summary <- intron_summary(all.introns)
cluster_summary <- cluster_summary(all.clusters) 
introns_to_plot <- leafcutter:::get_intron_meta(rownames(counts))
cluster_ids <- introns_to_plot$clu 

# to speed up the visualisation interface it would be useful to make all the plots in advance.


# plotFolder <- paste0(resultsFolder, "/plots")
# if( ! dir.exists(plotFolder)){ dir.create(plotFolder)}

# save all the objects needed by Leafcutter viz into single Rdata file
# include the mode variable 

save( introns, 
      clusters, 
      counts, 
      meta, 
      exon_table, 
      species, 
      pca, 
      intron_summary, 
      cluster_summary, 
      introns_to_plot,
      cluster_ids,
      mode,
      file = paste0( resultsFolder, "/results.Rdata")
)


quit()

# for all signiicant clusters, make a png of the visualisation
# for( i in 1:nrow(sig.annotated) ){

#   cluster_to_plot <- paste( sig.annotated$chr[i], sig.annotated$clusterID[i], sep = ":")
#   cluster_gene <- sig.annotated$gene[i]
#   if( cluster_gene == "."){ cluster_gene <- "Null"}

#   out_png <- paste0( plotFolder, "/",sig.annotated$clusterID[i], ".png")

#   png(out_png,width = 1024, height = 1024, units = "px", type = "cairo-png", pointsize = 40)
#   y=t(counts[ cluster_ids==cluster_to_plot, ])
#   make_differential_splicing_plot(y, meta$group, exons_table=exon_table)
#   dev.off()
# }


