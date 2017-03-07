# wrangle leafcutter output into something interpretable
library(data.table)
library(stringr)
library(dplyr)
library(optparse)

FDR_limit <- 0.05

# testing
setwd("/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/leafcutter/F210I_norm")
species <- "mouse"
outFolder <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/leafcutter/F210I_norm"


outFolder 
species


resultsFolder <- paste0(outFolder,"/results")

if( !dir.exists(resultsFolder)){
  dir.create(resultsFolder)
}

effectSizes <- fread(paste0(outFolder,"/leafcutter_ds_effect_sizes.txt"), stringsAsFactors = F)
effectSizesSplit <-  as.data.frame(str_split_fixed(effectSizes$intron, ":", 4), stringsAsFactors = F)
names(effectSizesSplit) <- c("chr","start","end","clusterID")
effectSizes <- cbind( effectSizes, effectSizesSplit)
effectSizes$cluster <- paste(effectSizesSplit$chr, effectSizesSplit$clusterID, sep = ":")

results <- fread(paste0(outFolder, "/leafcutter_ds_cluster_significance.txt"),stringsAsFactors = F)
results$FDR <- p.adjust( results$p, method = "fdr")

all <- merge(x = results, y = effectSizes, by = "cluster")
all <- all[ order(all$FDR),]

all <- subset( all, FDR <= FDR_limit )
all$start <- as.numeric(all$start)
all$end <- as.numeric(all$end)
# using AWK I've created lists of introns from each transcript in GENCODE comprehensive. I can match the introns in each cluster to these introns and work out if they are:
# 1) annotated: where the intron is found
# 1a) annotated as constitutive or alternate: how to define this? Present in all "basic" transcripts is probably safest
# 2) cryptic anchored: where either the start or end matches to an annotated intron
# 3) cryptic unanchored: where neither end matches to an annotated intron
# 4) unannotated alternate/skiptic - where both the start and the end are annotated splice sites but the combination has not been seen before
if( species == "mouse"){ intron.database.file <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/misc/gencode_mouse" }

#intron.database <- fread(intron.database.file)
#names(intron.database) <- c("chr","start","end","gene","transcript","intron.number", "tag")


# all.bed <- as.data.frame(select(all, chr, start, end, clusterID))
# all.bed$chr <- as.character(all.bed$chr)
# all.bed$start <- as.numeric( as.character(  all.bed$start ) ) # fucking R converts factors into their indexes so first convert to character
# all.bed$end <- as.numeric(  as.character( all.bed$end ) ) 

# all.bed <- bedr.sort.region(all.bed)

# 
# annotate_intron <- function( intron.start, intron.end, intron.strand, intron.database ){
  # first perform bedtools intersect between all and the intron database
  #all.bed.file <- paste0(outFolder, "/all.introns.bed")
  #write.table( all.bed, all.bed.file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t" )

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

 
  threeprime.cmd <- paste0( "bedtools intersect -a ", all.threeprime.file, " -b ", intron.database.file,"_threeprime.bed", " -wa -wb -loj -f 1" )
  
  threeprime_intersect <- fread(threeprime.cmd)

  fiveprime.cmd <- paste0( "bedtools intersect -a ", all.fiveprime.file, " -b ", intron.database.file,"_fiveprime.bed", " -wa -wb -loj -f 1" )
  
  fiveprime_intersect <- fread(fiveprime.cmd)

  # now I have two lists of splice site annotation
  # using these, can I annotate each intron in each cluster?
  cluster <- all[ all$clusterID == "clu_4879" , ]



verdict.list <- list()
coord.list <- list()
gene.list <- list()

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

  verdict <- c()
  coord <- c()
  gene <- c() 
  
  for( intron in 1:nrow(cluster) ){
    coord[intron] <- paste(cluster[intron]$chr,cluster[intron]$start, cluster[intron]$end )
    tgene <- names(sort(table( tprime[[intron]]$V8 ), decreasing = TRUE)[1])
    fgene <- names(sort(table( fprime[[intron]]$V8 ), decreasing = TRUE)[1])

    gene[intron] <- ifelse( tgene == ".",  no = fgene, yes = tgene )
    
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
  }
}

# verdict.list and coord.list are not equal lengths. why?

# need to match ideally
all$verdict <- unlist(verdict.list)[ match( paste( all$chr, all$start, all$end ), unlist(coord.list))]

all$gene <- unlist(gene.list)[ match( paste( all$chr, all$start, all$end ), unlist(coord.list))]

# coord.length <- lapply( coord.list, FUN = length )
# verdict.length <- lapply( verdict.list, FUN = length )

# #boxplots of the direction of change in each class of intron
# p <- ggplot( all, aes( x = verdict, y = deltapsi, group = cluster )) + geom_point() + geom_line()
# print(p);dev.off()


