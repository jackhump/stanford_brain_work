library(data.table)
library(stringr)
library(dplyr)
# lift barres mouse gene names to human 


mouse_human_genes <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/CIBERSORT/featureCounts_pipeline/mouse_human_genes.txt"

barres_mouse <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/CIBERSORT/Barres_mouse/Barres_mouse_all_FPKMs.txt"

out.table.pure <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/CIBERSORT/Barres_mouse/Barres_mouse_human_ensemblID_FPKMs_pure_cell_types.txt"

out.table.wholebrain <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/CIBERSORT/Barres_mouse/Barres_mouse_human_ensemblID_FPKMs_whole_brain.txt"

mh <- as.data.frame(fread(mouse_human_genes, na.strings = "."))
names(mh) <- c("mouse_ensemblID", "human_HUGO", "human_ensemblID")

d <- as.data.frame(fread(barres_mouse) )

d$ensemblID <- str_split_fixed( d$ensemblID, "\\.", 2)[,1]

# Do for EnsemblIDs

d$ensemblID <- mh$human_ensemblID[ match( d$ensemblID, mh$mouse_ensemblID)]

d <- d[ d$ensemblID != "", ]
# remove duplicate ensembl IDs. Table is sorted by median expression so should keep the highest expressed orthologue
d <- d[ !duplicated(d$ensemblID),]


# remove gene name and ID columns
d <- d[, c(-2)]

# split off whole brain - use in CIBERSORT as mixture file
wholebrain <- d[, grepl( "whole_brain", names(d)) ]

pure <- d[, !grepl( "whole_brain", names(d)) ]

write.table(pure, out.table.pure, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t" )
write.table(wholebrain, out.table.wholebrain,col.names = FALSE, row.names = TRUE, quote = FALSE, sep = "\t"  )

# do for HUGO IDs

d <- as.data.frame(fread(barres_mouse) )
d$ensemblID <- str_split_fixed( d$ensemblID, "\\.", 2)[,1]
d$HUGO <- mh$human_HUGO[ match( d$ensemblID, mh$mouse_ensemblID)]
d <- d[ d$HUGO != "", ]
d <- d[ !duplicated(d$HUGO),]
# reorder columns
d <- select(d, -external_gene_id, -ensemblID)
d <- select(d, HUGO, contains("_") )

wholebrain <- select(d, HUGO, contains("whole_brain"))
pure <- select(d, -contains("whole_brain"))

out.table.pure.hugo <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/CIBERSORT/Barres_mouse/Barres_mouse_human_hugo_FPKMs_pure_cell_types.txt"

out.table.wholebrain.hugo <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/CIBERSORT/Barres_mouse/Barres_mouse_human_hugo_FPKMs_whole_brain.txt"

write.table(pure, out.table.pure.hugo, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t" )
write.table(wholebrain, out.table.wholebrain.hugo,col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t"  )


