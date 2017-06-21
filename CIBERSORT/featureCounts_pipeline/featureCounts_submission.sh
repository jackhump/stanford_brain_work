# featureCounts FPKM work
PIPELINE="/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/CIBERSORT/featureCounts_pipeline/"

# should go in separate support
CODE=Prudencio_CBL
STRAND=reverse
PAIRED=yes
SPECIES=human
OUTFOLDER="/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/CIBERSORT/${CODE}"
SUPPORT=${OUTFOLDER}/${CODE}_support.tab
COUNT_FEATURES=yes
CREATE_FPKMS=yes
SUBMISSION=no


sh ${PIPELINE}/featureCounts_pipeline.sh --CODE $CODE \
					 --STRAND $STRAND \
					 --PAIRED $PAIRED \
					 --SPECIES $SPECIES \
					 --OUTFOLDER $OUTFOLDER \
					 --SUPPORT $SUPPORT \
					 --COUNT_FEATURES $COUNT_FEATURES \
					 --CREATE_FPKMS $CREATE_FPKMS \
					 --SUBMISSION $SUBMISSION 




