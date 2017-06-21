# Jack Humphrey 2017
# pipeline to use featureCounts from the Rsubread package to count overlapping fragments for each gene in the GENCODE GTF.
# second step then collates all the samples to produce an FPKM/RPKM table using the DESeq2 sizeFactors to normalise with.

set -x
set -euo pipefail

# featureCounts FPKM work 
R=/share/apps/R/bin/R 
FEATURECOUNTS="/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/CIBERSORT/featureCounts_pipeline/featureCounts_script.R" 
FPKMCREATOR="/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/CIBERSORT/featureCounts_pipeline/featureCounts_prepare.R" 
 
for i in $FEATURECOUNTS $FPKMCREATOR;do
	if [ ! -e $i ];then
		echo $file does not exist
		exit 1
	fi
done

# case statement
until [ -z "$1" ]
do
    # use a case statement to test vars. we always test $1 and shift at the end of the for block.
    case $1 in
        --SUPPORT)
            shift
            SUPPORT=$1;;
        --OUTFOLDER)
            shift
           	OUTFOLDER=$1;;
        --CODE) 
            shift 
            CODE=$1;; 
        --SPECIES) 
            shift 
            SPECIES=$1;;
        --STRAND) 
            shift
            STRAND=$1;;
        --PAIRED)
			shift
			PAIRED=$1;;
		--COUNT_FEATURES)
			shift
			COUNT_FEATURES=$1;;
		--CREATE_FPKMS)
			shift
			CREATE_FPKMS=$1;;
		--SUBMISSION)
			shift
			SUBMISSION=$1;;
        -* )
            stop "Unrecognized option: $1"
    esac
    shift
    if [ "$#" = "0" ]; then break; fi
    echo $1
done



if [ ! -e $SUPPORT ];then
	echo $SUPPORT does not exist
	exit 1 
fi


if [[ "$SPECIES" == "human" ]];then
	GTF=/SAN/vyplab/HuRNASeq/GENCODE/gencode.v25.annotation.gtf
elif [[ "$SPECIES" == "mouse" ]];then
	GTF=/SAN/vyplab/HuRNASeq/GENCODE/gencode.vM12.annotation.gtf
fi

SUBMISSIONFOLDER=${OUTFOLDER}
for folder in submission error out;do
	mkdir -p ${OUTFOLDER}/cluster/$folder
done

COUNTSFOLDER=${OUTFOLDER}/counts
if [ ! -e $COUNTSFOLDER ]; then
	mkdir $COUNTSFOLDER
fi

COUNTSLIST=${OUTFOLDER}/${CODE}_counts_list.txt

function featureCounts {

	if [ -e $COUNTSLIST ];then
		rm $COUNTSLIST
	fi

	SUBMISSIONLIST=${OUTFOLDER}/cluster/submission/feature_counts_submissions.txt
	if [ -e $SUBMISSIONLIST ];then
		rm $SUBMISSIONLIST
	fi

	echo "scripts" > $SUBMISSIONLIST


	cat $SUPPORT | while read BAM SAMPLE;do

		OUTFILE=${COUNTSFOLDER}/${SAMPLE}_featureCounts.txt

		echo $OUTFILE >> $COUNTSLIST

		JOBSCRIPT=${OUTFOLDER}/cluster/submission/${SAMPLE}_featureCounts.sh
		echo $JOBSCRIPT >> $SUBMISSIONLIST

		echo "
	${R}script $FEATURECOUNTS --bamFile $BAM \
						   --paired $PAIRED \
						   --countStrand $STRAND \
						   --GTF $GTF \
						   --outFile $OUTFILE
	" > $JOBSCRIPT

	done

	NJOBS=`wc -l $SUPPORT | awk '{print $1 + 1 }' `

	MASTERSCRIPT=${OUTFOLDER}/cluster/submission/feature_counts_master.sh

	echo "
#$ -S /bin/bash
#$ -l h_vmem=4G
#$ -l tmem=4G
#$ -l h_rt=24:00:00
#$ -pe smp 1
#$ -R y
#$ -o ${OUTFOLDER}/cluster/out
#$ -e ${OUTFOLDER}/cluster/error
#$ -N featureCounts_${CODE}
#$ -wd ${OUTFOLDER}
#$ -t 2-${NJOBS}
#$ -tc 20
echo \$HOSTNAME >&2
script=\`awk '{if (NR == '\$SGE_TASK_ID') print}' $SUBMISSIONLIST\`
sh \$script
	" > $MASTERSCRIPT
}

function createFPKM {
	FPKMSCRIPT=${OUTFOLDER}/cluster/submission/createFPKM.sh

	echo "
#$ -S /bin/bash
#$ -l h_vmem=4G
#$ -l tmem=4G
#$ -l h_rt=24:00:00
#$ -pe smp 1
#$ -R y
#$ -o ${OUTFOLDER}/cluster/out
#$ -e ${OUTFOLDER}/cluster/error
#$ -N FPKM_create_${CODE}
#$ -wd ${OUTFOLDER}
	
${R}script $FPKMCREATOR --support.frame $SUPPORT --featureCounts.list $COUNTSLIST --outFolder $OUTFOLDER --GTF $GTF --code $CODE
	" > $FPKMSCRIPT


}

# submission
QHOLD=""
if [[ "$COUNT_FEATURES" == "yes" ]];then
	featureCounts

	QHOLD="-hold_jid featureCounts_${CODE}"

	if [[ "$SUBMISSION" == "yes" ]];then
		qsub $MASTERSCRIPT
	fi
fi

if [[ "$CREATE_FPKMS" == "yes" ]];then
	createFPKM

	if [[ "$SUBMISSION" == "yes" ]];then
		qsub $QHOLD $FPKMSCRIPT
	fi
fi




