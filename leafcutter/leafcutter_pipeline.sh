#!/bin/bash 
# leafcutter pipeline
# example code to create a support file
# find `pwd` -name [^s]*_FC*unique.bam |  awk -F '/' 'BEGIN{OFS="\t"} { split( $NF, a, "_" ); print $0, a[1]' }
# script exits if return value of a command is not zero
# for debugging prints out every line before executing it
set -x
set -euo pipefail


# for testing
CODE="Prudencio_FC_MND_C9"
OUTFOLDER="/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/leafcutter/${CODE}"
SUPPORT=${OUTFOLDER}/${CODE}_support.tab
SPECIES=human
LEAFCUTTER="/SAN/vyplab/HuRNASeq/leafcutter" # where you keep the leafcutter repo
STEP1=yes
STEP2=yes
STEP3=yes


until [ -z "$1" ]
do
    # use a case statement to test vars. we always test $1 and shift at the end of the for block.
    case $1 in
        --support)
            shift
            SUPPORT=$1;;
        --outFolder)
            shift
           	OUTFOLDER=$1;;
        --code) 
            shift 
            CODE=$1;; 
        --species) 
            shift 
            SPECIES=$1;;
        --leafcutter) 
            shift
            LEAFCUTTER=$1;;
        --step1)
			shift
			STEP1=$1;;
		--step2)
			shift
			STEP2=$1;;
        -* )
            stop "Unrecognized option: $1"
    esac
    shift
    if [ "$#" = "0" ]; then break; fi
    echo $1
done

if [ ! -e $SUPPORT ];then
	echo "$SUPPORT not found" 1>&2
	exit 1
fi


# the base of where I'll output everything
if [ ! -e $OUTFOLDER ];then
	mkdir -p $OUTFOLDER
fi

for FOLDER in ${OUTFOLDER}/cluster/submission ${OUTFOLDER}/cluster/output ${OUTFOLDER}/cluster/error; do 
	if [ ! -e $FOLDER ];then
		mkdir -p $FOLDER
	fi
done




# provide a support file which is the the sample name,  and the condition for leafcutter to compare differential splicing
JUNCTIONDIR=${OUTFOLDER}/junctions # where I'll put the junction files
if [ ! -e $JUNCTIONDIR ];then  
	mkdir -p $JUNCTIONDIR
fi

JUNCTIONLIST=${JUNCTIONDIR}/junction_list.txt # whatever you call your list of junction files
# remove old junction list
if [ -e $JUNCTIONLIST ];then
	rm $JUNCTIONLIST
fi

# step 1 : junction counting
function step1 {
STEP1_MASTER=${OUTFOLDER}/cluster/submission/step1_submission.sh

STEP1_MASTER_TABLE=${OUTFOLDER}/cluster/submission/step1_submission.tab
  if [ -e $STEP1_MASTER_TABLE ]; then rm $STEP1_MASTER_TABLE;fi
  echo "scripts" > $STEP1_MASTER_TABLE


# Convert bam to junction files in batches of four parallel jobs
for BAMFILE in `cut -f1 $SUPPORT `;do
	
	if [ ! -e $BAMFILE ];then
		echo $BAMFILE does not exist 
		exit 1
	fi
	
	BAMNAME=`basename $BAMFILE`
	SCRIPTNAME=`echo $BAMNAME | sed 's/_unique.bam//' ` 
    STEP1_SCRIPT=${OUTFOLDER}/cluster/submission/step1_${SCRIPTNAME}

    echo "
    cd ${LEAFCUTTER}/leafcutter

    sh ${LEAFCUTTER}/scripts/bam2junc.sh $BAMFILE ${JUNCTIONDIR}/${BAMNAME}.junc

    # check if successful
	# if not then try running again
	if [ ! -e ${JUNCTIONDIR}/${BAMNAME}.junc ];then
		sh ${LEAFCUTTER}/scripts/bam2junc.sh $BAMFILE ${JUNCTIONDIR}/${BAMNAME}.junc
	fi


    " > $STEP1_SCRIPT

    echo ${JUNCTIONDIR}/${BAMNAME}.junc >> $JUNCTIONLIST 
    echo $STEP1_SCRIPT >> $STEP1_MASTER_TABLE
done

NJOBS_STEP1=`wc -l $SUPPORT | awk '{print $1 + 1 }' `

echo "
#$ -S /bin/bash
#$ -l h_vmem=4G
#$ -l tmem=4G
#$ -l h_rt=24:00:00
#$ -pe smp 1
#$ -R y
#$ -o ${OUTFOLDER}/cluster/out
#$ -e ${OUTFOLDER}/cluster/error
#$ -N leafcutter_step1_${CODE}
#$ -wd ${OUTFOLDER}
#$ -t 2-${NJOBS_STEP1}
#$ -tc 20
echo \$HOSTNAME >&2
script=\`awk '{if (NR == '\$SGE_TASK_ID') print}' $STEP1_MASTER_TABLE\`
sh \$script
" > $STEP1_MASTER
}
# steps 2 and 3 

function step2 {
STEP2_3_SCRIPT=${OUTFOLDER}/cluster/submission/step2_3_submission.sh
echo "
#$ -S /bin/bash
#$ -l h_vmem=3G
#$ -l tmem=3G
#$ -l h_rt=24:00:00
#$ -pe smp 4
#$ -R y
#$ -o ${OUTFOLDER}/cluster/out
#$ -e ${OUTFOLDER}/cluster/error
#$ -N leafcutter_step2_3_${CODE}
#$ -wd ${LEAFCUTTER}/leafcutter
echo \$HOSTNAME >&2

# change working directory to leafcutter - all the scripts rely on this for some reason
cd ${LEAFCUTTER}/leafcutter

# intron clustering

# it's insisting on running within the LEAFCUTTER/leafcutter directory
python ${LEAFCUTTER}/clustering/leafcutter_cluster.py -j $JUNCTIONLIST -m 50 -o ${CODE} -l 500000
" > $STEP2_3_SCRIPT

# Differential intron excision analysis
# take the support file and only retain the file names to create a new support
DS_SUPPORT=${OUTFOLDER}/${CODE}_ds_support.tab
awk -F'/' '{print $NF}' $SUPPORT > $DS_SUPPORT

echo "
# run differential expression analysis
${LEAFCUTTER}/scripts/leafcutter_ds.R -i 4 --num_threads 4 ${CODE}_perind_numers.counts.gz $DS_SUPPORT
" >> $STEP2_3_SCRIPT

# Visualise the significant splice events
# the GENCODE exons are for hg19 - need to recreate for hg38

if [[ "$SPECIES" == "human" ]];then
	EXONLIST=${LEAFCUTTER}/leafcutter/data/gencode38_exons.txt.gz
	# create the hg38 exon list if need be
	if [ ! -e ${LEAFCUTTER}/leafcutter/data/gencode38_exons.txt.gz ];then
		# remove the first line
		zcat ${LEAFCUTTER}/leafcutter/data/gencode19_exons.txt.gz | awk 'NR > 1'  > ${LEAFCUTTER}/leafcutter/data/gencode19_temp 
		liftOver ${LEAFCUTTER}/leafcutter/data/gencode19_temp \
			/cluster/project8/vyp/vincent/Software/liftover/hg19ToHg38.over.chain.gz \
			${LEAFCUTTER}/leafcutter/data/gencode38_exons.txt \
			${LEAFCUTTER}/leafcutter/data/gencode38_exons.txt.unmapped
		# add the first line back in
		echo "chr start end strand gene_name" >  ${LEAFCUTTER}/leafcutter/data/gencode38_header
		cat ${LEAFCUTTER}/leafcutter/data/gencode38_header ${LEAFCUTTER}/leafcutter/data/gencode38_exons.txt | 
			gzip > ${LEAFCUTTER}/leafcutter/data/gencode38_exons.txt.gz
	fi

elif [[ "$SPECIES" == "mouse" ]];then
	EXONLIST=${LEAFCUTTER}/leafcutter/data/gencode_mm10_exons.txt.gz
	if [ ! -e $EXONLIST ];then
		# gencode mouse gtf
		GENCODE_MOUSE_GTF="/SAN/vyplab/HuRNASeq/GENCODE/gencode.vM12.annotation.gtf.gz"

		zcat $GENCODE_MOUSE_GTF |
		awk -F '\t' 'BEGIN{
				print "chr start end strand gene_name"
			} $3 == "exon" {
				split( $NF, a, "; " );
				split( a[4], b, "\"" );
				print $1,$4,$5,$7,b[2]
			}' | gzip > $EXONLIST

	fi

fi


echo "

# splicing plots
${LEAFCUTTER}/scripts/ds_plots.R \
	 -e ${EXONLIST} \
	 ${CODE}_perind_numers.counts.gz \
	 $DS_SUPPORT \
	 leafcutter_ds_cluster_significance.txt \
	 -f 0.05 \
	 -m 100 

# move all the results to the outfolder#
mv -t $OUTFOLDER leafcutter_ds_cluster_significance.txt leafcutter_ds_effect_sizes.txt ds_plots.pdf ${CODE}_perind_numers.counts.gz
# clean up
rm *.bam.junc.${CODE}.sorted.gz ${CODE}_pooled ${CODE}_refined ${CODE}_sortedlibs ${CODE}_perind.counts.gz

" >> $STEP2_3_SCRIPT

}

# STEP 3 - wrangle the results, annotate them

function step3 {
	WRANGLE_SCRIPT="/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/leafcutter/wrangle_results.R"

	STEP3_SCRIPT=${OUTFOLDER}/cluster/submission/step3_submission.sh

}




# SUBMITTING

if [[ "$STEP1" == "yes" ]]; then
	QHOLD="-hold_jid leafcutter_step1_${CODE}"

	step1
	qsub $STEP1_MASTER
fi

if [[ "$STEP2" == "yes" ]]; then
# as leafcutter insists on running in its own directory, I can't have multiple step2 jobs running at the same time. 
# Therefore any step2 job should also wait for any pre-existing step2 jobs before executing
	if [[ "$STEP1" == "no" ]];then
		# run step1 anyway to get the junction list updated
		step1
		# check that step1 successfully ran
		for JUNCTIONFILE in `cat $JUNCTIONLIST`;do
			if [ ! -e $JUNCTIONFILE ];then
				echo $JUNCTIONFILE does not exist
				exit 1
			fi
		done
		QHOLD=""
	fi
	
	function qlong { # prints full job ids
        qstat -xml | 
        tr "\n" " " | 
        sed "s#<job_list[^>]*>#\n#g" | 
        sed "s#<[^>]*>##g" | 
        sed 's/\(-[0-9]*\)T\([0-9]*:\)/\1\t\2/g' | 
        grep " " | 
        sed 's/\(20[0-9]*\)-\([0-9]*\)-\([0-9]*\)/\3\/\2\/\1/g' | 
        column -t
	}


	# finds all currently queueing or running step2 jobs
	QHOLDNEW=` qlong | 
	awk 'BEGIN{ORS=","} $3 ~ /leafcutter_step2_3/ {print $3}' | 
	sed 's/,$//' `

	echo $QHOLDNEW

	echo "qholdnew is $QHOLDNEW"


	# if queuing jobs are found then add them to QHOLD
	if [[ "$QHOLDNEW" != "" && "$QHOLD" != "" ]];then
		QHOLD="${QHOLD},${QHOLDNEW}"
	fi
	if [[ "$QHOLDNEW" != "" && "$QHOLD" == "" ]];then
		QHOLD="-hold_jid $QHOLDNEW"
	fi
	step2
	qsub $QHOLD $STEP2_3_SCRIPT
fi


