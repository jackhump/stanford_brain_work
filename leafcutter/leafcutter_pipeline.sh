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
REGTOOLS="/SAN/vyplab/HuRNASeq/regtools/build/regtools"
STEP1=yes
STEP2=yes
STEP3=yes

# experimental speed up of step1
STEP1_REGTOOLS=yes

R=/share/apps/R/bin/R
STAR_TO_LEAFCUTTER="/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/leafcutter/star_to_leafcutter.sh"
WRANGLE_SCRIPT="/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/leafcutter/wrangle_results.R"




until [ -z "$1" ]
do
    # use a case statement to test vars. we always test $1 and shift at the end of the for block.
    case $1 in
    	--submit)
			shift
			submit=$1;;
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
		--step3)
			shift
			STEP3=$1;;
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
echo "outFolder is"
echo $OUTFOLDER
# the base of where I'll output everything
if [ ! -e $OUTFOLDER ];then
	mkdir -p $OUTFOLDER
fi

for FOLDER in ${OUTFOLDER}/cluster/submission ${OUTFOLDER}/cluster/output ${OUTFOLDER}/cluster/error; do 
	if [ ! -e $FOLDER ];then
		mkdir -p $FOLDER
	fi
done


if [[ "$SPECIES" == "human" ]];then
	ANNOTATION_CODE="/SAN/vyplab/HuRNASeq/leafcutter/leafcutter/data/gencode_hg38"
elif [[ "$SPECIES" == "mouse" ]];then
	ANNOTATION_CODE="/SAN/vyplab/HuRNASeq/leafcutter/leafcutter/data/gencode_mm10"
fi





# provide a support file which is the the sample name,  and the condition for leafcutter to compare differential splicing
JUNCTIONDIR=${OUTFOLDER}/junctions # where I'll put the junction files
if [ ! -e $JUNCTIONDIR ];then  
	mkdir -p $JUNCTIONDIR
fi

if [[ "$STEP1_REGTOOLS" == "yes" ]];then
    	JUNCTIONDIR=${JUNCTIONDIR}_regtools
    	if [ ! -e $JUNCTIONDIR ];then
    		mkdir $JUNCTIONDIR
    	fi
fi

JUNCTIONLIST=${JUNCTIONDIR}/junction_list.txt # whatever you call your list of junction files
# remove old junction list
if [ -e $JUNCTIONLIST ];then
	rm $JUNCTIONLIST
fi

# Differential intron excision analysis
# take the support file and only retain the file names to create a new support
DS_SUPPORT=${OUTFOLDER}/${CODE}_ds_support.tab
awk -F'/' '{print $NF}' $SUPPORT > $DS_SUPPORT


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

    if [[ $STEP1_REGTOOLS == "yes" ]];then

    	echo "
$REGTOOLS  junctions extract -a 8 -i 50 -I 500000 -o ${JUNCTIONDIR}/${BAMNAME}.junc $BAMFILE 
    	" > $STEP1_SCRIPT 
    fi

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



# step 2

function step2 {
STEP2_SCRIPT=${OUTFOLDER}/cluster/submission/step2_submission.sh
echo "
#$ -S /bin/bash
#$ -l h_vmem=3G
#$ -l tmem=3G
#$ -l h_rt=24:00:00
#$ -pe smp 4
#$ -R y
#$ -o ${OUTFOLDER}/cluster/out
#$ -e ${OUTFOLDER}/cluster/error
#$ -N leafcutter_step2_${CODE}
#$ -wd ${LEAFCUTTER}/leafcutter
echo \$HOSTNAME >&2

# change working directory to leafcutter - all the scripts rely on this for some reason
cd ${LEAFCUTTER}/leafcutter


ls -lhtr 
echo intron clustering

# -j JUNCFILES, --juncfiles=JUNCFILES
#                       text file with all junction files to be processed
# -o OUTPREFIX, --outprefix=OUTPREFIX
#                       output prefix (default leafcutter)
# -q, --quiet           don't print status messages to stdout
# -r RUNDIR, --rundir=RUNDIR
#                       write to directory (default ./)
# -l MAXINTRONLEN, --maxintronlen=MAXINTRONLEN
#                       maximum intron length in bp (default 100,000bp)
# -m MINCLUREADS, --minclureads=MINCLUREADS
#                       minimum reads in a cluster (default 30 reads)
# -p MINCLURATIO, --mincluratio=MINCLURATIO
#                       minimum fraction of reads in a cluster that support a
#                       junction (default 0.001)



# it's insisting on running within the LEAFCUTTER/leafcutter directory

python ${LEAFCUTTER}/clustering/leafcutter_cluster.py -j $JUNCTIONLIST -m 50 -o ${CODE} -l 500000
" > $STEP2_SCRIPT



echo "
#-i 1  # min samples per intron 
#	    -g 1 # min samples per group
#	    -c 1 # min read coverage that the min samples must have 
	    #--output_prefix ${CODE}_ds \
#--exon_file ${ANNOTATION_CODE}_all_exons.txt.gz \

echo differential expression analysis
${LEAFCUTTER}/scripts/leafcutter_ds.R \
		--output_prefix ${CODE}_ds \
		--num_threads 4 \
	    --min_samples_per_intron 1 \
	    --min_samples_per_group 1 \
	    --min_coverage 5 \
		${CODE}_perind_numers.counts.gz \
		$DS_SUPPORT

ls -lhtr 
" >> $STEP2_SCRIPT

# Visualise the significant splice events
# the GENCODE exons are for hg19 - need to recreate for hg38


echo "

echo splicing plots
${LEAFCUTTER}/scripts/ds_plots.R \
	 --output ${CODE}_ds_plots.pdf \
	 -e ${ANNOTATION_CODE}_all_exons.txt.gz \
	 ${CODE}_perind_numers.counts.gz \
	 $DS_SUPPORT \
	 ${CODE}_ds_cluster_significance.txt \
	 -f 0.05 \
	 -m 10 

ls -lhtr 

# move all the results to the outfolder#
mv -t $OUTFOLDER ${CODE}_ds_cluster_significance.txt ${CODE}_refined ${CODE}_ds_effect_sizes.txt ${CODE}_ds_plots.pdf ${CODE}_perind_numers.counts.gz
# clean up
rm *.bam.junc.${CODE}.sorted.gz ${CODE}_pooled  ${CODE}_sortedlibs ${CODE}_perind.counts.gz

" >> $STEP2_SCRIPT

}

# STEP 3 - wrangle the results, annotate them

function step3 {

	STEP3_SCRIPT=${OUTFOLDER}/cluster/submission/step3_submission.sh

	echo "
#$ -S /bin/bash
#$ -l h_vmem=4G
#$ -l tmem=4G
#$ -l h_rt=24:00:00
#$ -pe smp 1
#$ -R y
#$ -o ${OUTFOLDER}/cluster/out
#$ -e ${OUTFOLDER}/cluster/error
#$ -N leafcutter_step3_${CODE}
#$ -wd ${LEAFCUTTER}/leafcutter
echo \$HOSTNAME >&2

${R}script $WRANGLE_SCRIPT --outFolder $OUTFOLDER \
						   --species $SPECIES \
						   --groups_file $DS_SUPPORT \
						   --counts_file ${OUTFOLDER}/${CODE}_perind_numers.counts.gz \
						   --annotation_code $ANNOTATION_CODE \
						   --code $CODE
	" > $STEP3_SCRIPT

}




# SUBMITTING


QHOLD=""
if [[ "$STEP1" == "yes" ]]; then
	QHOLD="-hold_jid leafcutter_step1_${CODE}"

	step1
	if [[ "$submit" == "yes" ]];then
		qsub $STEP1_MASTER
	fi
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
	awk 'BEGIN{ORS=","} $3 ~ /leafcutter_step2/ {print $3}' | 
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
	if [[ "$submit" == "yes" ]];then
		qsub $QHOLD $STEP2_SCRIPT
	fi
	QHOLD="-hold_jid leafcutter_step2_${CODE}"
fi


if [[ "$STEP3" == "yes" ]];then
	step3
	if [[ "$submit" == "yes" ]];then
		qsub $QHOLD $STEP3_SCRIPT
	fi
fi

