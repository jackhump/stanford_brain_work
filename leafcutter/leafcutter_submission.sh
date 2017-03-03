

leafcutter=/SAN/vyplab/HuRNASeq/leafcutter # where you keep the leafcutter repo
leafcutterScript=/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/leafcutter/leafcutter_submission.sh
code=F210I_norm
species=mouse
outFolder=/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/brain_work_stanford/leafcutter/$code
support=${outFolder}/F210I_norm_support.tab
step1=yes
step2=yes

sh -x $leafcutterScript \
	--step1 $step1 \
	--step2 $step2 \
	--code $code \
	--species $species \
	--outFolder $outFolder \
	--support $support \
	--leafcutter $leafcutter
