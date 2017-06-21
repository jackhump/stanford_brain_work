if [[ $1 == "" || $1 == "--help" ]];then
	echo "
intronfinder.sh
-------------

finds introns in a GTF file

usage 
sh intron_finder.sh <GTF> <outfile_prefix> <species>

GTF can be gzipped
outfile_prefix can be file name or full path
species is either human or mouse
"
exit 0 
fi

if [[ $1 =~ gz$ ]];then
	command=zcat
else
	command=cat
fi

# $1 is the file
# $2 is the prefix code
CODE=$2
# $3 is the species , either mouse or human
SPECIES=$3
if [[ "$SPECIES" == "mouse" ]];then
	SPECIES=3
elif [[ "$SPECIES" == "human" ]];then
	SPECIES=4
else
	echo "\$3 must be either mouse or human"
	exit 1
fi
# take a list of exons in a GTF and return the introns
echo finding introns

$command $1 | awk -F '\t' \
		-v species=$SPECIES  '
	BEGIN{ OFS = "\t" } 
	$3 == "gene" {
		split( $9, a, ";" )
		split( a[species], b, " ")
		GENE=b[2]
		split( a[1], b, " ")
		GENEID=b[2] 
		#print GENE > "/dev/stderr" # print to stderr
	}

	$3 == "transcript" {
		split( $9, a, ";" )
		split( a[2], b, " ")
		TRANSCRIPT=b[2] # this is the transcript name
		split( a[10], c, " ")
		TAG=c[2]
		if ( TAG == ""){
			TAG="NA"
		}
		split( a[3], d, " ")
		TRANSCRIPTTYPE=d[2]
		INTRONCOUNT=0
		INTRONSTART=0; INTRONEND=0 # reset the coordinates 
	}
	$3 == "exon" {
		if ( INTRONCOUNT != 0 ){ # cannot happen on the first exon
			if ( $7 == "+" ){ positive strand
				INTRONEND=$4
				print $1, INTRONSTART, INTRONEND, GENE, GENEID, TRANSCRIPT, TRANSCRIPTTYPE, INTRONCOUNT,TAG
			}
			if ( $7 == "-" ){ # negative strand
				INTRONSTART=$5
				print $1, INTRONSTART, INTRONEND, GENE, GENEID, TRANSCRIPT, INTRONCOUNT, TRANSCRIPTTYPE, TAG
			}
		}
		INTRONCOUNT+=1
		if ( $7 == "+" ){
			INTRONSTART=$5 # begins on the first exon
		}
		if ( $7 == "-" ){
			INTRONEND=$4 # begins on the first exon
		}
	}' > $CODE"_all_introns.bed"


# Take a list of exons and return all the 5' splice sites to one file and all the 3' splice sites to another.
echo finding splice sites

$command $1 | awk -F '\t' \
	-v CODE=$2  \
	-v species=$SPECIES '
	BEGIN{ OFS = "\t" } 
	$3 == "gene" {
		split( $9, a, ";" )
		split( a[species], b, " ")
		GENE=b[2]
		split( a[1], b, " ")
		GENEID=b[2]  
	}

	$3 == "transcript" {
		split( $9, a, ";" )
		split( a[2], b, " ")
		TRANSCRIPT=b[2] # this is the transcript name
		split( a[10], c, " ")
		TAG=c[2] #
		if ( TAG == ""){
			TAG="NA"
		}
		split( a[3], d, " ")
		TRANSCRIPTTYPE=d[2]
		INTRONCOUNT=0
		INTRONSTART=0; INTRONEND=0 # reset the coordinates 
	}
	$3 == "exon" {
		if ( INTRONCOUNT != 0 ){ # cannot happen on the first exon
			if ( $7 == "+" ){ positive strand
				INTRONEND=$4				
			}
			if ( $7 == "-" ){ # negative strand
				INTRONSTART=$5
			}
				print $1, INTRONSTART, INTRONSTART + 1, GENE, GENEID, TRANSCRIPT, INTRONCOUNT, TRANSCRIPTTYPE, TAG > CODE"_fiveprime.bed"
				print $1, INTRONEND, INTRONEND + 1, GENE, GENEID, TRANSCRIPT, INTRONCOUNT, TRANSCRIPTTYPE, TAG > CODE"_threeprime.bed"
		}
		INTRONCOUNT+=1
		if ( $7 == "+" ){
			INTRONSTART=$5 # begins on the first exon
		}
		if ( $7 == "-" ){
			INTRONEND=$4 # begins on the first exon
		}
	}' 

echo creating exon list
# TODO: also create the exon lists for the graphing step.
$command $1 | awk -F '\t'   \
		-v species=$SPECIES ' # exons have the gene name in a different column to the gene entries. Fuck GTF files so hard.
BEGIN{
	print "chr start end strand gene_name"
			} $3 == "exon" {
				split( $NF, a, "; " );
				split( a[species + 1], b, "\"" );
				print $1,$4,$5,$7,b[2]
			}' | gzip > ${CODE}_all_exons.txt.gz
