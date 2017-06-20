if [[ $1 =~ gz$ ]];then
	command=zcat
else
	command=cat
fi

# $1 is the file
# $2 is the prefix code
CODE=$2
# take a list of exons in a GTF and return the introns
echo finding introns in each transcript

$command $1 | awk -F '\t'  '
	BEGIN{ OFS = "\t" } 
	$3 == "gene" {
		split( $9, a, ";" )
		split( a[3], b, " ")
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
				#print $1, INTRONSTART, INTRONEND, GENE, GENEID, TRANSCRIPT, TRANSCRIPTTYPE, INTRONCOUNT,TAG
				INTRONLENGTH = INTRONEND - INTRONSTART
			}
			if ( $7 == "-" ){ # negative strand
				INTRONSTART=$5
				#print $1, INTRONSTART, INTRONEND, GENE, GENEID, TRANSCRIPT, INTRONCOUNT, TRANSCRIPTTYPE, TAG

			}
			INTRONLENGTH[TRANSCRIPT] = INTRONEND - INTRONSTART
		}
		INTRONCOUNT+=1
		if ( $7 == "+" ){
			INTRONSTART=$5 # begins on the first exon
		}
		if ( $7 == "-" ){
			INTRONEND=$4 # begins on the first exon
		}
	}' > $CODE"_all_introns.bed"

# each line is an intron in a transcript of a gene
cat ${CODE}"_all_introns.bed" | awk -F '/t' '
	BEGIN{ OFS = "\t" }


'
