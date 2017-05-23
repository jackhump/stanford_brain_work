INFILE=$1  # the STAR SJ.out.tab file
OUTFILE=$2  # the new Leafcutter format junctions

awk 'BEGIN{
    OFS="\t"}
    {
    start=$2-1;
    strand=$4;
    if( strand == 1)
      strand="+";
    else if( strand == 2)
      strand="-";
    else
      strand=".";
    
    counts=$7+$8;

    print( $1, start, $3, ".", counts, strand ) 
}' $INFILE > $OUTFILE