#!/usr/bin/env bash

# Set the parameters
OUTDIR=$SCRATCH/IMG_soil/selected_scaffolds
TAXON_ID_LIST=$SCRATCH/IMG_soil/empty.txt
MIN_LENGTH=2000
MIN_CIRCULAR_LENGTH=2000

# Define a progress bar function
function progress_bar {
    # $1 is current state, #2 is total
    _bar_length=60
    _progress=$(((${1}*100/${2}*100)/100))
    _done=$(((${_progress}*_bar_length)/100))
    _left=$((_bar_length-$_done))
    _done=$(printf "%${_done}s")
    _left=$(printf "%${_left}s")
    printf "\rProgress: ｜${_done// /▮}${_left// /⋅}｜ $1/$2 (${_progress}%%)"
}

# Define a function to query scaffolds in a given genome
get_scaffolds () {
    OUTDIR=$1
    TAXON_ID=$2
    MIN_LENGTH=$3
    MIN_CIRCULAR_LENGTH=$4
    OUTFILE=${1}/${2}.fna
    DBDIR=/global/dna/projectdirs/microbial/img_web_data_merfs/${TAXON_ID}/assembled/fna
    if [ -d "$DBDIR" ]
    then
        echo -n > $OUTFILE
        for SQLDB in $(ls $DBDIR)
        do
            sqlite3 ${DBDIR}/${SQLDB} \
            "SELECT '>' || scaffold_oid || x'0a' || fna FROM scaffold_fna" |
            ./filter_seq.py $MIN_LENGTH $MIN_CIRCULAR_LENGTH |
            sed "s/^>/>${TAXON_ID}\|/" >> \
            $OUTFILE
        done
    fi
}

# Start the counter and the progress bar
COUNTER=0
TOTAL=$(wc -l $TAXON_ID_LIST | awk '{print $1}')
progress_bar $COUNTER $TOTAL

# Create the output file and loop through the list of Taxon IDs
mkdir -p $OUTDIR
cat $TAXON_ID_LIST | while read TAXON_ID
do
    get_scaffolds $OUTDIR $TAXON_ID $MIN_LENGTH $MIN_CIRCULAR_LENGTH
    COUNTER=$(($COUNTER+1))
    progress_bar $COUNTER $TOTAL
done
echo
