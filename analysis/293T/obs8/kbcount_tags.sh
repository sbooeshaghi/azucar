#!/bin/bash

# FQ=../../../../data/293T/rep3/lane1/tags/ && BASE=../../../../reference/293T/tags && ./kbcount_tags.sh -x 10xv3 -o tags_out/ -i $BASE/index.idx -g $BASE/t2g.txt -f $FQ

usage () {
    echo "Usage: $0 [options]
    
    Options:
    -o, --output            output folder
    -i, --index             pseudoalignment index
    -g, --genemap           transcripts to genes map
    -x, --technology        single-cell tech
    -f, --fastqdir          folder containing fastqs
    "
    exit 1
}

while getopts ":o:i:w:g:x:f:" opt; do
    case $opt in
        o|--output)
            OUTDIR=$OPTARG
            ;;
        i|--index)
            INDEX=$OPTARG
            ;;
        g|--genemap)
            T2G=$OPTARG
            ;;
        x|--technology)
            TECH=$OPTARG
            ;;
        f|--fastqdir)
            FASTQDIR=$OPTARG
            ;;
        h)
            usage
            ;;
        \?)
            echo "Invalid argument"
            usage
            ;;
        :)
            echo "Add arguments"
            usage
            ;;
    esac
done

# check options        
if [ -z "$OUTDIR" -o -z "$INDEX" -o -z "$T2G" -o -z "$TECH" -o -z "$FASTQDIR" ]
then
    echo "Error"
    usage
fi


# begin workflow
kb count \
--workflow kite:10xFB \
--h5ad \
--filter \
-t 16 \
-m 16G \
-o $OUTDIR \
-i $INDEX \
-g $T2G \
-x $TECH \
$(paste -d" " \
  <(ls $FASTQDIR | awk -v p=$FASTQDIR '{print p$0}' | grep R1) \
  <(ls $FASTQDIR | awk -v p=$FASTQDIR '{print p$0}' | grep R2))
