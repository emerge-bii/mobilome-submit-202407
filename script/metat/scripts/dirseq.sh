#!/bin/bash

########################################
# run dirseq for gene activity analysis
#########################################
### see dirseq details here: https://github.com/wwood/dirseq

Wkdir=~/mge/metat/nova-only.wkdir  # working directory
R1=nova-only/714E11014metaG_FD_JGI.1.fq.gz  # RNAseq paired-end read1 file


Sample=$(basename $R1 .1.fq.gz)
Output_dir=$Wkdir/$Sample

if [ -s $Output_dir/SS.dirseq.tsv ]; then
    echo "$Sample dirseq is already finished.."
else
    (cd $Output_dir && dirseq --measure-type count --bam final_bam/final.bam --gff annotate/combined_reference.gff > SS.dirseq.tsv)
fi
