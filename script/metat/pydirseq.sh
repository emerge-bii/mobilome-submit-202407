#!/bin/bash

### see installation instruction here: https://github.com/jiarong/pydirseq

Wkdir=~/mge/metat/nova-only.wkdir
R1=nova-only/714E11014metaG_FD_JGI.1.fq.gz
Sample=$(basename $R1 .1.fq.gz)
Output_dir=$Wkdir/$Sample

if [ -s $Output_dir/SS.final.tsv ]; then
    echo "$Sample pydirseq is already finished.."
else
    (cd $Output_dir && pydirseq run --bam final_bam/final.bam --gff annotate/combined_reference.gff all)
fi
