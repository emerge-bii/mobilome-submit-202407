#!/bin/bash

Wkdir=~/mge/metat/nova-only.wkdir
Contigdir=~/mge/contig-dram/jgi-jgi
R1=nova-only/714E11014metaG_FD_JGI.1.fq.gz
Cpu=8

Dname=$(dirname $R1)
Sample=$(basename $R1 .1.fq.gz)

P1=$Dname/$Sample.1.fq.gz
P2=$Dname/$Sample.2.fq.gz
Ref=$Contigdir/$Sample/scaffolds.fna
Gff=$Contigdir/$Sample/genes.gff
Output_dir=$Wkdir/$Sample
Memory=100  # in GB

date
mkdir -p $Wkdir

if [ -f $Output_dir/final_bam/final.bam ]; then
    echo "$Sample transcriptm is already finished.."
else
    rm -rf $Output_dir
    transcriptm count \
        -1 $P1 \
        -2 $P2 \
        --ref $Ref \
        --gff $Gff \
        --other-db ~/lab/jiarong/db/transcriptm/bowtie/rRNA_tRNA_tmRNA_PhiX \
        -o $Output_dir \
        -n $Cpu \
        -m $Memory
fi

date


