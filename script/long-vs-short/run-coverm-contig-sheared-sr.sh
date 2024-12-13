#!/bin/bash
#SBATCH --nodes=1 --ntasks=8 --time=12:00:00

set -e

READ1=contig-min1500bp-sheared/cmr2.MainAutochamber.201907_E_1_1to5.fasta
REF=mag/cmr6.MA.201907_E_1_1to5.fna
WKDIR=coverm-contig-sheared-sr
CPU=8
MEM=40g

BNAME=$(basename $REF .fna)
READ1=$(readlink -f $READ1)
REF=$(readlink -f $REF)

date
echo $BNAME contig
mkdir -p $WKDIR
cd $WKDIR
 
coverm contig -t $CPU -p minimap2-sr -m trimmed_mean --single $READ1 -r $REF --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 --min-covered-fraction 0.7 --bam-file-cache-directory bam_files --discard-unmapped > $BNAME.trimmed_mean.tsv
samtools calmd -b -@ $CPU --reference $REF bam_files/$(basename $REF).$(basename $READ1).bam > $BNAME.bam
samtools index $BNAME.bam -o $BNAME.bam.bai
date
