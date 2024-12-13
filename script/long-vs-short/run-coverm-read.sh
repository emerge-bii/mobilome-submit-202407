#!/bin/bash
#SBATCH --nodes=1 --ntasks=8 --time=12:00:00

set -e

READ1=read-trimmed/MainAutochamber.201907_E_1_1to5.1.trim_paired.fq.gz
READ2=read-trimmed/MainAutochamber.201907_E_1_1to5.2.trim_paired.fq.gz
REF=mag/cmr6.MA.201907_E_1_1to5.fna
WKDIR=coverm-read
CPU=8
MEM=40g

BNAME=$(basename $REF .fna)
READ1=$(readlink -f $READ1)
READ2=$(readlink -f $READ2)
REF=$(readlink -f $REF)

date
echo $BNAME read
mkdir -p $WKDIR
cd $WKDIR
coverm contig -t $CPU -p minimap2-sr -m trimmed_mean --coupled $READ1 $READ2 -r $REF --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 --min-covered-fraction 0.7 --bam-file-cache-directory bam_files --discard-unmapped > $BNAME.trimmed_mean.tsv

samtools calmd -b -@ $CPU --reference $REF bam_files/$(basename $REF).$(basename $READ1).bam > $BNAME.bam
samtools index $BNAME.bam -o $BNAME.bam.bai
date
