#!/bin/bash

Wkdir=~/mge/recombinase/partial-carriage/mapping/jgi/wkdir
Readdir=~/mge/recombinase/partial-carriage/mapping/jgi/reads
Contigfile=~/data/isogenie-novaseq-contig/PFST0metaG_FD_JGI.fasta
Scriptdir=~/scripts
Cpu=28

Contigfile=$(readlink -f $Contigfile)
Bname=$(basename $(dirname $Contigfile))
Arr=(${Bname//_FD/ })
Label=${Arr[0]}
Readfile=$Readdir/$Label.fastq.gz
echo $Bname
echo $Label

set -e
mkdir -p $Wkdir
cd $Wkdir

date
if [ ! -f $Readfile ]; then
    echo "*** WARNING: $Readfile does not exist"
    touch $Bname.covm.tsv
    exit 0
fi

if [ ! -s $Bname.covm_tmean.tsv ]; then
  rm -f bam.outdir/$(basename $Contigfile).$(basename $Readfile).bam
  coverm contig -p bwa-mem --bam-file-cache-directory bam.outdir --interleaved $Readfile --reference $Contigfile --min-read-percent-identity 95 --min-read-aligned-percent 75 --min-covered-fraction 70 -m trimmed_mean --discard-unmapped -t $Cpu > $Bname.covm_tmean.tsv
else
  echo "*** $Bname.covm_tmean.tsv already exists and not empty"
fi
date
