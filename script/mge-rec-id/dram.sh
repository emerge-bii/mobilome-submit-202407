#!/bin/bash

Wkdir=~/mge/contig-dram
Seqfile=~/mge/data/isogenie-hiseq-contigs/20170700_S35.contigs.fna
Cpu=18
Minlen=1500

mkdir -p $Wkdir
Outdir=$Wkdir/$(basename $Seqfile .fna)

DRAM.py annotate -i $Seqfile -o $Outdir --threads $Cpu --min_contig_size $Minlen
date
