#!/bin/bash

########################################
# DRAM annotation of contigs to produce 
#   ORFs and KEGG and PFAM annotations
########################################

Wkdir=~/mge/contig-dram # working directory
Seqfile=~/mge/data/isogenie-hiseq-contigs/20170700_S35.contigs.fna  # contigs
Cpu=18  # threads
Minlen=1500  # minimal contig length required

mkdir -p $Wkdir
Outdir=$Wkdir/$(basename $Seqfile .fna)

DRAM.py annotate -i $Seqfile -o $Outdir --threads $Cpu --min_contig_size $Minlen
date
