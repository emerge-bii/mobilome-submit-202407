#!/bin/bash


Seqfile=~/mge/recombinase/refine/Integron.contig.fna
Wkdir=~/mge/recombinase/refine/Integron.contig.wkdir
Cpu=4

mkdir -p $Wkdir
cd $Wkdir

date
integron_finder --cpu $Cpu $Seqfile
cp Results_Integron_Finder_integron_contig/integron_contig.integrons integron.if.raw.tsv
awk '$11=="complete"' integron.if.raw.tsv > integron.if.4curate.tsv
date

