#!/bin/bash

####################################################
# run integron-finder for integron boundary curation
####################################################

Seqfile=~/mge/recombinase/refine/Integron.contig.fna  # contigs
Wkdir=~/mge/recombinase/refine/Integron.contig.wkdir  # working directory
Cpu=4  # threads

mkdir -p $Wkdir
cd $Wkdir

date
integron_finder --cpu $Cpu $Seqfile
cp Results_Integron_Finder_integron_contig/integron_contig.integrons integron.if.raw.tsv
awk '$11=="complete"' integron.if.raw.tsv > integron.if.4curate.tsv
date

