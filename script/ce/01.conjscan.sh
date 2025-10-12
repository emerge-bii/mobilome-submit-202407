#!/bin/bash

########################################
# run conjscan for CE boundary curation
########################################
set -e

date
cd ~/mge/recombinase/refine  # working directory
Seqfile=~/mge/recombinase/refine/CE.contig.genes.faa  # ORFs from contigs (protein)
Cpu=24  # threads

### run plasmid models
Wkdir=~/mge/recombinase/refine/CE.contig.genes.wkdir1
macsyfinder --models CONJScan/Plasmids all \
	    --sequence-db "$Seqfile" \
	    -w $Cpu \
	    --db-type gembase \
	    -o $Wkdir >> /dev/null # stdout is already reported in output files from macsyfinder (macsyfinder.out)

### see annotations in "best_solution.tsv"
date

### run integrated CE models
Wkdir=~/mge/recombinase/refine/CE.contig.genes.wkdir2
macsyfinder --models CONJScan/Chromosome all \
	    --sequence-db "$Seqfile" \
	    -w $Cpu \
	    --db-type gembase \
	    -o $Wkdir >> /dev/null # stdout is already reported in output files from macsyfinder (macsyfinder.out)

### see annotations in "best_solution.tsv"
date

