#!/bin/bash

set -e

date
cd ~/mge/recombinase/refine
Seqfile=~/mge/recombinase/refine/CE.contig.genes.faa
Cpu=24

Wkdir=~/mge/recombinase/refine/CE.contig.genes.wkdir1
macsyfinder --models CONJScan/Plasmids all \
	    --sequence-db "$Seqfile" \
	    -w $Cpu \
	    --db-type gembase \
	    -o $Wkdir >> /dev/null # stdout is already reported in output files from macsyfinder (macsyfinder.out)

### see annotations in "best_solution.tsv"
date

Wkdir=~/mge/recombinase/refine/CE.contig.genes.wkdir2
macsyfinder --models CONJScan/Chromosome all \
	    --sequence-db "$Seqfile" \
	    -w $Cpu \
	    --db-type gembase \
	    -o $Wkdir >> /dev/null # stdout is already reported in output files from macsyfinder (macsyfinder.out)

### see annotations in "best_solution.tsv"
date

