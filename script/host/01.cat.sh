#!/bin/bash

#################################
# host taxon assignment with CAT
#################################

# CAT DB for GTDB notes here: https://github.com/MGXlab/gtdb2cat
module use /fs/project/PAS1117/modulefiles
module load CAT/5.2.3

Seqfile=~/mge/recombinase/recombinase.contig.fna  # contigs
Wkdir=~/mge/recombinase/host/contig-taxa  # working directory
DB=~/mge/db/cat/CAT_database  # DB
TAX=~/mge/db/cat/CAT_taxonomy # DB
CPU=28 # threads

Wkdir=$Seqfile.catout
mkdir -p $Wkdir
cd $Wkdir

date
CAT contigs -c $Seqfile -d $DB -t $TAX --sensitive -n $CPU

CAT add_names -i out.CAT.contig2classification.txt -o out.CAT.contig2classification.tax.txt -t $TAX
CAT add_names --only_official -i out.CAT.contig2classification.txt -o out.CAT.contig2classification.tax_only_official.txt -t $TAX
CAT summarise -c $Seqfile -i out.CAT.contig2classification.tax_only_official.txt -o out.CAT.contig2classification.tax_only_official.summary
date
