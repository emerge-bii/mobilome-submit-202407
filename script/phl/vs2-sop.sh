#!/bin/bash

Seqfile=~/mge/recombinase/refine/phage/Phage.certain.contig.fna
Wkdir=~/mge/recombinase/refine/phage
Label=0d5
Cpu=28

Scriptdir=~/scripts

set -e

Dir=$(cd $(dirname $Seqfile) && pwd)
Seqfile=$Dir/$(basename $Seqfile)
Bname=$(basename $Seqfile .fna)

date
mkdir -p $Wkdir/$Bname
cd $Wkdir/$Bname

virsorter run -i $Seqfile --keep-original-seq -w vs2-pass1 --min-score 0.5 -l $Label -j $Cpu --include-groups dsDNAphage,ssDNA --min-length 5000 all

checkv end_to_end vs2-pass1/${Label}-final-viral-combined.fa checkv -t $Cpu -d /fs/project/PAS1117/jiarong/db/checkv-db-v1.0

cat checkv/proviruses.fna checkv/viruses.fna > checkv/combined.fna

python $Scriptdir/merge-vs2-checkv-tab.py vs2-pass1/${Label}-final-viral-score.tsv checkv/contamination.tsv vs2-checkv-merged.tsv

virsorter run --seqname-suffix-off --viral-gene-enrich-off --provirus-off --prep-for-dramv -i checkv/combined.w_int__viral_gene_ge1.fna -w vs2-pass2 --include-groups dsDNAphage,ssDNA --min-length 0 --min-score 0.5 -j $Cpu all

DRAM-v.py annotate -i vs2-pass2/for-dramv/final-viral-combined-for-dramv.fa -v vs2-pass2/for-dramv/viral-affi-contigs-for-dramv.tab -o dramv-annotate --skip_trnascan --threads $Cpu --min_contig_size 1000
##step 2 summarize these anntotations
DRAM-v.py distill -i dramv-annotate/annotations.tsv -o dramv-distill


