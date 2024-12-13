#!/bin/bash

Wkdir=~/mge/recombinase
Seqfile=~/mge/contig-dram/20170700_S35.contigs/genes.faa
Hmm=~/mge/db/recombinase-hmm/bork.68subfam.hmm
Cpu=4

Bname=$(basename $(dirname $Seqfile) .contigs)

mkdir -p $Wkdir
cd $Wkdir

date
### get top hit
hmmsearch --noali --notextw -o /dev/null --cut_ga --cpu $Cpu --tblout $Bname.hmmtblout $Hmm $Seqfile
grep -v '^#'  $Bname.hmmtblout | tr -s ' ' | sort -t ' ' -k1,1 -k6,6rn -k9,9rn | cut -d ' ' -f1-18 --output-delimiter=$'\t' | tee $Bname.hmmtblout.tsv | awk '!seen[$1]++' >  $Bname.hmmtblout.tophit
### get protein seqs
cut -f1 $Bname.hmmtblout.tophit | python ~/scripts/pick-seq-in-list.py $Seqfile - > $Bname.hmmtblout.tophit.faa
date
