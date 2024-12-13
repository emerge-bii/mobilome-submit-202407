#!/bin/bash

cd ~/mge/recombinase/is-split-gene
Cpu=8

date
awk '$5=="IS_Tn"' ../recombinase.allinfo.tsv | cut -f1 | python get-gene-neighborhood-exclude-gene-exclude-intergenic.py --sep-ends --min-oneside-len 50 --min-twoside-len 500 --radius 4000 ../recombinase.dram_up2_down2_w_anchor.tsv - ../recombinase.contig.fna > is_tn.gene_neighbor.fna 2> is_tn.gene_neighbor.log

minimap2 -N 10 -t $Cpu ~/mge/pc/mmseq/SS.95_rep_seq.fna is_tn.gene_neighbor.fna > is_tn.gene_neighbor.ali2pc_rep_cds.paf
date


