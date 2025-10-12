#!/bin/bash

######################################################################
# Eggnog-mapper annotation of recombinases from hmmsearch for curation
########################################################################

cd ~/mge/recombinase # working directory
FAA=gene.prefilter.hmmtblout.tophit.faa # recombinase protein seqs
DB=~/mge/db/emapper # emapper DB
CPU=48

emapper.py --cut_ga --report_no_hits --hmm_maxhits 0 --cpu $CPU --pfam_realign denovo -i $FAA -o gene.prefilter.hmmtblout.tophit.emapper --data_dir $DB --override
date

