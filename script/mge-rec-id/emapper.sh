#!/bin/bash

cd ~/mge/recombinase # working directory
Cpu=48

emapper.py --cut_ga --report_no_hits --hmm_maxhits 0 --cpu $Cpu --pfam_realign denovo -i gene.prefilter.hmmtblout.tophit.faa -o gene.prefilter.hmmtblout.tophit.emapper --data_dir ~/mge/db/emapper --override
date

