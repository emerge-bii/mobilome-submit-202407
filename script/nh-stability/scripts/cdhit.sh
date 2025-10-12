#!/bin/bash

####################################
# run CDHIT to cluster recombinases
###################################

cd ~/mge/recombinase/stability  # working directory
Faa=../recombinase.faa  # recombinase sequence file (protein)
Cpu=24  # threads

Bname=$(basename $Faa .faa)

date
cd-hit -i $Faa -c 1 -g 0 -G 1 -l 30 -d 0 -M 0 -T $Cpu -o $Bname.cdhit1d0.faa
date
