#!/bin/bash

cd ~/mge/recombinase/stability
Faa=../recombinase.faa
Bname=$(basename $Faa .faa)
Cpu=24

date
cd-hit -i $Faa -c 1 -g 0 -G 1 -l 30 -d 0 -M 0 -T $Cpu -o $Bname.cdhit1d0.faa
date
