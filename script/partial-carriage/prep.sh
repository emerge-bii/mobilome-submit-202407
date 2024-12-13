#!/bin/bash

set -e

cd ~/mge/recombinase/partial-carriage/all_recombinase.coord.min_up_down_20000.list.splitdir

REC_ALLINFO=~/mge/recombinase/recombinase.allinfo.tsv
REC_LIST=all_recombinase.coord.min_up_down_20000.list

date

cat $REC_LIST | sed -E 's/^((.+)_[0-9]+)$/\1\t\2/' | perl -pe 's/^(SS.fna.[0-9]+_(.+?)(_FD_JGI_|_contigs_NODE_).+?\t)/\1\2\t/' | tee all_recombinase.coord.min_up_down_20000.recombinase.tsv | sed '1i recombinase\tsample\tcontig' | csvtk join -L -t --na NA -f "recombinase" - $REC_ALLINFO | (sed -u 1q && sort -nr -k9,9) | tee all_recombinase.coord.min_up_down_20000.recombinase.add_allinfo.sorted_by_contig_length.tsv | cut -f1-3 | sed '1d' > all_recombinase.coord.min_up_down_20000.recombinase.add_allinfo.sorted_by_contig_length.filt.tsv
