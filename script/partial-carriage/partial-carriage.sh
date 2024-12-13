#!/bin/bash

set -e

cd ~/mge/recombinase/partial-carriage/all_recombinase.coord.min_up_down_20000.list.splitdir

REC_ALLINFO=~/mge/recombinase/recombinase.allinfo.tsv
REC_DRAM=~/mge/recombinase/recombinase.allinfo.dram.tsv
LABEL=all_recombinase.coord.min_up_down_20000
SCRIPTDIR=~/scripts

date

cut -f3 all_recombinase.coord.min_up_down_20000.recombinase.add_allinfo.sorted_by_contig_length.filt.tsv | sort | uniq | while read line; do 
    CONTIG=$line
    if [ ! -s $LABEL.wkdir/$CONTIG.sam ]; then
        echo "[INFO] $LABEL.wkdir/$CONTIG.sam does not exist, ie. mapping data is not available, skipping.."
	continue
    fi
    mkdir -p $LABEL.wkdir/$CONTIG.outdir
    grep $CONTIG $REC_DRAM | awk -v OFS=$'\t' '{print $1, $2"_"$3, $5, $6}' | sed '1i recombinase\tcontig\tstart\tend' | csvtk join -L -t --na NA -f "recombinase" - $REC_ALLINFO | tee $LABEL.wkdir/$CONTIG.outdir/tmp.recombinase.coord.more.tsv | csvtk cut -t -f"contig,start,end,origin" > $LABEL.wkdir/$CONTIG.outdir/tmp.recombinase.coord.tsv
    samtools view -h $LABEL.wkdir/$CONTIG.sam 2>/dev/null | python $SCRIPTDIR/filtersam.py -m identity -c 95 - | python $SCRIPTDIR/filtersam.py -m aligned -c 90 - | samtools depth -a - | sed '1i recombinase\tposition\tdepth_cov' > $LABEL.wkdir/$CONTIG.outdir/tmp.sam.depth
    python $SCRIPTDIR/rec-vs-contig-depth-cov.py $LABEL.wkdir/$CONTIG.outdir/tmp.sam.depth $LABEL.wkdir/$CONTIG.outdir/tmp.recombinase.coord.more.tsv > $LABEL.wkdir/$CONTIG.outdir/tmp.rec_vs_contig_depth_cov.tsv
done

printf "%s\n" $LABEL.wkdir/*.outdir/tmp.rec_vs_contig_depth_cov.tsv | xargs awk 'NR==1 || FNR!=1' > $LABEL.rec_vs_contig_depth_cov.tsv

date

