#!/bin/bash

cd ~/mge/recombinase/stability
Cpu=4

Scriptdir=~/scripts
Clusttab=gene.prefilter.hmmtblout.tophit.cdhit1d0.faa.clstr
Rec_allinfo_f=../recombinase.allinfo.tsv
Rec_dram_neighborhood_f=../gene.prefilter.hmmtblout.tophit.dram_up2_down2_w_anchor.tsv
Rec_contig_f=../gene.prefilter.hmmtblout.tophit.contig.fna
Label=all

date
set -e

list_descendants () {
    local children=$(ps -o pid= --ppid "$1")
    for pid in $children
    do
        list_descendants "$pid"
    done
    echo "$children"
}

clean_up() {
    #pkill -P $$
    kill $(list_descendants $$)
}

trap clean_up SIGINT SIGTERM

python ~/scripts/parse-cdhit-aa.py $Clusttab > $Clusttab.tsv

csvtk join -t -L --na "NA" -f "mem;seqname" $Clusttab.tsv $Rec_allinfo_f > $Clusttab.add_allinfo.tsv

# check OTU in prevalence across years
awk -v FS=$'\t' '$19!="NA" && $20!="Collapsed Palsa"' $Clusttab.add_allinfo.tsv | python $Scriptdir/check-prevalence.py - Year > $Clusttab.add_allinfo.pick_otu.Year.tsv


# filter number of years detected (col 3) and OTU size (col 5)
awk -v FS=$'\t' '$19!="NA" && $20!="Collapsed Palsa"' $Clusttab.add_allinfo.tsv | python ~/scripts/pick-tabrow-from-list.py - <(awk '$3>=0 && $5>=2' $Clusttab.add_allinfo.pick_otu.Year.tsv | cut -f1) > $Clusttab.add_allinfo.filt_sample.$Label.tsv


# exclude mobile gene and include inter-genic regions up- and down-stream in its neighborhood
cat $Clusttab.add_allinfo.filt_sample.$Label.tsv | cut -f2 | python $Scriptdir/get-gene-neighborhood-exclude-gene-exclude-intergenic.py --fixed-length --sep-ends --radius 300 $Rec_dram_neighborhood_f - $Rec_contig_f > $Label.gene_neighbor.fna 2> $Label.gene_neighbor.log

date
mkdir -p $Label.splitdir
cut -f1 $Clusttab.add_allinfo.filt_sample.$Label.tsv | tail -n+2 | uniq | split -d -l 1000 --suffix-length 4 --additional-suffix .split - $Label.splitdir/x


find $Label.splitdir -type f -name "x*.split" | head -n 1 | xargs -n1 -P $Cpu --process-slot-var=IDX -I {} bash -c '
Split_f={}
mkdir -p $Split_f.outdir
cat $Split_f | while read OTU; do
  sleep 0.2
  Tmpfile=$(mktemp); trap "rm -f $Tmpfile" 0 2 3 15
  grep --no-group-separator -A 1 -f <(awk -v OTU=$OTU '\''$1==OTU {print $2}'\'' $Clusttab.add_allinfo.filt_sample.$Label.tsv) $Label.gene_neighbor.fna > $Tmpfile || :
  if [ ! -s $Tmpfile ]; then
    echo "[INFO] no seqs found for $OTU after neighborhood screening, skipping.." && continue
  fi
  bash ~/scripts/blastn-ani-clust.sh $Tmpfile $OTU $Split_f.outdir 1 blastn
  python $Scriptdir/rm-dup-pairwise-ani.py $Split_f.outdir/$OTU.blast.anicalc.tsv | python $Scriptdir/add-unaligned-to-pairwise-ani.py - $Tmpfile > $Split_f.outdir/$OTU.blast.anicalc.dedup.add_unaligned.tsv
  if [ ! -f $Split_f.outdir/$OTU.blast.anicalc.dedup.add_unaligned.tsv ]; then
    continue
  fi
  # NOTE: there is one extra line (header) counted in .unstable_cnt.tsv; all "NR" col needs to decrease by 1
  cat $Split_f.outdir/$OTU.blast.anicalc.dedup.add_unaligned.tsv | python $Scriptdir/sort-up-down-stream-pairwise.py - | awk -v OFS=$'\''\t'\'' -v OTU=$OTU '\''NR==1 {print $0, "OTU"; next} {print $0, OTU}'\'' | tee $Split_f.outdir/$OTU.blast.anicalc.dedup.add_unaligned.sep_ends_sum.tsv | awk -v OTU=$OTU -v OFS=$'\''\t'\'' '\''BEGIN{Cnt=0} $3<1 {Cnt+=1} END{if  (NR>1) {print OTU, Cnt, NR-1, Cnt/(NR-1)}}'\'' > $Split_f.outdir/$OTU.blast.anicalc.dedup.add_unaligned.sep_ends_sum.unstable_cnt.tsv
  cat $Split_f.outdir/$OTU.blast.anicalc.dedup.add_unaligned.sep_ends_sum.tsv | cut -f1-3 | python $Scriptdir/aniclust-mod.py --one-mem-per-line --min-ani 90 - | awk -v OFS=$'\''\t'\'' -v OTU=$OTU '\''{print $0, OTU}'\'' > $Split_f.outdir/$OTU.blast.anicalc.dedup.add_unaligned.sep_ends_sum.ani_clust.tsv
  sed '1i Q1\tQ2\tNH\tOTU' $Split_f.outdir/$OTU.blast.anicalc.dedup.add_unaligned.sep_ends_sum.ani_clust.tsv | cut -f3,4 |  python ~/scripts/groupby-tab-col2.py - OTU > $Split_f.outdir/$OTU.blast.anicalc.dedup.add_unaligned.sep_ends_sum.ani_clust.nh_per_otu.tsv
  echo "[INFO] $OTU"
done
'

printf "%s\n" $Label.splitdir/x*.split.outdir/*.blast.anicalc.dedup.add_unaligned.tsv | xargs cat | awk 'NR==1||FNR!=1'  > $Label.gene_neighbor.blastn.pw_iden_within_cluster.tsv

printf "%s\n" $Label.splitdir/x*.split.outdir/*.blast.anicalc.dedup.add_unaligned.sep_ends_sum.tsv | xargs awk 'NR==1||FNR!=1'  > $Label.gene_neighbor.blastn.pw_iden_within_cluster.sep_ends_sum.tsv

printf "%s\n" $Label.splitdir/x*.split.outdir/*.unstable_cnt.tsv | xargs cat > $Label.gene_neighbor.blastn.pw_iden_within_cluster.sep_ends_sum.unstable_cnt.tsv

printf "%s\n" $Label.splitdir/x*.split.outdir/*.blast.anicalc.dedup.add_unaligned.sep_ends_sum.ani_clust.tsv | xargs cat > $Label.gene_neighbor.blastn.pw_iden_within_cluster.sep_ends_sum.ani_clust.tsv

printf "%s\n" $Label.splitdir/x*.split.outdir/*.blast.anicalc.dedup.add_unaligned.sep_ends_sum.ani_clust.nh_per_otu.tsv | xargs awk 'NR==1||FNR!=1' > $Label.gene_neighbor.blastn.pw_iden_within_cluster.sep_ends_sum.ani_clust.nh_per_otu.tsv

csvtk join -t -L --na NA -f"OTU" <(sed '1i OTU\tunstable_pair\ttotal_pair\tratio' $Label.gene_neighbor.blastn.pw_iden_within_cluster.sep_ends_sum.unstable_cnt.tsv) $Label.gene_neighbor.blastn.pw_iden_within_cluster.sep_ends_sum.ani_clust.nh_per_otu.tsv > $Label.gene_neighbor.blastn.pw_iden_within_cluster.sep_ends_sum.unstable_cnt.add_nh_per_otu.tsv

date



