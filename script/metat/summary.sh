#!/bin/bash

set -e

date

### metaT setting
Pydirseqtab=~/mge/metat/nova-only.wkdir/714E11014metaG_FD_JGI/SS.final.tsv
Rectab=~/mge/recombinase/recombinase.allinfo.tsv
Lab=cov_0d9
Mincov=0.9

### paired metaG setting
Metag_min_depth_cov=0
Lab2=depth_cov_0
Metag_cov_dir=~/mge/metat/paired-metag/all.wkdir

Wkdir=$(dirname $Pydirseqtab)
Sample=$(basename $Wkdir)
Sample=${Sample//_FD_JGI/}
Sample=${Sample//_contigs/}
Pydirseqtab=$(readlink -f $Pydirseqtab)

echo $Sample

# locate metaG depth cov file
Metag_cov_f=$(ls $Metag_cov_dir/"$Sample"*coverm_trimmed_mean.tsv | head -n1)
Metag_cov_f=$(readlink -f $Metag_cov_f)

# locate anno file
if [[ $Sample =~ ^.+metaG(_2){0,1}$ ]]; then
  # 714E11014metaG
  Anno_f=~/mge/contig-dram/jgi-jgi/${Sample}_FD_JGI/annotations.tsv
elif [[ $Sample =~ ^201....._[PSE]..$ ]]; then
  # 20170700_S35
  Anno_f=~/mge/contig-dram/field/field/${Sample}_contigs/annotations.tsv
fi

cd $Wkdir

### start metaT process
### recombinase only; conservative with coverage ($10) >=0.9;
Arr=($(awk '{fbp+=$7*($5-$4); rbp+=$11*($5-$4)} END{total=fbp+rbp; printf "\t%.3f\t%.3f\t%d\n", fbp/total, rbp/total, total}' $Pydirseqtab))
if (( $(echo "${Arr[0]} > ${Arr[1]}" | bc -l) )); then
  python ~/scripts/pick-tabrow-from-list.py $Rectab <(awk -v mincov=$Mincov '$10>=mincov' $Pydirseqtab | tail -n+2 | cut -f1) | cut -f1 | tail -n+2 > $Lab.rec.active.list
  ### all genes
  awk -v mincov=$Mincov '$10>=mincov' $Pydirseqtab | tail -n+2 | cut -f1 > $Lab.all.active.list
else
  python ~/scripts/pick-tabrow-from-list.py $Rectab <(awk -v mincov=$Mincov '$14>=mincov' $Pydirseqtab | tail -n+2 | cut -f1) | cut -f1 | tail -n+2 > $Lab.rec.active.list
  ### all genes
  awk -v mincov=$Mincov '$14>=mincov' $Pydirseqtab | tail -n+2 | cut -f1 > $Lab.all.active.list
fi

python ~/scripts/pick-tabrow-from-list.py $Anno_f $Lab.all.active.list > $Lab.all.active.dram.tsv

python ~/scripts/pick-tabrow-not-from-list.py $Lab.all.active.dram.tsv $Lab.rec.active.list > $Lab.nonrec.active.dram.tsv

Cnt_rec=$(wc -l $Lab.rec.active.list | cut -d' ' -f1) && cut -f9 $Lab.nonrec.active.dram.tsv | tail -n+2 | sort | uniq -c | tr -s ' ' | sed 's/^ //' | sed '1s/$/Noko/' | sed "1i $Cnt_rec\trecombinase" | sort -k1,1 -rn | awk -v OFS=$'\t' '{print $2, $1}' | tee $Lab.all.active.dram.ko_cnt.tsv | csvtk join -t -H -L --na NA -f"1;1" - <(cut -f4,5 ~/mge/db/kegg/ko00001.tsv | tail -n+2 | sort | uniq) > $Lab.all.active.dram.ko_cnt.add_anno.tsv


### metaG
# filt dram table by depth cov
tail -n+2 $Metag_cov_f | awk -v metag_min_depth_cov=$Metag_min_depth_cov '$2 >=metag_min_depth_cov {print $1}' | python ~/scripts/pick-dram-row-from-contig-list.py $Anno_f - > $Lab2.metag.all.filt_depth.dram.tsv

# get rec gene list with enough detph cov
python ~/scripts/pick-tabrow-from-list.py $Rectab <(tail -n+2 $Lab2.metag.all.filt_depth.dram.tsv | cut -f1) | cut -f1 | tail -n+2 > $Lab2.metag.rec.filt_depth.list

# get nonrec gene dram table
python ~/scripts/pick-tabrow-not-from-list.py $Lab2.metag.all.filt_depth.dram.tsv $Lab2.metag.rec.filt_depth.list > $Lab2.metag.nonrec.filt_depth.dram.tsv

#  
Cnt_rec=$(wc -l $Lab2.metag.rec.filt_depth.list | cut -d' ' -f1) && cut -f9 $Lab2.metag.nonrec.filt_depth.dram.tsv | tail -n+2 | sort | uniq -c | tr -s ' ' | sed 's/^ //' | sed '1s/$/Noko/' | sed "1i $Cnt_rec\trecombinase" | sort -k1,1 -rn | awk -v OFS=$'\t' '{print $2, $1}' | tee $Lab2.metag.all.filt_depth.dram.ko_cnt.tsv | csvtk join -t -H -L --na NA -f"1;1" - <(cut -f4,5 ~/mge/db/kegg/ko00001.tsv | tail -n+2 | sort | uniq) > $Lab2.metag.all.filt_depth.dram.ko_cnt.add_anno.tsv

csvtk join -t -H -L --na NA -f"1;1" $Lab.all.active.dram.ko_cnt.add_anno.tsv $Lab2.metag.all.filt_depth.dram.ko_cnt.tsv > $Lab.all.active.dram.ko_cnt.add_anno.add_metag_cnt_w_${Lab2}.tsv

### optional: remove "Noko" row
sed -i '/^Noko/d' $Lab.all.active.dram.ko_cnt.add_anno.add_metag_cnt_w_${Lab2}.tsv

awk -v FS=$'\t' -v OFS=$'\t' '{print $0,$2/$4}' $Lab.all.active.dram.ko_cnt.add_anno.add_metag_cnt_w_${Lab2}.tsv | sort -t $'\t' -k5,5 -nr | awk -v OFS=$'\t' -v FS=$'\t' '{print $0,NR}' | sort -t $'\t' -k4,4 -nr | awk -v OFS=$'\t' -v FS=$'\t' '{print $0,NR}' | sort -t $'\t' -k2,2 -nr | awk -v OFS=$'\t' -v FS=$'\t' '{print $0,NR}' | awk -v sample=$Sample -v OFS=$'\t' -v FS=$'\t' '{print $8,$7,$6,$1,$2,$4,$5,sample,$3}' | sed '1i rank_metat\trank_metag\trank_ratio\tko\tcnt_metat\tcnt_metag\tratio\tsample\tanno' > $Lab.all.active.dram.ko_cnt.add_anno.add_metag_cnt_w_${Lab2}.add_ratio.add_rank.tsv

date
