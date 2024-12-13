#!/bin/bash
#SBATCH --nodes=1 --ntasks=28 --time=96:00:00
#SBATCH -A PAS1117

if [ $# -ne 5 ]; then
  echo "*** bash $(basename $0) <seqfile.fna> <basename> <outdir> cpu_n task"
  echo "***  task options: 'blastn', 'blastn-short', 'dc-megablast', 'megablast', 'rmblastn'"
  exit 1
fi

#Seqfile=vs2_sop-len_5k.fa
Seqfile=$1
Bname=$2
Outdir=$3
Cpu=$4
Task=$5


set -e
mkdir -p $Outdir


echo
echo "[INFO] $Bname"
echo "[INFO] Starting blastn: $(date)"
makeblastdb -in $Seqfile -dbtype nucl
##blastn -query $Seqfile -db $Seqfile -outfmt '6 std qlen slen' -max_target_seqs 10000 -perc_identity 90 -out $Outdir/$Bname.blast.tmp -num_threads $Cpu
blastn -task $Task -query $Seqfile -db $Seqfile -outfmt '6 std qlen slen' -max_target_seqs 10000 -out $Outdir/$Bname.blast.tsv -num_threads $Cpu
awk '$12>=50 && $1!=$2' $Outdir/$Bname.blast.tsv > $Outdir/$Bname.blast.filt.tsv

echo "[INFO] Starting anicalc.py: $(date)"
python ~/lab/jiarong/tools/checkv/scripts/anicalc.py -i $Outdir/$Bname.blast.filt.tsv -o $Outdir/$Bname.blast.anicalc.tsv

#echo "[INFO] Starting aniclust.py: $(date)"
#python ~/lab/jiarong/tools/checkv/scripts/aniclust.py --fna $Seqfile --ani $Outdir/$Bname.blast.anicalc.tsv --out $Outdir/$Bname.blast.anicalc.cluster.tsv --min_ani 95 --min_qcov 0 --min_tcov 80
#python ~/lab/jiarong/tools/checkv/scripts/aniclust.py --fna $Seqfile --ani $Outdir/$Bname.blast.anicalc.tsv --out $Outdir/$Bname.blast.anicalc.cluster.tsv --min_ani 80 --min_qcov 0 --min_tcov 90

