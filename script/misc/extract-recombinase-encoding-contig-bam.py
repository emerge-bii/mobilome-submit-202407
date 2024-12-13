#!/usr/bin/env python
import sys
import os
import glob
import subprocess

list_f = 'all_recombinase.coord.min_up_down_20000.recombinase.add_allinfo.sorted_by_contig_length.filt.tsv'
bam_dir = 'mapping/jgi/wkdir/bam.outdir'
outdir = 'all_recombinase.coord.min_up_down_20000.wkdir'

os.makedirs(outdir, exist_ok=True)
with open(list_f) as fh:
    for line in fh:
        rec, sample, contig = line.rstrip().split('\t')
        # assume there is only one matching
        bams = glob.glob(f'{bam_dir}/*{sample}*.bam')
        bams = sorted(bams)
        if len(bams) == 0:
            mes = f'[WARN] {rec} from sample {sample} not find in {bam_dir}, skipping\n'
            sys.stderr.write(mes)
            continue
        elif len(bams) > 1:
            mes = f'[WARN] muliple bam files matched {sample} and the first one is selected: {bams}\n'
            sys.stderr.write(mes)
        bam = bams[0]
        command = f'samtools view -h -o {outdir}/{contig}.tmpsam {bam} {contig}'
        #sys.stderr.write(f'{command}\n')
        command_lst = command.split()
        # require python >=3.7

        res = subprocess.run(command_lst, capture_output=True, text=True)
        try:
            res.check_returncode()
        except subprocess.CalledProcessError as e:
            if 'Could not retrieve index file' in res.stderr:
                index_command = f'samtools index {bam}'
                index_command_lst = index_command.split()
                index_res = subprocess.run(index_command_lst, capture_output=True, text=True)
                index_res.check_returncode()
                res = subprocess.run(command_lst, capture_output=True, text=True)
                res.check_returncode()
            else:
                raise
        #https://www.digitalocean.com/community/tutorials/how-to-use-subprocess-to-run-external-programs-in-python-3
        if 'specifies an invalid region or unknown reference' in res.stderr:
            mes = f'[Error] {contig} not found in {bam}\n'
            sys.stderr.write(mes)
            sys.exit(1)

        command = ['awk', f'($1=="@SQ" && $2=="SN:{contig}") || $1!="@SQ"', f'{outdir}/{contig}.tmpsam']
        with open(f'{outdir}/{contig}.sam', 'w') as fw:
            res = subprocess.run(command, check=True, stdout=fw)
        os.remove(f'{outdir}/{contig}.tmpsam')

