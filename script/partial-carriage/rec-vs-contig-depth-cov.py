#!/usr/bin/env python

import sys
import os
import numpy as np
import pandas as pd

mask_start = 150
mask_end = 150
min_depth_cov = 1
min_breadth_cov = 0.9

def main():
    if len(sys.argv) != 3:
        mes = '[INFO] Usage: python {} <contig.sam.depth> <rec2start2end.tsv>\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    depth_f = sys.argv[1]
    coord_f = sys.argv[2]

    df_depth = pd.read_csv(depth_f, sep='\t', header=0)
    df_coord = pd.read_csv(coord_f, sep='\t', header=0)

    contig_len = len(df_depth)

    df_depth = df_depth.iloc[mask_start:-1*(mask_end), :]
    contig_mean = np.mean(df_depth['depth_cov'])
    contig_std = np.std(df_depth['depth_cov'])

    # header: rec | cohen_d | carriage_ratio | rec_len | rec_mean | rec_std | contig_len | contig_mean | contig_std
    mes = '\t'.join(['recombinase', 'cohen_d', 'carriage_ratio', 'rec_len', 'rec_mean', 'rec_std', 'contig_len', 'contig_mean', 'contig_std'])
    sys.stdout.write(f'{mes}\n')

    for i in range(len(df_coord)):
        row = df_coord.iloc[i, :]
        rec = row.loc['recombinase']
        start = row.loc['start']
        end = row.loc['end']
        if start <= mask_start or end >= contig_len - mask_end:
            mes = f'[WARN] recombinase {rec} ({start} to {end} of {contig_len} total) is too close to edge, results may be not reliable..\n'
            sys.stderr.write(mes)

        start = start - mask_start - 1
        end = end - mask_start - 1

        rec_len = end - start + 1

        sel = (df_depth['position'] >= start) & (df_depth['position'] <= end)
        df_rec = df_depth.loc[sel, :]

        rec_mean = np.mean(df_rec['depth_cov'])
        rec_std = np.std(df_rec['depth_cov'])
        #pooled_std = contig_std  # treat overall std as pooled std

        pooled_std = ((rec_std**2 + contig_std**2)/2)**0.5
        cohen_d = abs((rec_mean - contig_mean)/pooled_std)

        carriage_ratio = rec_mean / contig_mean
        mes = f'{rec}\t{cohen_d}\t{carriage_ratio}\t{rec_len}\t{rec_mean}\t{rec_std}\t{contig_len}\t{contig_mean}\t{contig_std}\n'
        sys.stdout.write(mes)

if __name__ == '__main__':
    main()
