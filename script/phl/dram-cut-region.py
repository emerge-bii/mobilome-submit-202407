#!/usr/bin/env python

import sys
import os
import argparse
import pandas as pd

CONTIG_COLNAME_IN_DRAM_DICT = {'bp': 'scaffold', 'orf': 'contig'}

def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--bp', action='store_true', help='use "start" and "end" column for bp positions')
    group.add_argument('--orf', action='store_true', help='use "start_orf" and "end_orf" column for gene positions')
    parser.add_argument('boundary', help='table with boundary of elements, ie. "start" and "end" or "start_orf" and "end_orf" columns')
    parser.add_argument('dram', help='dram annotations')
    args = parser.parse_args()
    args_d = vars(args)
	
    bdy_f = args_d['boundary']
    dram_f = args_d['dram']

    if bdy_f == '-' and dram_f == '-':
        mes = '[ERROR] <vs2sop_boundary.tsv> and <contig.dram.tsv> cannot both be from stdin\n'
    elif bdy_f == '-':
        bdy_f = '/dev/stdin'
    elif dram_f == '-':
        dram_f = '/dev/stdin'

    df_dram = pd.read_csv(dram_f, sep='\t', header=0, index_col=0)
    df_bdy = pd.read_csv(bdy_f, sep='\t', header=0)

    df_dram['contig'] = df_dram.index.map(lambda x: x.rsplit('_', 1)[0])

    if args_d['bp']:
        col = CONTIG_COLNAME_IN_DRAM_DICT['bp']
    elif args_d['orf']:
        col = CONTIG_COLNAME_IN_DRAM_DICT['orf']
    # dram scaffold <-> bdy contig
    overlap_scaffold_st = set(df_dram[col]) & set(df_bdy['contig'])
    filt = df_bdy['contig'].isin(overlap_scaffold_st)
    df_bdy = df_bdy.loc[filt,:]

    filt = df_dram[col].isin(overlap_scaffold_st)
    df_dram = df_dram.loc[filt, :]

    lis = []
    for index, row in df_bdy.iterrows():
        contig = row.loc['contig']
        filt = (df_dram[col] == contig)
        df_dram_sub = df_dram.loc[filt, :]
        total_orf_num = len(df_dram_sub)

        if args_d['bp']:
            phage_start_bp = row.loc['start']
            phage_end_bp = row.loc['end']

            # get gene coord from dram
            # check if the above coord fall within range in bdy

            # prophage
            phage_start_orf = sum(df_dram_sub['start_position'] <= phage_start_bp)
            # the start boudary in bdy_f may not match exactly at orf start if its the first orf
            if phage_start_orf == 0:
                phage_start_orf = 1
            phage_end_orf = sum(df_dram_sub['end_position'] <= phage_end_bp)

        elif args_d['orf']:
            phage_start_orf = row.loc['start_orf']
            phage_end_orf = row.loc['end_orf']

        df_dram_sub_prophage = df_dram_sub.iloc[(phage_start_orf-1):phage_end_orf]
        lis.append(df_dram_sub_prophage)

    df_dram_prophage_merged = pd.concat(lis)
    df_dram_prophage_merged.to_csv('/dev/stdout', sep='\t', index=True, header=True, index_label=None)

if __name__ == '__main__':
    main()
