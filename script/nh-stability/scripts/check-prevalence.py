#!/usr/bin/env python
import sys
import os
import pandas as pd
from collections import Counter

def main():
    if len(sys.argv) < 2:
        mes = '[INFO] Usage: python {} <cdhit.add_allinfo.tsv> col1 col2..\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    tab = sys.argv[1]
    cols = sys.argv[2:]

    if tab == '-':
        tab = '/dev/stdin'
    cols_str = '___'.join(cols)
    n_cols_str = f'n_{cols_str}'
    df = pd.read_csv(tab, sep='\t', header=0)
    if 'Year' in df.columns:
        df['Year'] = df['Year'].map(lambda x: "NA" if pd.isna(x) else int(x))

    df[cols] = df[cols].apply(lambda x: x.astype(str), axis=1)
    gb = df.groupby('OTU')
    for otu, df_sub in gb:
        df_sub[cols_str] = df_sub[cols].agg('___'.join, axis=1)
        ser = df_sub[cols_str].drop_duplicates()
        lst = [i for i in ser if i != 'p__Other']  # remove 'p__Other' for phylum col
        n = len(lst)
        s = ';'.join(ser)
        d_cnt = Counter(df_sub[cols_str])
        s_more = ';'.join([f'{key}|{d_cnt[key]}' for key in d_cnt])
        n_mem = len(df_sub)
        mes = f'{otu}\t{s}\t{n}\t{s_more}\t{n_mem}\n'
        sys.stdout.write(mes)

if __name__ == '__main__':
    main()
