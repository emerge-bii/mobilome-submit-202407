#!/usr/bin/env python

import os
import sys
import pandas as pd

def main():
    if len(sys.argv) != 2:
        mes = '[INFO] Usage: python {} <table.blast.anicalc.dedup.add_unaligned.tsv>\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    tabfile = sys.argv[1]
    if tabfile == '-':
        tabfile = '/dev/stdin'

    st = set()

    df = pd.read_csv(tabfile, sep='\t', header=0)

    # need to fix the order of qname and tname in case that they are not sorted
    names = df.loc[:, ['qname', 'tname']].agg(sorted, axis=1) # array of list
    q_lst, t_lst = zip(*names.values)
    ani_lst = df.loc[:, ['qcov', 'tcov']].agg(min, axis=1) * df['pid'] / 100

    # make new df
    df_new = pd.DataFrame({'qname': q_lst, 'tname': t_lst, 'ANI': ani_lst})

    # filter out comparison from the same anchor gene
    df_new['qgene'] = df_new['qname'].map(lambda x: x.rsplit('_', 1)[0])
    df_new['tgene'] = df_new['tname'].map(lambda x: x.rsplit('_', 1)[0])
    filt = (df_new['qgene'] == df_new['tgene'])
    df_new = df_new.loc[~filt, :]

    df_new['qlocation'] = df_new['qname'].map(lambda x: x.rsplit('_', 1)[-1])
    df_new['tlocation'] = df_new['tname'].map(lambda x: x.rsplit('_', 1)[-1])

    mes = 'qname\ttname\tlow_ani\thigh_ani\n'
    sys.stdout.write(mes)
    gb = df_new.groupby(['qgene', 'tgene'])
    for (qname, tname), df_sub in gb:
        #print(qname, tname)
        assert len(df_sub) == 4, len(df_sub)
        filt1 = (df_sub['qlocation'] == 'upstream') & (df_sub['tlocation'] == 'upstream')
        filt2 = (df_sub['qlocation'] == 'downstream') & (df_sub['tlocation'] == 'downstream')
        filt3 = (df_sub['qlocation'] == 'upstream') & (df_sub['tlocation'] == 'downstream')
        filt4 = (df_sub['qlocation'] == 'downstream') & (df_sub['tlocation'] == 'upstream')
        ani1 = df_sub.loc[filt1, 'ANI'].iloc[0]
        ani2 = df_sub.loc[filt2, 'ANI'].iloc[0]
        ani3 = df_sub.loc[filt3, 'ANI'].iloc[0]
        ani4 = df_sub.loc[filt4, 'ANI'].iloc[0]
        ### if disregard orientation, ie. mobile gene is inserted in either strand 
        #if max([ani1, ani2]) >= max([ani3, ani4]):
        #    lo, hi = sorted([ani1, ani2])
        #else:
        #    lo, hi = sorted([ani3, ani4])
        ### if there is only one orientation, ie. mobile gene can only be inserted in one strand
        lo, hi = sorted([ani1, ani2])

        mes = f'{qname}\t{tname}\t{lo}\t{hi}\n'
        sys.stdout.write(mes)
    
if __name__ == '__main__':
    main()
