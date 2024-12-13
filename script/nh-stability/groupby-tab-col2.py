#!/usr/bin/env python

import sys
import os
import pandas as pd

#IGNORE = []
IGNORE = ['d__Other', 'p__Other', 'c__Other', 'o__Other', 'f__Other', 'g__Other', 's__Other']
IGNORE = set(IGNORE)

def main():
    if len(sys.argv) < 2:
        mes = '[INFO] Usage: python {} <table.tsv> col1 col2 ..\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    tab = sys.argv[1]
    cols = sys.argv[2:]

    if tab == '-':
        tab = '/dev/stdin'

    df = pd.read_csv(tab, sep='\t', header=0)
    cols_all = df.columns
    for i in cols:
        assert i in cols_all, f'{i} is not a column name'
    cols_other = [i for i in list(cols_all) if not i in set(cols)]
    mes = '\t'.join(cols + ['total_cnt', 'uniq_cnt'] + cols_other)
    mes = f'{mes}\n'
    sys.stdout.write(mes)
    
    df = df.fillna('NA')
    gb = df.groupby(cols)
    for item, sub_df in gb:
        if len(cols) == 1:
            lst = [item]
        else:
            lst = list(item)

        item_cnt = len(sub_df)
        lst.append(str(item_cnt))
        for n in cols_other:
            cnt_ser = sub_df[n].value_counts(dropna=False)
            _l = [f'{m}:{cnt_ser.loc[m]}' for m in cnt_ser.index if not m in IGNORE]
            uniq_item_cnt = len(_l)
            lst.append(str(uniq_item_cnt))
            tmp_str = ','.join(_l)
            lst.append(tmp_str)

        mes = '\t'.join(lst)
        mes = f'{mes}\n'
        sys.stdout.write(mes)

if __name__ == '__main__':
    main()
