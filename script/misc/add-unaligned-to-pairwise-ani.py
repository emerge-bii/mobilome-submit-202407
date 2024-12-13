#!/usr/bin/env python

import os
import sys
import screed

def main():
    if len(sys.argv) != 3:
        mes = '[INFO] Usage: python {} <seqs.blastn.anicalc.tsv> <seqs.fna>\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    tabfile = sys.argv[1]
    seqfile = sys.argv[2]
    if tabfile == '-' and seqfile == '-':
        mes = '[ERROR] <seqs.blastn.anicalc.tsv> and <seqs.fna> can not be both from stdin\n'
        sys.stderr.write(mes)
    elif tabfile == '-':
        tabfile = '/dev/stdin'
    elif seqfile == '-':
        tabfile = '/dev/stdin'

    st = set()
    with open(tabfile) as fh1:
        for n, line in enumerate(fh1):
            sys.stdout.write(line)
            _line = line.rstrip()
            name1, name2, _ = _line.split('\t', 2)
            tu = tuple(sorted([name1, name2]))
            if tu in st:
                continue
            else:
                st.add(tu)

    cnt = 0
    for rec1 in screed.open(seqfile):
        name1 = rec1.name.split(None, 1)[0]
        for rec2 in screed.open(seqfile):
            name2 = rec2.name.split(None, 1)[0]
            if name1 == name2:
                continue
            tu = tuple(sorted([name1, name2]))
            if tu in st:
                continue
            else:
                st.add(tu)
                cnt += 1
                mes = f'{name1}\t{name2}\t0\t0\t0\t0\n'
                sys.stdout.write(mes)

    sys.stderr.write(f'[INFO] unaligned pairs added from {os.path.basename(seqfile)}: {cnt}\n')

if __name__ == '__main__':
    main()
