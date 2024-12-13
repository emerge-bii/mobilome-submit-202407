#!/usr/bin/env python

import sys
import os
import screed


def main():
    '''
    pick seqs that are in list file
    '''
    if len(sys.argv) != 3:
        mes = ('*** Usage: python {} <seqfile.fa> <prot-hits.list>\n')
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    seqfile = sys.argv[1]
    list_f = sys.argv[2]
    if seqfile == '-' and list_f == '-':
        mes = ('*** Error: only one of <seqfile.fa> and <file.list> '
                'can be "-"\n')
        sys.stderr.write(mes)
        sys.exit(1)

    if seqfile == '-':
        seqfile = '/dev/stdin'

    if list_f == '-':
        list_f = '/dev/stdin'

    st = set()
    with open(list_f) as fp:
        for line in fp:
            if line.startswith('#'):
                continue
            line = line.rstrip()
            st.add(line)

    with screed.open(seqfile) as sp:
        for rec in sp:
            name = rec.name.split(None, 1)[0]
            if not name in st:
                continue
            seq = rec.sequence
            sys.stdout.write('>{}\n{}\n'.format(rec.name, seq))

if __name__ == '__main__':
    main()
