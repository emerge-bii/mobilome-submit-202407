#!/usr/bin/env python

import sys
import os
import re


def parse_cdhit_aa(f):
    with open(f) as fh:
        mes = 'OTU\tmem\tlength\tiden\tis_repr\n'
        sys.stdout.write(mes)
        trigger = False
        clust_id = None
        for line in fh:
            if line.startswith('#'):
                continue
            line = line.rstrip()
            if not line:
                continue
            if line.startswith('>'):
                clust_id = line.split('>Cluster')[1].strip()
                trigger = True
                continue
            lst = line.split()
            if len(lst) == 5:
                _, length, mem, _, iden = lst
            elif len(lst) == 4:
                _, length, mem, iden = lst
                #print(lst)
            length = length.strip().strip('a,')
            mem = mem[1:] # remve '>'
            mem = mem.rsplit('...', 1)[0].strip() # remove '...' at the end
            is_repr = False if iden != '*' else True
            iden = iden.strip('%') if iden != '*' else 100

            assert clust_id != None
            mes = f'OTU{clust_id}\t{mem}\t{length}\t{iden}\t{is_repr}\n'
            sys.stdout.write(mes)

def main():
    if len(sys.argv) != 2:
        mes = '***Usage: python {} <file.cdhit.clstr>\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    clsr_f = sys.argv[1]
    if clsr_f == '-':
        clsr_f = '/dev/stdin'

    parse_cdhit_aa(clsr_f)

if __name__ == '__main__':
    main()
