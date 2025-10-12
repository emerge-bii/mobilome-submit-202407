#!/usr/bin/env python

import os
import sys

def main():
    if len(sys.argv) != 2:
        mes = '[INFO] Usage: python {} <table.blast.anicalc.tsv>\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    tabfile = sys.argv[1]
    if tabfile == '-':
        tabfile = '/dev/stdin'

    st = set()
    with open(tabfile) as fh1:
        for n, line in enumerate(fh1):
            if n == 0:
                sys.stdout.write(line)
                continue
            _line = line.rstrip()
            name1, name2, _ = _line.split('\t', 2)
            tu = tuple(sorted([name1, name2]))
            if tu in st:
                continue
            else:
                st.add(tu)
                sys.stdout.write(line)
    
if __name__ == '__main__':
    main()
