#!/usr/bin/env python

import os
import sys

def main():
    if len(sys.argv) != 3:
        mes = '*** Usage: python {} <table.tsv> <name.list>\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    tabfile = sys.argv[1]
    listfile = sys.argv[2]
    if tabfile == '-' and listfile == '-':
        mes = '[Error] <table.tsv> and <name.list> can not be both from stdin\n'
        sys.stderr.write(mes)
        sys.exit(1)
    elif tabfile == '-':
        tabfile = '/dev/stdin'
    elif listfile == '-':
        listfile = '/dev/stdin'
    else:
        pass

    st = set()
    with open(tabfile) as fh1, open(listfile) as fh2:
        for line in fh2:
            if line.startswith('#'):
                continue
            line = line.strip()
            st.add(line)

        for n, line in enumerate(fh1):
            if n == 0:
                sys.stdout.write(line)
                continue
            _line = line.rstrip()
            name = _line.split('\t', 1)[0]
            if not name in st:
                sys.stdout.write(line)
 
if __name__ == '__main__':
    main()
