#!/usr/bin/env python

import sys
import os
import screed


def main():
    if len(sys.argv) != 3:
        mes = '[INFO] Usage: python {} <seqfile> size\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    seqfile = sys.argv[1]
    size = int(sys.argv[2])

    with screed.open(seqfile) as sh:
        for rec in sh:
            _l = rec.name.split(None, 1)
            name = _l[0]
            if len(_l) > 1:
                desc = _l[1]
            else:
                desc = ''
            seq = rec.sequence
            contig_length = len(seq)
            for i in range(0, contig_length, size):
                start = i
                end = i+size
                if end > contig_length:
                    end = contig_length
                subseq = seq[start:end]
                mes = f'>{name}__{start+1}to{end+1}  {desc}\n{subseq}\n'
                sys.stdout.write(mes)

if __name__ == '__main__':
    main()
