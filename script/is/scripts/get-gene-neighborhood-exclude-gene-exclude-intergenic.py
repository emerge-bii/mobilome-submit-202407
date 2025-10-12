#!/usr/bin/env/python

import os
import sys
import platform
import resource
import logging
import argparse
import screed
import pandas as pd

def get_mem_usage():
    max_mem_self = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    max_mem_child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    max_mem = max_mem_self + max_mem_child
    if platform.system() == 'Linux':
        max_mem_gb = max_mem / float(1e6)
    else:
        max_mem_gb = max_mem / float(1e9)
    return max_mem_gb

def get_logger(quiet=False):
    logger = logging.getLogger(__name__)
    if not quiet:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)
    formatter = logging.Formatter(fmt="[INFO] %(asctime)s %(message)s", datefmt='%Y-%m-%d %H:%M:%S')
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.handlers.clear()
    logger.addHandler(stream_handler)
    return logger

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--include-anchor-gene', dest='include_anchor_gene', action='store_true',
            help='exclude the anchor gene sequence in neighborhood')
    parser.add_argument('--exclude-intergenic', dest='exclude_intergenic',
            action='store_true', 
            help="exclude the intergenic region near anchor gene")
    parser.add_argument('--exclude-extra-gene', dest='exclude_extra_gene',
            type=int,
            help="exclude extra gene upstream and downstream of anchor gene")
    parser.add_argument('--exclude-extra-bp', dest='exclude_extra_bp',
            type=int,
            help="exclude extra bp upstream and downstream of anchor gene, could be used on top of --exclude-intergenic and --exclude-extra-gene")
    parser.add_argument('--radius', dest='radius',
            default=2000,
            type=int,
            help="take radius upstream and downstream, could be used on top of --exclude-intergenic and --exclude-extra-bp")
    parser.add_argument('--min-oneside-len', dest='min_oneside_len',
            default=100,
            type=int,
            help="minimal length required for upstream or downstream seq")
    parser.add_argument('--min-twoside-len', dest='min_twoside_len',
            default=400,
            type=int,
            help="minimal length required for the sum of upstream and downstream seq")
    parser.add_argument('--fixed-length', dest='fixed_length',
            action='store_true',
            help="only take fixed radius bp upstream and downstream, discard if shorter and trim if longer, could be used on top of --exclude-intergenic and --exclude-extra-bp")
    parser.add_argument('--sep-ends', dest='sep_ends',
            action='store_true',
            help="separate two ends as two fasta records, instead of concatenating them")
    parser.add_argument('dramtab')
    parser.add_argument('namelist')
    parser.add_argument('contig')
    return vars(parser.parse_args())

def main():
    args = parse_arguments()
    if args['sep_ends'] and args['include_anchor_gene']:
        mes = '[ERROR] --include-anchor-gene and --sep-ends can not be both used\n'
        sys.stderr.write(mes)
        sys.exit(1)

    tabfile = args['dramtab']
    listfile = args['namelist']
    contigfile = args['contig']
    radius = args['radius']

    # minimal length required for neighborhood at one side of mobile gene
    #   or else this side is discarded
    MIN_ONESIDE_LEN = args['min_oneside_len']
    # minimal length required for neighborhood combined from two sides
    MIN_NH_LEN = args['min_twoside_len']

    logger = get_logger(quiet=False)

    if tabfile == '-' and listfile == '-':
        mes = '[Error] <dram.tsv> and <name.list> can not be both from stdin\n'
        sys.stderr.write(mes)
        sys.exit(1)
    elif tabfile == '-':
        tabfile = '/dev/stdin'
    elif listfile == '-':
        listfile = '/dev/stdin'
    else:
        pass

    st = set()
    d = {}
    with open(listfile) as fh2:
        for line in fh2:
            if line.startswith('#'):
                continue
            line = line.strip()
            st.add(line)

    df1 = pd.read_csv(tabfile, sep='\t', header=0, index_col=0, dtype={'pfam_hits': str, 'cazy_hits': str})
    gb = df1.groupby(['fasta', 'scaffold'])
    for _, df in gb:
        if args['exclude_intergenic']:
            for i, name in enumerate(df.index):
                if not name in st:
                    continue
                if i == 0:
                    #start = df.iloc[i, 'start_position']
                    # start extend to 0 to make sure make sure possible intergenic region 
                    #  at the beginning is excluded
                    start = 1
                else:
                    start = df['end_position'].iloc[i-1] + 1

                last_cds = False
                if i == len(df) - 1:
                    end = df['end_position'].iloc[i]
                    last_cds = True
                else:
                    end = df['start_position'].iloc[i+1] - 1

                contig_name, ind = name.rsplit('_', 1)
                d.setdefault(contig_name, []).append((ind, start, end, last_cds))

        elif args['exclude_extra_gene']:
            for i, name in enumerate(df.index):
                if not name in st:
                    continue
                start_gene_ind =  i - args['exclude_extra_gene']
                if start_gene_ind < 0:
                    #start = df.iloc[i, 'start_position']
                    # start extend to 0 to make sure make sure possible intergenic region 
                    #  at the beginning is excluded
                    start = 1
                else:
                    start = df['start_position'].iloc[start_gene_ind]

                last_cds = False
                end_gene_ind = i + args['exclude_extra_gene']
                if end_gene_ind > len(df) - 1:
                    end = df['end_position'].iloc[len(df)-1]
                    last_cds = True
                else:
                    end = df['end_position'].iloc[end_gene_ind]

                contig_name, ind = name.rsplit('_', 1)
                d.setdefault(contig_name, []).append((ind, start, end, last_cds))
        else:
            for i, name in enumerate(df.index):
                if not name in st:
                    continue
                start = df['start_position'].iloc[i]
                end = df['end_position'].iloc[i]
                contig_name, ind = name.rsplit('_', 1)
                d.setdefault(contig_name, []).append((ind, start, end, False))

    mes = 'Anchor gene position loaded in dict'
    logger.info(mes)
    max_mem = get_mem_usage()
    logger.info(f'Max mem usage: {max_mem}GB')

    with screed.open(contigfile) as sh:
        cnt = 0
        total = 0
        for rec in sh:
            header = rec.name
            name = header.split(None, 1)[0]
            if name in d:
                for ind, start, end, last_cds in d[name]:
                    total += 1
                    # convert to python index (starting from 0)
                    start = start - 1 
                    end = end

                    seq = rec.sequence
                    length = len(seq)

                    if args['exclude_intergenic'] or args['exclude_extra_gene']:
                        if last_cds:
                            end = length

                    if args['exclude_extra_bp']:
                        start = start - args['exclude_extra_bp']
                        end = end + args['exclude_extra_bp']

                    # expand to neighborhood
                    nh_start = start - radius
                    nh_end = end + radius
                    if nh_start < 0:
                        nh_start = 0

                    if nh_end > length:
                        nh_end = length 

                    upstream = seq[nh_start:start] # if start is 0, nh_start should 0 too 
                    downstream = seq[end:nh_end]   # if end is length, nh_end should length too

                    if args['fixed_length']:
                        if len(upstream) < radius:
                            err_mes = f'[INFO] {name} has upstream neighborhood < {radius}bp, discarded..\n'
                            sys.stderr.write(err_mes)
                            cnt += 1
                            continue

                        if len(downstream) < radius:
                            err_mes = f'[INFO] {name} has downstream neighborhood < {radius}bp, discarded..\n'
                            sys.stderr.write(err_mes)
                            cnt += 1
                            continue

                        nh = f'{upstream}NNN{downstream}'

                    else:
                        if len(upstream) < MIN_ONESIDE_LEN:
                            upstream = ''
                        if len(downstream) < MIN_ONESIDE_LEN:
                            downstream = ''

                        nh = f'{upstream}NNN{downstream}'
                        if len(nh) - 3 < MIN_NH_LEN:
                            cnt += 1
                            err_mes = f'[INFO] {name} has neighborhood < {MIN_NH_LEN}bp, discarded..\n'
                            sys.stderr.write(err_mes)
                            continue

                    nh = f'{upstream}NNN{downstream}'
                    if args['sep_ends']:
                        mes = f'>{name}_{ind}_upstream  {nh_start+1}|{start}\n{upstream}\n>{name}_{ind}_downstream  {end+1}|{nh_end}\n{downstream}\n'
                    else:
                        mes = f'>{name}_{ind}  {nh_start}|{start}|{end}|{nh_end}\n{nh}\n'
                    sys.stdout.write(mes)

        if args['fixed_length']:
            err_mes = ('[INFO] total sequences discarded due to '
                    f'up/down-stream < {radius} bp: {cnt}')
            err_mes2 = f'[INFO] total sequences screended for neighborhood: {total}'
        else:
            err_mes = (f'Total sequences discarded due to < {MIN_NH_LEN} '
                        f'neighorhood length: {cnt}')
            err_mes2 = f'[INFO] total sequences screended for neighborhood: {total}'
        logger.info(err_mes)
        logger.info(err_mes2)

    
    mes = 'Neighborhood extraction finished'
    logger.info(mes)
    max_mem = get_mem_usage()
    logger.info(f'Max mem usage: {max_mem}GB')

if __name__ == '__main__':
    main()
