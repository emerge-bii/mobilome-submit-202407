#!/usr/bin/env python

import time, resource, platform, sys, argparse, gzip

def parse_seqs(path):
    handle = gzip.open(path) if path.split('.')[-1] == 'gz' else open(path)
    id = next(handle).split()[0][1:]
    seq = ''
    for line in handle:
        if line[0] == '>':
            yield id, seq
            id = line.split()[0][1:]
            seq = ''
        else:
            seq += line.rstrip()
    yield id, seq
    handle.close()

def log_time(start):
    current_time = time.time()
    program_time = round(current_time - start, 2)
    peak_ram = round(max_mem_usage(), 2)
    sys.stderr.write("time: %s seconds, peak RAM: %s GB\n" % (program_time, peak_ram))

def max_mem_usage():
    """ Return max mem usage (Gb) of self and child processes """
    max_mem_self = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    max_mem_child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if platform.system() == 'Linux':
        return round((max_mem_self + max_mem_child)/float(1e6), 2)
    else:
        return round((max_mem_self + max_mem_child)/float(1e9), 2)

def parse_arguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        usage=argparse.SUPPRESS,
        description="Centroid based sequence clustering"
    )
    parser.add_argument('ani', type=str, metavar='ANI',
        help="""Path to tab-delimited file with fields: [qname, tname, ani]""")
    parser.add_argument('--exclude', type=str, metavar='PATH',
        help="""Path to list of sequence ids to exclude from clustering""")
    parser.add_argument('--keep', type=str, metavar='PATH',
        help="""Path to list of sequence ids to keep from clustering""")
    parser.add_argument('--min-ani', type=float, metavar='FLOAT', default=95,
        help="""Minimum average nucleotide identity (0...100, default=95)""")
    parser.add_argument('--one-mem-per-line', action='store_true',
        help="""Each member of cluster on one line with cols: mem, repr, clu_idx""")
    return vars(parser.parse_args())

# main
start = time.time()

# args
args = parse_arguments()

if args['ani'] == '-':
    args['ani'] = '/dev/stdin'

exclude = set([_.rstrip() for _ in open(args['exclude'])]) if args['exclude'] else None
keep = set([_.rstrip() for _ in open(args['keep'])]) if args['keep'] else None

# store edges
sys.stderr.write("storing edges...\n")
num_edges = 0
edges = {}
names = set()
handle = gzip.open(args['ani']) if args['ani'].split('.')[-1] == 'gz' else open(args['ani'])
for index, line in enumerate(handle):
    if line.startswith('qname\ttname\t'):
        continue
    qname, tname, ani = line.split()
    num_edges += 1
    if qname == tname:
        continue
    elif exclude and (qname in exclude or tname in exclude):
        continue
    elif keep and ((not qname in keep) or (not tname in keep)):
        continue

    names.add(qname)
    names.add(tname)
    if float(ani) < args['min_ani']:
        continue
    edges.setdefault(qname, []).append(tname)
    edges.setdefault(tname, []).append(qname)

#   if not num_edges % 1000000:
#       current_time = time.time()
#       program_time = round(current_time - start, 2)
#       peak_ram = round(max_mem_usage(), 2)
#       print("num edges: %s, seconds: %s, GB: %s" % (num_edges, program_time, peak_ram))
handle.close()
sys.stderr.write("%s edges (non-directional) retained from ani\n" % num_edges)
sys.stderr.write("%s edges (directional) currently stored\n" % sum([len(_) for _ in edges.values()]))
log_time(start)

# cluster
sys.stderr.write("clustering...\n")
clust_to_seqs = {}
seq_to_clust = {}
# loop over seqs in sorted order
for seq_id in sorted(names):
    # seq already been assigned; cant be centroid
    if seq_id in seq_to_clust:
        continue
    # seq is centroid for new cluster
    else:
        # add self to cluster
        clust_to_seqs[seq_id] = [seq_id]
        seq_to_clust[seq_id] = seq_id
        # update with cluster members
        for mem_id in edges.get(seq_id, []):
            if mem_id not in seq_to_clust:
                clust_to_seqs[seq_id].append(mem_id)
                seq_to_clust[mem_id] = seq_id
sys.stderr.write("%s total clusters\n" % len(clust_to_seqs))
log_time(start)

# write
sys.stderr.write("writing clusters...\n")
clu_id = 0
for seq_id, mem_ids in clust_to_seqs.items():
    if args['one_mem_per_line']:
        for mem_id in  mem_ids:
            mes = '{}\t{}\tnh{}\n'.format(mem_id, seq_id, clu_id)
            sys.stdout.write(mes)

    else:
        mes = '{}\t{}\tnh{}\n'.format(seq_id, ','.join(mem_ids), clu_id)
        sys.stdout.write(mes)
        
    clu_id += 1

log_time(start)

