#!/usr/bin/env python
import sys
import os
import resource
import platform
import pandas as pd

COLNAMES = ['qname', 'qlen', 'qstart', 'qend', 'strand', 'tname', 'tlen', 'tstart', 'tend', 'match_len', 'ali_len', 'mapq']
MIN_COV_FRAC = 0.5
MAX_OL_FRAC = 0.1  # over the whole target gene
MAX_GAP_FRAC = 0.1   # the gap between the alignments of up- dand down-stream to target gene should not be too big, ideally they should be continuous parts of the target gene
MAX_OL_FRAC_OF_ALI = 0.1  # over the shorter ali of up- and down-stream

def get_peak_mem():
    mem_self = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    mem_children = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    peak_mem = mem_self + mem_children
    if platform.system() == 'Linux':
        peak_mem_gb = peak_mem / float(1e6)
    else:
        peak_mem_gb = peak_mem / float(1e9)
    return peak_mem_gb

def get_nh_stats(up_row_ser, down_row_ser):

    up_qlen = up_row_ser.loc['qlen']
    up_qstrand = up_row_ser.loc['strand']
    up_qstart = up_row_ser.loc['qstart']
    up_qend = up_row_ser.loc['qend']

    down_qlen = down_row_ser.loc['qlen']
    down_qstrand = down_row_ser.loc['strand']
    down_qstart = down_row_ser.loc['qstart']
    down_qend = down_row_ser.loc['qend']

    tlen = up_row_ser.loc['tlen']
    up_tstart = up_row_ser.loc['tstart']
    up_tend = up_row_ser.loc['tend']
    down_tstart = down_row_ser.loc['tstart']
    down_tend = down_row_ser.loc['tend']

    up_match_len = up_row_ser.loc['match_len']
    up_ali_len = up_row_ser.loc['ali_len']
    down_match_len = down_row_ser.loc['match_len']
    down_ali_len = down_row_ser.loc['ali_len']

    st_up = set(range(up_tstart, up_tend))
    st_down = set(range(down_tstart, down_tend))
    covered = len(st_up | st_down)
    overlap = len(st_up & st_down)
    if overlap == 0:
        _, gap_start, gap_end, _ = sorted([up_tstart, up_tend, down_tstart, down_tend])
        # gap_start is 1 bp before the actual gap starting point
        # gap_end is 1 bp after the actual gap ending point
        gap = gap_end - 1 - gap_start
    else:
        gap = 0

    # filter
    cov_frac = covered / tlen
    ol_frac = overlap / tlen
    ol_frac_of_ali = overlap / min(up_ali_len, down_ali_len)
    gap_frac = gap / tlen
    identity = (up_match_len + down_match_len) / (up_ali_len + down_ali_len)

    return cov_frac, ol_frac, identity, gap_frac, ol_frac_of_ali, up_qstrand, tlen, up_tstart, up_tend, down_tstart, down_tend, up_qlen, up_qstart, up_qend, up_match_len, up_ali_len, down_qlen, down_qstart, down_qend, down_match_len, down_ali_len


tab_f = 'is_tn.gene_neighbor.ali2pc_rep_cds.paf'

df = pd.read_csv(tab_f, sep='\t', header=None, usecols=range(0,12))
df.columns = COLNAMES

tmp_df = df['qname'].str.rsplit('_', 1, expand=True)
tmp_df.columns = ['anchor_gene', 'side']
df = pd.concat([df, tmp_df], axis=1)

# headers
mes = 'anchor_gene\ttname\ttlen\tcov_frac\tol_frac\tidentity\tup_tstart\tup_tend\tdown_tstart\tdown_tend\tstrand\tup_qlen\tup_qstart\tup_qend\tup_match_len\tup_ali_len\tdown_qlen\tdown_qstart\tdown_qend\tdown_match_len\tdown_ali_len\n'
sys.stdout.write(mes)

for (anchor_gene, tname), sub_df in df.groupby(['anchor_gene', 'tname']):
    if len(set(sub_df['side'])) < 2:
        # not sure if minimap2 output secondary alignment for the same
        #   qname and tname pair; in case that it does,
        #   skip if one end is not aligned
        continue
    sel = sub_df['side'] == 'upstream'
    sub_df_up = sub_df.loc[sel,:]
    sel = sub_df['side'] == 'downstream'
    sub_df_down = sub_df.loc[sel,:]
    tlen = sub_df['tlen'].iloc[0]
    # simple case
    if len(sub_df_up) == 1 and len(sub_df_down) == 1:
        up_row_ser = sub_df_up.iloc[0, :]
        down_row_ser = sub_df_down.iloc[0, :]
        if up_row_ser.loc['strand'] != down_row_ser.loc['strand']:
            # make sure up and down match the same strand of target
            continue
        if up_row_ser.loc['strand'] == '-' and up_row_ser.loc['tstart'] <= down_row_ser.loc['tstart']:
            continue
        if up_row_ser.loc['strand'] == '+' and up_row_ser.loc['tstart'] >= down_row_ser.loc['tstart']:
            continue
        cov_frac, ol_frac, identity, gap_frac, ol_frac_of_ali, up_qstrand, tlen, up_tstart, up_tend, down_tstart, down_tend, up_qlen, up_qstart, up_qend, up_match_len, up_ali_len, down_qlen, down_qstart, down_qend, down_match_len, down_ali_len =  get_nh_stats(up_row_ser, down_row_ser)

    else:
        # there are multiple alignment in up or downstream
        #   need to check all combinations and find the one 
        #   w/ highest cov_frac and lowest ol_frac
        max_cov_frac = 0
        tmp_ol_frac = 9999999
        tmp_lst = None
        for i in range(len(sub_df_up)):
            for j in range(len(sub_df_down)):
                up_row_ser = sub_df_up.iloc[i,:]
                down_row_ser = sub_df_down.iloc[j,:]
                if up_row_ser.loc['strand'] != down_row_ser.loc['strand']:
                    # make sure up and down match the same strand of target
                    continue
                if up_row_ser.loc['strand'] == '-' and up_row_ser.loc['tstart'] <= down_row_ser.loc['tstart']:
                    continue
                if up_row_ser.loc['strand'] == '+' and up_row_ser.loc['tstart'] >= down_row_ser.loc['tstart']:
                    continue
                cov_frac, ol_frac, identity, gap_frac, ol_frac_of_ali, up_qstrand, tlen, up_tstart, up_tend, down_tstart, down_tend, up_qlen, up_qstart, up_qend, up_match_len, up_ali_len, down_qlen, down_qstart, down_qend, down_match_len, down_ali_len =  get_nh_stats(up_row_ser, down_row_ser)


                if cov_frac > max_cov_frac:
                    max_cov_frac = cov_frac
                    tmp_ol_frac = ol_frac
                    tmp_lst = [cov_frac, ol_frac, identity, gap_frac, ol_frac_of_ali, up_qstrand, tlen, up_tstart, up_tend, down_tstart, down_tend, up_qlen, up_qstart, up_qend, up_match_len, up_ali_len, down_qlen, down_qstart, down_qend, down_match_len, down_ali_len]
                elif cov_frac == max_cov_frac and ol_frac < tmp_ol_frac:
                    max_cov_frac = cov_frac
                    tmp_ol_frac = ol_frac
                    tmp_lst = [cov_frac, ol_frac, identity, gap_frac, ol_frac_of_ali, up_qstrand, tlen, up_tstart, up_tend, down_tstart, down_tend, up_qlen, up_qstart, up_qend, up_match_len, up_ali_len, down_qlen, down_qstart, down_qend, down_match_len, down_ali_len]
                else:
                    pass

        if tmp_lst != None:
            cov_frac, ol_frac, identity, gap_frac, ol_frac_of_ali, up_qstrand, tlen, up_tstart, up_tend, down_tstart, down_tend, up_qlen, up_qstart, up_qend, up_match_len, up_ali_len, down_qlen, down_qstart, down_qend, down_match_len, down_ali_len = tmp_lst
        else:
            continue

    # filter
    if cov_frac < MIN_COV_FRAC:
        continue
    if ol_frac > MAX_OL_FRAC:
        continue
    if gap_frac > MAX_GAP_FRAC:
        continue
    if ol_frac_of_ali > MAX_OL_FRAC_OF_ALI:
        continue

    mes = f'{anchor_gene}\t{tname}\t{tlen}\t{cov_frac}\t{ol_frac}\t{identity}\t{up_tstart}\t{up_tend}\t{down_tstart}\t{down_tend}\t{up_qstrand}\t{up_qlen}\t{up_qstart}\t{up_qend}\t{up_match_len}\t{up_ali_len}\t{down_qlen}\t{down_qstart}\t{down_qend}\t{down_match_len}\t{down_ali_len}\n'
    sys.stdout.write(mes)

peak_mem_gb = get_peak_mem()
bname= os.path.basename(sys.argv[0])
mes = f'[INFO] Peak mem for {bname}: {peak_mem_gb}GB\n'
sys.stderr.write(mes)
