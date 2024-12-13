#!/usr/bin/env python

'''
NOTE if vs2 circular, the first incomplete gene is cut and moved to the end in vs2 boundary
NOTE if vs2 circular and partial, the vs2 seq output is trimmed for provirus for sure
NOTE if vs2 (not partial), the vs2 seq output is NOT trimmed to remove overhangs in SOP (--keep-original-seq), but trimmed by default

Rule:
if circular and vs2 provirus, two copies of original seq need 
  to be concatenated for the trim_start_bp and trim_end_bp to apply

if vs2 partial and checkv provirus, use vs2 boundary and checkv boundary;
elif vs2 partial and not checkv provirus, use vs2 boundary;
otherwise, use checkv boundary only;
'''

import sys
import os
import screed
import pandas as pd

pd.options.display.max_rows = 1000
pd.options.display.max_columns = 10

PARTIAL = 'partial'
FULL = 'full'
LT2GENE = 'lt2gene'

python get-vs2sop-checkv-trimmed-gene-feature.py vs2-pass1/0d5-final-viral-boundary.tsv checkv/contamination.tsv checkv/tmp/gene_features.tsv $Phage_int_coord_f > vs2sop.add_gene_feature_w_int.tsv

vs2_bdy_f = 'vs2-pass1/0d5-final-viral-boundary.tsv'
contam_f = 'checkv/contamination.tsv'
gene_ftr_f = 'checkv/tmp/gene_features.tsv'
int_coord_f = '~/mge/recombinase/recombinase.dram.tsv'  # gene_id |start_bp |end_bp

df = pd.read_csv(vs2_bdy_f, sep='\t', header=0)
df.index = df['seqname_new']

df_contam = pd.read_csv(contam_f, sep='\t', header=0)
df_gene_ftr = pd.read_csv(gene_ftr_f, sep='\t', header=0)
df_int = pd.read_csv(int_coord_f, sep='\t', header=0)

df_int.columns = ['gene_id', 'start_bp', 'end_bp']
df_int['contig_id'] = df_int['gene_id'].map(lambda x: x.rsplit('_', 1)[0])
df_int['gene_idx'] = df_int['gene_id'].map(lambda x: int(x.rsplit('_', 1)[1]))

# TODO: 1) parse df_contam to find prophage or phage; 2) adjust to original boundary based on vs2_bdy; 3) check if phage int are within the boundary, and change int gene "type" in df_contam to value of 2; 4) define a conservative boundary using type value 1 and 2 within the checkv boundary

# add contig_id_ori (original contig id before vs2 adding suffix) for matching w/ names in df_int
df_contam['contig_id_ori'] = df_contam['contig_id'].map(lambda x: x.rsplit('||', 1)[0])

lst = []
for i in range(len(df_contam)):
    row =  df_contam.iloc[i, :]
    contig_id = row.loc['contig_id']
    contig_id_ori = row.loc['contig_id_ori']
    contig_len = row.loc['contig_length']
    total_gene_cnt = row.loc['total_genes']
    viral_gene_cnt = row.loc['viral_genes']
    host_gene_cnt = row.loc['host_genes']
    has_provirus = row.loc['provirus']

    # no viral gene; can not form a region even after adding int
    #if viral_gene_cnt == 0:
    #    continue

    sel = df_gene_ftr['contig_id'] == contig_id
    df_gene_ftr_sub = df_gene_ftr.loc[sel, :].copy()

    # scan phage int
    # make sure "SS.fna.0" are consistent for df_int and df_contam
    sel = df_int['contig_id'] == contig_id_ori
    df_int_sub = df_int.loc[sel, :].copy()

    coord_lst = []
    if has_provirus == 'No':
        viral_region_bp = (1, contig_len)
        viral_region_orf = (1, total_gene_cnt)
        cnt = 1  # viral region index in contig
        coord_lst.append((cnt, viral_region_bp, viral_region_orf))
    else:
        region_type_lst = row.loc['region_types'].split(',')
        region_bp_lst = row.loc['region_coords_bp'].split(',')
        region_orf_lst = row.loc['region_coords_genes'].split(',')
        
        cnt = 0
        for type, range_bp, range_orf in zip(region_type_lst, region_bp_lst, region_orf_lst):
            if type == 'host':
                continue
            cnt += 1
            viral_region_bp = list(map(int, range_bp.split('-')))
            viral_region_orf = list(map(int, range_orf.split('-')))
            coord_lst.append((cnt, viral_region_bp, viral_region_orf))

        assert cnt != 0, f'[ERROR] no viral region found in {contig_id}'


    ### get prophage start postion in vs2 boundary
    name = contig_id
    ori_name = name.rsplit('||', 1)[0]

    if PARTIAL in name:
        # need shape if circular and vs2 provirus; 
        #   concatenated seq of two copy
        shape = df.loc[name, 'shape'] 
        vs2_start = df.loc[name, 'trim_bp_start']
        vs2_end = df.loc[name, 'trim_bp_end']
        vs2_orf_start = df.loc[name, 'trim_orf_index_start']
        vs2_orf_end = df.loc[name, 'trim_orf_index_end']
        # need to get original contig length
        #  the ones in boundary files are two dup concatenated for circular shape
        #  need to add it in the vs2 boundary output; skip for now
    elif FULL in name:
        shape = df.loc[name, 'shape'] 
        vs2_start = 1
        vs2_end = 'NA'
        vs2_orf_start = 1
        vs2_orf_end = 'NA'
        # length here equals to original contig length since
        #  not trimmed is done with --keep-original-seq
    elif LT2GENE in name:
        shape = 'lt2gene'
        vs2_start = 1
        vs2_end = 'NA'
        vs2_orf_start = 1
        vs2_orf_end = 'NA'
        # length here equals to original contig length since
        #  not trimmed is done with --keep-original-seq
    else:
        raise f'[ERROR] {name} does NOT look like seqname of checkv step in VS2 SOP..\n'

    # change phage integrase "hmm_cat" to 2
    df_int_sub['checkv_gene_idx'] = df_int_sub['gene_idx'] - vs2_orf_start + 1 
    sel = df_gene_ftr_sub['gene_num'].isin(set(df_int_sub['checkv_gene_idx']))
    df_gene_ftr_sub.loc[sel, 'hmm_cat'] = 2


    #print(contig_id)
    #print(coord_lst)
    for _l in coord_lst:
        cnt = _l[0]
        viral_region_bp = _l[1]
        viral_region_orf = _l[2]
        #print(_l)
        #print(viral_region_bp)
        #print('aaa')
        checkv_bp_start, checkv_bp_end = viral_region_bp
        checkv_orf_start, checkv_orf_end = viral_region_orf

        ### check viral region position on original seq
        checkv_bp_start_oriseq = vs2_start + checkv_bp_start - 1
        checkv_bp_end_oriseq = vs2_start + checkv_bp_end - 1
        checkv_orf_start_oriseq = vs2_orf_start + checkv_orf_start - 1
        checkv_orf_end_oriseq = vs2_orf_start + checkv_orf_end - 1

        ### pick the viral region defined by checkv
        sel1 = df_gene_ftr_sub['gene_num'] >= checkv_orf_start
        sel2 = df_gene_ftr_sub['gene_num'] <= checkv_orf_end
        sel = (sel1 & sel2)
        df_gene_ftr_sub_checkv = df_gene_ftr_sub.loc[sel, :].copy()

        df_gene_ftr_sub_checkv['vs2_trimmed_orf_idx'] = df_gene_ftr_sub_checkv['gene_num']
        df_gene_ftr_sub_checkv['checkv_trimmed_orf_idx'] = df_gene_ftr_sub_checkv['gene_num'] - checkv_orf_start + 1
        df_gene_ftr_sub_checkv['original_orf_idx'] = df_gene_ftr_sub_checkv['gene_num'] + vs2_orf_start - 1

        sel_col_lst = ['contig_id', 'original_orf_idx', 'vs2_trimmed_orf_idx', 'checkv_trimmed_orf_idx', 'hmm_cat', 'hmm_db', 'hmm_name']
        lst.append(df_gene_ftr_sub_checkv.loc[:, sel_col_lst])


    ## end for the loop for df_contam
    #break

df_merged = pd.concat(lst, axis=0)
df_merged.to_csv('vs2sop.add_gene_feature_w_int.tsv', sep='\t', header=True, index=False, na_rep='NA')

