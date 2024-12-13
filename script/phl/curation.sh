#!/bin/bash

# get conservative boundary

perl -pe 's/^((.+?)-cat_[0-9]+_[0-9]+\t)/\2\t\1/' dramv-annotate/annotations.tsv | sed '1s/^/contig_id\tgene_id/' | csvtk join -L -t -f"contig_id,gene_position;contig_id,checkv_trimmed_orf_idx" - <(sed 's/||/__/' vs2sop.w_int.gene_feature_w_int.tsv | cut -f1,4-) > vs2_checkv.dramv.add_gene_feature.4curate.tsv

### Manual step
# from `vs2_checkv.dramv.add_gene_feature.4curate.tsv`, use the phage integrase (hmm_cat==2) and checkv viral gene (hmm_cat==1) as conservative boundaries to curate a prophage boundary file `vs2sop.conservative_boundary.tsv` w/ three columns: 1) contig, 2) start_orf, 3) end_orf
# those viral contigs w/o phage integrase or checkv viral gene should be removed
#

# cut the conservative region
python dram-cut-region.py --orf --trim-dramv-seqname <(sed 's/||/__/' vs2sop.conservative_boundary.tsv) dramv-annotate/annotations.tsv > vs2sop.w_int.conservative_boundary.dramv.tsv
date
