import sys
import os
import pandas
import csv


def main():
    if len(sys.argv) != 4:
        mes = (
                '*** Usage: python {} '
                '<vs2-score.tsv> <checkv-contamination.tsv> <out.tsv>\n'
        )
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    diamond_f = sys.argv[1] # vs2
    ncbi_f = sys.argv[2]    # checkv
    out_f = sys.argv[3]

    vs2_id_col = 'seqname'
    checkv_id_col = 'contig_id'
    df1 = pandas.read_table(diamond_f, header=0)
    df2 = pandas.read_table(ncbi_f, header=0)

    df1 = df1[pandas.notnull(df1.loc[:,vs2_id_col])]
    df2 = df2[pandas.notnull(df2.loc[:,checkv_id_col])]

    df1.loc[:,vs2_id_col] = [str(i).strip() for i in df1.loc[:,vs2_id_col].tolist()]
    df2.loc[:,checkv_id_col] = [str(i).strip() for i in df2.loc[:,checkv_id_col].tolist()]

    ### check accession # in two files
    ncbi_accs = df2.loc[:,checkv_id_col]
    diamond_accs = df1.loc[:,vs2_id_col]
    ncbi_st = set(ncbi_accs)
    diamond_st = set(diamond_accs)

    if diamond_st != ncbi_st:
        sys.stderr.write('*** Unique accessions in {}: {}\n'.format(os.path.basename(diamond_f), repr(diamond_st.difference(ncbi_st))))
        sys.stderr.write('*** Unique accessions in {}: {}\n'.format(os.path.basename(ncbi_f), repr(ncbi_st.difference(diamond_st))))

    ### merge
    ncbi_other_cols = list(df2.columns)
    ncbi_other_cols.remove(checkv_id_col)
    df_list = [df1,]
    for col in ncbi_other_cols:
        l = df2.loc[:,col].tolist()
        d = dict(zip(ncbi_accs, l))
        ser = diamond_accs.map(d)
        ser.rename(col, inplace=True)
        df_list.append(ser)

    df_merged = pandas.concat(df_list, axis=1)
    df_merged.to_csv(out_f, header=True, index=False, sep='\t', na_rep='NA', quoting=csv.QUOTE_NONE)

if __name__ == '__main__':
    main()
