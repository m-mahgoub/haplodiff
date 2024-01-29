#!/usr/bin/env python
# This script processes RNA and methylation report data.

import argparse
import pandas as pd

def main(args):
    rna_report = args.rna_report
    meth_report = args.meth_report

    rna_table_output = args.rna_table_output
    meth_table_output = args.meth_table_output
    genes_table_output = args.genes_table_output

    rna_vaf_threshold = args.rna_vaf_threshold
    meth_pval_threshold = args.meth_pval_threshold
    rna_df = pd.read_csv(rna_report,sep='\t')
    rna_df = rna_df.drop_duplicates(subset=['chr', 'start', 'end'])
    rna_df = rna_df[rna_df['wgbs_vaf'] >= 0.3]
    rna_df = rna_df[((rna_df['alt_vaf'] >= 0.5 + rna_vaf_threshold) \
    | (rna_df['alt_vaf'] <= rna_vaf_threshold)) \
    & (rna_df['gene_id_exon'] != '.') \
    & (rna_df['chr'] != 'chrX') \
    & (rna_df['chr'] != 'chrY') ]
    rna_df = rna_df.copy()
    rna_df['pos'] = rna_df.apply(lambda x: f"{x['chr']}:{x['start']}-{x['end']}", axis=1)
    rna_biased_genes = rna_df[['pos','gene_id_exon', 'symbol_exon']]
    rna_biased_genes = rna_biased_genes.copy()
    rna_biased_genes['gene_id_exon'] = rna_biased_genes['gene_id_exon'].astype(float).astype(int)
    rna_biased_genes = rna_biased_genes.drop_duplicates().groupby(['gene_id_exon', 'symbol_exon'])['pos'].agg(';'.join).reset_index()
    rna_biased_genes.columns = ['gene_id', 'gene_symbol', 'rna_snp_position']
    meth_df = pd.read_csv(meth_report,sep='\t')
    meth_df  = meth_df .drop_duplicates(subset=['chr', 'start', 'end'])
    meth_df = meth_df[((meth_df['FEp+'] < meth_pval_threshold) \
    | (meth_df['FEp-'] < meth_pval_threshold)) \
    & (meth_df['annot.gene_id'] != '.') \
    & (meth_df['chr'] != 'chrX') \
    & (meth_df['chr'] != 'chrY') ]
    meth_df = meth_df.copy()
    meth_df['pos'] = meth_df.apply(lambda x: f"{x['chr']}:{x['start']}-{x['end']}", axis=1)
    meth_biased_genes = meth_df[['pos','annot.gene_id', 'annot.symbol']]
    meth_biased_genes = meth_biased_genes.copy()
    meth_biased_genes['annot.gene_id'] = meth_biased_genes['annot.gene_id'].astype(float).astype(int)
    meth_biased_genes = meth_biased_genes.drop_duplicates().groupby(['annot.gene_id', 'annot.symbol'])['pos'].agg(';'.join).reset_index()
    meth_biased_genes.columns = ['gene_id', 'gene_symbol', 'meth_snp_position']
    df = meth_biased_genes.merge(rna_biased_genes,how='inner')
    rna_df['gene_id_exon'] = rna_df['gene_id_exon'].astype(float).astype(int)
    meth_df['annot.gene_id'] = meth_df['annot.gene_id'].astype(float).astype(int)
    rna_df.to_csv(rna_table_output,sep='\t',header=True, index=False)
    meth_df.to_csv(meth_table_output,sep='\t',header=True, index=False)
    df.to_csv(genes_table_output,sep='\t',header=True, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process RNA and methylation report data.')
    parser.add_argument('--rna_report', required=True, help='Path to RNA report file')
    parser.add_argument('--meth_report', required=True, help='Path to methylation report file')
    parser.add_argument('--rna_table_output', required=True, help='Output file path for RNA table')
    parser.add_argument('--meth_table_output', required=True, help='Output file path for methylation table')
    parser.add_argument('--genes_table_output', required=True, help='Output file path for genes table')
    parser.add_argument('--rna_vaf_threshold', type=float, default=0.2, help='RNA VAF threshold (default: 0.2)')
    parser.add_argument('--meth_pval_threshold', type=float, default=0.01, help='Methylation p-value threshold (default: 0.01)')

    args = parser.parse_args()
    main(args)
