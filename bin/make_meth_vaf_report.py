#!/usr/bin/env python
# This script processes epiallele methylation data.

import argparse
import pandas as pd

def main(args):
    epiallele_meth = args.epiallele_meth
    epiallele_annotation = args.epiallele_annotation
    output = args.output

    meth_snps = pd.read_csv(epiallele_meth,sep='\t')
    meth_snps['start'] = meth_snps['range'] - 1
    meth_snps = meth_snps[['seqnames', 'start', 'range', 'name'] + list(meth_snps.columns[3:-1])]
    meth_snps.columns =['chr', 'start', 'end'] + list(meth_snps.columns[3:])
    meth_snps_annotations = pd.read_csv(epiallele_annotation,sep='\t')
    meth_snps_annotations = meth_snps_annotations[['seqnames', 'start', 'end', 'name', 'width', 'strand'] + list(meth_snps_annotations.columns[6:-1])]
    meth_snps_annotations.columns =['chr', 'start', 'end'] + list(meth_snps_annotations.columns[3:])
    meth_snps_annotations = meth_snps_annotations.dropna(subset='annot.gene_id')
    meth_snps_annotations = meth_snps_annotations.copy()
    meth_snps_annotations['start'] = meth_snps_annotations['start'] - 1
    df_meth = meth_snps.merge(meth_snps_annotations,how='right')
    # Columns to fill NA with 0
    meth_count_columns = [
    'M+Ref',
    'U+Ref',
    'M-Ref',
    'U-Ref',
    'M+Alt',
    'U+Alt',
    'M-Alt',
    'U-Alt'
    ]
    df_meth = df_meth.dropna(subset=meth_count_columns,how='all')
    df_meth.to_csv(output,sep='\t',header=True, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process epiallele methylation data.')
    parser.add_argument('--epiallele_meth', required=True, help='Path to epiallele methylation data file')
    parser.add_argument('--epiallele_annotation', required=True, help='Path to epiallele annotation file')
    parser.add_argument('--output', required=True, help='Output file path')

    args = parser.parse_args()
    main(args)
