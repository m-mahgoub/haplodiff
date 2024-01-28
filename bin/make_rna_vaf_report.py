#!/usr/bin/env python
# This script processes VCF and BED files for RNA VAF analysis.

import argparse
import pandas as pd
from pybedtools import BedTool

def main(args):
    output = args.output
    vaf_tag = args.vaf_tag
    depth_tag = args.depth_tag
    rna_snpdep_txt = args.rna_snpdep_txt
    exons_bed_txt = args.exons_bed_txt
    dmrs_bed_txt = args.dmrs_bed_txt

    exon_bed = BedTool.from_dataframe(pd.read_csv(exons_bed_txt,sep='\t'))
    exon_col = list(pd.read_csv(exons_bed_txt,sep='\t').columns)
    rna_snps = pd.read_csv(rna_snpdep_txt,sep='\t',header=None)
    rna_snps[10] = rna_snps.apply(lambda x: float(x[9].split(':')[x[8].split(':').index(vaf_tag)]), axis =1)
    rna_snps[11] = rna_snps.apply(lambda x: int(x[9].split(':')[x[8].split(':').index(depth_tag)]), axis =1)
    rna_snps[12] = rna_snps.apply(lambda x: int(x[9].split(':')[x[8].split(':').index('RNA')].split(',')[0]), axis =1)
    rna_snps[13] = rna_snps.apply(lambda x: int(x[9].split(':')[x[8].split(':').index('RNA')].split(',')[1]), axis =1)
    rna_snps['start'] = rna_snps[1] - 1
    rna_snps = rna_snps[[0,'start',1,3,4,11,10,12,13]]
    rna_snps.columns = ['chr', 'start', 'end', 'ref', 'alt', 'wgbs_depth','wgbs_vaf', 'ref_rna_count', 'alt_rna_count']
    rna_snps = rna_snps.copy()
    rna_snps['alt_vaf'] = rna_snps['alt_rna_count']/(rna_snps['ref_rna_count']+rna_snps['alt_rna_count'])
    rna_snps_bed = BedTool.from_dataframe(rna_snps)
    rna_snps_bed_with_genes = rna_snps_bed.intersect(exon_bed, wo=True)
    df_rna = rna_snps_bed_with_genes.to_dataframe(disable_auto_names=True, header=None)
    df_rna.drop(columns=20,inplace=True)
    df_rna.columns = list(rna_snps.columns) +  [f'{col}_exon' for col in exon_col]
    df_rna.to_csv(output,sep='\t',header=True, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process VCF and BED files for RNA VAF analysis.')
    parser.add_argument('--output', required=True, help='Output file path')
    parser.add_argument('--rna_snpdep_txt', required=True, help='Path to RNA VCF SNPDEP file')
    parser.add_argument('--exons_bed_txt', required=True, help='Path to exons BED text file')
    parser.add_argument('--dmrs_bed_txt', required=True, help='Path to DMRs BED text file')
    parser.add_argument('--vaf_tag', default='AF1', help='VAF tag (default: AF1)')
    parser.add_argument('--depth_tag', default='DP', help='Depth tag (default: DP)')
    
    args = parser.parse_args()
    main(args)
