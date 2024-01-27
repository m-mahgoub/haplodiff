#!/usr/bin/env python

import argparse
from pysam import VariantFile

# Set up argument parser
parser = argparse.ArgumentParser(description='Modify VCF file header.')
parser.add_argument('--in_vcf', help='Input VCF file path', required=True)
parser.add_argument('--out_vcf', help='Output VCF file path', required=True)
args = parser.parse_args()

# Use arguments
vcf_in = VariantFile(args.in_vcf)
header_str = str(vcf_in.header)

chr_list = [f'chr{num}' for num in range(1,23)] + ['chrX', 'chrY']
contig_assembled = {}
for chr in chr_list:
    contig_assembled[chr] = None

for line in header_str.split('\n'):
    if line.startswith('##contig='):
        contig = line.split('ID=')[1].split(':')[0].split(',')[0]
        if contig in chr_list:
            contig_assembled[contig] = line

header_str_new = ''
contig_added = False
for line in header_str.split('\n'):
    if line.startswith('##contig='):
        if contig_added == False:
            header_str_new += '\n'.join(list(contig_assembled.values())) + '\n'
            contig_added = True
    if line.startswith('##contig=') == False:
        if line.startswith('#CHROM'):
            header_str_new += line
        else:
            header_str_new += line + '\n'

            
with open(args.out_vcf, 'w') as f:
    f.write(header_str_new)

# Append the non-header lines from the input file to the output file
# ! bcftools view -H {args.in_vcf} >> {args.out_vcf}
