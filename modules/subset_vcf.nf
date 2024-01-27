process subset_vcf {
    container 'mdivr/basic-tools:v2'
    cpus 6
    memory '60 GB'

    input:
        tuple val(meta), path (filtered_vcf)
        path (exons_bed)
        path (promoters_bed)
        path (dmrs_bed)

    output:
        publishDir "${params.outdir}/${meta.id}/subset_vcf"
        tuple val(meta), path ("${meta.id}.meth_roi.vcf"), emit: meth_roi_vcf
        tuple val(meta), path ("${meta.id}.rna_roi.vcf"), emit: rna_roi_vcf

    script:
        """
        bedtools intersect -header -a ${meta.id}.filtered.vcf.gz -b ${exons_bed} -wa -u > ${meta.id}.meth_roi.vcf
        bedtools intersect -header -a ${meta.id}.filtered.vcf.gz -b ${promoters_bed} ${dmrs_bed} -wa -u > ${meta.id}.rna_roi.vcf
        """

}