process make_rna_report {
    container 'mdivr/basic-tools:v2'
    cpus 4
    memory '64 GB'

    input:
        tuple val(meta), path (rna_snpdep_txt)
    output:
        publishDir "${params.outdir}/${meta.id}/snpdep"
        tuple val(meta), path ("${meta.id}.rna_vaf_report.txt"), emit: rna_vaf_report

    script:
        """
        make_rna_vaf_report.py \\
        --output ${meta.id}.rna_vaf_report.txt \\
        --rna_snpdep_txt ${rna_snpdep_txt} \\
        --exons_bed_txt ${params.exons_bed_txt} \\
        --dmrs_bed_txt ${params.dmrs_bed_txt}
        """
}