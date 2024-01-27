process make_meth_report {
    container 'mdivr/basic-tools:v2'
    cpus 4
    memory '64 GB'

    input:
        tuple val(meta), path (epiallele_txt), path (annotated_epiallele_txt)
    output:
        publishDir "${params.outdir}/${meta.id}/epialleleR"
        tuple val(meta), path ("${meta.id}.meth_vaf_report.txt"), emit: meth_vaf_report

    script:
        """
        make_meth_vaf_report.py \\
        --epiallele_meth ${epiallele_txt} \\
        --epiallele_annotation ${annotated_epiallele_txt} \\
        --output ${meta.id}.meth_vaf_report.txt
        """
}