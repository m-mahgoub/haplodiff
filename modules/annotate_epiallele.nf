process annotate_epiallele {
    container 'mdivr/quarto:v2'
    cpus 4
    memory '64 GB'

    input:
        tuple val(meta), path (epiallele_txt)
    output:
        publishDir "${params.outdir}/${meta.id}/epialleleR", mode: 'copy'
        tuple val(meta), path (epiallele_txt), path ("${meta.id}.hapmeth.annotated.txt"), emit: annotated_epiallele_txt

    script:
        """
        annotate_epiallele.R \\
        --input ${epiallele_txt} \\
        --output ${meta.id}.hapmeth.annotated.txt
        """
}