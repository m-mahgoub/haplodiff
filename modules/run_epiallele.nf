process run_epiallele {
    container 'mdivr/quarto:v2'
    cpus 4
    memory '120 GB'

    input:
        tuple val(meta), path (meth_sorted_bam), path (meth_roi_vcf)
    output:
        publishDir "${params.outdir}/${meta.id}/epialleleR", mode: 'copy'
        tuple val(meta), path ("${meta.id}.hapmeth.txt"), emit: epiallele_txt

    script:
        """
        epiallele.R \\
        --input ${meth_sorted_bam} \\
        --genome ${params.meth_fasta} \\
        --vcf ${meth_roi_vcf} \\
        --output ${meta.id}.hapmeth.txt
        """

}