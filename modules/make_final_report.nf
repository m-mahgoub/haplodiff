process make_final_report {
    container 'mdivr/basic-tools:v2'
    cpus 4
    memory '64 GB'

    input:
        tuple val(meta), path (rna_vaf_report), path (meth_vaf_report)
    output:
        publishDir "${params.outdir}/${meta.id}", mode: 'copy'
        tuple val(meta), path ("${meta.id}*_bias.txt"), emit: final_reports

    script:
        """
        make_final_report.py \\
        --rna_vaf_threshold ${params.rna_vaf_threshold} \\
        --meth_pval_threshold ${params.meth_pval_threshold} \\
        --rna_report ${rna_vaf_report} \\
        --meth_report ${meth_vaf_report} \\
        --rna_table_output ${meta.id}.rna_bias.txt \\
        --meth_table_output ${meta.id}.meth_bias.txt \\
        --genes_table_output ${meta.id}.genes_with_meth_rna_bias.txt
        """
}