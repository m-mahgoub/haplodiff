process run_snpdep {
    container 'mdivr/snpdep:v3'
    cpus 8
    memory '64 GB'

    input:
        tuple val(meta), path (rna_bam), path (rna_bam_index), path (rna_roi_vcf)
    output:
        publishDir "${params.outdir}/${meta.id}/snpdep"
        tuple val(meta), path ("${meta.id}.rna_roi.annotated.vcf"), emit: snpdep_vcf
        tuple val(meta), path ("${meta.id}.rna_roi.annotated.txt"), emit: snpdep_txt

    script:
        def bam_format_args = rna_bam.getName().endsWith('.cram') ? 
            "--reads-format cram --reference ${params.rna_fasta}" : 
            ""
        """
        snpdep \\
        --format-field-id RNA \\
        --format-field-name RNAseq \\
        --output ${meta.id}.rna_roi.annotated.vcf \\
        --threads ${task.cpus} \\
        --min-count ${params.snpdep_min_count} \\
        ${bam_format_args} \\
        ${rna_roi_vcf} \\
        ${rna_bam}

        bcftools view -H ${meta.id}.rna_roi.annotated.vcf | grep ":RNA" > ${meta.id}.rna_roi.annotated.txt
        """

}