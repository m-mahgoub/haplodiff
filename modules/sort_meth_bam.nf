process sort_meth_bam {
    container 'mdivr/nextflow:20230925'
    cpus 6
    memory '60 GB'

    input:
        tuple val(meta), path (meth_bam), path (meth_bam_index)

    output:
        publishDir "${params.outdir}/${meta.id}/sorted_meth_bam"
        tuple val(meta), path ("${meta.id}.nsorted.bam"), emit: meth_sorted_bam


    script:
        """
        mkdir -p tmp
        java -jar /storage1/fs1/dspencer/Active/spencerlab/mohamed/apps/picard/picard.jar SortSam \\
        I=${meth_bam} \\
        O=${meta.id}.nsorted.bam \\
        R=${params.meth_fasta} \\
        SORT_ORDER=queryname \\
        TMP_DIR=tmp/
        """

}