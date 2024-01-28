process filter_vcf {
    container 'mdivr/basic-tools:v2'
    cpus 12
    memory '100 GB'

    input:
    tuple val(meta), path (vcf)
    // path (fastq1)
    // path (fastq2)

    output:
        publishDir "${params.outdir}/${meta.id}/filtered_vcf"
        tuple val(meta), path ("${meta.id}.filtered.*"), emit: filter_vcf

    script:
        """
        tabix --csi -p vcf ${vcf}
        modify_vcf_header.py --in_vcf ${vcf} --out_vcf ${meta.id}.header_modified.vcf
        bcftools view --threads ${task.cpus} -H ${vcf} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY >> ${meta.id}.header_modified.vcf
        bcftools sort -o ${meta.id}.sorted.for_pop.vcf ${meta.id}.header_modified.vcf
        bedtools intersect -sorted -header -a ${meta.id}.sorted.for_pop.vcf -b ${params.population_snps_vcf} -g ${params.vcf_genome} -wa -u > ${meta.id}.only_pop.vcf
        bgzip -c ${meta.id}.only_pop.vcf > ${meta.id}.only_pop.vcf.gz
        tabix -p vcf ${meta.id}.only_pop.vcf.gz
        bcftools view --threads ${task.cpus} -f 'PASS' -g 'het' --exclude 'ALT="N"' --exclude-types indels ${meta.id}.only_pop.vcf.gz chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY | bedtools intersect -header -a - -b ${params.cpgs_bed} -wa -v > ${meta.id}.filtered.vcf
        bgzip -c ${meta.id}.filtered.vcf > ${meta.id}.filtered.vcf.gz
        tabix -p vcf ${meta.id}.filtered.vcf.gz
        """

        stub:
        """
        cp ${params.stub_path}/filtered_vcf/${meta.id}.filtered.* .
        sleep 10
        """
}
