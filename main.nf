

include { filter_vcf            } from './modules/filter_vcf'
include { subset_vcf            } from './modules/subset_vcf'
include { sort_meth_bam         } from './modules/sort_meth_bam'
include { run_epiallele         } from './modules/run_epiallele'
include { run_snpdep            } from './modules/run_snpdep'
include { make_rna_report       } from './modules/make_rna_report'
include { annotate_epiallele    } from './modules/annotate_epiallele'
include { make_meth_report      } from './modules/make_meth_report'
include { make_final_report     } from './modules/make_final_report'

Channel.fromPath(params.input)
    .splitCsv(header:true, sep:',')
    .map { row ->
        // Determine the index file extension based on the meth_bam file extension
        def index_extension = row.meth_bam.endsWith('.bam') ? '.bai' : row.meth_bam.endsWith('.cram') ? '.crai' : ''
        def meth_bam_index = file(row.meth_bam + index_extension, checkIfExists: true)

        [ ['id'       : row.sample],
          file(row.meth_bam, checkIfExists: true),
          meth_bam_index
        ]
    }
    // .dump()
    .set { ch_meth_bam }


Channel.fromPath(params.input)
    .splitCsv(header: true, sep: ',')
    .map { row ->
        // Determine the index file extension based on the rna_bam file extension
        def index_extension = row.rna_bam.endsWith('.bam') ? '.bai' : row.rna_bam.endsWith('.cram') ? '.crai' : ''
        def rna_bam_index = file(row.rna_bam + index_extension, checkIfExists: true)

        [ ['id'      : row.sample],
          file(row.rna_bam, checkIfExists: true),
          rna_bam_index
        ] 
    }
    // .dump()
    .set { ch_rna_bam }


Channel.fromPath(params.input)
    .splitCsv(header:true, sep:',')
    .map { 
        row -> [ ['id'      : row.sample],
                 file(row.vcf, checkIfExists: true)
                  ] 
        }
    // .dump()
    .set { ch_vcf }


workflow {


    filter_vcf (
        ch_vcf
    )

    subset_vcf (
        filter_vcf.out.filter_vcf,
        params.exons_bed,
        params.promoters_bed,
        params.dmrs_bed
    )

    sort_meth_bam (
        ch_meth_bam
    )
    // epiallele
    sort_meth_bam.out.meth_sorted_bam
        .join(subset_vcf.out.meth_roi_vcf)
        .set { ch_meth_bam_vcf }

    run_epiallele (
        ch_meth_bam_vcf
    )

    annotate_epiallele (
        run_epiallele.out.epiallele_txt
    )

    make_meth_report (
        annotate_epiallele.out.annotated_epiallele_txt
    )
    // snpdep
    ch_rna_bam
        .join(subset_vcf.out.rna_roi_vcf)
        .set { ch_rna_bam_vcf }

    run_snpdep (
        ch_rna_bam_vcf
    )
   
    make_rna_report (
        run_snpdep.out.snpdep_txt
    )

    make_rna_report.out.rna_vaf_report
        .join(make_meth_report.out.meth_vaf_report)
        .set { ch_rna_meth_reports }
    
    make_final_report (
        ch_rna_meth_reports
    )
}

