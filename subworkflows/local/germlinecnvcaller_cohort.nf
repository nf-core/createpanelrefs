include { GATK4_COLLECTREADCOUNTS         } from '../../modules/nf-core/gatk4/collectreadcounts/main'
include { GATK4_PREPROCESSINTERVALS       } from '../../modules/nf-core/gatk4/preprocessintervals/main'
include { PICARD_CREATESEQUENCEDICTIONARY } from '../../modules/nf-core/picard/createsequencedictionary/main'
include { SAMTOOLS_FAIDX                  } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_INDEX                  } from '../../modules/nf-core/samtools/index/main'

workflow GERMLINECNVCALLER_COHORT {
    take:
        ch_bam    // channel: [mandatory] [ val(meta), [path(bam)] ]
        ch_fasta  // channel: [mandatory] [ val(meta), [path(fasta)] ]

    main:
        ch_versions = Channel.empty()

        ch_fai  = SAMTOOLS_FAIDX (ch_fasta, [[:],[]]).fai
        ch_dict = PICARD_CREATESEQUENCEDICTIONARY (ch_fasta).reference_dict
        ch_bai  = SAMTOOLS_INDEX (ch_bam)

        ch_bam_bai = ch_bam.join(SAMTOOLS_INDEX.out.bai)

        GATK4_PREPROCESSINTERVALS (ch_fasta, ch_fai, ch_dict, [[:],[]], [[:],[]])

        ch_bam_bai
            .combine(GATK4_PREPROCESSINTERVALS.out.interval_list.map{it -> it[1]})
            .set {ch_readcounts_in}

        GATK4_COLLECTREADCOUNTS(ch_readcounts_in, ch_fasta, ch_fai, ch_dict)

        ch_versions = ch_versions.mix(GATK4_PREPROCESSINTERVALS.out.versions)

    emit:
        versions = ch_versions
}
