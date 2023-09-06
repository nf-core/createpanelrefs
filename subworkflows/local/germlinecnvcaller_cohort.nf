include { GATK4_ANNOTATEINTERVALS             } from '../../modules/nf-core/gatk4/annotateintervals/main'
include { GATK4_COLLECTREADCOUNTS             } from '../../modules/nf-core/gatk4/collectreadcounts/main'
include { GATK4_DETERMINEGERMLINECONTIGPLOIDY } from '../../modules/nf-core/gatk4/determinegermlinecontigploidy/main'
include { GATK4_FILTERINTERVALS               } from '../../modules/nf-core/gatk4/filterintervals/main'
include { GATK4_GERMLINECNVCALLER             } from '../../modules/nf-core/gatk4/germlinecnvcaller/main'
include { GATK4_POSTPROCESSGERMLINECNVCALLS   } from '../../modules/nf-core/gatk4/postprocessgermlinecnvcalls/main'
include { GATK4_PREPROCESSINTERVALS           } from '../../modules/nf-core/gatk4/preprocessintervals/main'
include { PICARD_CREATESEQUENCEDICTIONARY     } from '../../modules/nf-core/picard/createsequencedictionary/main'
include { SAMTOOLS_FAIDX                      } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_INDEX                      } from '../../modules/nf-core/samtools/index/main'

workflow GERMLINECNVCALLER_COHORT {
    take:
        ch_bam           // channel: [mandatory] [ val(meta), [path(bam)] ]
        ch_fasta         // channel: [mandatory] [ val(meta), [path(fasta)] ]
        ch_ploidy_priors // channel: [mandatory] [ path(tsv) ]

    main:
        ch_versions = Channel.empty()

        ch_fai  = SAMTOOLS_FAIDX (ch_fasta, [[:],[]]).fai
        ch_dict = PICARD_CREATESEQUENCEDICTIONARY (ch_fasta).reference_dict
        ch_bai  = SAMTOOLS_INDEX (ch_bam)

        ch_bam_bai = ch_bam.join(SAMTOOLS_INDEX.out.bai)

        GATK4_PREPROCESSINTERVALS (ch_fasta,
                                   ch_fai,
                                   ch_dict,
                                   [[:],[]], [[:],[]])

        ch_bam_bai
            .combine(GATK4_PREPROCESSINTERVALS.out.interval_list.map{it -> it[1]})
            .set {ch_readcounts_in}

        GATK4_COLLECTREADCOUNTS (ch_readcounts_in,
                                 ch_fasta,
                                 ch_fai,
                                 ch_dict)
                                .tsv
                                .collect { it[1] }
                                .map {tsv -> [[id:'cohort'],tsv]}
                                .set { ch_readcounts_out }

        GATK4_ANNOTATEINTERVALS (GATK4_PREPROCESSINTERVALS.out.interval_list,
                                 ch_fasta,
                                 ch_fai,
                                 ch_dict,
                                 [[:],[]], [[:],[]], [[:],[]], [[:],[]])

        GATK4_FILTERINTERVALS (GATK4_PREPROCESSINTERVALS.out.interval_list,
                               ch_readcounts_out,
                               GATK4_ANNOTATEINTERVALS.out.annotated_intervals)

        ch_readcounts_out
            .combine(GATK4_FILTERINTERVALS.out.interval_list)
            .map{ meta, counts, meta2, il -> [meta,  counts, il, []] }
            .set {ch_contigploidy_in}

        GATK4_DETERMINEGERMLINECONTIGPLOIDY (ch_contigploidy_in,
                                             [[:],[]],
                                             ch_ploidy_priors)


        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        ch_versions = ch_versions.mix(PICARD_CREATESEQUENCEDICTIONARY.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
        ch_versions = ch_versions.mix(GATK4_PREPROCESSINTERVALS.out.versions)
        ch_versions = ch_versions.mix(GATK4_COLLECTREADCOUNTS.out.versions.first())
        ch_versions = ch_versions.mix(GATK4_ANNOTATEINTERVALS.out.versions)
        ch_versions = ch_versions.mix(GATK4_FILTERINTERVALS.out.versions)
        ch_versions = ch_versions.mix(GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.versions)

    emit:
        versions = ch_versions
}
