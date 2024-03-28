include { GATK4_COLLECTREADCOUNTS             } from '../../modules/nf-core/gatk4/collectreadcounts/main'
include { GATK4_CREATEREADCOUNTPANELOFNORMALS } from '../../modules/nf-core/gatk4/createreadcountpanelofnormals/main'
include { GATK4_PREPROCESSINTERVALS           } from '../../modules/nf-core/gatk4/preprocessintervals/main'
include { PICARD_CREATESEQUENCEDICTIONARY     } from '../../modules/nf-core/picard/createsequencedictionary/main'
include { SAMTOOLS_FAIDX                      } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_INDEX                      } from '../../modules/nf-core/samtools/index/main'

workflow GENS_PON {
    take:
        ch_user_dict     // channel: [mandatory] [ val(meta), path(dict) ]
        ch_user_fai      // channel: [mandatory] [ val(meta), path(fai) ]
        ch_fasta         // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_input         // channel: [mandatory] [ val(meta), path(bam/cram), path(bai/crai) ]

    main:
        ch_versions = Channel.empty()

        //
        //  Prepare references
        //
        SAMTOOLS_FAIDX ( ch_fasta, [[:],[]] )

        PICARD_CREATESEQUENCEDICTIONARY ( ch_fasta )

        ch_user_dict
            .mix(PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict)
            .collect()
            .set { ch_dict }

        ch_user_fai
            .mix(SAMTOOLS_FAIDX.out.fai)
            .collect()
            .set { ch_fai }

        GATK4_PREPROCESSINTERVALS ( ch_fasta, ch_fai, ch_dict, [[:],[]], [[:],[]] )

        //
        // Filter out files that lack indices, and generate them
        //
        ch_input
            .branch { meta, alignment, index ->
                alignment_with_index: index.size() > 0
                    return [meta, alignment, index]
                alignment_without_index: index.size() == 0
                    return [meta, alignment]
            }
            .set { ch_for_mix }

        SAMTOOLS_INDEX ( ch_for_mix.alignment_without_index )

        SAMTOOLS_INDEX.out.bai
            .mix(SAMTOOLS_INDEX.out.crai)
            .set { ch_index }

        //
        // Collect alignment files and their indices
        //
        ch_for_mix.alignment_without_index
            .join(ch_index)
            .mix(ch_for_mix.alignment_with_index)
            .combine(GATK4_PREPROCESSINTERVALS.out.interval_list.map{it -> it[1]})
            .set {ch_readcounts_in}

        //
        // Collect read counts, and generate models
        //
        GATK4_COLLECTREADCOUNTS ( ch_readcounts_in, ch_fasta, ch_fai, ch_dict )

        GATK4_COLLECTREADCOUNTS.out.tsv
            .mix(GATK4_COLLECTREADCOUNTS.out.hdf5)
            .collect { it[1] }
            .map { it ->
                    return [[id:"gens_pon"], it]
                }
            .set { ch_readcounts_out }

        GATK4_CREATEREADCOUNTPANELOFNORMALS ( ch_readcounts_out )

        ch_versions = ch_versions.mix(GATK4_COLLECTREADCOUNTS.out.versions.first())
        ch_versions = ch_versions.mix(GATK4_CREATEREADCOUNTPANELOFNORMALS.out.versions.first())
        ch_versions = ch_versions.mix(GATK4_PREPROCESSINTERVALS.out.versions)
        ch_versions = ch_versions.mix(PICARD_CREATESEQUENCEDICTIONARY.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
        genspon     = GATK4_CREATEREADCOUNTPANELOFNORMALS.out.pon
        readcounts  = ch_readcounts_out
        versions    = ch_versions
}
