include { GATK4_ANNOTATEINTERVALS             } from '../../modules/nf-core/gatk4/annotateintervals/main'
include { GATK4_COLLECTREADCOUNTS             } from '../../modules/nf-core/gatk4/collectreadcounts/main'
include { GATK4_DETERMINEGERMLINECONTIGPLOIDY } from '../../modules/nf-core/gatk4/determinegermlinecontigploidy/main'
include { GATK4_FILTERINTERVALS               } from '../../modules/nf-core/gatk4/filterintervals/main'
include { GATK4_GERMLINECNVCALLER             } from '../../modules/nf-core/gatk4/germlinecnvcaller/main'
include { GATK4_INTERVALLISTTOOLS             } from '../../modules/nf-core/gatk4/intervallisttools/main'
include { GATK4_PREPROCESSINTERVALS           } from '../../modules/nf-core/gatk4/preprocessintervals/main'
include { PICARD_CREATESEQUENCEDICTIONARY     } from '../../modules/nf-core/picard/createsequencedictionary/main'
include { SAMTOOLS_FAIDX                      } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_INDEX                      } from '../../modules/nf-core/samtools/index/main'

workflow GERMLINECNVCALLER_COHORT {
    take:
        ch_user_dict     // channel: [mandatory] [ val(meta), path(dict) ]
        ch_user_fai      // channel: [mandatory] [ val(meta), path(fai) ]
        ch_fasta         // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_input         // channel: [mandatory] [ val(meta), path(bam/cram), path(bai/crai) ]
        ch_ploidy_priors // channel: [mandatory] [ path(tsv) ]

    main:
        ch_versions = Channel.empty()

        //
        //  Prepare references
        //
        SAMTOOLS_FAIDX (ch_fasta, [[:],[]])

        PICARD_CREATESEQUENCEDICTIONARY (ch_fasta)

        ch_user_dict
            .mix(PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict)
            .collect()
            .set { ch_dict }

        ch_user_fai
            .mix(SAMTOOLS_FAIDX.out.fai)
            .collect()
            .set { ch_fai }

        GATK4_PREPROCESSINTERVALS (ch_fasta,
                                   ch_fai,
                                   ch_dict,
                                   [[:],[]], [[:],[]])

        GATK4_ANNOTATEINTERVALS (GATK4_PREPROCESSINTERVALS.out.interval_list,
                                 ch_fasta,
                                 ch_fai,
                                 ch_dict,
                                 [[:],[]], [[:],[]], [[:],[]], [[:],[]])

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

        SAMTOOLS_INDEX (ch_for_mix.alignment_without_index)

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
        GATK4_COLLECTREADCOUNTS (ch_readcounts_in,
                                 ch_fasta,
                                 ch_fai,
                                 ch_dict)

        GATK4_COLLECTREADCOUNTS.out.tsv
                                .mix(GATK4_COLLECTREADCOUNTS.out.hdf5)
                                .collect { it[1] }
                                .map {tsv -> [[id:'cohort'],tsv]}
                                .set { ch_readcounts_out }


        GATK4_FILTERINTERVALS (GATK4_PREPROCESSINTERVALS.out.interval_list,
                               ch_readcounts_out,
                               GATK4_ANNOTATEINTERVALS.out.annotated_intervals)

        GATK4_INTERVALLISTTOOLS(GATK4_FILTERINTERVALS.out.interval_list)
                                .interval_list
                                .map {meta, it -> it}
                                .flatten()
                                .set { ch_intervallist_out }

        ch_readcounts_out
                .combine(GATK4_FILTERINTERVALS.out.interval_list)
                .map{ meta, counts, meta2, il -> [meta, counts, il, []] }
                .set {ch_contigploidy_in}

        GATK4_DETERMINEGERMLINECONTIGPLOIDY (ch_contigploidy_in,
                                             [[:],[]],
                                             ch_ploidy_priors)

        ch_readcounts_out
                .combine(ch_intervallist_out)
                .combine(GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.calls)
                .map{ meta, counts, il, meta2, calls -> [meta + [id:il.baseName],  counts, il, calls, []] }
                .set {ch_cnvcaller_in}

        GATK4_GERMLINECNVCALLER (ch_cnvcaller_in)

        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        ch_versions = ch_versions.mix(PICARD_CREATESEQUENCEDICTIONARY.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())
        ch_versions = ch_versions.mix(GATK4_PREPROCESSINTERVALS.out.versions)
        ch_versions = ch_versions.mix(GATK4_COLLECTREADCOUNTS.out.versions.first())
        ch_versions = ch_versions.mix(GATK4_ANNOTATEINTERVALS.out.versions)
        ch_versions = ch_versions.mix(GATK4_FILTERINTERVALS.out.versions)
        ch_versions = ch_versions.mix(GATK4_INTERVALLISTTOOLS.out.versions)
        ch_versions = ch_versions.mix(GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.versions)
        ch_versions = ch_versions.mix(GATK4_GERMLINECNVCALLER.out.versions.first())

    emit:
        cnvmodel    = GATK4_GERMLINECNVCALLER.out.cohortmodel
        ploidymodel = GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.model
        readcounts  = ch_readcounts_out
        versions    = ch_versions
}
