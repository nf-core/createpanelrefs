include { GATK4_ANNOTATEINTERVALS                                      } from '../../../modules/nf-core/gatk4/annotateintervals'
include { GATK4_BEDTOINTERVALLIST as GATK4_BEDTOINTERVALLIST_TARGETS   } from '../../../modules/nf-core/gatk4/bedtointervallist'
include { GATK4_BEDTOINTERVALLIST as GATK4_BEDTOINTERVALLIST_EXCLUDE   } from '../../../modules/nf-core/gatk4/bedtointervallist'
include { GATK4_COLLECTREADCOUNTS                                      } from '../../../modules/nf-core/gatk4/collectreadcounts'
include { GATK4_DETERMINEGERMLINECONTIGPLOIDY                          } from '../../../modules/nf-core/gatk4/determinegermlinecontigploidy'
include { GATK4_FILTERINTERVALS                                        } from '../../../modules/nf-core/gatk4/filterintervals'
include { GATK4_GERMLINECNVCALLER                                      } from '../../../modules/nf-core/gatk4/germlinecnvcaller'
include { GATK4_INDEXFEATUREFILE as GATK4_INDEXFEATUREFILE_MAPPABILITY } from '../../../modules/nf-core/gatk4/indexfeaturefile'
include { GATK4_INDEXFEATUREFILE as GATK4_INDEXFEATUREFILE_SEGDUP      } from '../../../modules/nf-core/gatk4/indexfeaturefile'
include { GATK4_INTERVALLISTTOOLS                                      } from '../../../modules/nf-core/gatk4/intervallisttools'
include { GATK4_PREPROCESSINTERVALS                                    } from '../../../modules/nf-core/gatk4/preprocessintervals'
include { SAMTOOLS_INDEX                                               } from '../../../modules/nf-core/samtools/index'

workflow GERMLINECNVCALLER_COHORT {
    take:
    ch_input                      // channel: [mandatory] [ val(meta), path(bam/cram), path(bai/crai) ]
    val_pon_name                  //  string: [optional] name for panel of normals
    ch_dict                       // channel: [optional] [ val(meta), path(dict) ]
    ch_fai                        // channel: [optional] [ val(meta), path(fai) ]
    ch_fasta                      // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_exclude_bed                // channel: [optional] [ val(meta), path(bed) ]
    ch_user_exclude_interval_list // channel: [optional] [ val(meta), path(intervals) ]
    ch_mappable_regions           // channel: [optional] [ val(meta), path(bed) ]
    ch_ploidy_priors              // channel: [mandatory] [ path(tsv) ]
    ch_segmental_duplications     // channel: [optional] [ val(meta), path(bed) ]
    ch_target_bed                 // channel: [optional] [ val(meta), path(bed) ]
    ch_user_target_interval_list  // channel: [optional] [ val(meta), path(intervals) ]

    main:
    versions = Channel.empty()

    //  Prepare references
    GATK4_INDEXFEATUREFILE_MAPPABILITY(ch_mappable_regions)
    GATK4_INDEXFEATUREFILE_SEGDUP(ch_segmental_duplications)

    //Runs for wes analysis, when target_bed file is provided instead of target_interval_list
    GATK4_BEDTOINTERVALLIST_TARGETS(ch_target_bed, ch_dict)

    //Runs for wes analysis, when exclude_bed file is provided instead of target_interval_list
    GATK4_BEDTOINTERVALLIST_EXCLUDE(ch_exclude_bed, ch_dict)

    ch_user_target_interval_list
        .combine(GATK4_BEDTOINTERVALLIST_TARGETS.out.interval_list.ifEmpty(null))
        .branch { it ->
            intervallistfrompath: it[2].equals(null)
            return [it[0], it[1]]
            intervallistfrombed: !it[2].equals(null)
            return [it[2], it[3]]
        }
        .set { ch_targets_for_mix }

    ch_targets_for_mix.intervallistfrompath
        .mix(ch_targets_for_mix.intervallistfrombed)
        .collect()
        .set { ch_target_interval_list }

    ch_user_exclude_interval_list
        .combine(GATK4_BEDTOINTERVALLIST_EXCLUDE.out.interval_list.ifEmpty(null))
        .branch { it ->
            intervallistfrompath: it[2].equals(null)
            return [it[0], it[1]]
            intervallistfrombed: !it[2].equals(null)
            return [it[2], it[3]]
        }
        .set { ch_exclude_for_mix }

    ch_exclude_for_mix.intervallistfrompath
        .mix(ch_exclude_for_mix.intervallistfrombed)
        .collect()
        .set { ch_exclude_interval_list }

    GATK4_PREPROCESSINTERVALS(
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_target_interval_list,
        ch_exclude_interval_list,
    )

    GATK4_ANNOTATEINTERVALS(
        GATK4_PREPROCESSINTERVALS.out.interval_list,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_mappable_regions,
        GATK4_INDEXFEATUREFILE_MAPPABILITY.out.index.ifEmpty([[:], []]),
        ch_segmental_duplications,
        GATK4_INDEXFEATUREFILE_SEGDUP.out.index.ifEmpty([[:], []]),
    )

    // Filter out files that lack indices, and generate them
    ch_input
        .branch { meta, alignment, index ->
            alignment_with_index: index.size() > 0
            return [meta, alignment, index]
            alignment_without_index: index.size() == 0
            return [meta, alignment]
        }
        .set { ch_for_mix }

    SAMTOOLS_INDEX(ch_for_mix.alignment_without_index)

    SAMTOOLS_INDEX.out.bai
        .mix(SAMTOOLS_INDEX.out.crai)
        .set { ch_index }

    // Collect alignment files and their indices
    ch_for_mix.alignment_without_index
        .join(ch_index)
        .mix(ch_for_mix.alignment_with_index)
        .combine(GATK4_PREPROCESSINTERVALS.out.interval_list.map { it -> it[1] })
        .set { ch_readcounts_in }

    // Collect read counts, and generate models
    GATK4_COLLECTREADCOUNTS(
        ch_readcounts_in,
        ch_fasta,
        ch_fai,
        ch_dict,
    )

    GATK4_COLLECTREADCOUNTS.out.tsv
        .mix(GATK4_COLLECTREADCOUNTS.out.hdf5)
        .collect { it[1] }
        .map { tsv -> [[id: val_pon_name], tsv] }
        .set { ch_readcounts_out }


    GATK4_FILTERINTERVALS(
        GATK4_PREPROCESSINTERVALS.out.interval_list,
        ch_readcounts_out,
        GATK4_ANNOTATEINTERVALS.out.annotated_intervals,
    )

    GATK4_INTERVALLISTTOOLS(GATK4_FILTERINTERVALS.out.interval_list).interval_list.map { meta, it -> it }.flatten().set { ch_intervallist_out }

    ch_readcounts_out
        .combine(GATK4_FILTERINTERVALS.out.interval_list)
        .map { meta, counts, meta2, il -> [meta, counts, il, []] }
        .set { ch_contigploidy_in }

    GATK4_DETERMINEGERMLINECONTIGPLOIDY(
        ch_contigploidy_in,
        [[:], []],
        ch_ploidy_priors,
    )

    ch_readcounts_out
        .combine(ch_intervallist_out)
        .combine(GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.calls)
        .map { meta, counts, il, meta2, calls -> [meta + [id: il.baseName], counts, il, calls, []] }
        .set { ch_cnvcaller_in }

    GATK4_GERMLINECNVCALLER(ch_cnvcaller_in)

    versions = versions.mix(SAMTOOLS_INDEX.out.versions)
    versions = versions.mix(GATK4_PREPROCESSINTERVALS.out.versions)
    versions = versions.mix(GATK4_BEDTOINTERVALLIST_TARGETS.out.versions)
    versions = versions.mix(GATK4_BEDTOINTERVALLIST_EXCLUDE.out.versions)
    versions = versions.mix(GATK4_COLLECTREADCOUNTS.out.versions)
    versions = versions.mix(GATK4_ANNOTATEINTERVALS.out.versions)
    versions = versions.mix(GATK4_FILTERINTERVALS.out.versions)
    versions = versions.mix(GATK4_INDEXFEATUREFILE_MAPPABILITY.out.versions)
    versions = versions.mix(GATK4_INDEXFEATUREFILE_SEGDUP.out.versions)
    versions = versions.mix(GATK4_INTERVALLISTTOOLS.out.versions)
    versions = versions.mix(GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.versions)
    versions = versions.mix(GATK4_GERMLINECNVCALLER.out.versions)

    emit:
    cnvmodel    = GATK4_GERMLINECNVCALLER.out.cohortmodel
    ploidymodel = GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.model
    readcounts  = ch_readcounts_out
    versions
}
