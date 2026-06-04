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

workflow GERMLINECNVCALLER_COHORT {
    take:
    ch_input // channel: [mandatory] [ val(meta), path(bam/cram), path(bai/crai) ]
    val_pon_name //  string: [optional] name for panel of normals
    val_analysis_type // string: [mandatory] type of analysis ('wes' or 'wgs')
    ch_dict // channel: [optional] [ val(meta), path(dict) ]
    ch_fai // channel: [optional] [ val(meta), path(fai) ]
    ch_fasta // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_exclude_bed // channel: [optional] [ val(meta), path(bed) ]
    ch_user_exclude_interval_list // channel: [optional] [ val(meta), path(intervals) ]
    ch_mappable_regions // channel: [optional] [ val(meta), path(bed) ]
    ch_ploidy_priors // channel: [mandatory] [ path(tsv) ]
    ch_segmental_duplications // channel: [optional] [ val(meta), path(bed) ]
    ch_target_bed // channel: [optional] [ val(meta), path(bed) ]
    ch_user_target_interval_list // channel: [optional] [ val(meta), path(intervals) ]

    main:
    //  Index feature files — only when a real file is provided
    GATK4_INDEXFEATUREFILE_MAPPABILITY(ch_mappable_regions.filter { _meta, regions -> !(regions instanceof List) })
    GATK4_INDEXFEATUREFILE_SEGDUP(ch_segmental_duplications.filter { _meta, segdup -> !(segdup instanceof List) })

    // Bed to interval list conversion — only for WES when bed is provided and no interval list given
    ch_target_bed_for_conversion = channel.of(val_analysis_type)
        .combine(ch_user_target_interval_list.filter { _meta, interval -> interval instanceof List })
        .combine(ch_target_bed.filter { _meta, bed -> !(bed instanceof List) })
        .filter { analysis_type, _meta_interval, _interval, _meta_bed, _bed -> analysis_type == "wes" }
        .map { _analysis_type, _meta_interval, _interval, meta_bed, bed -> [meta_bed, bed] }

    GATK4_BEDTOINTERVALLIST_TARGETS(ch_target_bed_for_conversion, ch_dict)

    ch_exclude_bed_for_conversion = channel.of(val_analysis_type)
        .combine(ch_user_exclude_interval_list.filter { _meta, interval -> interval instanceof List })
        .combine(ch_exclude_bed.filter { _meta, bed -> !(bed instanceof List) })
        .filter { analysis_type, _meta_interval, _interval, _meta_bed, _bed -> analysis_type == "wes" }
        .map { _analysis_type, _meta_interval, _interval, meta_bed, bed -> [meta_bed, bed] }

    GATK4_BEDTOINTERVALLIST_EXCLUDE(ch_exclude_bed_for_conversion, ch_dict)

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

    ch_input
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
        .collect { _meta, file -> [file] }
        .map { tsv -> [[id: val_pon_name], tsv] }
        .set { ch_readcounts_out }


    GATK4_FILTERINTERVALS(
        GATK4_PREPROCESSINTERVALS.out.interval_list,
        ch_readcounts_out,
        GATK4_ANNOTATEINTERVALS.out.annotated_intervals,
    )

    GATK4_INTERVALLISTTOOLS(GATK4_FILTERINTERVALS.out.interval_list).interval_list.map { _meta, it -> it }.flatten().set { ch_intervallist_out }

    ch_readcounts_out
        .combine(GATK4_FILTERINTERVALS.out.interval_list)
        .map { meta, counts, _meta2, il -> [meta, counts, il, []] }
        .set { ch_contigploidy_in }

    GATK4_DETERMINEGERMLINECONTIGPLOIDY(
        ch_contigploidy_in,
        [[:], []],
        ch_ploidy_priors,
    )

    ch_readcounts_out
        .combine(ch_intervallist_out)
        .combine(GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.calls)
        .map { meta, counts, il, _meta2, calls -> [meta + [id: il.baseName], counts, il, calls, []] }
        .set { ch_cnvcaller_in }

    GATK4_GERMLINECNVCALLER(ch_cnvcaller_in)

    emit:
    cnvmodel    = GATK4_GERMLINECNVCALLER.out.cohortmodel
    ploidymodel = GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.model
    readcounts  = ch_readcounts_out
}
