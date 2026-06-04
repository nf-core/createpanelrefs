include { GATK4_CREATESEQUENCEDICTIONARY                              } from '../../../modules/nf-core/gatk4/createsequencedictionary'
include { GATK4_PREPROCESSINTERVALS as GATK4_PREPROCESSINTERVALS_GENS } from '../../../modules/nf-core/gatk4/preprocessintervals'
include { GATK4_SPLITINTERVALS                                        } from '../../../modules/nf-core/gatk4/splitintervals'
include { GAWK as BUILD_INTERVALS                                     } from '../../../modules/nf-core/gawk'
include { SAMTOOLS_FAIDX                                              } from '../../../modules/nf-core/samtools/faidx'

workflow PREPARE_GENOME {
    take:
    ch_fasta // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_user_dict // channel: [optional]  [ val(meta), path(dict) ]
    ch_user_fai // channel: [optional]  [ val(meta), path(fai) ]
    ch_user_gens_interval_list // channel: [optional]  [ val(meta), path(gens_interval_list) ]
    ch_user_mutect2_target_bed // channel: [optional]  [ val(meta), path(mutect2_target_bed) ]
    mutect2_intervals_num //   value: [optional] number of intervals for mutect2 scatter
    tools //   list: [mandatory] tools to run

    main:
    ch_dict = channel.empty()
    ch_fai = channel.empty()
    ch_gens_interval_list = channel.empty()
    intervals_num = channel.empty()
    ch_mutect2_target_bed = channel.empty()

    // If a user_dict is provided, no fasta will be used to generate a dict
    // Otherwise, GATK4_CREATESEQUENCEDICTIONARY will be run to generate a dict
    ch_fasta_for_dict = ch_fasta
        .mix(ch_user_dict)
        .groupTuple()
        .filter { _meta, files -> !files[1] }

    GATK4_CREATESEQUENCEDICTIONARY(ch_fasta_for_dict)

    ch_dict = ch_user_dict.mix(GATK4_CREATESEQUENCEDICTIONARY.out.dict).collect()

    // If a user_fai is provided, no fasta will be used to generate a fai
    // Otherwise, SAMTOOLS_FAIDX will be run to generate a fai
    ch_fasta_for_fai = ch_fasta
        .mix(ch_user_fai)
        .groupTuple()
        .filter { _meta, files -> !files[1] }
        .map { meta, fasta_ -> [meta, fasta_, []] }

    SAMTOOLS_FAIDX(ch_fasta_for_fai, false)

    ch_fai = ch_user_fai.mix(SAMTOOLS_FAIDX.out.fai).collect()

    // If a user_gens_interval_list is provided or if gens is not a specified tools, no fasta will be used to generate an interval list
    // Otherwise, GATK4_PREPROCESSINTERVALS_GENS will be run to generate an interval list
    ch_fasta_for_interval_list = ch_fasta
        .mix(ch_user_gens_interval_list)
        .groupTuple()
        .filter { _meta, files -> ('gens' in tools && !files[1]) }

    GATK4_PREPROCESSINTERVALS_GENS(ch_fasta_for_interval_list, ch_fai.collect(), ch_dict.collect(), [[:], []], [[:], []])

    ch_gens_interval_list = ch_user_gens_interval_list.mix(GATK4_PREPROCESSINTERVALS_GENS.out.interval_list).collect()

    // If a user_mutect2_target_bed is provided or if mutect2 is not a specified tools, no fai will be used to generate a target bed
    // Otherwise, BUILD_INTERVALS will be run to generate a target bed
    ch_fai_for_intervals = ch_fai
        .mix(ch_user_mutect2_target_bed)
        .groupTuple()
        .filter { _meta, files -> ('mutect2' in tools && !files[1]) }

    BUILD_INTERVALS(ch_fai_for_intervals, [], false)

    ch_mutect2_target_bed = ch_user_mutect2_target_bed.mix(BUILD_INTERVALS.out.output).collect()

    // If mutect2 is in tools and mutect2_intervals_num > 1, split intervals for scatter/gather strategy
    // intervals_num: [ path(intervals), val(num_intervals) ]
    // num_intervals > 1 triggers per-interval scatter and merge of outputs
    if ('mutect2' in tools && mutect2_intervals_num > 1) {
        GATK4_SPLITINTERVALS(
            ch_mutect2_target_bed,
            ch_fasta,
            ch_fai,
            ch_dict,
        )

        intervals_num = GATK4_SPLITINTERVALS.out.split_intervals.flatMap { _meta, intervals -> intervals.collect { interval -> [interval, intervals.size()] } }
    }
    else {
        intervals_num = channel.of([[], 1])
    }

    emit:
    dict = ch_dict // channel: [ val(meta), path(dict) ]
    fai = ch_fai // channel: [ val(meta), path(fai) ]
    gens_interval_list = ch_gens_interval_list // channel: [ val(meta), path(gens_interval_list) ]
    intervals_num // channel: [ path(intervals), val(num_intervals) ]
    mutect2_target_bed = ch_mutect2_target_bed // channel: [ val(meta), path(mutect2_target_bed) ]
}
