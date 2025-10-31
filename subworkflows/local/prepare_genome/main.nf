include { GATK4_CREATESEQUENCEDICTIONARY                              } from '../../../modules/nf-core/gatk4/createsequencedictionary'
include { GATK4_PREPROCESSINTERVALS as GATK4_PREPROCESSINTERVALS_GENS } from '../../../modules/nf-core/gatk4/preprocessintervals'
include { GAWK as BUILD_INTERVALS                                     } from '../../../modules/nf-core/gawk'
include { SAMTOOLS_FAIDX                                              } from '../../../modules/nf-core/samtools/faidx'

//  Prepare references
workflow PREPARE_GENOME {
    take:
    fasta                   // channel: [mandatory] [ val(meta), path(fasta) ]
    user_dict               // channel: [optional]  [ val(meta), path(dict) ]
    user_fai                // channel: [optional]  [ val(meta), path(fai) ]
    user_gens_interval_list // channel: [optional]  [ val(meta), path(gens_interval_list) ]
    user_mutect2_target_bed // channel: [optional]  [ val(meta), path(mutect2_target_bed) ]
    tools                   //   array: [mandatory] [ tools ]

    main:
    dict = Channel.empty()
    fai = Channel.empty()
    gens_interval_list = Channel.empty()
    mutect2_target_bed = Channel.empty()
    versions = Channel.empty()

    // If a user_dict is provided, no fasta will be used to generate a dict
    // Otherwise, GATK4_CREATESEQUENCEDICTIONARY will be run to generate a dict
    fasta_for_dict = fasta
        .mix(user_dict)
        .groupTuple()
        .filter { _meta, files -> !files[1] }

    GATK4_CREATESEQUENCEDICTIONARY(fasta_for_dict)

    dict = user_dict.mix(GATK4_CREATESEQUENCEDICTIONARY.out.dict).collect()

    // If a user_fai is provided, no fasta will be used to generate a fai
    // Otherwise, SAMTOOLS_FAIDX will be run to generate a fai
    fasta_for_fai = fasta
        .mix(user_fai)
        .groupTuple()
        .filter { _meta, files -> !files[1] }

    SAMTOOLS_FAIDX(fasta_for_fai, [[:], []], false)

    fai = user_fai.mix(SAMTOOLS_FAIDX.out.fai).collect()

    // If a user_gens_interval_list is provided or if gens is not a specified tools, no fasta will be used to generate an interval list
    // Otherwise, GATK4_PREPROCESSINTERVALS_GENS will be run to generate an interval list
    fasta_for_interval_list = fasta
        .mix(user_gens_interval_list)
        .groupTuple()
        .filter { _meta, files -> (tools.split(',').contains('gens') && !files[1]) }

    GATK4_PREPROCESSINTERVALS_GENS(fasta_for_interval_list, fai.collect(), dict.collect(), [[:], []], [[:], []])

    gens_interval_list = user_gens_interval_list.mix(GATK4_PREPROCESSINTERVALS_GENS.out.interval_list).collect()

    // If a user_mutect2_target_bed is provided or if mutect2 is not a specified tools, no fai will be used to generate a target bed
    // Otherwise, BUILD_INTERVALS will be run to generate a target bed
    fai_for_intervals = fai
        .mix(user_mutect2_target_bed)
        .groupTuple()
        .filter { _meta, files -> (tools.split(',').contains('mutect2') && !files[1]) }

    BUILD_INTERVALS(fai_for_intervals, [], false)

    mutect2_target_bed = user_mutect2_target_bed.mix(BUILD_INTERVALS.out.output).collect()

    versions = versions.mix(BUILD_INTERVALS.out.versions)
    versions = versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    versions = versions.mix(GATK4_PREPROCESSINTERVALS_GENS.out.versions)
    versions = versions.mix(SAMTOOLS_FAIDX.out.versions)

    emit:
    dict               // channel: [ val(meta), path(dict) ]
    fai                // channel: [ val(meta), path(fai) ]
    gens_interval_list // channel: [ val(meta), path(gens_interval_list) ]
    mutect2_target_bed // channel: [ val(meta), path(mutect2_target_bed) ]
    versions           // channel: [ path(versions.yml)]
}
