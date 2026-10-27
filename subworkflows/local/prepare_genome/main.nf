include { GATK4_CREATESEQUENCEDICTIONARY                              } from '../../../modules/nf-core/gatk4/createsequencedictionary'
include { GATK4_PREPROCESSINTERVALS as GATK4_PREPROCESSINTERVALS_GENS } from '../../../modules/nf-core/gatk4/preprocessintervals'
include { GATK4_SPLITINTERVALS                                        } from '../../../modules/nf-core/gatk4/splitintervals'
include { GAWK as BUILD_INTERVALS                                     } from '../../../modules/nf-core/gawk'
include { SAMTOOLS_FAIDX                                              } from '../../../modules/nf-core/samtools/faidx'

workflow PREPARE_GENOME {
    take:
    tools //   list: [mandatory] tools to run
    params_cnvkit_targets // [optional] path(cnvkit_targets)
    params_dict // [optional]  path(dict)
    params_fai // [optional]  path(fai)
    params_fasta // [mandatory] path(fasta)
    params_gcnv_exclude_bed // [optional] path(gcnv_exclude_bed)
    params_gcnv_exclude_interval_list // [optional] path(gcnv_exclude_interval_list)
    params_gcnv_mappable_regions // [optional] path(gcnv_mappable_regions)
    params_gcnv_ploidy_priors // [optional] path(gcnv_ploidy_priors)
    params_gcnv_segmental_duplications // [optional] path(gcnv_segmental_duplications)
    params_gcnv_target_bed // [optional] path(gcnv_target_bed)
    params_gcnv_target_interval_list // [optional] path(gcnv_target_interval_list)
    params_gens_interval_list // [optional]  path(gens_interval_list)
    params_mutect2_intervals_num //   value: [optional] number of intervals for mutect2 scatter
    params_mutect2_target_bed // [optional]  path(mutect2_target_bed)

    main:
    def ch_fasta = channel.fromPath(params_fasta).map { fasta -> [[id: 'genome'], fasta] }.collect()

    // Initialize cnvkit specific parameters
    def is_cnvkit_targets = params_cnvkit_targets
    def ch_cnvkit_targets = is_cnvkit_targets
        ? channel.fromPath(params_cnvkit_targets).map { targets -> [[id: 'genome'], targets] }.collect()
        : channel.value([[id: 'genome'], []])

    // Initialize germlinecnvcaller specific parameters
    def is_gcnv_exclude_bed = params_gcnv_exclude_bed
    def ch_gcnv_exclude_bed = is_gcnv_exclude_bed
        ? channel.fromPath(params_gcnv_exclude_bed).map { exclude -> [[id: 'genome'], exclude] }.collect()
        : channel.value([[id: 'genome'], []])

    def is_gcnv_exclude_interval_list = params_gcnv_exclude_interval_list
    def ch_gcnv_exclude_interval_list = is_gcnv_exclude_interval_list
        ? channel.fromPath(params_gcnv_exclude_interval_list).map { exclude -> [[id: 'genome'], exclude] }.collect()
        : channel.value([[id: 'genome'], []])

    def is_gcnv_mappable_regions = params_gcnv_mappable_regions
    def ch_gcnv_mappable_regions = is_gcnv_mappable_regions
        ? channel.fromPath(params_gcnv_mappable_regions).collect()
        : channel.value([[id: 'genome'], []])

    def is_gcnv_ploidy_priors = params_gcnv_ploidy_priors
    def ch_gcnv_ploidy_priors = is_gcnv_ploidy_priors
        ? channel.fromPath(params_gcnv_ploidy_priors).collect()
        : channel.empty()

    def is_gcnv_target_bed = params_gcnv_target_bed
    def ch_gcnv_target_bed = is_gcnv_target_bed
        ? channel.fromPath(params_gcnv_target_bed).map { targets -> [[id: 'genome'], targets] }.collect()
        : channel.value([[id: 'genome'], []])

    def is_gcnv_target_interval_list = params_gcnv_target_interval_list
    def ch_gcnv_target_interval_list = is_gcnv_target_interval_list
        ? channel.fromPath(params_gcnv_target_interval_list).map { targets -> [[id: 'genome'], targets] }.collect()
        : channel.value([[id: 'genome'], []])

    def is_gcnv_segmental_duplications = params_gcnv_segmental_duplications
    def ch_gcnv_segmental_duplications = is_gcnv_segmental_duplications
        ? channel.fromPath(params_gcnv_segmental_duplications).collect()
        : channel.value([[id: 'genome'], []])



    def run_createsequencedictionary = !params_dict

    GATK4_CREATESEQUENCEDICTIONARY(ch_fasta.filter { run_createsequencedictionary })

    def ch_dict = run_createsequencedictionary
        ? GATK4_CREATESEQUENCEDICTIONARY.out.dict.collect()
        : channel.fromPath(params_dict).map { dict -> [[id: 'genome'], dict] }.collect()

    def run_faidx = !params_fai

    SAMTOOLS_FAIDX(ch_fasta.filter { run_faidx }, [])

    def ch_fai = run_faidx
        ? SAMTOOLS_FAIDX.out.fai.collect()
        : channel.fromPath(params_fai).map { fai -> [[id: 'genome'], fai] }.collect()

    def run_preprocessintervals_gens = !params_gens_interval_list

    GATK4_PREPROCESSINTERVALS_GENS(
        ch_fasta.filter { run_preprocessintervals_gens && ('gens' in tools) },
        ch_fai.collect(),
        ch_dict.collect(),
        [[:], []],
        [[:], []],
    )

    def ch_gens_interval_list = run_preprocessintervals_gens
        ? GATK4_PREPROCESSINTERVALS_GENS.out.interval_list.collect()
        : channel.fromPath(params_gens_interval_list).map { gens_interval_list -> [[id: 'genome'], gens_interval_list] }.collect()

    def run_build_intervals = !params_mutect2_target_bed

    BUILD_INTERVALS(
        ch_fai.filter { run_build_intervals && ('mutect2' in tools) },
        [],
        false,
    )

    def ch_mutect2_target_bed = run_build_intervals
        ? BUILD_INTERVALS.out.output.collect()
        : channel.fromPath(params_mutect2_target_bed).map { targets -> [[id: 'genome'], targets] }.collect()

    // If mutect2 is in tools and params_mutect2_intervals_num > 1, split intervals for scatter/gather strategy
    // intervals_num: [ path(intervals), val(num_intervals) ]
    // num_intervals > 1 triggers per-interval scatter and merge of outputs
    if ('mutect2' in tools && params_mutect2_intervals_num > 1) {
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
    cnvkit_targets              = ch_cnvkit_targets
    dict                        = ch_dict // channel: [ val(meta), path(dict) ]
    fai                         = ch_fai // channel: [ val(meta), path(fai) ]
    fasta                       = ch_fasta
    gcnv_exclude_bed            = ch_gcnv_exclude_bed
    gcnv_exclude_interval_list  = ch_gcnv_exclude_interval_list
    gcnv_mappable_regions       = ch_gcnv_mappable_regions
    gcnv_ploidy_priors          = ch_gcnv_ploidy_priors
    gcnv_segmental_duplications = ch_gcnv_segmental_duplications
    gcnv_target_bed             = ch_gcnv_target_bed
    gcnv_target_interval_list   = ch_gcnv_target_interval_list
    gens_interval_list          = ch_gens_interval_list // channel: [ val(meta), path(gens_interval_list) ]
    intervals_num // channel: [ path(intervals), val(num_intervals) ]
    mutect2_target_bed          = ch_mutect2_target_bed // channel: [ val(meta), path(mutect2_target_bed) ]
}
