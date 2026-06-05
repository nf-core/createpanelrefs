/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BAM_CREATE_SOM_PON_GATK  } from '../subworkflows/nf-core/bam_create_som_pon_gatk'
include { CNVKIT_PON               } from '../subworkflows/local/cnvkit_pon'
include { GENS_PON                 } from '../subworkflows/local/gens_pon'
include { GERMLINECNVCALLER_COHORT } from '../subworkflows/local/germlinecnvcaller_cohort'
include { PREPARE_ALIGNMENT        } from '../subworkflows/local/prepare_alignment'

workflow CREATEPANELREFS {
    take:
    samplesheet // channel: samplesheet read in from --input
    tools // list: tools to run
    cnvkit_pon_name // string: name of cnvkit pon
    gcnv_model_name // string: name of gcnv model
    gcnv_analysis_type // string: type of analysis for germlinecnvcaller ('wes' or 'wgs')
    gens_analysis_type // string: type of analysis for gens pon ('lrs' or 'srs')
    gens_pon_name // string: name of gens pon
    mutect2_pon_name // string: name of mutect2 pon
    cnvkit_targets // channel: [meta, cnvkit_targets]
    dict // channel: [meta, dict]
    fai // channel: [meta, fai]
    fasta // channel: [meta, fasta]
    gcnv_exclude_bed // channel: [meta, gcnv_exclude_bed]
    gcnv_exclude_interval_list // channel: [meta, gcnv_exclude_interval_list]
    gcnv_mappable_regions // channel: [meta, gcnv_mappable_regions]
    gcnv_ploidy_priors // channel: [meta, gcnv_ploidy_priors]
    gcnv_segmental_duplications // channel: [meta, gcnv_segmental_duplications]
    gcnv_target_bed // channel: [meta, gcnv_target_bed]
    gcnv_target_interval_list // channel: [meta, gcnv_target_interval_list]
    gens_interval_list // channel: [meta, gens_interval_list]
    intervals_num // channel: [ path(intervals), val(num_intervals) ]
    mutect2_target_bed // channel: [meta, mutect2_target_bed]

    main:
    ch_gens_bed = channel.empty()
    ch_gens_pon = channel.empty()
    ch_gens_read_counts = channel.empty()
    ch_germlinecnvcaller_cnv = channel.empty()
    ch_germlinecnvcaller_ploidy_model = channel.empty()
    ch_germlinecnvcaller_read_counts = channel.empty()
    ch_gatk4_pon = channel.empty()
    ch_gatk4_mutect2 = channel.empty()
    ch_gatk4_genomicsdb = channel.empty()

    // Auto-index alignment files if indexes are missing from the samplesheet
    PREPARE_ALIGNMENT(samplesheet, tools)

    //CNVKIT
    CNVKIT_PON(
        PREPARE_ALIGNMENT.out.reads_index.filter { 'cnvkit' in tools },
        fasta,
        cnvkit_targets,
        cnvkit_pon_name,
    )

    // GENS
    if ('gens' in tools) {
        GENS_PON(
            PREPARE_ALIGNMENT.out.reads_index,
            gens_analysis_type,
            gens_pon_name,
            dict,
            fai,
            fasta,
            gens_interval_list,
        )

        ch_gens_bed = GENS_PON.out.bed
        ch_gens_pon = GENS_PON.out.pon
        ch_gens_read_counts = GENS_PON.out.read_counts
    }

    // GERMLINECNVCALLER
    if ('germlinecnvcaller' in tools) {
        GERMLINECNVCALLER_COHORT(
            PREPARE_ALIGNMENT.out.reads_index,
            gcnv_model_name,
            gcnv_analysis_type,
            dict,
            fai,
            fasta,
            gcnv_exclude_bed,
            gcnv_exclude_interval_list,
            gcnv_mappable_regions,
            gcnv_ploidy_priors,
            gcnv_segmental_duplications,
            gcnv_target_bed,
            gcnv_target_interval_list,
        )

        ch_germlinecnvcaller_cnv = GERMLINECNVCALLER_COHORT.out.cnv_calls.mix(GERMLINECNVCALLER_COHORT.out.cnv_model)
        ch_germlinecnvcaller_ploidy_model = GERMLINECNVCALLER_COHORT.out.ploidy_model
        ch_germlinecnvcaller_read_counts = GERMLINECNVCALLER_COHORT.out.read_counts
    }

    // MUTECT2
    if ('mutect2' in tools) {
        BAM_CREATE_SOM_PON_GATK(
            PREPARE_ALIGNMENT.out.reads_index,
            fasta,
            fai.map { meta, fai_ -> [meta, fai_, []] },
            dict,
            mutect2_pon_name,
            mutect2_target_bed.map { _meta, target -> [target] },
            intervals_num,
        )
        ch_gatk4_genomicsdb = BAM_CREATE_SOM_PON_GATK.out.genomicsdb
        ch_gatk4_mutect2 = BAM_CREATE_SOM_PON_GATK.out.mutect2_vcf.mix(BAM_CREATE_SOM_PON_GATK.out.mutect2_index, BAM_CREATE_SOM_PON_GATK.out.mutect2_stats)
        ch_gatk4_pon = BAM_CREATE_SOM_PON_GATK.out.pon_vcf.mix(BAM_CREATE_SOM_PON_GATK.out.pon_index)
    }

    emit:
    cnvkit_bed                     = CNVKIT_PON.out.bed
    cnvkit_out                     = CNVKIT_PON.out.cnn.mix(CNVKIT_PON.out.cnr)
    gens_bed                       = ch_gens_bed
    gens_pon                       = ch_gens_pon
    gens_read_counts               = ch_gens_read_counts
    germlinecnvcaller_cnv          = ch_germlinecnvcaller_cnv
    germlinecnvcaller_ploidy_model = ch_germlinecnvcaller_ploidy_model
    germlinecnvcaller_read_counts  = ch_germlinecnvcaller_read_counts
    gatk4_genomicsdb               = ch_gatk4_genomicsdb
    gatk4_mutect2                  = ch_gatk4_mutect2
    gatk4_pon                      = ch_gatk4_pon
    reads_index                    = PREPARE_ALIGNMENT.out.reads_index
}
