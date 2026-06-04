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
    gcnv_model_name // string: name of gcnv model
    gens_analysis_type // string: type of analysis for gens pon ('lrs' or 'srs')
    gens_pon_name // string: name of gens pon
    mutect2_pon_name // string: name of mutect2 pon
    fasta // channel: [meta, fasta]
    dict // channel: [meta, dict]
    fai // channel: [meta, fai]
    cnvkit_targets // channel: [meta, cnvkit_targets]
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
    // Auto-index alignment files if indexes are missing from the samplesheet
    PREPARE_ALIGNMENT(samplesheet, tools)

    //CNVKIT
    CNVKIT_PON(
        PREPARE_ALIGNMENT.out.reads_index.filter { 'cnvkit' in tools },
        fasta,
        cnvkit_targets,
    )

    // GENS
    GENS_PON(
        PREPARE_ALIGNMENT.out.reads_index.filter { 'gens' in tools },
        gens_analysis_type,
        gens_pon_name,
        dict,
        fai,
        fasta,
        gens_interval_list.filter { 'gens' in tools },
    )

    // GERMLINECNVCALLER
    GERMLINECNVCALLER_COHORT(
        PREPARE_ALIGNMENT.out.reads_index.filter { 'germlinecnvcaller' in tools },
        gcnv_model_name,
        dict,
        fai,
        fasta,
        gcnv_exclude_bed.filter { 'germlinecnvcaller' in tools },
        gcnv_exclude_interval_list.filter { 'germlinecnvcaller' in tools },
        gcnv_mappable_regions.filter { 'germlinecnvcaller' in tools },
        gcnv_ploidy_priors,
        gcnv_segmental_duplications.filter { 'germlinecnvcaller' in tools },
        gcnv_target_bed.filter { 'germlinecnvcaller' in tools },
        gcnv_target_interval_list.filter { 'germlinecnvcaller' in tools },
    )

    // MUTECT2
    BAM_CREATE_SOM_PON_GATK(
        PREPARE_ALIGNMENT.out.reads_index.filter { 'mutect2' in tools },
        fasta,
        fai.map { meta, fai_ -> [meta, fai_, []] },
        dict.filter { 'mutect2' in tools },
        mutect2_pon_name,
        mutect2_target_bed.filter { 'mutect2' in tools }.map { _meta, target -> [target] },
        intervals_num,
    )
}
