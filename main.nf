#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/createpanelrefs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/createpanelrefs
    Website: https://nf-co.re/createpanelrefs
    Slack  : https://nfcore.slack.com/channels/createpanelrefs
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.dict                        = getGenomeAttribute('dict')
params.fai                         = getGenomeAttribute('fai')
params.fasta                       = getGenomeAttribute('fasta')
params.gcnv_exclude_bed            = getGenomeAttribute('gcnv_exclude_bed')
params.gcnv_exclude_interval_list  = getGenomeAttribute('gcnv_exclude_interval_list')
params.gcnv_mappable_regions       = getGenomeAttribute('gcnv_mappable_regions')
params.gcnv_ploidy_priors          = getGenomeAttribute('gcnv_ploidy_priors')
params.gcnv_segmental_duplications = getGenomeAttribute('gcnv_segmental_duplications')
params.gcnv_target_bed             = getGenomeAttribute('gcnv_target_bed')
params.gcnv_target_interval_list   = getGenomeAttribute('gcnv_target_interval_list')
params.gens_interval_list          = getGenomeAttribute('gens_interval_list')
params.mutect2_target_bed          = getGenomeAttribute('mutect2_target_bed')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CREATEPANELREFS         } from './workflows/createpanelrefs'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_createpanelrefs_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_createpanelrefs_pipeline'
include { PREPARE_GENOME          } from './subworkflows/local/prepare_genome'
include { MULTIQC                 } from './modules/nf-core/multiqc'
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc    } from './subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML  } from './subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText  } from './subworkflows/local/utils_nfcore_createpanelrefs_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    versions = Channel.empty()
    multiqc_files = Channel.empty()

    // Initialize file channels based on params, defined in the params.genomes[params.genome] scope
    user_dict = params.dict
        ? Channel.fromPath(params.dict).map { dict -> [[id: 'genome'], dict] }.collect()
        : Channel.empty()

    user_fai = params.fai
        ? Channel.fromPath(params.fai).map { fai -> [[id: 'genome'], fai] }.collect()
        : Channel.empty()

    fasta = params.fasta
        ? Channel.fromPath(params.fasta).map { fasta -> [[id: 'genome'], fasta] }.collect()
        : Channel.empty()

    // Initialize cnvkit specific parameters
    cnvkit_targets = params.cnvkit_targets
        ? Channel.fromPath(params.cnvkit_targets).map { targets -> [[id: 'genome'], targets] }.collect()
        : Channel.value([[id: 'genome'], []])

    // Initialize gens interval list specific parameters
    user_gens_interval_list = params.gens_interval_list
        ? Channel.fromPath(params.gens_interval_list).map { gens_interval_list -> [[id: 'genome'], gens_interval_list] }.collect()
        : Channel.empty()

    // Initialize germlinecnvcaller specific parameters
    gcnv_exclude_bed = params.gcnv_exclude_bed
        ? Channel.fromPath(params.gcnv_exclude_bed).map { exclude -> [[id: 'genome'], exclude] }.collect()
        : Channel.value([[id: 'genome'], []])
    gcnv_exclude_interval_list = params.gcnv_exclude_interval_list
        ? Channel.fromPath(params.gcnv_exclude_interval_list).map { exclude -> [[id: 'genome'], exclude] }.collect()
        : Channel.value([[id: 'genome'], []])
    gcnv_mappable_regions = params.gcnv_mappable_regions
        ? Channel.fromPath(params.gcnv_mappable_regions).collect()
        : Channel.value([[id: 'genome'], []])
    gcnv_ploidy_priors = params.gcnv_ploidy_priors
        ? Channel.fromPath(params.gcnv_ploidy_priors).collect()
        : Channel.empty()
    gcnv_target_bed = params.gcnv_target_bed
        ? Channel.fromPath(params.gcnv_target_bed).map { targets -> [[id: 'genome'], targets] }.collect()
        : Channel.value([[id: 'genome'], []])
    gcnv_target_interval_list = params.gcnv_target_interval_list
        ? Channel.fromPath(params.gcnv_target_interval_list).map { targets -> [[id: 'genome'], targets] }.collect()
        : Channel.value([[id: 'genome'], []])
    gcnv_segmental_duplications = params.gcnv_segmental_duplications
        ? Channel.fromPath(params.gcnv_segmental_duplications).collect()
        : Channel.value([[id: 'genome'], []])

    // Initialize mutect2 specific parameters
    user_mutect2_target_bed = params.mutect2_target_bed
        ? Channel.fromPath(params.mutect2_target_bed).map { targets -> [[id: 'genome'], targets] }.collect()
        : Channel.empty()

    // SUBWORKFLOW: Run initialisation tasks
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
    )

    PREPARE_GENOME(fasta, user_dict, user_fai, user_gens_interval_list, user_mutect2_target_bed, params.tools ?: "no_tools")

    dict = PREPARE_GENOME.out.dict
    fai = PREPARE_GENOME.out.fai
    gens_interval_list = PREPARE_GENOME.out.gens_interval_list
    mutect2_target_bed = PREPARE_GENOME.out.mutect2_target_bed


    multiqc_config = Channel.fromPath("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)
    multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    multiqc_logo = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)

    versions = versions.mix(PREPARE_GENOME.out.versions)

    // WORKFLOW: Run main workflow
    NFCORE_CREATEPANELREFS(
        PIPELINE_INITIALISATION.out.samplesheet,
        params.tools ?: "no_tools",
        params.gcnv_model_name,
        params.gens_pon_name,
        params.mutect2_pon_name,
        fasta,
        dict,
        fai,
        cnvkit_targets,
        gcnv_exclude_bed,
        gcnv_exclude_interval_list,
        gcnv_mappable_regions,
        gcnv_ploidy_priors,
        gcnv_segmental_duplications,
        gcnv_target_bed,
        gcnv_target_interval_list,
        gens_interval_list,
        mutect2_target_bed,
    )

    versions = versions.mix(NFCORE_CREATEPANELREFS.out.versions)

    // Collate and save software versions
    collated_versions = softwareVersionsToYAML(versions).collectFile(
        storeDir: "${params.outdir}/pipeline_info",
        name: 'nf_core_createpanelrefs_software_mqc_versions.yml',
        sort: true,
        newLine: true,
    )

    // MODULE: MultiQC
    multiqc_config = Channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )
    multiqc_custom_config = params.multiqc_config
        ? Channel.fromPath(params.multiqc_config, checkIfExists: true)
        : Channel.empty()
    multiqc_logo = params.multiqc_logo
        ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : Channel.empty()

    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    multiqc_files = multiqc_files.mix(
        workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )
    multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    methods_description = Channel.value(
        methodsDescriptionText(multiqc_custom_methods_description)
    )

    multiqc_files = multiqc_files.mix(collated_versions)
    multiqc_files = multiqc_files.mix(
        methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )

    MULTIQC(
        multiqc_files.collect(),
        multiqc_config.toList(),
        multiqc_custom_config.toList(),
        multiqc_logo.toList(),
        [],
        [],
    )

    // SUBWORKFLOW: Run completion tasks
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        MULTIQC.out.report.toList(),
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// WORKFLOW: Run main analysis pipeline depending on type of input
workflow NFCORE_CREATEPANELREFS {
    take:
    samplesheet                 // channel: samplesheet read in from --input
    tools                       // string: comma separated list of tools to run
    gcnv_model_name             // string: name of gcnv model
    gens_pon_name               // string: name of gens pon
    mutect2_pon_name            // string: name of mutect2 pon
    fasta                       // channel: [meta, fasta]
    dict                        // channel: [meta, dict]
    fai                         // channel: [meta, fai]
    cnvkit_targets              // channel: [meta, cnvkit_targets]
    gcnv_exclude_bed            // channel: [meta, gcnv_exclude_bed]
    gcnv_exclude_interval_list  // channel: [meta, gcnv_exclude_interval_list]
    gcnv_mappable_regions       // channel: [meta, gcnv_mappable_regions]
    gcnv_ploidy_priors          // channel: [meta, gcnv_ploidy_priors]
    gcnv_segmental_duplications // channel: [meta, gcnv_segmental_duplications]
    gcnv_target_bed             // channel: [meta, gcnv_target_bed]
    gcnv_target_interval_list   // channel: [meta, gcnv_target_interval_list]
    gens_interval_list          // channel: [meta, gens_interval_list]
    mutect2_target_bed          // channel: [meta, mutect2_target_bed]

    main:
    // WORKFLOW: Run pipeline
    CREATEPANELREFS(samplesheet, tools, gcnv_model_name, gens_pon_name, mutect2_pon_name, fasta, dict, fai, cnvkit_targets, gcnv_exclude_bed, gcnv_exclude_interval_list, gcnv_mappable_regions, gcnv_ploidy_priors, gcnv_segmental_duplications, gcnv_target_bed, gcnv_target_interval_list, gens_interval_list, mutect2_target_bed)

    emit:
    versions = CREATEPANELREFS.out.versions // channel: versions.yml
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEFINE FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Get attribute from genome config file e.g. fasta
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[params.genome].containsKey(attribute)) {
            return params.genomes[params.genome][attribute]
        }
    }
    return null
}
