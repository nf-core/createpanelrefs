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
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CREATEPANELREFS         } from './workflows/createpanelrefs'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_createpanelrefs_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_createpanelrefs_pipeline'
include { PREPARE_GENOME          } from './subworkflows/local/prepare_genome'
include { MULTIQC                 } from './modules/nf-core/multiqc'
include { defineToolsList         } from './subworkflows/local/utils_nfcore_createpanelrefs_pipeline'
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc    } from './subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText  } from './subworkflows/local/utils_nfcore_createpanelrefs_pipeline'
include { getGenomeAttribute      } from 'plugin/nf-core-utils'
include { softwareVersionsToYAML  } from 'plugin/nf-core-utils'

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

    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    // Define list of tools to run
    def tools = defineToolsList(params.tools)

    // Initialize file channels based on params, defined in the params.genomes[params.genome] scope
    user_dict = params.dict
        ? channel.fromPath(params.dict).map { dict -> [[id: 'genome'], dict] }.collect()
        : channel.empty()

    user_fai = params.fai
        ? channel.fromPath(params.fai).map { fai -> [[id: 'genome'], fai] }.collect()
        : channel.empty()

    fasta = params.fasta
        ? channel.fromPath(params.fasta).map { fasta -> [[id: 'genome'], fasta] }.collect()
        : channel.empty()

    // Initialize cnvkit specific parameters
    cnvkit_targets = params.cnvkit_targets
        ? channel.fromPath(params.cnvkit_targets).map { targets -> [[id: 'genome'], targets] }.collect()
        : channel.value([[id: 'genome'], []])

    // Initialize gens interval list specific parameters
    user_gens_interval_list = params.gens_interval_list
        ? channel.fromPath(params.gens_interval_list).map { gens_interval_list -> [[id: 'genome'], gens_interval_list] }.collect()
        : channel.empty()

    // Initialize germlinecnvcaller specific parameters
    gcnv_exclude_bed = params.gcnv_exclude_bed
        ? channel.fromPath(params.gcnv_exclude_bed).map { exclude -> [[id: 'genome'], exclude] }.collect()
        : channel.value([[id: 'genome'], []])
    gcnv_exclude_interval_list = params.gcnv_exclude_interval_list
        ? channel.fromPath(params.gcnv_exclude_interval_list).map { exclude -> [[id: 'genome'], exclude] }.collect()
        : channel.value([[id: 'genome'], []])
    gcnv_mappable_regions = params.gcnv_mappable_regions
        ? channel.fromPath(params.gcnv_mappable_regions).collect()
        : channel.value([[id: 'genome'], []])
    gcnv_ploidy_priors = params.gcnv_ploidy_priors
        ? channel.fromPath(params.gcnv_ploidy_priors).collect()
        : channel.empty()
    gcnv_target_bed = params.gcnv_target_bed
        ? channel.fromPath(params.gcnv_target_bed).map { targets -> [[id: 'genome'], targets] }.collect()
        : channel.value([[id: 'genome'], []])
    gcnv_target_interval_list = params.gcnv_target_interval_list
        ? channel.fromPath(params.gcnv_target_interval_list).map { targets -> [[id: 'genome'], targets] }.collect()
        : channel.value([[id: 'genome'], []])
    gcnv_segmental_duplications = params.gcnv_segmental_duplications
        ? channel.fromPath(params.gcnv_segmental_duplications).collect()
        : channel.value([[id: 'genome'], []])

    // Initialize mutect2 specific parameters
    user_mutect2_target_bed = params.mutect2_target_bed
        ? channel.fromPath(params.mutect2_target_bed).map { targets -> [[id: 'genome'], targets] }.collect()
        : channel.empty()

    // SUBWORKFLOW: Run initialisation tasks
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden,
        tools,
        params.mutect2_pon_name,
        params.genome,
        params.genomes,
    )

    PREPARE_GENOME(
        fasta,
        user_dict,
        user_fai,
        user_gens_interval_list,
        user_mutect2_target_bed,
        params.mutect2_intervals_num,
        tools,
    )

    // WORKFLOW: Run main workflow
    NFCORE_CREATEPANELREFS(
        PIPELINE_INITIALISATION.out.samplesheet,
        tools,
        params.gcnv_model_name,
        params.gcnv_analysis_type,
        params.gens_analysis_type,
        params.gens_pon_name,
        params.mutect2_pon_name,
        fasta,
        PREPARE_GENOME.out.dict,
        PREPARE_GENOME.out.fai,
        cnvkit_targets,
        gcnv_exclude_bed,
        gcnv_exclude_interval_list,
        gcnv_mappable_regions,
        gcnv_ploidy_priors,
        gcnv_segmental_duplications,
        gcnv_target_bed,
        gcnv_target_interval_list,
        PREPARE_GENOME.out.gens_interval_list,
        PREPARE_GENOME.out.intervals_num,
        PREPARE_GENOME.out.mutect2_target_bed,
    )

    def collated_versions = softwareVersionsToYAML(
        softwareVersions: channel.topic("versions"),
        nextflowVersion: workflow.nextflow.version,
    ).collectFile(
        storeDir: "${params.outdir}/pipeline_info",
        name: 'nf_core_' + 'createpanelrefs_software_' + 'mqc_' + 'versions.yml',
        sort: true,
        newLine: true,
    )

    def collated_reports = channel.topic("multiqc_files")
        .map { _meta, _process, _tool, reports -> reports }

    // MODULE: MultiQC
    // Present summary of reads, alignment, duplicates, BSQR stats for all samples as well as workflow summary/parameters as single report
    def multiqc_report = channel.empty()

    // MULTIQC
    def multiqc_files = channel.empty()

    multiqc_files = multiqc_files.mix(collated_versions)
    multiqc_files = multiqc_files.mix(collated_reports)

    def summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    def multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def methods_description = channel.value(methodsDescriptionText(multiqc_custom_methods_description))

    multiqc_files = multiqc_files.mix(workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    multiqc_files = multiqc_files.mix(methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

    MULTIQC(
        multiqc_files.flatten().collect().map { files ->
            [
                [id: 'createpanelrefs'],
                files,
                params.multiqc_config
                    ? file(params.multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                params.multiqc_logo ? file(params.multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList()

    // SUBWORKFLOW: Run completion tasks
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
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
    samplesheet // channel: samplesheet read in from --input
    tools // list: tools to run
    gcnv_model_name // string: name of gcnv model
    gcnv_analysis_type // string: type of analysis for germlinecnvcaller ('wes' or 'wgs')
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
    mutect2_intervals_num
    mutect2_target_bed // channel: [meta, mutect2_target_bed]

    main:
    // WORKFLOW: Run pipeline
    CREATEPANELREFS(
        samplesheet,
        tools,
        gcnv_model_name,
        gcnv_analysis_type,
        gens_analysis_type,
        gens_pon_name,
        mutect2_pon_name,
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
        mutect2_intervals_num,
        mutect2_target_bed,
    )

    emit:
    cnvkit_bed                     = CREATEPANELREFS.out.cnvkit_bed
    cnvkit_cnn                     = CREATEPANELREFS.out.cnvkit_cnn
    cnvkit_cnr                     = CREATEPANELREFS.out.cnvkit_cnr
    gens_pon                       = CREATEPANELREFS.out.gens_pon
    gens_read_counts               = CREATEPANELREFS.out.gens_read_counts
    germlinecnvcaller_cnv_model    = CREATEPANELREFS.out.germlinecnvcaller_cnv_model
    germlinecnvcaller_ploidy_model = CREATEPANELREFS.out.germlinecnvcaller_ploidy_model
    germlinecnvcaller_read_counts  = CREATEPANELREFS.out.germlinecnvcaller_read_counts
    som_pon_gatk_genomicsdb        = CREATEPANELREFS.out.som_pon_gatk_genomicsdb
    som_pon_gatk_index             = CREATEPANELREFS.out.som_pon_gatk_index
    som_pon_gatk_mutect2_stats     = CREATEPANELREFS.out.som_pon_gatk_mutect2_stats
    som_pon_gatk_vcf               = CREATEPANELREFS.out.som_pon_gatk_vcf
}
