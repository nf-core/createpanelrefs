/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowCreatepanelrefs.initialise(params, log)

// Check input path parameters to see if they exist

def checkPathParamList = [
    params.fasta
]

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHECK MANDATORY PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

for (param in checkPathParamList) if (param) file(param, checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MANAGE SAMPLESHEET
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_from_samplesheet = Channel.fromSamplesheet("input")

ch_input = ch_from_samplesheet.map{meta, bam, bai, cram, crai ->
    if (bam)  return [ meta + [data_type:"bam"], bam, bai ]
    if (cram) return [ meta + [data_type:"cram"], cram, crai ]
}

// Initialize file channels based on params, defined in the params.genomes[params.genome] scope
ch_cnvkit_targets        = params.cnvkit_targets        ? Channel.fromPath(params.cnvkit_targets).map { targets -> [[id:targets.baseName],targets]}.collect()
                                                    : Channel.value([[:],[]])
ch_dict                  = params.dict                  ? Channel.fromPath(params.dict).map { dict -> [[id:dict.baseName],dict]}.collect()
                                                    : Channel.empty()
ch_exclude_bed           = params.exclude_bed           ? Channel.fromPath(params.exclude_bed).map { exclude -> [[id:exclude.baseName],exclude]}.collect()
                                                    : Channel.value([[:],[]])
ch_exclude_interval_list = params.exclude_interval_list ? Channel.fromPath(params.exclude_interval_list).map { exclude -> [[id:exclude.baseName],exclude]}.collect()
                                                    : Channel.value([[:],[]])
ch_fai                   = params.fai                   ? Channel.fromPath(params.fai).map { fai -> [[id:fai.baseName],fai]}.collect()
                                                    : Channel.empty()
ch_fasta                 = params.fasta                 ? Channel.fromPath(params.fasta).map { fasta -> [[id:fasta.baseName],fasta]}.collect()
                                                    : Channel.empty()
ch_ploidy_priors         = params.ploidy_priors         ? Channel.fromPath(params.ploidy_priors).collect()
                                                    : Channel.empty()
ch_target_bed            = params.target_bed            ? Channel.fromPath(params.target_bed).map { targets -> [[id:targets.baseName],targets]}.collect()
                                                    : Channel.value([[:],[]])
ch_target_interval_list  = params.target_interval_list  ? Channel.fromPath(params.target_interval_list).map { targets -> [[id:targets.baseName],targets]}.collect()
                                                    : Channel.value([[:],[]])

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { GERMLINECNVCALLER_COHORT    } from '../subworkflows/local/germlinecnvcaller_cohort'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { CNVKIT_BATCH                } from '../modules/nf-core/cnvkit/batch/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow CREATEPANELREFS {

    ch_versions = Channel.empty()

    if (params.tools && params.tools.split(',').contains('cnvkit')) {

        ch_input
        .map{ meta, align, index ->
            new_meta = meta + [id:"panel"]
            [new_meta, align]
        }
        .groupTuple()
        .branch{
            bam:  it[0].data_type == "bam"
        }
        .bam
        .map {meta, bam -> [ meta, [], bam ]}
        .set { ch_cnvkit_input }

        CNVKIT_BATCH ( ch_cnvkit_input, ch_fasta, [[:],[]], ch_cnvkit_targets, [[:],[]], true )
        ch_versions = ch_versions.mix(CNVKIT_BATCH.out.versions)
    }

    if (params.tools && params.tools.split(',').contains('germlinecnvcaller')) {

        GERMLINECNVCALLER_COHORT(ch_dict,
                                 ch_fai,
                                 ch_fasta,
                                 ch_input,
                                 ch_ploidy_priors,
                                 ch_target_bed,
                                 ch_target_interval_list,
                                 ch_exclude_bed,
                                 ch_exclude_interval_list)

        ch_versions = ch_versions.mix(GERMLINECNVCALLER_COHORT.out.versions)
    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    // MULTIQC
    workflow_summary    = WorkflowCreatepanelrefs.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowCreatepanelrefs.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
