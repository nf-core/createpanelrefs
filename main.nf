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

    main:
    // Define list of tools to run
    def tools = defineToolsList(params.tools)

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
        params.cnvkit_pon_name,
        params.gcnv_model_name,
        params.gens_pon_name,
        params.mutect2_pon_name,
        params.genome,
        params.genomes,
    )

    PREPARE_GENOME(
        tools,
        params.cnvkit_targets,
        params.dict,
        params.fai,
        params.fasta,
        params.gcnv_exclude_bed,
        params.gcnv_exclude_interval_list,
        params.gcnv_mappable_regions,
        params.gcnv_ploidy_priors,
        params.gcnv_segmental_duplications,
        params.gcnv_target_bed,
        params.gcnv_target_interval_list,
        params.gens_interval_list,
        params.mutect2_intervals_num,
        params.mutect2_target_bed,
    )

    // WORKFLOW: Run main workflow
    NFCORE_CREATEPANELREFS(
        PIPELINE_INITIALISATION.out.samplesheet,
        tools,
        params.cnvkit_pon_name,
        params.gcnv_model_name,
        params.gcnv_analysis_type,
        params.gens_analysis_type,
        params.gens_pon_name,
        params.mutect2_pon_name,
        PREPARE_GENOME.out.cnvkit_targets,
        PREPARE_GENOME.out.dict,
        PREPARE_GENOME.out.fai,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.gcnv_exclude_bed,
        PREPARE_GENOME.out.gcnv_exclude_interval_list,
        PREPARE_GENOME.out.gcnv_mappable_regions,
        PREPARE_GENOME.out.gcnv_ploidy_priors,
        PREPARE_GENOME.out.gcnv_segmental_duplications,
        PREPARE_GENOME.out.gcnv_target_bed,
        PREPARE_GENOME.out.gcnv_target_interval_list,
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

    publish:
    multiqc                       = MULTIQC.out.data.mix(MULTIQC.out.plots, MULTIQC.out.report)
    cnvkit_bed                    = NFCORE_CREATEPANELREFS.out.cnvkit_bed
    cnvkit_out                    = NFCORE_CREATEPANELREFS.out.cnvkit_out
    fasta_refs                    = PREPARE_GENOME.out.dict.mix(PREPARE_GENOME.out.fai)
    gatk4_genomicsdb              = NFCORE_CREATEPANELREFS.out.gatk4_genomicsdb
    gatk4_mutect2                 = NFCORE_CREATEPANELREFS.out.gatk4_mutect2
    gatk4_mutect2_bed             = PREPARE_GENOME.out.mutect2_target_bed
    gatk4_pon                     = NFCORE_CREATEPANELREFS.out.gatk4_pon
    gens_intervals                = NFCORE_CREATEPANELREFS.out.gens_bed.mix(PREPARE_GENOME.out.gens_interval_list)
    gens_pon                      = NFCORE_CREATEPANELREFS.out.gens_pon
    gens_read_counts              = NFCORE_CREATEPANELREFS.out.gens_read_counts
    germlinecnvcaller_cnv         = NFCORE_CREATEPANELREFS.out.germlinecnvcaller_cnv
    germlinecnvcaller_ploidy      = NFCORE_CREATEPANELREFS.out.germlinecnvcaller_ploidy_model
    germlinecnvcaller_read_counts = NFCORE_CREATEPANELREFS.out.germlinecnvcaller_read_counts
}

output {
    multiqc {
        path "reports/multiqc"
    }
    cnvkit_bed {
        path { _meta, file ->
            file >> "references/cnvkit/"
        }
    }
    cnvkit_out {
        path { _meta, file ->
            file >> "cnvkit/"
        }
    }
    fasta_refs {
        path "references/"
    }
    gatk4_genomicsdb {
        path { _meta, file ->
            file >> "gatk4/genomicsdb/"
        }
    }
    gatk4_mutect2 {
        path { _meta, file ->
            file >> "gatk4/mutect2/"
        }
    }
    gatk4_mutect2_bed {
        path "references/mutect2/"
    }
    gatk4_pon {
        path { _meta, file ->
            file >> "gatk4/createsomaticpanelofnormals/"
        }
    }
    gens_intervals {
        path "references/gens/"
    }
    gens_pon {
        path { _meta, file ->
            file >> "gens/createreadcountpanelofnormals/"
        }
    }
    gens_read_counts {
        path { _meta, file ->
            file >> "gens/readcounts/"
        }
    }
    germlinecnvcaller_cnv {
        path { _meta, file ->
            file >> "germlinecnvcaller/germlinecnvcaller/"
        }
    }
    germlinecnvcaller_ploidy {
        path { _meta, file ->
            file >> "germlinecnvcaller/determinecontigploidy/"
        }
    }
    germlinecnvcaller_read_counts {
        path { _meta, file ->
            file >> "germlinecnvcaller/readcounts/"
        }
    }
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
    intervals_num
    mutect2_target_bed // channel: [meta, mutect2_target_bed]

    main:
    // WORKFLOW: Run pipeline
    CREATEPANELREFS(
        samplesheet,
        tools,
        cnvkit_pon_name,
        gcnv_model_name,
        gcnv_analysis_type,
        gens_analysis_type,
        gens_pon_name,
        mutect2_pon_name,
        cnvkit_targets,
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
        gens_interval_list,
        intervals_num,
        mutect2_target_bed,
    )

    emit:
    cnvkit_bed                     = CREATEPANELREFS.out.cnvkit_bed
    cnvkit_out                     = CREATEPANELREFS.out.cnvkit_out
    gatk4_genomicsdb               = CREATEPANELREFS.out.gatk4_genomicsdb
    gatk4_mutect2                  = CREATEPANELREFS.out.gatk4_mutect2
    gatk4_pon                      = CREATEPANELREFS.out.gatk4_pon
    gens_bed                       = CREATEPANELREFS.out.gens_bed
    gens_pon                       = CREATEPANELREFS.out.gens_pon
    gens_read_counts               = CREATEPANELREFS.out.gens_read_counts
    germlinecnvcaller_cnv          = CREATEPANELREFS.out.germlinecnvcaller_cnv
    germlinecnvcaller_ploidy_model = CREATEPANELREFS.out.germlinecnvcaller_ploidy_model
    germlinecnvcaller_read_counts  = CREATEPANELREFS.out.germlinecnvcaller_read_counts
}

// Get workflow summary for MultiQC
def paramsSummaryMultiqc(summary_params) {
    def summary_section = ''
    summary_params
        .keySet()
        .each { group ->
            def group_params = summary_params.get(group)
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>${group}</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                group_params
                    .keySet()
                    .sort()
                    .each { param ->
                        summary_section += "        <dt>${param}</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                    }
                summary_section += "    </dl>\n"
            }
        }

    def yaml_file_text = "id: '${workflow.manifest.name.replace('/', '-')}-summary'\n" as String
    yaml_file_text     += "description: ' - this information is collected when the pipeline is started.'\n"
    yaml_file_text     += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
    yaml_file_text     += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
    yaml_file_text     += "plot_type: 'html'\n"
    yaml_file_text     += "data: |\n"
    yaml_file_text     += "${summary_section}"

    return yaml_file_text
}
