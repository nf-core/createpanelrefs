//
// Subworkflow with functionality specific to the nf-core/createpanelrefs pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

include { checkCondaChannels   } from 'plugin/nf-core-utils'
include { checkConfigProvided  } from 'plugin/nf-core-utils'
include { checkProfileProvided } from 'plugin/nf-core-utils'
include { completionEmail      } from 'plugin/nf-core-utils'
include { completionSummary    } from 'plugin/nf-core-utils'
include { dumpParametersToJSON } from 'plugin/nf-core-utils'
include { getReferences        } from 'plugin/nf-core-utils'
include { getReferencePath     } from 'plugin/nf-core-utils'
include { getWorkflowVersion   } from 'plugin/nf-core-utils'
include { paramsHelp           } from 'plugin/nf-schema'
include { paramsSummaryLog     } from 'plugin/nf-schema'
include { paramsSummaryMap     } from 'plugin/nf-schema'
include { samplesheetToList    } from 'plugin/nf-schema'
include { validateParameters   } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

workflow PIPELINE_INITIALISATION {
    take:
    version // boolean: Display version and exit
    validate_params // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir //  string: The output directory where the results will be saved
    input //  string: Path to input samplesheet
    help // boolean or string: Display help message (for one parameter if a string) and exit
    help_full // boolean: Show the full help message
    show_hidden // boolean: Show hidden parameters in the help message
    tools
    cnvkit_annotate
    cnvkit_method
    cnvkit_pon_name
    cnvkit_targets
    gcnv_model_name
    gens_pon_name
    mutect2_pon_name
    genome // string: iGenomes key, or datasheet name with igenomes_ignore
    genomes // map: iGenomes catalogue
    igenomes_ignore // boolean: Read --genome as a datasheet name
    references_datasheet // string: Path to a references datasheet
    references_datasheets_path // string: Location of the datasheets used with igenomes_ignore
    references_schema // string: Path to the references datasheet schema
    references_params // map: Reference params: null reads the datasheet or genome config, false turns one off
    assay_datasheet // string: Path to an assay datasheet
    assay_schema // string: Path to the assay datasheet schema
    assay_params // map: Assay params: null reads the assay datasheet, false turns one off

    main:

    //
    // Print workflow version and exit on --version
    //
    if (version) {
        log.info("${workflow.manifest.name} ${getWorkflowVersion()}")
        System.exit(0)
    }

    //
    // Dump pipeline parameters to a JSON file
    //
    if (outdir) {
        dumpParametersToJSON(outdir, params)
    }

    //
    // When running with Conda, warn if channels have not been set-up appropriately
    //
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        checkCondaChannels()
    }

    checkConfigProvided()

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/createpanelrefs ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/', '')}" }.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/createpanelrefs/blob/main/CITATIONS.md
"""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    if (help || help_full) {
        def help_options = [
            beforeText: before_text,
            afterText: after_text,
            command: command,
            showHidden: show_hidden,
            fullHelp: help_full,
            parameter: help instanceof String && help != "true" ? help : "",
        ]
        log.info(paramsHelp(help_options))
        exit(0)
    }

    checkProfileProvided(nextflow_cli_args)

    //
    // Print parameter summary to stdout. This will display the parameters
    // that differ from the default given in the JSON schema
    //

    log.info(before_text)
    log.info(paramsSummaryLog())
    log.info(after_text)

    extra_text = """
\033[1;37mExtra informations\033[0m
\033[0;34m  Tools selected to be run  :\033[0;32m ${tools.join(",")} \033[0m
-\033[2m----------------------------------------------------\033[0m-
"""

    if (monochrome_logs) {
        extra_text = extra_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    log.info(extra_text)

    //
    // Validate the parameters using nextflow_schema.json or the schema
    // given via the validation.parametersSchema configuration option
    //
    if (validate_params) {
        validateParameters([:])
    }

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters(
        cnvkit_annotate,
        cnvkit_method,
        cnvkit_pon_name,
        cnvkit_targets,
        gcnv_model_name,
        gens_pon_name,
        mutect2_pon_name,
        tools,
    )

    //
    // References come from --genome (iGenomes) or from a references datasheet; a non-null parameter wins, false turns it off
    //
    def references_path = getReferencePath(
        references: references_datasheet,
        genome: genome,
        igenomesIgnore: igenomes_ignore,
        datasheetPath: references_datasheets_path,
    )
    def references_rows = references_path ? samplesheetToList(references_path, references_schema, unrecognisedHeaders: 'skip') : null
    def references_resolved = getReferences(
        schema: references_schema,
        rows: references_rows,
        params: references_params,
        genome: genome,
        genomes: genomes,
    )

    //
    // Assay target files come from an assay datasheet built for the same genome; a non-null parameter wins, false turns it off
    //
    def assay_rows = assay_datasheet ? samplesheetToList(assay_datasheet, assay_schema, unrecognisedHeaders: 'skip') : null
    if (assay_rows) {
        assayGenomeMatchError(datasheetGenome(assay_rows), references_rows ? datasheetGenome(references_rows) : genome)
    }
    def assay_resolved = getReferences(schema: assay_schema, rows: assay_rows, params: assay_params)

    //
    // Create channel from input file provided through input
    //
    ch_samplesheet = channel.fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))

    emit:
    references = references_resolved + assay_resolved
    samplesheet = ch_samplesheet
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

workflow PIPELINE_COMPLETION {
    take:
    email //  string: email address
    email_on_fail //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.isBound() ? multiqc_reports.getVal() : [],
            )
        }

        completionSummary(monochrome_logs)
    }

    workflow.onError {
        log.error("Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

//
// Check and validate pipeline parameters
//
def validateInputParameters(cnvkit_annotate, cnvkit_method, cnvkit_pon_name, cnvkit_targets, gcnv_model_name, gens_pon_name, mutect2_pon_name, tools) {
    if ('cnvkit' in tools && cnvkit_annotate && !cnvkit_targets) {
        error("--cnvkit_annotate requires --cnvkit_targets.")
    }
    if ('cnvkit' in tools && cnvkit_method != 'wgs' && !cnvkit_targets) {
        error("--cnvkit_method ${cnvkit_method} requires --cnvkit_targets.")
    }
    checkPonName(
        tools,
        [
            cnvkit_pon_name: cnvkit_pon_name,
            gcnv_model_name: gcnv_model_name,
            gens_pon_name: gens_pon_name,
            mutect2_pon_name: mutect2_pon_name,
        ],
    )
}

//
// Genome of the first row of a datasheet whose genome column maps to meta.genome
//
def datasheetGenome(rows) {
    return rows[0][0].genome
}

//
// Exit pipeline if the references have no genome name, or the assay datasheet was built for another genome
//
def assayGenomeMatchError(assay_genome, reference_genome) {
    if (!reference_genome) {
        error("--assay needs the references from a datasheet or iGenomes: provide --references or --genome.")
    }
    if (assay_genome != reference_genome) {
        error("The assay datasheet is for genome '${assay_genome}', but the references are for genome '${reference_genome}'.")
    }
}

//
// Check PON names: error on null/empty, warn if default for a selected tool
//
def checkPonName(tools, pon_names) {
    def defaults = [cnvkit_pon_name: 'cnvkit', mutect2_pon_name: 'mutect2', gens_pon_name: 'gens', gcnv_model_name: 'germlinecnvcaller']
    pon_names.each { param, value ->
        if (!value) {
            error("--${param} is not set. Please specify a name for the panel of normals.")
        }
        if (defaults[param] in tools && value == defaults[param]) {
            log.warn("--${param} is set to the default value '${defaults[param]}'.")
        }
    }
}

//
// Define list of tools to run
//
def defineToolsList(input_tools) {
    def tools_list = input_tools ? input_tools.tokenize(',') : []
    return tools_list.sort().unique()
}
