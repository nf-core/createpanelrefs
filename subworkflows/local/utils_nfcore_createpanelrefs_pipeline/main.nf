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
    help // boolean: Display help message and exit
    help_full // boolean: Show the full help message
    show_hidden // boolean: Show hidden parameters in the help message
    tools
    cnvkit_pon_name
    gcnv_model_name
    gens_pon_name
    mutect2_pon_name
    genome
    genomes

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
        ]
        log.info(
            paramsHelp(
                help_options,
                params.help instanceof String && params.help != "true" ? params.help : "",
            )
        )
        exit(0)
    }

    checkProfileProvided(nextflow_cli_args)

    //
    // Print parameter summary to stdout. This will display the parameters
    // that differ from the default given in the JSON schema
    //

    log.info(before_text)
    log.info(paramsSummaryLog([:], workflow))
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
        genome,
        genomes,
        cnvkit_pon_name,
        gcnv_model_name,
        gens_pon_name,
        mutect2_pon_name,
        tools,
    )

    //
    // Create channel from input file provided through input
    //
    ch_samplesheet = channel.fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))

    emit:
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
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
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
                multiqc_reports.getVal(),
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
def validateInputParameters(genome, genomes, cnvkit_pon_name, gcnv_model_name, gens_pon_name, mutect2_pon_name, tools) {
    genomeExistsError(genome, genomes)
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
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError(genome, genomes) {
    if (genomes && genome && !genomes.containsKey(genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" + "  Genome '${genome}' not found in any config files provided to the pipeline.\n" + "  Currently, the available genome keys are:\n" + "  ${genomes.keySet().join(", ")}\n" + "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
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

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
        "Tools used in the workflow included:",
        "MultiQC (Ewels et al. 2016)",
        ".",
    ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = ["<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    }
    else {
        meta["doi_text"] = ""
    }
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine = new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
