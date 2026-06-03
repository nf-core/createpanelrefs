/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BAM_CREATE_SOM_PON_GATK  } from '../subworkflows/nf-core/bam_create_som_pon_gatk'
include { CNVKIT_BATCH             } from '../modules/nf-core/cnvkit/batch'
include { GENS_PON                 } from '../subworkflows/local/gens_pon'
include { GERMLINECNVCALLER_COHORT } from '../subworkflows/local/germlinecnvcaller_cohort'
include { PREPARE_ALIGNMENT        } from '../subworkflows/local/prepare_alignment'
include { SAMTOOLS_VIEW            } from '../modules/nf-core/samtools/view'

workflow CREATEPANELREFS {
    take:
    samplesheet // channel: samplesheet read in from --input
    tools // array: tools to run, or no_tools if none (it's actually comma separated values string, but close enough)
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
    mutect2_target_bed // channel: [meta, mutect2_target_bed]

    main:
    // Build alignment channels from samplesheet, allowing missing indexes
    ch_bam = samplesheet
        .filter { meta, bam, bai, cram, crai -> bam }
        .map { meta, bam, bai, cram, crai -> [meta, bam, bai] }

    ch_cram = samplesheet
        .filter { meta, bam, bai, cram, crai -> cram }
        .map { meta, bam, bai, cram, crai -> [meta, cram, crai] }

    // Auto-index alignment files if indexes are missing from the samplesheet
    PREPARE_ALIGNMENT(ch_bam, ch_cram, tools)

    if (tools.split(',').contains('cnvkit')) {

        input_by_fmt = PREPARE_ALIGNMENT.out.reads_index.branch { meta, alignment, index ->
            bam: alignment.extension == "bam"
            return [meta, alignment]
            cram: alignment.extension == "cram"
            return [meta, alignment, index]
        }

        cnvkit_input = SAMTOOLS_VIEW(
            input_by_fmt.cram,
            fasta.map { meta, fasta_ -> [meta, fasta_, []] },
            [[:], []],
            [[:], []],
            false,
        ).bam.mix(input_by_fmt.bam).map { meta, bam ->
            [meta + [id: 'panel'], bam]
        }.groupTuple().map { meta, bam ->
            [meta, [], [], bam, []]
        }

        CNVKIT_BATCH(
            cnvkit_input,
            fasta.map { meta, fasta_ -> [meta, fasta_, []] },
            cnvkit_targets,
            [[:], []],
            true,
        )
    }

    if (tools.split(',').contains('germlinecnvcaller')) {

        germlinecnvcaller_input = PREPARE_ALIGNMENT.out.reads_index
            .map { meta, alignment, index ->
                if (alignment.extension == "bam") {
                    return [meta + [data_type: 'bam'], alignment, index]
                }
                return [meta + [data_type: 'cram'], alignment, index]
            }

        GERMLINECNVCALLER_COHORT(
            germlinecnvcaller_input,
            gcnv_model_name,
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
    }

    if (tools.split(',').contains('mutect2')) {

        mutect2_input = PREPARE_ALIGNMENT.out.reads_index
            .map { meta, alignment, index ->
                if (alignment.extension == "bam") {
                    return [meta + [data_type: 'bam'], alignment, index, []]
                }
                return [meta + [data_type: 'cram'], alignment, index, []]
            }

        BAM_CREATE_SOM_PON_GATK(
            mutect2_input,
            fasta,
            fai.map { meta, fai_ -> [meta, fai_, []] },
            dict,
            mutect2_pon_name,
            mutect2_target_bed.map { _meta, target -> [target] },
        )
    }

    if (tools.split(',').contains('gens')) {

        gens_input = PREPARE_ALIGNMENT.out.reads_index
            .map { meta, alignment, index ->
                if (alignment.extension == "bam") {
                    return [meta + [data_type: 'bam'], alignment, index]
                }
                return [meta + [data_type: 'cram'], alignment, index]
            }

        GENS_PON(
            gens_input,
            gens_analysis_type,
            gens_pon_name,
            dict,
            fai,
            fasta,
            gens_interval_list,
        )
    }
}
