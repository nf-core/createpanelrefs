/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BAM_CREATE_SOM_PON_GATK  } from '../subworkflows/nf-core/bam_create_som_pon_gatk'
include { CNVKIT_BATCH             } from '../modules/nf-core/cnvkit/batch'
include { GENS_PON                 } from '../subworkflows/local/gens_pon'
include { GERMLINECNVCALLER_COHORT } from '../subworkflows/local/germlinecnvcaller_cohort'
include { SAMTOOLS_VIEW            } from '../modules/nf-core/samtools/view'

workflow CREATEPANELREFS {
    take:
    samplesheet                 // channel: samplesheet read in from --input
    tools                       // array: tools to run, or no_tools if none (it's actually comma separated values string, but close enough)
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
    versions = Channel.empty()

    if (tools.split(',').contains('cnvkit')) {

        input_by_fmt = samplesheet.branch { meta, bam, _bai, cram, crai ->
            bam: bam
            return [meta, bam]
            cram: cram
            return [meta, cram, crai]
        }

        cnvkit_input = SAMTOOLS_VIEW(input_by_fmt.cram, fasta, [], "").bam
            .mix(input_by_fmt.bam)
            .map { meta, bam ->
                return [meta + [id: 'panel'], bam]
            }
            .groupTuple()
            .map { meta, bam ->
                return [meta, [], bam]
            }

        CNVKIT_BATCH(cnvkit_input, fasta, [[:], []], cnvkit_targets, [[:], []], true)

        versions = versions.mix(CNVKIT_BATCH.out.versions)
    }

    if (tools.split(',').contains('germlinecnvcaller')) {

        germlinecnvcaller_input = samplesheet.map { meta, bam, bai, cram, crai ->
            if (bam) {
                return [meta + [data_type: 'bam'], bam, bai]
            }
            if (cram) {
                return [meta + [data_type: 'cram'], cram, crai]
            }
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

        versions = versions.mix(GERMLINECNVCALLER_COHORT.out.versions)
    }

    if (tools.split(',').contains('mutect2')) {

        mutect2_input = samplesheet.map { meta, bam, bai, cram, crai ->
            if (bam) {
                return [meta + [data_type: 'bam'], bam, bai, []]
            }
            if (cram) {
                return [meta + [data_type: 'cram'], cram, crai, []]
            }
        }

        BAM_CREATE_SOM_PON_GATK(
            mutect2_input,
            fasta,
            fai,
            dict,
            mutect2_pon_name,
            mutect2_target_bed.map { _meta, target -> [target] },
        )

        versions = versions.mix(BAM_CREATE_SOM_PON_GATK.out.versions)
    }

    if (tools.split(',').contains('gens')) {

        gens_input = samplesheet.map { meta, bam, bai, cram, crai ->
            if (bam) {
                return [meta + [data_type: 'bam'], bam, bai]
            }
            if (cram) {
                return [meta + [data_type: 'cram'], cram, crai]
            }
        }

        GENS_PON(
            gens_input,
            gens_pon_name,
            dict,
            fai,
            fasta,
            gens_interval_list,
        )

        versions = versions.mix(GENS_PON.out.versions)
    }

    emit:
    versions // channel: [ path(versions.yml) ]
}
