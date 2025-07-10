include { GATK4_COLLECTREADCOUNTS             } from '../../../modules/nf-core/gatk4/collectreadcounts'
include { GATK4_CREATEREADCOUNTPANELOFNORMALS } from '../../../modules/nf-core/gatk4/createreadcountpanelofnormals'
include { SAMTOOLS_INDEX                      } from '../../../modules/nf-core/samtools/index'

workflow GENS_PON {
    take:
    ch_input         // channel: [mandatory] [ val(meta), path(bam/cram), path(bai/crai) ]
    val_pon_name     //  string: [optional] name for panel of normals
    ch_dict          // channel: [optional] [ val(meta), path(dict) ]
    ch_fai           // channel: [optional] [ val(meta), path(fai) ]
    ch_fasta         // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_interval_list // channel: [mandatory] [ val(meta), path(interval_list) ]

    main:
    versions = Channel.empty()

    // Filter out files that lack indices, and generate them
    ch_input
        .branch { meta, alignment, index ->
            alignment_with_index: index.size() > 0
            return [meta, alignment, index]
            alignment_without_index: index.size() == 0
            return [meta, alignment]
        }
        .set { ch_for_mix }

    SAMTOOLS_INDEX(ch_for_mix.alignment_without_index)

    SAMTOOLS_INDEX.out.bai
        .mix(SAMTOOLS_INDEX.out.crai)
        .set { ch_index }

    // Collect alignment files and their indices
    ch_for_mix.alignment_without_index
        .join(ch_index)
        .mix(ch_for_mix.alignment_with_index)
        .combine(ch_interval_list.map { it -> it[1] })
        .set { ch_readcounts_in }

    // Collect read counts, and generate models
    GATK4_COLLECTREADCOUNTS(ch_readcounts_in, ch_fasta, ch_fai, ch_dict)

    GATK4_COLLECTREADCOUNTS.out.tsv
        .mix(GATK4_COLLECTREADCOUNTS.out.hdf5)
        .collect { it[1] }
        .map { it ->
            return [[id: val_pon_name], it]
        }
        .set { ch_readcounts_out }

    GATK4_CREATEREADCOUNTPANELOFNORMALS(ch_readcounts_out)

    versions = versions.mix(GATK4_COLLECTREADCOUNTS.out.versions)
    versions = versions.mix(GATK4_CREATEREADCOUNTPANELOFNORMALS.out.versions)
    versions = versions.mix(SAMTOOLS_INDEX.out.versions)

    emit:
    genspon    = GATK4_CREATEREADCOUNTPANELOFNORMALS.out.pon
    readcounts = ch_readcounts_out
    versions
}
