include { GATK4_COLLECTREADCOUNTS             } from '../../../modules/nf-core/gatk4/collectreadcounts'
include { GATK4_CREATEREADCOUNTPANELOFNORMALS } from '../../../modules/nf-core/gatk4/createreadcountpanelofnormals'
include { SAMTOOLS_INDEX                      } from '../../../modules/nf-core/samtools/index'
include { SAMTOOLS_VIEW                       } from '../../../modules/nf-core/samtools/view'
include { MOSDEPTH                            } from '../../../modules/nf-core/mosdepth'
include { GAWK as MOSDEPTH_GATK_HEADER        } from '../../../modules/nf-core/gawk'
include { GAWK as MOSDEPTH_GATK_FORMAT        } from '../../../modules/nf-core/gawk'
include { GAWK as INTERVAL_LIST_TO_BED        } from '../../../modules/nf-core/gawk'
include { CAT_CAT                             } from '../../../modules/nf-core/cat/cat'

workflow GENS_PON {
    take:
    ch_input // channel: [mandatory] [ val(meta), path(bam/cram), path(bai/crai) ]
    val_analysis_type // string: [mandatory] type of analysis ('lrs' or 'srs')
    val_pon_name //  string: [optional] name for panel of normals
    ch_dict // channel: [optional] [ val(meta), path(dict) ]
    ch_fai // channel: [optional] [ val(meta), path(fai) ]
    ch_fasta // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_interval_list // channel: [mandatory] [ val(meta), path(interval_list) ]

    main:
    versions = channel.empty()
    ch_readcounts_out = channel.empty()

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
    versions = versions.mix(SAMTOOLS_INDEX.out.versions)

    SAMTOOLS_INDEX.out.bai
        .mix(SAMTOOLS_INDEX.out.crai)
        .set { ch_index }

    // Collect alignment files and their indices
    ch_for_mix.alignment_without_index
        .join(ch_index)
        .mix(ch_for_mix.alignment_with_index)
        .set { ch_bam_bai }

    if (val_analysis_type == 'srs') {
        ch_bam_bai
            .combine(ch_interval_list.map { _meta, interval_list -> interval_list })
            .set { ch_readcounts_in }

        // Collect read counts, and generate models
        GATK4_COLLECTREADCOUNTS(ch_readcounts_in, ch_fasta, ch_fai, ch_dict)
        versions = versions.mix(GATK4_COLLECTREADCOUNTS.out.versions)

        GATK4_COLLECTREADCOUNTS.out.tsv
            .mix(GATK4_COLLECTREADCOUNTS.out.hdf5)
            .set { ch_readcounts }

    } else if (val_analysis_type == 'lrs') {

        INTERVAL_LIST_TO_BED(
            ch_interval_list, [], []
        )
        versions = versions.mix(INTERVAL_LIST_TO_BED.out.versions)

        ch_bam_bai
            .combine(INTERVAL_LIST_TO_BED.out.output)
            .map { meta, bam, bai, _bins_meta, bins ->
            [meta, bam, bai, bins]
        }
        .set { ch_mosdepth_in }

        // Prepare the body
        MOSDEPTH(
            ch_mosdepth_in,
            [[],[]]
        )

        // Prepare the header
        SAMTOOLS_VIEW(
            ch_bam_bai,
            [[],[]],
            [],
            false
        )
        versions = versions.mix(SAMTOOLS_VIEW.out.versions)

        MOSDEPTH_GATK_HEADER(
            SAMTOOLS_VIEW.out.sam,
            [],
            false
        )
        versions = versions.mix(MOSDEPTH_GATK_HEADER.out.versions)


        MOSDEPTH_GATK_FORMAT(
            MOSDEPTH.out.regions_bed,
            [],
            false
        )
        versions = versions.mix(MOSDEPTH_GATK_FORMAT.out.versions)

        // Prepare GATK inputs
        MOSDEPTH_GATK_HEADER.out.output.view()
            .join(MOSDEPTH_GATK_FORMAT.out.output)
            .map { meta, header, body -> [meta, [header, body]] }
            .set { ch_cat_in }

        CAT_CAT(ch_cat_in)

        CAT_CAT.out.file_out
            .map { meta, gatk_input ->
                return [meta, gatk_input]
            }
            .set { ch_readcounts }

    }

    ch_readcounts
        .collect { _meta, readcounts -> readcounts }
        .map { readcounts -> [[id: val_pon_name], readcounts] }
        .set { ch_create_pon_in }

    GATK4_CREATEREADCOUNTPANELOFNORMALS(ch_create_pon_in)

    versions = versions.mix(GATK4_CREATEREADCOUNTPANELOFNORMALS.out.versions)

    emit:
    genspon    = GATK4_CREATEREADCOUNTPANELOFNORMALS.out.pon
    readcounts = ch_readcounts_out
    versions
}
