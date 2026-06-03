include { FIND_CONCATENATE                    } from '../../../modules/nf-core/find/concatenate'
include { GATK4_COLLECTREADCOUNTS             } from '../../../modules/nf-core/gatk4/collectreadcounts'
include { GATK4_CREATEREADCOUNTPANELOFNORMALS } from '../../../modules/nf-core/gatk4/createreadcountpanelofnormals'
include { GAWK as INTERVAL_LIST_TO_BED        } from '../../../modules/nf-core/gawk'
include { GAWK as MOSDEPTH_GATK_FORMAT        } from '../../../modules/nf-core/gawk'
include { GAWK as MOSDEPTH_GATK_HEADER        } from '../../../modules/nf-core/gawk'
include { MOSDEPTH                            } from '../../../modules/nf-core/mosdepth'
include { SAMTOOLS_VIEW                       } from '../../../modules/nf-core/samtools/view'

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
    ch_readcounts_out = channel.empty()

    ch_input.set { ch_bam_bai }

    if (val_analysis_type == 'srs') {
        ch_bam_bai
            .combine(ch_interval_list.map { _meta, interval_list -> interval_list })
            .set { ch_readcounts_in }

        // Collect read counts, and generate models
        GATK4_COLLECTREADCOUNTS(ch_readcounts_in, ch_fasta, ch_fai, ch_dict)

        GATK4_COLLECTREADCOUNTS.out.tsv
            .mix(GATK4_COLLECTREADCOUNTS.out.hdf5)
            .set { ch_readcounts }
    }
    else if (val_analysis_type == 'lrs') {

        INTERVAL_LIST_TO_BED(
            ch_interval_list,
            [],
            [],
        )

        ch_bam_bai
            .combine(INTERVAL_LIST_TO_BED.out.output)
            .map { meta, bam, bai, _bins_meta, bins ->
                [meta, bam, bai, bins]
            }
            .set { ch_mosdepth_in }

        // Prepare the body
        MOSDEPTH(
            ch_mosdepth_in,
            [[], []],
            false,
        )

        // Prepare the header
        SAMTOOLS_VIEW(
            ch_bam_bai,
            [[:], [], []],
            [[:], []],
            [[:], []],
            false,
        )

        MOSDEPTH_GATK_HEADER(
            SAMTOOLS_VIEW.out.sam,
            [],
            false,
        )

        MOSDEPTH_GATK_FORMAT(
            MOSDEPTH.out.regions_bed,
            [],
            false,
        )
        // Prepare GATK inputs
        MOSDEPTH_GATK_HEADER.out.output
            .join(MOSDEPTH_GATK_FORMAT.out.output)
            .map { meta, header, body -> [meta, [header, body]] }
            .set { ch_cat_in }

        FIND_CONCATENATE(ch_cat_in)

        FIND_CONCATENATE.out.file_out
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

    emit:
    genspon    = GATK4_CREATEREADCOUNTPANELOFNORMALS.out.pon
    readcounts = ch_readcounts_out
}
