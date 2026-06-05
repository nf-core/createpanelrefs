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
    ch_reads_index // channel: [mandatory] [ val(meta), path(bam/cram), path(bai/crai) ]
    val_analysis_type // string: [mandatory] type of analysis ('lrs' or 'srs')
    val_pon_name //  string: [optional] name for panel of normals
    ch_dict // channel: [optional] [ val(meta), path(dict) ]
    ch_fai // channel: [optional] [ val(meta), path(fai) ]
    ch_fasta // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_interval_list // channel: [mandatory] [ val(meta), path(interval_list) ]

    main:
    ch_readcounts = channel.empty()
    ch_bed = channel.empty()

    if (val_analysis_type == 'srs') {
        // Collect read counts, and generate models
        GATK4_COLLECTREADCOUNTS(
            ch_reads_index.combine(ch_interval_list.map { _meta, interval_list -> interval_list }),
            ch_fasta,
            ch_fai,
            ch_dict,
        )

        ch_readcounts = GATK4_COLLECTREADCOUNTS.out.tsv.mix(GATK4_COLLECTREADCOUNTS.out.hdf5)
    }
    else if (val_analysis_type == 'lrs') {

        INTERVAL_LIST_TO_BED(
            ch_interval_list,
            [],
            [],
        )

        ch_bed = INTERVAL_LIST_TO_BED.out.output

        // Prepare the body
        MOSDEPTH(
            ch_reads_index.combine(ch_bed).map { meta, bam, bai, _bins_meta, bins -> [meta, bam, bai, bins] },
            [[], []],
            false,
        )

        // Prepare the header
        SAMTOOLS_VIEW(
            ch_reads_index,
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
        FIND_CONCATENATE(MOSDEPTH_GATK_HEADER.out.output.join(MOSDEPTH_GATK_FORMAT.out.output).map { meta, header, body -> [meta, [header, body]] })

        ch_readcounts = FIND_CONCATENATE.out.file_out.map { meta, gatk_input ->
            return [meta, gatk_input]
        }
    }

    GATK4_CREATEREADCOUNTPANELOFNORMALS(ch_readcounts.collect { _meta, readcounts -> readcounts }.map { readcounts -> [[id: val_pon_name], readcounts] })

    emit:
    bed         = ch_bed
    pon         = GATK4_CREATEREADCOUNTPANELOFNORMALS.out.pon
    read_counts = ch_readcounts
}
