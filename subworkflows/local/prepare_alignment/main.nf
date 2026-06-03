//
// Prepare input alignment files
//

include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index'

workflow PREPARE_ALIGNMENT {
    take:
    samplesheet // [ val(meta), path(bam), path(bai), path(cram), path(crai) ]
    tools // list: tools to run

    main:
    def index_bams = tools.any { tool -> tool in ['germlinecnvcaller', 'gens', 'mutect2'] }
    def index_crams = tools.any { tool -> tool in ['germlinecnvcaller', 'gens', 'mutect2'] }

    ch_bam = samplesheet
        .filter { _meta, bam, _bai, _cram, _crai -> bam }
        .map { meta, bam, bai, _cram, _crai -> [meta, bam, bai] }

    ch_cram = samplesheet
        .filter { _meta, _bam, _bai, cram, _crai -> cram }
        .map { meta, _bam, _bai, cram, crai -> [meta, cram, crai] }

    def input_reads = ch_bam
        .mix(ch_cram)
        .branch { meta, reads, index ->
            indexed: index
            return [meta, reads, index]
            not_indexed: !index && reads
            return [meta, reads]
        }

    def to_index = input_reads.not_indexed.filter { _meta, reads -> (reads.extension == "bam" && index_bams) || (reads.extension == "cram" && index_crams) }

    SAMTOOLS_INDEX(to_index)

    ch_reads_index = input_reads.indexed.mix(input_reads.not_indexed.filter { _meta, reads -> !((reads.extension == "bam" && index_bams) || (reads.extension == "cram" && index_crams)) }.map { meta, reads -> [meta, reads, []] }).mix(to_index.join(SAMTOOLS_INDEX.out.index, failOnMismatch: true, failOnDuplicate: true)) // [ val(meta), path(bam|cram), path(bai|crai) ]

    emit:
    reads_index = ch_reads_index // [ val(meta), path(bam|cram), path(bai|crai) ]
    bam_index = ch_reads_index.filter { _meta, reads, _index -> reads.extension == "bam" } // [ val(meta), path(bam), path(bai) ]
    cram_index = ch_reads_index.filter { _meta, reads, _index -> reads.extension == "cram" } // [ val(meta), path(cram), path(crai) ]
}
