//
// Prepare input alignment files
//

include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index'

workflow PREPARE_ALIGNMENT {
    take:
    bam // [ val(meta), path(bam), path(bai) ]
    cram // [ val(meta), path(cram), path(crai) ]
    tools // string: comma-separated list of tools

    main:
    def index_bams = tools.split(',').find { it in ['germlinecnvcaller', 'gens', 'mutect2'] } != null
    def index_crams = tools.split(',').find { it in ['germlinecnvcaller', 'gens', 'mutect2'] } != null

    def input_reads = bam
        .mix(cram)
        .branch { meta, reads, index ->
            indexed: index
            return [meta, reads, index]
            not_indexed: !index && reads
            return [meta, reads]
        }

    def to_index = input_reads.not_indexed.filter { meta, reads ->
        (reads.extension == "bam" && index_bams) || (reads.extension == "cram" && index_crams)
    }

    SAMTOOLS_INDEX(to_index)

    emit:
    reads_index = input_reads.indexed
        .mix(input_reads.not_indexed.filter { meta, reads ->
            !((reads.extension == "bam" && index_bams) || (reads.extension == "cram" && index_crams))
        }.map { meta, reads -> [meta, reads, []] })
        .mix(to_index.join(SAMTOOLS_INDEX.out.index, failOnMismatch: true, failOnDuplicate: true)) // [ val(meta), path(bam|cram), path(bai|crai) ]
}
