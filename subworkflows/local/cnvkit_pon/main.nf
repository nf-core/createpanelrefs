include { CNVKIT_BATCH  } from '../../../modules/nf-core/cnvkit/batch'
include { SAMTOOLS_VIEW } from '../../../modules/nf-core/samtools/view'

workflow CNVKIT_PON {
    take:
    ch_reads_index // channel: [mandatory] [ val(meta), path(bam|cram), path(bai|crai) ]
    ch_fasta // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_cnvkit_targets // channel: [mandatory] [ val(meta), path(targets) ]
    pon_name // string: name of the cnvkit pon

    main:
    ch_reads = ch_reads_index.branch { meta, reads, index ->
        bam: meta.data_type == "bam"
        return [meta, reads]
        cram: meta.data_type == "cram"
        return [meta, reads, index]
    }

    SAMTOOLS_VIEW(
        ch_reads.cram,
        ch_fasta.map { meta, fasta_ -> [meta, fasta_, []] },
        [[:], []],
        [[:], []],
        false,
    )

    CNVKIT_BATCH(
        ch_reads.bam.mix(SAMTOOLS_VIEW.out.bam).map { _meta, bam -> [[id: pon_name], bam] }.groupTuple().map { meta, bam -> [meta, [], [], bam, []] },
        ch_fasta.map { meta, fasta_ -> [meta, fasta_, []] },
        ch_cnvkit_targets,
        [[:], []],
        true,
    )

    emit:
    cnn = CNVKIT_BATCH.out.cnn // channel: [ val(meta), path(cnn) ]
    bed = CNVKIT_BATCH.out.bed // channel: [ val(meta), path(bed) ]
    cnr = CNVKIT_BATCH.out.cnr // channel: [ val(meta), path(cnr) ]
    cns = CNVKIT_BATCH.out.cns // channel: [ val(meta), path(cns) ]
}
