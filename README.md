<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-createpanelrefs_logo_dark.png">
    <img alt="nf-core/createpanelrefs" src="docs/images/nf-core-createpanelrefs_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/createpanelrefs/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/createpanelrefs/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/createpanelrefs/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/createpanelrefs/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/createpanelrefs/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.10.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/createpanelrefs)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23createpanelrefs-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/createpanelrefs)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)
[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)
[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/createpanelrefs** is a bioinformatics helper pipeline that will help in creating panel of normals and other models.

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Build Panel of Normals for [`CNVKIT`](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004873)
3. Build ploidy and cnv calling models for [`GATK's germlinecnvcaller workflow`](https://genome.cshlp.org/content/20/9/1297)
4. Build Panel of Normals for [`GENS`](https://github.com/Clinical-Genomics-Lund/gens)
5. Build Panel of Normals for [`Mutect2`](https://genome.cshlp.org/content/20/9/1297)
6. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,bam,bai,cram,crai
sample1,sample1.bam,sample1.bai,,
sample2,sample2.bam,,,
sample3,sample3.bam,sample3.bai,,
sample4,sample4.bam,,,
```

Each row in the samplesheet represents an alignment file, and it is important that you provide the files in the right format for the analysis you want to run.

| Tool              | Alignment format             |
| ----------------- | ---------------------------- |
| cnvkit            | bam                          |
| germlinecnvcaller | bam or cram or a mix of both |

Now, you can run the pipeline using:

```bash
nextflow run nf-core/createpanelrefs \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --tools <cnvkit/germlinecnvcaller> \
   --genome GATK.GRCh38 \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/createpanelrefs/usage) and the [parameter documentation](https://nf-co.re/createpanelrefs/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/createpanelrefs/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/createpanelrefs/output).

## Credits

nf-core/createpanelrefs was originally written by @maxulysse.
@marrip contributed in the idea that started it all.
@matthdsm and @FriederikeHanssen contributed in the actual design.
@ramprasadn's interest was the final push that led to the creation.

We thank the following people for their extensive assistance in the development of this pipeline:

- @jfy133
- @JoseEspinosa

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#createpanelrefs` channel](https://nfcore.slack.com/channels/createpanelrefs) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/createpanelrefs for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
