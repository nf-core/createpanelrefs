# nf-core/createpanelrefs: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [CNVKit](#cnvkit) - Create reference files for copy number variant detection from sequencing data.
- [GATK's germlinecnvcaller](#germlinecnvcaller) - Publish read counts, ploidy and cnvcalling models that can be used to call cnv's in the case mode.
- [GATK's Mutect2](#gatk-mutect2) - Create panel of normals for somatic variant calling.
- [GENS](#gens) - Create panel of normals for read-count denoising.
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### CNVKit

<details markdown="1">
<summary>Output files</summary>

- `results/reference/cnvkit/`
  - "panel.cnn": Panel reference file containing coverage information for copy number.
  - "<SAMPLE>.antitargetcoverage.cnn": Antitarget coverage file for each sample.
  - "<SAMPLE>.targetcoverage.cnn": Target coverage file for each sample.

</details>

[CNVKit](https://cnvkit.readthedocs.io/en/stable/index.html) is a Python library and command-line software toolkit to infer and visualize copy number from high-throughput DNA sequencing data.
In this pipeline, CNVKit creates reference files that can be used for copy number variant detection.
The workflow processes normal samples to generate a reference CNN file that captures the baseline coverage patterns, which can then be used for tumor-only or tumor-normal CNV analysis in downstream applications.
The reference file contains coverage information normalized across the cohort and is essential for accurate copy number calling.

### GATK germlinecnvcaller

<details markdown="1">
<summary>Output files</summary>

- `results/germlinecnvcaller/`
  - `determinecontigploidy`
    - `cohort-model`: Contig ploidy model.
  - `germlinecnvcaller`
    - `*_model`: CNV caller model for each scattered shard.
  - `readcounts`
    - `*.hdf5|.tsv`: Read count statistics for each sample.
  - `references`
    - `*.dict`: Sequence dictionary file. This file is not published if user supplies this file to the pipeline using the `--dict` parameter.
    - `*.fai`: Fasta index file. This file is not published if user supplies this file to the pipeline using the `--fai` parameter.

</details>

[GATK](https://github.com/broadinstitute/gatk) is a toolkit which offers a wide variety of tools with a primary focus on variant discovery and genotyping. In this pipeline we have implemented GATK's germlinecnvcalling workflow for analysing a cohort of samples. The output files generated from this analysis can be used for analysing samples in case mode. For more information about the workflow and output files, see GATK's documentation [here.](https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-common-and-rare-germline-copy-number-variants)

### GATK Mutect2

<details markdown="1">
<summary>Output files</summary>

- `results/gatk4/`
  - `mutect2/`
    - `*.vcf.gz`: Compressed VCF files containing somatic variant calls for each sample.
    - `*.vcf.gz.tbi`: Tabix index files for the VCF files.
    - `*.vcf.gz.stats`: Statistics files containing detailed metrics for each sample.
    - `*.f1r2.tar.gz`: Files containing information for LearnReadOrientationModel (only output when tumor-normal pair mode is run).
  - `genomicsdb/`
    - `{pon_name}/`: GenomicsDB workspace containing all sample VCFs combined.
  - `createsomaticpanelofnormals/`
    - `{pon_name}.vcf.gz`: Final panel of normals VCF file.
    - `{pon_name}.vcf.gz.tbi`: Tabix index for the panel of normals VCF.

</details>

[GATK Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Mutect2) creates a panel of normals from normal samples for somatic variant calling. The workflow: (1) calls variants in each normal sample using Mutect2 in panel of normals mode, (2) imports all VCFs into a GenomicsDB workspace, and (3) creates a final panel of normals VCF file. This panel can be used with Mutect2 in case mode via the `--panel-of-normals` parameter to filter out common germline variants and sequencing artifacts.

### GENS

<details markdown="1">
<summary>Output files</summary>

- `results/gatk4/`
  - `collectreadcounts/`
    - `*.hdf5`: Read count data in HDF5 format for each sample.
    - `*.tsv`: Read count data in TSV format for each sample.
  - `createreadcountpanelofnormals/`
    - `{pon_name}.hdf5`: Final panel of normals file in HDF5 format.
  - `references/intervals/gens_pon/`
    - `*.interval_list`: Interval list file used for read count collection.
    - `*.bed`: BED versions of interval list file used for read count collection for long-reads.

</details>

[GENS](https://github.com/Clinical-Genomics-Lund/gens) creates a panel of normals for read-count denoising to improve somatic variant detection. The workflow: (1) indexes BAM/CRAM files if needed, (2) collects read counts at specified intervals using GATK's CollectReadCounts, and (3) creates a panel of normals using GATK's CreateReadCountPanelOfNormals. This panel can be used with GENS for somatic variant calling to reduce technical noise and improve variant detection sensitivity.

When `gens_analysis_type` is set to 'lrs', a modified version of the workflow above is run where coverage calculated by mosdepth is used instead of read counts.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
