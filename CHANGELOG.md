# nf-core/createpanelrefs: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0](https://github.com/nf-core/createpanelrefs/releases/tag/1.0.0) - Hell's Gate

Hell's Gate National Park is a national park situated near Lake Naivasha in Kenya.
Initial release of nf-core/createpanelrefs, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- [#5](https://github.com/nf-core/createpanelrefs/pull/5) - `CNVKIT` can be used to create a PON
- [#5](https://github.com/nf-core/createpanelrefs/pull/5) - Usage of nf-validation
- [#5](https://github.com/nf-core/createpanelrefs/pull/5) - Usage of nf-test
- [#8](https://github.com/nf-core/createpanelrefs/pull/8) - `Mutect2` can be used to create a PON
- [#10](https://github.com/nf-core/createpanelrefs/pull/10) - `GATK germlinecnvcaller` can be used to create a PON
- [#17](https://github.com/nf-core/createpanelrefs/pull/17) - `GENS` can be used to create a PON
- [#50](https://github.com/nf-core/createpanelrefs/pull/50) - Add auto creation of interval_list file from gens, and bed file for mutect2
- [#62](https://github.com/nf-core/createpanelrefs/pull/62) - Add megatests
- [#78](https://github.com/nf-core/createpanelrefs/pull/78) - Add `mutect2_intervals_num` params
- [#83](https://github.com/nf-core/createpanelrefs/pull/83) - Created a local `cnvkit_pon` subworkflow encapsulating CRAM-to-BAM conversion + PON generation

### `Changed`

- [#19](https://github.com/nf-core/createpanelrefs/pull/19) - Updates germlinecnvcaller subworkflow to handle exome samples
- [#24](https://github.com/nf-core/createpanelrefs/pull/24) - Updates germlinecnvcaller and gens subworkflows to use custom names for panel of normals.
- [#24](https://github.com/nf-core/createpanelrefs/pull/24) - Updates germlinecnvcaller subworkflow to use mappability and segmental duplications track
- [#28](https://github.com/nf-core/createpanelrefs/pull/28) - Updates default args for gens subworkflow and made the parameters available from the command line.
- [#31](https://github.com/nf-core/createpanelrefs/pull/31) - Publish interval_list file from gens subworkflow by default.
- [#35](https://github.com/nf-core/createpanelrefs/pull/35) - Improve pipeline level tests
- [#35](https://github.com/nf-core/createpanelrefs/pull/35) - Template update for nf-core/tools v3.0.2
- [#48](https://github.com/nf-core/createpanelrefs/pull/48) - Improve CI (early failure + automatic nf-test shards + [RunsOn](https://runs-on.com/))
- [#49](https://github.com/nf-core/createpanelrefs/pull/49) - Improve CI (Test Mutect2 with CRAM + better usage of test references)
- [#49](https://github.com/nf-core/createpanelrefs/pull/49) - Move all parameters in the schema that are references in the references section
- [#50](https://github.com/nf-core/createpanelrefs/pull/50) - Heavy refactoring of the pipeline
- [#50](https://github.com/nf-core/createpanelrefs/pull/50) - Improve references related files handling
- [#52](https://github.com/nf-core/createpanelrefs/pull/52) - Template update for nf-core/tools v3.2.1
- [#54](https://github.com/nf-core/createpanelrefs/pull/54) - Template update for nf-core/tools v3.3.1
- [#54](https://github.com/nf-core/createpanelrefs/pull/54) - Update nft-utils to 0.0.4
- [#55](https://github.com/nf-core/createpanelrefs/pull/55) - Prepare relase 1.0.0
- [#63](https://github.com/nf-core/createpanelrefs/pull/63) - Template update for nf-core/tools v3.5.0dev
- [#66](https://github.com/nf-core/createpanelrefs/pull/66) - Update `GENS` to allow for creating a long-read PON
- [#69](https://github.com/nf-core/createpanelrefs/pull/69) - Replace `CAT_CAT` by `FIND_CONCATENATE`
- [#69](https://github.com/nf-core/createpanelrefs/pull/69) - Update all dependencies (modules, subworfklows and plugins)
- [#74](https://github.com/nf-core/createpanelrefs/pull/74) - Update all modules to work with singularity and apptainer
- [#76](https://github.com/nf-core/createpanelrefs/pull/76) - Template update for nf-core/tools v4.0.1
- [#78](https://github.com/nf-core/createpanelrefs/pull/78) - Add intervals for Mutect2 in PON creation
- [#78](https://github.com/nf-core/createpanelrefs/pull/78) - Update MultiQC
- [#79](https://github.com/nf-core/createpanelrefs/pull/79) - Extract alignment indexing into a dedicated `prepare_alignment` subworkflow used across all tools
- [#80](https://github.com/nf-core/createpanelrefs/pull/80) - Update all modules/subworkflows to latest
- [#81](https://github.com/nf-core/createpanelrefs/pull/81) - Refactored tools handling to use a centralized `defineToolsList()` function, aligned with the pattern used in nf-core/rnavar and nf-core/seqinspector
- [#83](https://github.com/nf-core/createpanelrefs/pull/83) - Fixed `gens_pon` bug where `ch_readcounts_out` was initialized but never assigned
- [#83](https://github.com/nf-core/createpanelrefs/pull/83) - Normalized `take:` parameter names to `ch_` prefix and emit names to snake_case across all local subworkflows
- [#83](https://github.com/nf-core/createpanelrefs/pull/83) - Removed all `.set` operators
- [#83](https://github.com/nf-core/createpanelrefs/pull/83) - Removed time requirements for Mutect2 as intervals usage is now possible via [#78](https://github.com/nf-core/createpanelrefs/pull/78)
- [#83](https://github.com/nf-core/createpanelrefs/pull/83) - Simplified `prepare_alignment` subworkflow usage
- [#84](https://github.com/nf-core/createpanelrefs/pull/84) - Removed `ext.when` from `germlinecnvcaller_cohort` config and moved conditional logic into the subworkflow via channel filtering
- [#85](https://github.com/nf-core/createpanelrefs/pull/85) - Refactored `prepare_genome` subworkflow to consolidate reference file channel creation
- [#85](https://github.com/nf-core/createpanelrefs/pull/85) - Update documentation
- [#85](https://github.com/nf-core/createpanelrefs/pull/85) - Update minimum Nextflow version from 25.10.4 to 26.04.0
- [#85](https://github.com/nf-core/createpanelrefs/pull/85) - Workflow outputs for all files
- [#88](https://github.com/nf-core/createpanelrefs/pull/88) - Added default values for `cnvkit_pon_name`, `gens_pon_name`, `gcnv_model_name` and `mutect2_pon_name`, warn when they are unchanged, err when null.
- [#88](https://github.com/nf-core/createpanelrefs/pull/88) - Remove all nf-core utils subworkflows in favor of the nf-core-utils plugin
- [#88](https://github.com/nf-core/createpanelrefs/pull/88) - Update nf-core-utils to 0.5.0 and nf-schema to 2.7.2

### `Fixed`

- [#50](https://github.com/nf-core/createpanelrefs/pull/50) - Fix mutect2 that wasn't working without a bed file
- [#53](https://github.com/nf-core/createpanelrefs/pull/53) - Minor syntax fixes due to [#50](https://github.com/nf-core/createpanelrefs/pull/50)
- [#54](https://github.com/nf-core/createpanelrefs/pull/54) - Fix name for `_mqc_versions.yml` file
- [#56](https://github.com/nf-core/createpanelrefs/pull/56) - Fix gcnv interval list
- [#57](https://github.com/nf-core/createpanelrefs/pull/57) - Fix missing documentation for GATK Mutect2 and GENS
- [#57](https://github.com/nf-core/createpanelrefs/pull/57) - Improve syntax in `assets/schema_input.json` file, from @nvnieuwk in [#46](https://github.com/nf-core/createpanelrefs/pull/46)
- [#70](https://github.com/nf-core/createpanelrefs/pull/70) - Fix CI issues with conda
- [#71](https://github.com/nf-core/createpanelrefs/pull/71) - Fix time resource requirement
- [#72](https://github.com/nf-core/createpanelrefs/pull/72) - Adjust time resource requirement for MUTECT2
- [#73](https://github.com/nf-core/createpanelrefs/pull/73) - More time for Mutect2
- [#82](https://github.com/nf-core/createpanelrefs/pull/82) - Fix schema requiring `mutect2_pon_name` even when not running mutect2; replaced with runtime validation
- [#83](https://github.com/nf-core/createpanelrefs/pull/83) - Fixed `gens_pon` bug where `ch_readcounts_out` was initialized but never assigned

### `Dependencies` - modules

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| cnvkit     |             | 0.9.12      |
| - htslib   |             | 1.21        |
| - samtools |             | 1.21        |
| gatk4      |             | 4.6.2.0     |
| gawk       |             | 5.3.1       |
| mosdepth   |             | 0.3.14      |
| - htslib   |             | 1.23.1      |
| multiqc    |             | 1.35        |
| samtools   |             | 1.23.1      |
| - htslib   |             | 1.23.1      |

### `Dependencies` - Nextflow plugins

| Dependency    | Old version | New version |
| ------------- | ----------- | ----------- |
| nf-core-utils |             | 0.5.0       |
| nf-schema     |             | 2.7.2       |

### `Deprecated`
