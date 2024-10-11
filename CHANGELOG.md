# nf-core/createpanelrefs: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [date]

Initial release of nf-core/createpanelrefs, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- [#5](https://github.com/nf-core/createpanelrefs/pull/5) - `CNVKIT` can be used to create a PON
- [#5](https://github.com/nf-core/createpanelrefs/pull/5) - Usage of nf-validation
- [#5](https://github.com/nf-core/createpanelrefs/pull/5) - Usage of nf-test
- [#8](https://github.com/nf-core/createpanelrefs/pull/8) - `Mutect2` can be used to create a PON
- [#10](https://github.com/nf-core/createpanelrefs/pull/10) - `GATK germlinecnvcaller` can be used to create a PON
- [#17](https://github.com/nf-core/createpanelrefs/pull/17) - `GENS` can be used to create a PON

### `Updated`

- [#19](https://github.com/nf-core/createpanelrefs/pull/19) - Updates germlinecnvcaller subworkflow to handle exome samples
- [#24](https://github.com/nf-core/createpanelrefs/pull/24) - Updates germlinecnvcaller subworkflow to use mappability and segmental duplications track
- [#24](https://github.com/nf-core/createpanelrefs/pull/24) - Updates germlinecnvcaller and gens subworkflows to use custom names for panel of normals.
- [#28](https://github.com/nf-core/createpanelrefs/pull/28) - Updates default args for gens subworkflow and made the parameters available from the command line.
- [#31](https://github.com/nf-core/createpanelrefs/pull/31) - Publish interval_list file from gens subworkflow by default.
- [#35](https://github.com/nf-core/createpanelrefs/pull/35) - Template update for nf-core/tools v3.0.2

### `Fixed`

### `Dependencies`

### `Deprecated`
