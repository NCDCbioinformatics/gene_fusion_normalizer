# gene_fusion_normalizer

Direction-preserving fusion-gene normalization component of the CURE-NGS
panel harmonization framework.

> **Supported deployment:** use the unified
> [CURE-NGS Docker/OCI distribution](https://github.com/NCDCbioinformatics/cure-ngs-panel-harmonization-framework).
> This repository is retained as component provenance. The supported container
> package is published only from the umbrella repository; **No packages
> published** here is therefore expected.

## Role in the unified project

| Item | Value |
| --- | --- |
| Historical responsibility | Parse and standardize fusion partners using GTF and HGNC |
| Supported command | `cure-ngs normalize-fusion` |
| Latest audited release | `gene_fusion_normalizer` / release name `gene_fusion_normalizer_0.2.1` |
| Required data | GTF containing `gene_id`/`gene_name` and HGNC complete-set TSV |

## Install the supported Docker distribution

The supported unified distribution is 0.2.6. This component can run in the
core image; no host Python installation is required. The
[Docker-only quickstart](https://github.com/NCDCbioinformatics/cure-ngs-panel-harmonization-framework/blob/v0.2.6/docs/DOCKER_ONLY_QUICKSTART.md)
explains image-contained self-tests, local output mounts, and user-selected
external resource paths. Historical component releases remain unchanged.

Install [Docker Desktop](https://docs.docker.com/desktop/) or
[Docker Engine](https://docs.docker.com/engine/install/), then pull the public
core image without a GitHub login:

```bash
docker pull ghcr.io/ncdcbioinformatics/cure-ngs-harmonizer:0.2.6-core
```

To build the identical `v0.2.6` release source instead:

```bash
git clone --branch v0.2.6 --depth 1 https://github.com/NCDCbioinformatics/cure-ngs-panel-harmonization-framework.git
cd cure-ngs-panel-harmonization-framework
docker build --file docker/Dockerfile.core --tag cure-ngs-harmonizer:0.2.6-core .
```

The supported container is the umbrella repository's audited
[`v0.2.6` distribution](https://github.com/NCDCbioinformatics/cure-ngs-panel-harmonization-framework/releases/tag/v0.2.6).

## Verify and run this capability

The reviewer walkthrough verifies directional normalization of `EML4-ALK`:

```bash
git clone --branch v0.2.6 --depth 1 https://github.com/NCDCbioinformatics/cure-ngs-panel-harmonization-framework.git
cd cure-ngs-panel-harmonization-framework
bash scripts/run_reviewer_demo.sh
```

Direct component command with bundled synthetic resources:

```bash
docker run --rm \
  --volume "$PWD/examples:/examples:ro" \
  ghcr.io/ncdcbioinformatics/cure-ngs-harmonizer:0.2.6-core \
  normalize-fusion EML4-ALK \
  --gtf /examples/synthetic/genes.gtf \
  --hgnc /examples/synthetic/hgnc.tsv
```

Ambiguous partners are reported instead of guessed, and partner direction is
retained in the normalized output.

## Historical standalone package

The `gene_fusion_normalizer_0.2.1` release remains available for provenance and
supports tabular inputs, automatic fusion-column detection, and exploded
outputs. The release asset is labelled 0.2.1 while its archived internal
metadata reports 0.2.0; the umbrella release lock preserves both facts.

## Documentation and test data

- [Project structure](https://github.com/NCDCbioinformatics/cure-ngs-panel-harmonization-framework/blob/main/docs/PROJECT_STRUCTURE.md)
- [Gene/fusion commands](https://github.com/NCDCbioinformatics/cure-ngs-panel-harmonization-framework/blob/main/docs/COMMAND_REFERENCE.md#gene-and-fusion-normalization)
- [GTF and HGNC setup](https://github.com/NCDCbioinformatics/cure-ngs-panel-harmonization-framework/blob/main/docs/REFERENCE_DATA.md#4-install-gtf-and-hgnc-resources)
- [Synthetic GTF/HGNC fixtures](https://github.com/NCDCbioinformatics/cure-ngs-panel-harmonization-framework/tree/main/examples/synthetic)
- [Clean public-image validation](https://github.com/NCDCbioinformatics/cure-ngs-panel-harmonization-framework/actions/runs/33350796468)

License: MIT. No CURE-NGS patient-level data are distributed here.
