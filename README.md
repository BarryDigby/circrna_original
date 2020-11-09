<<<<<<< HEAD
# nf-circ
Nextflow pipeline to scan high throughput sequencing data for circRNAs, perform differential expression and generate circRNA-miRNA targets.

# Notes
The pipeline uses RNA-Seq data to scan for the presence of circRNAs using a combination of `CIRCexplorer2`, `CIRIquant`, `DCC`, `find_circ`, `circrna_finder`, `mapsplice`. There is an option to use all 6 tools `--tool combine` or to use a tool on its own `--tool CIRIquant`. Currently developing the tool using a small test dataset, which works for the combine flag however using just one tool on the toy dataset can result in very few circRNAs being called which is an issue for DESeq2.

***

Differential expression -- the contrast of interest must be labelled as `condition` as this is hardcoded into the automated DE analysis script. The reference level, i.e the wild type condition to compare against (tumor vs **normal**, treated vs **untreated**) **must be called normal**. 

For example:

| sample   	| condition 	|
|----------	|-----------	|
| rep1_trt 	| treated   	|
| rep2_trt 	| treated   	|
| rep3_trt 	| treated   	|
| rep1     	| normal    	|
| rep2     	| normal    	|
| rep3     	| normal    	|

***

# TO-DO
Facilitate single end data. 
Incorporate 'stranded' datasets.

# TOOLS
The container hosts the following circRNA tools:

- CIRIquant (CIRI2)
- CIRCexplorer2
- circrna_finder
- MapSplice
- find_circ
- DCC
- UROBORUS **still under testing

## Gotchas
Thank you Simone Coughlan for the catch :) 
```
*--params.star_index:* "/data/bdigby/grch38/index/star_index"

*--params.bwa_index:* "/data/bdigby/grch38/index/bwa"

*--params.hisat2_index:* "/data/bdigby/grch38/index/hisat2"

*--params.bowtie2_index:* "/data/bdigby/grch38/index/bowtie2/*" (need to point to files to collect for find_circ)

*--params.bowtie_index:* "/data/bdigby/grch38/index/bowtie/*" 
```

=======
# ![nf-core/circrna](docs/images/nf-core-circrna_logo.png)

**circRNA quantification, differential expression analysis and miRNA target prediction of RNA-Seq data**.

[![GitHub Actions CI Status](https://github.com/nf-core/circrna/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/circrna/actions)
[![GitHub Actions Linting Status](https://github.com/nf-core/circrna/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/circrna/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/circrna.svg)](https://hub.docker.com/r/nfcore/circrna)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23circrna-4A154B?logo=slack)](https://nfcore.slack.com/channels/circrna)

## Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`Podman`](https://podman.io/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```bash
    nextflow run nf-core/circrna -profile test,<docker/singularity/podman/conda/institute>
    ```

    > Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.

4. Start running your own analysis!

    <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

    ```bash
    nextflow run nf-core/circrna -profile <docker/singularity/podman/conda/institute> --input '*_R{1,2}.fastq.gz' --genome GRCh37
    ```

See [usage docs](https://nf-co.re/circrna/usage) for all of the available options when running the pipeline.

## Documentation

The nf-core/circrna pipeline comes with documentation about the pipeline: [usage](https://nf-co.re/circrna/usage) and [output](https://nf-co.re/circrna/output).

<!-- TODO nf-core: Add a brief overview of what the pipeline does and how it works -->

## Credits

nf-core/circrna was originally written by Barry Digby.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#circrna` channel](https://nfcore.slack.com/channels/circrna) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citation

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi. -->
<!-- If you use  nf-core/circrna for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
> ReadCube: [Full Access Link](https://rdcu.be/b1GjZ)
>>>>>>> TEMPLATE
