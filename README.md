<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-oncoseq_logo_dark.png">
    <img alt="nf-core/oncoseq" src="docs/images/nf-core-oncoseq_logo_light.png">
  </picture>
</h1>

![Status: Pre-Release](https://img.shields.io/badge/status-pre--release-important)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/oncoseq/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/oncoseq)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23oncoseq-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/oncoseq)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/oncoseq** is a bioinformatics pipeline that calls variants of interests for Oncology research from raw data to vcf files and digestible reports.

<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

![nf-core-oncoseq summary of full workflow](assets/nf-core-oncoseq_schema.jpg)

---

## Pipeline Overview

The general steps are as followed:

1. **Basecalling** – [Dorado](https://github.com/nanoporetech/dorado)
   *(Optional: skip using `--skip_basecalling` if FASTQ input is provided)*
   *(Optional: Demultiplexing is done if kit is included in samplesheet and demux_samplesheet is provided*
2. **Read QC** – [Seqkit](https://bioinf.shenwei.me/seqkit/)
3. **Alignment** – [minimap2](https://lh3.github.io/minimap2/minimap2.html)
   *(Optional: skip using `--skip_mapping` if bam input is provided)*
4. **Tumor Classification** - [Marlin](https://github.com/hovestadt/MARLIN), [Sturgeon](https://github.com/UMCUGenetics/sturgeon), [CrossNN](https://gitlab.com/euskirchen-lab/crossNN), [Tucan](https://github.com/UMCUGenetics/tucan) implemented in Rust, with liftovers for reference used.
5. **Alignment QC** – [Cramino](https://github.com/wdecoster/cramino)
6. **CNV/cnLOH calling** – [QDNAseq](https://www.bioconductor.org/packages/QDNAseq), [SubChrom](https://github.com/Shaohua-Lei/SubChrom) and [Delly](https://github.com/dellytools/delly)
7. **Variant calling** – [ClairS-TO](https://github.com/HKU-BAL/ClairS-TO), [Clair3](https://github.com/HKU-BAL/Clair3)
8. **Structural variants** – [Sniffles2](https://github.com/fritzsedlazeck/Sniffles)
9. **VCF Annotation** – [SnpEff](https://pcingola.github.io/SnpEff/) and [Ensemblvep](https://grch37.ensembl.org/info/docs/tools/vep/index.html)
10. **Phasing** – [WhatsHap](https://whatshap.readthedocs.io/) (only for `--adaptive` and `--wgs` modes)
11. **Reporting** - [Quarto](https://quarto.org/) report summarizing coverage, variant information and methylation based tumor classification

---

#### Adaptive specific:

When `--adaptive` is used, additional steps are included to show region specific coverage

1. **Coverage calculation** – [mosdepth](https://github.com/brentp/mosdepth)
2. **Visualization and reporting** – [R](https://www.r-project.org/)

---

#### Time Series Breakdown Analysis

For prospectively run the workflow on different timepoints, [Ontime](https://github.com/mbhall88/ontime) is run after the **Alignment** step.

---

#### Tumor classification

Runs on basecalled 5mC/5hmC reads.

Downsampling (`Ontime`) rules:

- `--adaptive` / `--wgs`: → downsample to **1h**
- `--cfdna`: → downsample to **8h** (if ≥10×)

Tools:

All classifiers are used, but their relevency is the following:

- **Leukemia** → *Marlin*
- **CNS** → *Sturgeon* and *CrossNN Caper*
- **Pan cancer** → *CrossNN Pancancer* and *Tucan*

---

#### cf-DNA specific:

- Filters reads **>700 bp** using *chopper* (if `filter = yes`)
- CNV/SV calling: *IchorCNA*, *QDNAseq*, *Delly*, *Sniffles2*
- SNP calling if ≥10×: *ClairS-TO*, *Clair3*
- SubChrom if ≥15×

If you want to run IchorCNA as part of the cfDNA mode, please include `-profile ichor_hg38` or `-profile ichor_hg19` depending on your reference genome.

![nf-core-oncoseq summary of workflow with `--cfdna`](assets/nf-core-oncoseq_cfdna.jpg)

---

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate): -->

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

| sample <br> *(required)* | project <br> *(only with `--demux`)* | input <br> *(required)* | ref <br> *(required if not provided as parameter)* | ref_path <br> *(required if not provided as parameter)* | kit <br> *(only with `--demux`)* | bed <br> *(only with `--adaptive` and if not provided as parameter)* | low_fidelity <br> *(only with `--adaptive` and if not provided as parameter)* | padding <br> *(only with `--adaptive` and if not provided as parameter)* | purity <br> *(only with `--cfdna`)* | filter <br> *(only with `--cfdna`)* |
|---------------------------|-------------------------|-----------------------|-----------------------------|----------------------------------|------------------------------------|--------------------------------------------------------------------------|--------------------------------------------|--------------------------------------------|---------------------------|-----------------------------|-------------------------------------------|
| sample1 | project ID for multiplexed input | path to <br> pod5, fastq or bam | hg38 | Path to reference genome, <br> index must exist in this path | Multiplexing kit used | Path to BED file containing <br> regions targeted by adaptive sampling | List of low-fidelity genes <br> in panel | Padding used around ROI <br> in bp | Estimated tumor purity <br> (between 0 and 1) | `yes` or `no` — whether to <br> filter out genomic DNA



Each row represents the pod5/fastq/bam directory for one sample, and the reference to map it to.

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run nf-core/oncoseq \
   -profile narval \
   --input samplesheet.csv \
   --outdir results \
   --reference_cache_dir /path/to/reference/cache \
   [--adaptive | --wgs | --cfdna] \
   [--realtime INT] \
   [--skip_basecalling] \
   [--skip_mapping]
```

By default, the pipeline will run in adaptive mode starting from basecalling, but the pipeline can also be run in WGS or cf-DNA mode using `--wgs` or `--cfdna` parameters respectively. Please see the pipeline output section to see which outputs are included with each mode. Please note that `--cfdna` mode is still in development.

If your input is already basecalled, use the parameter `--skip_basecalling` and provide the path to fastq files or folder containing fastq files in input. If your input is already mapped with minimap2, use the parameter `--skip_mapping` and provide the path to bam files or folder containing bam files as input.

The parameter `--vep_cache` indicates the path to the Ensembl vep database.

To run in real time while data is still sequencing, use the `--realtime [INT]` where you must provide the time of sequencing as an integer. **This is only available for the adaptive mode**

![nf-core-oncoseq when run in real time](assets/nf-core-oncoseq_schema_realtime.jpg)

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/oncoseq/usage) and the [parameter documentation](https://nf-co.re/oncoseq/parameters).

## Pipeline input and parameters

### General and/or Required Parameters

| Parameter              | Type     | Default | Description                                                                 |
|------------------------|----------|---------|-----------------------------------------------------------------------------|
| `--input`              | `path`   | _None_  | Path to input samplesheet (**required**).                                   |
| `--outdir`             | `path`   | _None_  | Directory where outputs will be published (**required**).                   |
| `--reference_cache_dir`| `path`   | workspace | Directory for staging reference assets and Dorado models.                   |
| `--basecall_model`     | `string` | _None_  | Dorado basecalling model name (e.g., `sup`, `hac`, `fast`). Supports version pinning with `@v`. |
| `--m_bases`            | `string` | _None_  | Basecalling modification model (e.g., `5mC_5hmC`, `5mC`).                   |
| `--ref`                | `path`   | _None_  | Path to reference genome (must include `.fai` index). Required if not in samplesheet. |
| `--ref_id`             | `string` | hg38    | Reference genome ID: `hg38` (GRCh38) or `hg19` (GRCh37).                    |
| `--skip_basecalling`   | `bool`   | `false` | Skip basecalling; input is FASTQ files.                                     |
| `--skip_mapping`       | `bool`   | `false` | Skip basecalling and mapping; input is BAM files.                           |
| `--demux`              | `bool`   | `false` | Enable demultiplexing via `--demux_samplesheet` (requires kit column).       |
| `--demux_samplesheet`  | `path`   | _None_  | Samplesheet with demux info: barcode, sample name (requires `--demux`).     |

---

### Pipeline Modes

Choose exactly one mode:

| Parameter        | Type   | Default | Description                                                   |
|------------------|--------|---------|---------------------------------------------------------------|
| `--adaptive`     | `bool` | `false` | Adaptive sampling mode (requires `--bed`, `--padding`).        |
| `--wgs`          | `bool` | `false` | Whole Genome Sequencing mode.                                 |
| `--cfdna`        | `bool` | `false` | cf-DNA sequencing mode (in development).                      |
| `--realtime`     | `int`  | _None_  | Real-time mode for ongoing sequencing (hours). Adaptive only. |
| `--time_series`  | `bool` | `false` | Enable time series analysis with `--time_points`.             |

---

### Adaptive Sampling Parameters

| Parameter               | Type   | Default | Description                                                                 |
|--------------------------|--------|---------|-----------------------------------------------------------------------------|
| `--bed`                  | `path` | _None_  | BED file with target regions (**required** for `--adaptive`).               |
| `--adaptive_samplesheet` | `path` | _None_  | Samplesheet specifying per-sample BED files (overrides `--bed`).            |
| `--padding`              | `int`  | `20000` | Padding (bp) around target regions (**required** for `--adaptive`).         |
| `--low_fidelity`         | `path` | default | Path to text file listing low-fidelity genes to exclude from coverage.     |

---

### Basecalling and Mapping Options

| Parameter         | Type     | Default | Description                                                                 |
|-------------------|----------|---------|-----------------------------------------------------------------------------|
| `--minqs`         | `int`    | `10`    | Minimum Phred quality score for read filtering.                             |
| `--device`        | `string` | _None_  | GPU devices for Dorado (e.g., `cuda:0,1`). Auto-detects if not specified.  |
| `--batch`         | `int`    | _None_  | Batch size for basecalling (auto-tuned if not specified).                  |
| `--mapping_small` | `bool`   | `true`  | Use reduced resources for mapping (recommended for adaptive/cfDNA).        |

---

### Variant Calling Options

| Parameter          | Type     | Default                            | Description                                                                 |
|--------------------|----------|------------------------------------|-----------------------------------------------------------------------------|
| `--clairsto_model` | `string` | `ont_r10_dorado_sup_5khz_ssrs`     | ClairS-TO model for somatic SNV/indel calling.                              |
| `--qdnaseq_binsize`| `int`    | `500`                              | Bin size (kb) for QDNAseq CNV calling.                                      |
| `--subchrom_binsize`| `int`   | `500`                              | Bin size (kb) for SubChrom CNV calling.                                     |
| `--delly_bin_size` | `int`    | `50000`                            | Bin size (bp) for Delly CNV calling.                                        |
| `--sv_targets`     | `path`   | `NOCSV`                            | CSV file for SV visualization: columns `GENE`, `pos` (Figeno output).       |
| `--sv_exclude`     | `path`   | default                            | Path to text file listing genes to exclude from SV analysis.                |

---

### VEP Annotation Options

| Parameter              | Type     | Default | Description                                                                 |
|------------------------|----------|---------|-----------------------------------------------------------------------------|
| `--vep_cache`          | `path`   | _None_  | Path to VEP cache directory (**required** for annotation).                  |
| `--vep_version`        | `int`    | `115`   | Ensembl VEP version matching your cache.                                    |
| `--filtervep_expression` | `string` | (see schema) | FilterVEP expression for SNV filtering (impacts, frequencies, etc.). |

---

### cf-DNA Mode Parameters

| Parameter               | Type     | Default | Description                                                                 |
|------------------------|----------|---------|-----------------------------------------------------------------------------|
| `--max_length`         | `int`    | `700`   | Maximum read length (bp) for cf-DNA filtering (removes longer reads).       |
| `--ichor_bin_size`     | `int`    | `500000` | Bin size (bp) for IchorCNA CNV calling.                                     |
| `--min_mapq_ichor`     | `int`    | `20`    | Minimum MAPQ threshold for IchorCNA.                                        |
| `--rca`                | `bool`   | `false` | Enable rolling circle amplification mode (runs TideHunter).                 |

---

### Time Series Analysis Parameters

| Parameter             | Type     | Default | Description                                                                 |
|------------------------|----------|---------|-----------------------------------------------------------------------------|
| `--time_series`       | `bool`   | `false` | Enable time series analysis (downsamples via Ontime).                      |
| `--time_points`       | `string` | `1,3,6,12,18,24,48,60,72` | Comma-separated downsampling timepoints (hours).  |
| `--include_full`      | `bool`   | `true`  | Include full-depth analysis alongside downsampled timepoints.              |

---

### Reference Cache and Resource Parameters

| Parameter             | Type     | Default | Description                                                                 |
|------------------------|----------|---------|-----------------------------------------------------------------------------|
| `--ubam`              | `path`   | `NOFILE` | Path to incomplete uBAM for basecalling resume.                             |
| `--ubam_samplesheet`  | `path`   | _None_  | Samplesheet pointing to uBAM files to restart basecalling from.             |
| `--max_memory`        | `string` | `16G`   | Maximum memory any process can consume.                                     |
| `--max_cpus`          | `int`    | `4`     | Maximum CPU cores any process can use.                                      |
| `--max_time`          | `string` | `4h`    | Maximum execution time any process can take.                                |
| `--huge`              | `bool`   | `false` | Scale resource allocations for very large datasets.                         |


## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/oncoseq/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/oncoseq/output).

| File Path             | Description | Condition        |
| --------------------- | ----------- | ---------------- |
| reads/{sample}_passed.fq.gz | Merged raw reads that have passed filter of average QS >= `--minqs` | If --skip_basecalling is not used |
| reads/{sample}_failed.fq.gz | Merged raw reads that have failed filter of average QS >= `--minqs` | If --skip_basecalling is not used |
| alignments/{sample}.sorted.bam<br>alignments/{sample}.sorted.bam.bai | Aligned and sorted bam file mapped to reference along with it's index | Always |
| variants/{sample}_snp_somatic_phased.vcf.gz<br>variants/{sample}_snp_somatic_vep.vcf.gz | VCF files of phased SNV and indels called by [ClairS-TO](https://github.com/HKU-BAL/ClairS-TO)  | If `--adaptive` or `--wgs` mode is used, or if `--cfdna` coverage is > 10X |
| variants/{sample}_snp_germline_phased.vcf.gz<br>variants/{sample}_snp_germline_vep.vcf.gz | VCF files of phased SNV and indels called by [Clair3](https://github.com/HKU-BAL/Clair3)  | If `--adaptive` or `--wgs` mode is used, or if `--cfdna` coverage is > 10X |
| variants/{sample}_sv.vcf.gz | VCF file of phased SV called by [Sniffles2](https://github.com/fritzsedlazeck/Sniffles) | Always |
| variants/qdnaseq/{sample}_cnv_calls.vcf | VCF file of CNV called by [QDNAseq](https://www.bioconductor.org/packages/release/bioc/html/QDNAseq.html) | Always |
| variants/delly/{sample}_cnv_calls.vcf | VCF file of CNV called by [Delly](https://github.com/dellytools/delly) | Always |
| variants/subchrom/{sample}_cnv_calls.vcf | VCF file of CNV called by [SubChrom](https://github.com/Shaohua-Lei/SubChrom) | If `--adaptive` or `--wgs` mode is used: Always, if `--cfdna` is used, only with samples >= 15X |
| phasing/{sample}_haplotagged.bam<br>phasing/{sample}_haplotagged.bam.bai | Aligned bam and index file including phasing HP tags added by [WhatsHap](https://whatshap.readthedocs.io/en/latest/index.html) | If `--adaptive` or `--wgs` mode is used |
| phasing/{sample}.haploblocks.gtf | Gtf files containing phase blocks | If `--adaptive` or `--wgs` mode is used |
| reports/{sample}_report_output | Directory containing Quarto files and report | Always |
| reports/tumor_classifiers | Results of tumor classifiers | Always |
| reports/seqkit/{sample}_pass.tsv<br>reports/{sample}_fail.tsv | Passed and failed fastq files stats | Always |
| reports/cramino/{sample}_cramino_stats.txt | Alignment summary stats | Always |
| reports/{sample}_coverage_mapq.pdf | Plot showing ROIs coverage | If `--adaptive` mode is used |


## Credits

nf-core/oncoseq was originally written by CHUSJ-MPGI.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#oncoseq` channel](https://nfcore.slack.com/channels/oncoseq) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/oncoseq for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
