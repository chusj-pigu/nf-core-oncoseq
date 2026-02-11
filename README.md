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
4. **Tumor Classification** - [Marlin](https://github.com/hovestadt/MARLIN) or [Sturgeon](https://github.com/UMCUGenetics/sturgeon) *if tumor type is leukemia or cns and basecalled with 5mCG-5hmCG*
5. **Alignment QC** – [Cramino](https://github.com/wdecoster/cramino)
6. **CNV/cnLOH calling** – [QDNAseq](https://www.bioconductor.org/packages/QDNAseq), [SubChrom](https://github.com/Shaohua-Lei/SubChrom) and [Delly](https://github.com/dellytools/delly)
7. **Variant calling** – [ClairS-TO](https://github.com/HKU-BAL/ClairS-TO), [Clair3](https://github.com/HKU-BAL/Clair3)
8. **Structural variants** – [Sniffles2](https://github.com/fritzsedlazeck/Sniffles)
9. **VCF Annotation** – [SnpEff](https://pcingola.github.io/SnpEff/)
10. **Phasing** – [WhatsHap](https://whatshap.readthedocs.io/)
11. **Reporting** - [Quarto](https://quarto.org/) report summarizing coverage and variant information

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

Runs on basecalled 5mC/5hmC reads when `--realtime INT` is used.

Downsampling (`Ontime`) rules:

- `--adaptive` / `--wgs`: → downsample to **1h**
- `--cfdna`: → downsample to **8h** (if ≥10×)

Tools:

- **Leukemia** → *Marlin*
- **CNS** → *Sturgeon*
  - Uses the **T2T CHM13v2.0** reference

---

#### cf-DNA specific:

- Filters reads **>700 bp** using *chopper* (if `filter = yes`)
- CNV/SV calling: *IchorCNA*, *QDNAseq*, *Delly*, *Sniffles2*
- SNP calling if ≥10×: *ClairS-TO*, *Clair3*
- SubChrom if ≥15×

If you want to run IchorCNA as part of the cfDNA mode, please include `-profile ichor_hg38` or `-profile ichor_hg19` depending on your reference genome.

---

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate): -->

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

`samplesheet.csv`:

| sample <br> *(required)* | project <br> *(only with `--demux`)* | input <br> *(required)* | ref <br> *(required if not provided as parameter)* | ref_path <br> *(required if not provided as parameter)* | kit <br> *(only with `--demux`)* | bed <br> *(only with `--adaptive` and if not provided as parameter)* | low_fidelity <br> *(only with `--adaptive` and if not provided as parameter)* | padding <br> *(only with `--adaptive` and if not provided as parameter)* | purity <br> *(only with `--cfdna`)* | filter <br> *(only with `--cfdna`)* | tumor_type <br> |
|---------------------------|-------------------------|-----------------------|-----------------------------|----------------------------------|------------------------------------|--------------------------------------------------------------------------|--------------------------------------------|--------------------------------------------|---------------------------|-----------------------------|-------------------------------------------|
| sample1 | project ID for multiplexed input | path to <br> pod5, fastq or bam | hg38 | Path to reference genome, <br> index must exist in this path | Multiplexing kit used | Path to BED file containing <br> regions targeted by adaptive sampling | List of low-fidelity genes <br> in panel | Padding used around ROI <br> in bp | Estimated tumor purity <br> (between 0 and 1) | `yes` or `no` — whether to <br> filter out genomic DNA | leukemia, cns or other — <br> will run tumor classifier |



Each row represents the pod5/fastq/bam directory for one sample, and the reference to map it to.

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run chusj-pigu/nf-core-oncoseq \
   -profile narval \
   --input samplesheet.csv \
   --outdir results \
   --clin_database /path/to/clinvar.vcf.gz \
   --model sup \
   [--adaptive | --wgs | --cfdna] \
   [--realtime INT] \
   [--skip_basecalling] \
   [--skip_mapping]
```

By default, the pipeline will run in adaptive mode starting from basecalling, but the pipeline can also be run in WGS or cf-DNA mode using `--wgs` or `--cfdna` parameters respectively. Please see the pipeline output section to see which outputs are included with each mode. Please note that `--cfdna` mode is still in development.

If your input is already basecalled, use the parameter `--skip_basecalling` and provide the path to fastq files or folder containing fastq files in input. If your input is already mapped with minimap2, use the parameter `--skip_mapping` and provide the path to bam files or folder containing bam files as input.

The parameter `--clin_database` indicates the path to the ClinVar database that [SnpEff](https://pcingola.github.io/SnpEff/) will use to annotate the SNP vcf files. They can be downloaded here with https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/ for hg38 or with https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/ for hg19.

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
| `--basecall_model`     | `string` | _None_  | Dorado basecalling model to use (**required** if not using `--basecall_model_path`). |
| `--m_bases`            | `string` | _None_  | Basecalling modification model name.                                        |
| `--basecall_model_path`| `path`   | _None_  | Path to local copy of Dorado model (**required** if no network connection). |
| `--ref`                | `path`   | _None_  | Path to reference genome (must include index) (**required** if not in samplesheet)                              |
| `--ref_id`             | `path`   | hg38    | Reference ID (either hg19/GRCh37 or hg39/GRCh38)                            |
| `--demux`              | `bool`   | `false` | Enable demultiplexing after basecalling.                                    |
| `--skip_basecalling`   | `bool`   | `false` | Skip basecalling (use FASTQ files as input).                                |
| `--skip_mapping`       | `bool`   | `false` | Skip basecalling + mapping (use BAM files as input).                        |
| `--clin_database`      | `path`   | _None_  | Path to clinical database for annotating VCF files (must include index) (**required**).          |

---

### Modes

| Parameter        | Type   | Default | Description                                                   |
|------------------|--------|---------|---------------------------------------------------------------|
| `--adaptive`     | `bool` | `false` | Enable adaptive sampling mode.                                |
| `--wgs`          | `bool` | `false` | Whole Genome Sequencing (WGS) mode.                           |
| `--cfdna`        | `bool` | `false` | cf-DNA mode (in development).                                 |
| `--realtime`     | `int`  | _None_  | Enable realtime mode (only valid with `--adaptive`).          |

---

### Adaptive Mode Parameters

| Parameter               | Type   | Default | Description                                                                 |
|--------------------------|--------|---------|-----------------------------------------------------------------------------|
| `--bed`                  | `path` | _None_  | Path to BED file for adaptive sampling (**required** if not in samplesheet). |
| `--low_fidelity`         | `path` | _None_  | List of low-fidelity genes for adaptive sampling.                           |
| `--padding`              | `int`  | _None_  | Padding (in bp) around target regions (**required** if not in samplesheet). |

---

### Basecalling Options

| Parameter         | Type     | Default | Description                                                                 |
|-------------------|----------|---------|-----------------------------------------------------------------------------|
| `--minqs`         | `int`    | `10`    | Minimum Phred quality score for filtering reads.                            |
| `--device`        | `string` | _None_  | GPUs to use for basecalling (e.g. `cuda:0,1`).                              |
| `--batch`         | `int`    | _None_  | Batch size for basecalling.                                                 |
| `--mapping_small` | `bool`   | `true`  | Use reduced compute resources for mapping (recommended for adaptive/cfDNA). |

---

### Variant Calling Options

| Parameter          | Type     | Default                            | Description                                                                 |
|--------------------|----------|------------------------------------|-----------------------------------------------------------------------------|
| `--clairsto_model` | `string` | `ont_r10_dorado_sup_5khz_ssrs`     | ClairS-TO model to use.                                                     |
| `--sv_targets`     | `path`   | _None_                             | Regions of interest for SV visualization (CSV with `GENE,pos` columns).     |
| `--qdnaseq_binsize`| `int`    | `500`                              | Bin size (kb) for QDNAseq calling.                                          |
| `--subchrom_binsize`| `int`   | `500`                              | Bin size (kb) for Subchrom calling.                                         |
| `--delly_binsize`  | `int`   | `50000`                              | Bin size (bp) for Delly CNV calling.                                         |


## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/oncoseq/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/oncoseq/output).

| File Path             | Description | Condition        |
| --------------------- | ----------- | ---------------- |
| reads/{sample}_passed.fq.gz | Merged raw reads that have passed filter of average QS >= `--minqs` | If --skip_basecalling is not used |
| reads/{sample}_failed.fq.gz | Merged raw reads that have failed filter of average QS >= `--minqs` | If --skip_basecalling is not used |
| alignments/{sample}.sorted.bam<br>alignments/{sample}.sorted.bam.bai | Aligned and sorted bam file mapped to reference along with it's index | Always |
| variants/{sample}_snp_somatic_phased.vcf.gz<br>variants/{sample}_snp_somatic_clinvar_phased.vcf.gz | VCF files of phased SNV and indels called by [ClairS-TO](https://github.com/HKU-BAL/ClairS-TO)  | If `--adaptive` or `--wgs` mode is used |
| variants/{sample}_snp_germline_phased.vcf.gz<br>variants/{sample}_snp_germline_clinvar_phased.vcf.gz | VCF files of phased SNV and indels called by [Clair3](https://github.com/HKU-BAL/Clair3)  | If `--adaptive` or `--wgs` mode is used |
| variants/{sample}_sv.vcf.gz | VCF file of phased SV called by [Sniffles2](https://github.com/fritzsedlazeck/Sniffles) | If `--adaptive` or `--wgs` mode is used |
| variants/qdnaseq/{sample}_cnv_calls.vcf | VCF file of CNV called by [QDNAseq](https://www.bioconductor.org/packages/release/bioc/html/QDNAseq.html) | Always |
| variants/delly/{sample}_cnv_calls.vcf | VCF file of CNV called by [Delly](https://github.com/dellytools/delly) | Always |
| variants/subchrom/{sample}_cnv_calls.vcf | VCF file of CNV called by [SubChrom](https://github.com/Shaohua-Lei/SubChrom) | If `--adaptive` or `--wgs` mode is used: Always, if `--cfdna` is used, only with samples >= 15X |
| phasing/{sample}_haplotagged.bam<br>phasing/{sample}_haplotagged.bam.bai | Aligned bam and index file including phasing HP tags added by [WhatsHap](https://whatshap.readthedocs.io/en/latest/index.html) | If `--adaptive` or `--wgs` mode is used |
| phasing/{sample}.haploblocks.gtf | Gtf files containing phase blocks | If `--adaptive` or `--wgs` mode is used |
| reports/{sample}_report_output | Directory containing Quarto files and report | Always |
| reports/tumor_classifiers | Results of tumor classifiers | Depends on tumor type |
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
