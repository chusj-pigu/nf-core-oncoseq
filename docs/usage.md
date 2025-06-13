# nf-core/oncoseq: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/oncoseq/usage](https://nf-co.re/oncoseq/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

nf-core/oncoseq is a comprehensive bioinformatics pipeline for analyzing Oxford Nanopore Technologies (ONT) long-read sequencing data in oncology contexts. The pipeline supports three distinct workflow modes optimized for different sequencing strategies and sample types:

- **Adaptive Sampling Mode** (`--adaptive`): Optimized for targeted sequencing with adaptive sampling, including time-series analysis capabilities
- **Cell-free DNA Mode** (`--cfdna`): Specialized for circulating tumor DNA (ctDNA) analysis from liquid biopsies
- **Whole Genome Sequencing Mode** (`--wgs`): Complete genome analysis for comprehensive variant detection

The pipeline can process data starting from raw ONT signals (POD5/FAST5 files) through basecalling, or begin with pre-basecalled FASTQ files. It performs comprehensive variant analysis including SNV/indel calling, structural variant detection, copy number variation analysis, and variant phasing.

## Pipeline Modes

### Adaptive Sampling Mode (`--adaptive`)

This mode is designed for targeted sequencing experiments using ONT's adaptive sampling technology. It includes:

- **Targeted Coverage Analysis**: Evaluates sequencing coverage across specified genomic regions
- **Time-series Evaluation**: Tracks sequencing progress over time for adaptive sampling optimization
- **Region-specific Variant Calling**: Focuses analysis on regions of interest defined by BED files
- **Coverage Separation**: Separates on-target and off-target reads for quality assessment

**Key Parameters:**
- `--bed`: BED file defining regions of interest for adaptive sampling
- `--adaptive_samplesheet`: Optional samplesheet specifying BED files per sample
- `--time_series`: Enable time-series analysis for tracking sequencing progress
- `--time_points`: Comma-separated list of time points for analysis
- `--padding`: Padding around regions of interest (default: based on BED file)
- `--low_fidelity`: File listing genes with known sequencing challenges

### Cell-free DNA Mode (`--cfdna`)

Optimized for liquid biopsy analysis focusing on:

- **Copy Number Variation Detection**: Specialized CNV calling optimized for ctDNA
- **Low-input Sample Processing**: Handles the unique challenges of cfDNA analysis
- **Multiplexed Sample Support**: Designed for high-throughput cfDNA screening

### Whole Genome Sequencing Mode (`--wgs`)

Comprehensive genome analysis including:

- **Complete Variant Calling**: SNVs, indels, structural variants, and CNVs
- **Variant Phasing**: Links variants to their parental chromosomes
- **Structural Variant Detection**: Large-scale genomic rearrangements
- **Comprehensive Annotation**: Clinical and functional variant annotation

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

The pipeline supports different input types depending on your starting data:

### FASTQ Input (Skip Basecalling)

For pre-basecalled data, use the `--skip_basecalling` flag:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
SAMPLE1,/path/to/sample1.fastq.gz,
SAMPLE2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz
```

### Raw Signal Input (With Basecalling)

For raw ONT data requiring basecalling:

```csv title="samplesheet.csv"
sample,pod5,ubam
SAMPLE1,/path/to/sample1_pod5_dir,/path/to/sample1.ubam
SAMPLE2,/path/to/sample2_pod5_dir,/path/to/sample2.ubam
```

**Additional Input Files:**

- `--demux_samplesheet`: For multiplexed samples requiring demultiplexing
- `--adaptive_samplesheet`: For adaptive sampling with sample-specific BED files
- `--ubam_samplesheet`: For unaligned BAM files from previous basecalling runs

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz
```

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 3 columns to match those defined in the table below.

A final samplesheet file consisting of both single- and paired-end data may look something like the one below. This is for 6 samples, where `TREATMENT_REP3` has been sequenced twice.

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L004_R1_001.fastq.gz,
```

| Column    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1` | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2` | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Running the pipeline

The typical command for running the pipeline varies depending on the workflow mode:

### Adaptive Sampling Mode

```bash
nextflow run nf-core/oncoseq \
    --input ./samplesheet.csv \
    --outdir ./results \
    --adaptive \
    --bed ./regions_of_interest.bed \
    --time_series \
    --time_points "1,2,4,8,12,24" \
    -profile docker
```

### Cell-free DNA Mode

```bash
nextflow run nf-core/oncoseq \
    --input ./samplesheet.csv \
    --outdir ./results \
    --cfdna \
    --demux_samplesheet ./demux.csv \
    -profile docker
```

### Whole Genome Sequencing Mode

```bash
nextflow run nf-core/oncoseq \
    --input ./samplesheet.csv \
    --outdir ./results \
    --wgs \
    --clairsto_model ont_r10_dorado_sup_5khz_ssrs \
    --basecall_model sup \
    -profile docker
```

### Skip Basecalling (Pre-basecalled FASTQ)

```bash
nextflow run nf-core/oncoseq \
    --input ./samplesheet.csv \
    --outdir ./results \
    --adaptive \
    --skip_basecalling \
    --bed ./regions_of_interest.bed \
    -profile docker
```

## Key Parameters

### Basecalling Parameters

- `--basecall_model`: Dorado basecalling model (`fast`, `hac`, `sup`)
- `--basecall_model_path`: Path to custom basecalling model
- `--device`: Computing device for basecalling (`cpu`, `cuda:0`)
- `--batch`: Batch size for basecalling processes
- `--m_bases`: Modified base detection model
- `--minqs`: Minimum quality score for read filtering (default: 10)

### Variant Calling Parameters

- `--clairsto_model`: ClairS-TO model for variant calling
- `--clin_database`: Clinical variant database for annotation (ClinVar)
- `--qdnaseq_binsize`: Bin size for CNV calling (default: 500)

### Adaptive Sampling Parameters

- `--bed`: BED file defining target regions
- `--adaptive_samplesheet`: Sample-specific BED file assignments
- `--padding`: Padding around target regions
- `--low_fidelity`: File listing problematic genomic regions
- `--time_series`: Enable time-series analysis
- `--time_points`: Time points for analysis (comma-separated)
- `--include_full`: Include full sequencing time in analysis

### Resource Parameters

- `--max_memory`: Maximum memory per process (default: 16G)
- `--max_cpus`: Maximum CPUs per process (default: 4)
- `--max_time`: Maximum time per process (default: 4h)

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

## Output Structure

The pipeline generates organized output directories based on the analysis performed:

```bash
results/
├── alignments/          # BAM files and alignment statistics
├── basecalling/         # Basecalled FASTQ files (if basecalling performed)
├── variants/            # Variant calling results (VCF files)
│   ├── clairs/         # ClairS-TO somatic variant calls
│   ├── clair3/         # Clair3 germline variant calls
│   ├── phasing/        # Phased variant files
│   └── structural/     # Structural variant calls
├── cnv/                # Copy number variation analysis
├── coverage/           # Coverage analysis and statistics
├── time_series/        # Time-series analysis results (adaptive mode)
├── reports/            # MultiQC and other analysis reports
└── pipeline_info/      # Pipeline execution information
```

## Test Profiles

The pipeline includes several test profiles for validation:

- `test`: Minimal test with small dataset
- `test_fastq`: Test with pre-basecalled FASTQ input
- `test_full`: Comprehensive test with realistic data
- `test_full_ts`: Full test with time-series analysis
- `test_drac`: Test configuration for DRAC cluster

Example test run:

```bash
nextflow run nf-core/oncoseq -profile test,docker --outdir test_results
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/oncoseq -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
adaptive: true
bed: './target_regions.bed'
time_series: true
time_points: '1,2,4,8,12,24'
clairsto_model: 'ont_r10_dorado_sup_5khz_ssrs'
basecall_model: 'sup'
max_memory: '32.GB'
max_cpus: 8
```

## Analysis Capabilities

### Variant Detection and Analysis

The pipeline provides comprehensive variant analysis capabilities:

**Small Variant Calling:**
- **ClairS-TO**: Somatic variant calling optimized for tumor samples
- **Clair3**: Germline variant calling for constitutional variants
- **Quality Filtering**: Configurable quality thresholds and filters
- **Clinical Annotation**: Integration with ClinVar and other databases

**Structural Variant Detection:**
- **Sniffles2**: Long-read optimized SV calling
- **Large Deletions, Insertions**: Detection of major genomic rearrangements
- **Copy Number Variations**: Both targeted and genome-wide CNV analysis
- **Breakpoint Resolution**: Single-nucleotide resolution for many SVs

**Variant Phasing:**
- **WhatsHap**: Statistical phasing of heterozygous variants
- **Haplotype-resolved Analysis**: Separates maternal and paternal variants
- **Long-range Phasing**: Leverages long-read advantages for extended phasing

### Coverage and Quality Analysis

**Adaptive Sampling Specific:**
- **On/Off-target Analysis**: Quantifies enrichment efficiency
- **Time-series Tracking**: Monitors coverage accumulation over time
- **Region-specific Metrics**: Per-gene and per-exon coverage statistics
- **Adaptive Efficiency**: Measures real-time selection performance

**General Coverage Metrics:**
- **Depth Distribution**: Coverage histograms and statistics
- **Uniformity Assessment**: Coefficient of variation across regions
- **Quality Metrics**: Read quality, mapping quality, and alignment statistics

### Specialized Analysis Modes

**Cell-free DNA Analysis:**
- **Low-input Optimization**: Adapted algorithms for cfDNA characteristics
- **Fragment Analysis**: Size distribution and fragmentation patterns
- **CNV Detection**: Optimized for detecting circulating tumor DNA
- **Quality Control**: cfDNA-specific QC metrics and thresholds

**Time-series Analysis:**
- **Progressive Coverage**: Tracks sequencing progress over defined time points
- **Adaptive Efficiency**: Measures on-target enrichment over time
- **Stopping Criteria**: Supports analysis of optimal sequencing duration
- **Comparative Analysis**: Multiple time point comparison and visualization

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/oncoseq
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/oncoseq releases page](https://github.com/nf-core/oncoseq/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Workflow Mode Selection

The pipeline supports exactly one workflow mode at a time. You must specify one of:

- `--adaptive`: Adaptive sampling analysis
- `--cfdna`: Cell-free DNA analysis  
- `--wgs`: Whole genome sequencing analysis

**Note:** These modes are mutually exclusive and cannot be combined.

## Advanced Usage

### Time-series Analysis (Adaptive Mode Only)

For adaptive sampling experiments, you can track sequencing progress over time:

```bash
nextflow run nf-core/oncoseq \
    --input ./samplesheet.csv \
    --adaptive \
    --time_series \
    --time_points "0.5,1,2,4,8,12,24" \
    --include_full \
    --outdir ./results \
    -profile docker
```

### Custom Basecalling Models

For specialized basecalling requirements:

```bash
nextflow run nf-core/oncoseq \
    --input ./samplesheet.csv \
    --basecall_model_path /path/to/custom/model \
    --m_bases 5mCG_5hmCG \
    --device cuda:0 \
    --outdir ./results \
    -profile docker
```

### Multi-sample Demultiplexing

For multiplexed runs requiring demultiplexing:

```bash
nextflow run nf-core/oncoseq \
    --input ./samplesheet.csv \
    --demux_samplesheet ./demux_barcodes.csv \
    --demux SQK-RBK004 \
    --outdir ./results \
    -profile docker
```

## Troubleshooting

### Common Issues

1. **Memory Errors**: Increase `--max_memory` for large datasets
2. **Time Limits**: Adjust `--max_time` for complex analyses  
3. **GPU Issues**: Ensure CUDA drivers are installed for GPU basecalling
4. **File Permissions**: Ensure input files are readable and output directory is writable

### Resource Optimization

For large datasets or resource-constrained environments:

```bash
# Increase resources
--max_memory 64GB --max_cpus 16 --max_time 24h

# Reduce parallel processes
process.maxForks = 2

# Use specific resource profiles
-profile drac  # For DRAC cluster
-profile mpgi_local  # For local MPGI systems
```

### Debugging

Enable detailed logging for troubleshooting:

```bash
nextflow run nf-core/oncoseq \
    --input ./samplesheet.csv \
    --outdir ./results \
    -profile docker \
    -with-trace \
    -with-report \
    -with-timeline \
    -with-dag flowchart.html
```
