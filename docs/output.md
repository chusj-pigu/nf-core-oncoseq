# nf-core/oncoseq: Output

## Introduction

This document describes the output produced by the nf-core/oncoseq pipeline. The pipeline generates comprehensive results for Oxford Nanopore Technologies (ONT) long-read sequencing analysis in oncology contexts, supporting three distinct workflow modes: adaptive sampling, cell-free DNA analysis, and whole genome sequencing.

Most plots and quality metrics are summarized in the MultiQC report, which provides an integrated view of all analysis results. The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data through several key analysis stages depending on the selected workflow mode:

**Common Processing Steps:**
- [Basecalling](#basecalling) - Converting raw ONT signals to sequences (optional)
- [Read Quality Control](#read-quality-control) - Quality assessment and filtering
- [Alignment](#alignment) - Mapping reads to reference genome
- [Variant Calling](#variant-calling) - SNV, indel, and structural variant detection
- [Variant Annotation](#variant-annotation) - Functional and clinical annotation
- [Coverage Analysis](#coverage-analysis) - Depth and uniformity assessment
- [MultiQC](#multiqc) - Comprehensive quality control report
- [Pipeline Information](#pipeline-information) - Execution metrics and provenance

**Mode-specific Steps:**
- [Adaptive Sampling Analysis](#adaptive-sampling-analysis) - On/off-target coverage separation
- [Time-series Analysis](#time-series-analysis) - Progressive sequencing evaluation
- [Copy Number Variation](#copy-number-variation) - CNV detection and visualization
- [Variant Phasing](#variant-phasing) - Haplotype-resolved variant analysis

## Basecalling

<details markdown="1">
<summary>Output files</summary>

- `basecalling/`
  - `*_unaligned.bam`: Raw basecalled sequences in unaligned BAM format
  - `*_unaligned_final.bam`: Final basecalled sequences (when resuming from existing uBAM)
  - `demux/`: Demultiplexed BAM files (multiplex mode only)
    - `barcode*/`: Individual barcode directories
    - `unclassified/`: Reads that couldn't be assigned to barcodes

</details>

Basecalling converts raw ONT signals (POD5/FAST5) into DNA sequences using [Dorado](https://github.com/nanoporetech/dorado). The pipeline supports both simplex and duplex basecalling modes, with optional modified base detection.

**Simplex Basecalling:** Single-pass basecalling for standard sequencing runs.
**Multiplex Basecalling:** Includes demultiplexing for barcoded samples using ONT kits.
**Modified Base Calling:** Optional detection of DNA modifications (5mC, 6mA, etc.).

## Read Quality Control

<details markdown="1">
<summary>Output files</summary>

- `reads/`
  - `*.fastq.gz`: Quality-filtered FASTQ files
  - `*_pass.fastq.gz`: Reads passing quality thresholds
  - `*_fail.fastq.gz`: Reads failing quality thresholds
- `reports/seqkit/`
  - `*.txt`: Read statistics (length, quality, N50, etc.)
- `reports/cramino/`
  - `*.txt`: Alignment and coverage statistics

</details>

Quality control involves filtering reads based on quality scores using [SAMtools](https://www.htslib.org/) and generating comprehensive statistics with [SeqKit](https://bioinf.shenwei.me/seqkit/) and [Cramino](https://github.com/wdecoster/cramino).

**Quality Filtering:** Removes reads below minimum quality threshold (default: Q10).
**Read Statistics:** Length distributions, quality scores, and sequence composition.
**Alignment QC:** Mapping statistics and coverage metrics.

## Alignment

<details markdown="1">
<summary>Output files</summary>

- `alignments/`
  - `*.bam`: Sorted and indexed BAM files
  - `*.bam.bai`: BAM index files
- `reports/cramino/`
  - `*_alignment_stats.txt`: Detailed alignment metrics

</details>

Read alignment to the reference genome using [minimap2](https://lh3.github.io/minimap2/) optimized for ONT long reads. Subsequent processing with [SAMtools](https://www.htslib.org/) for sorting and indexing.

**Alignment Features:**
- ONT-optimized alignment parameters (`-ax map-ont`)
- Automatic detection of primary/secondary alignments
- Comprehensive mapping quality assessment
- Efficient BAM compression and indexing

## Variant Calling

<details markdown="1">
<summary>Output files</summary>

- `variants/`
  - `clairs/`: Somatic variant calls from ClairS-TO
    - `*_somatic_snp.vcf.gz`: Somatic SNVs and indels
    - `*_somatic_snp_clinvar.vcf.gz`: ClinVar-annotated somatic variants
  - `clair3/`: Germline variant calls from Clair3
    - `*_germline_snp.vcf.gz`: Germline SNVs and indels
    - `*_germline_snp_clinvar.vcf.gz`: ClinVar-annotated germline variants
  - `structural/`: Structural variant calls
    - `*_sv.vcf.gz`: Structural variants from Sniffles2
  - `cnv/`: Copy number variations
    - `*_calls.vcf`: CNV calls in VCF format
    - `*_calls.bed`: CNV calls in BED format
    - `*_segs.seg`: Segmentation results

</details>

Comprehensive variant detection using multiple specialized callers:

**Small Variant Calling:**
- **[ClairS-TO](https://github.com/HKU-BAL/ClairS-TO)**: Somatic variant calling optimized for tumor samples
- **[Clair3](https://github.com/HKU-BAL/Clair3)**: Germline variant calling for constitutional variants

**Structural Variant Calling:**
- **[Sniffles2](https://github.com/fritzsedlazeck/Sniffles)**: Detection of large structural variants (deletions, insertions, inversions, duplications, translocations)

**Copy Number Variation:**
- **[QDNAseq](https://bioconductor.org/packages/QDNAseq/)**: Read-depth based CNV detection with visualization

## Variant Annotation

<details markdown="1">
<summary>Output files</summary>

- `variants/`
  - `*_annotated.vcf.gz`: Functionally annotated variants
  - `*_clinvar.vcf.gz`: Clinically annotated variants
- `reports/`
  - `annotation_summary.html`: Annotation statistics

</details>

Comprehensive variant annotation using multiple databases and tools:

**Functional Annotation:**
- **[SNPEff](https://pcingola.github.io/SnpEff/)**: Gene-based functional impact prediction
- **Effect Prediction**: Missense, nonsense, splice site impacts
- **Gene Annotation**: HGNC symbols, transcript IDs, protein changes

**Clinical Annotation:**
- **[SnpSift](https://pcingola.github.io/SnpEff/ss_introduction/)**: Clinical database integration
- **[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)**: Clinical significance and disease associations
- **Pathogenicity Assessment**: Benign, likely benign, VUS, likely pathogenic, pathogenic

## Variant Phasing

<details markdown="1">
<summary>Output files</summary>

- `phasing/`
  - `*_phased.vcf.gz`: Phased variant calls
  - `*_haplotagged.bam`: Haplotype-tagged BAM files
  - `*_haplotagged.bam.bai`: BAM index files
  - `*_phasing_stats.txt`: Phasing quality metrics

</details>

Variant phasing using [WhatsHap](https://whatshap.readthedocs.io/) to determine which variants are on the same chromosome copy:

**Phasing Features:**
- **Statistical Phasing**: Uses read linkage information for accurate phasing
- **Long-range Phasing**: Leverages ONT long reads for extended haplotype blocks
- **Haplotype Tagging**: Tags reads with haplotype information for downstream analysis
- **Quality Metrics**: Phasing statistics including block lengths and switch errors

## Copy Number Variation

<details markdown="1">
<summary>Output files</summary>

- `variants/qdnaseq/`
  - `*_calls.vcf`: CNV calls in VCF format
  - `*_calls.bed`: CNV calls in BED format
  - `*_segs.seg`: Segmentation results (IGV compatible)
  - `*_segs.bed`: Segmentation in BED format
  - `*_cov.png`: Coverage visualization plots
  - `*_noise_plot.png`: Noise assessment plots
  - `*_isobar_plot.png`: Read count distribution plots

</details>

Copy number variation analysis using QDNAseq with comprehensive visualization:

**CNV Analysis Features:**
- **Read-depth Based Detection**: Uses coverage patterns to identify CNVs
- **Binned Analysis**: Configurable bin sizes for resolution tuning
- **Statistical Modeling**: Noise correction and normalization
- **Multiple Output Formats**: VCF, BED, SEG for different downstream tools
- **Visualization**: Coverage plots, noise assessment, and quality metrics## Coverage Analysis

<details markdown="1">
<summary>Output files</summary>

- `reports/mosdepth/`
  - `*.mosdepth.global.dist.txt`: Global coverage distribution
  - `*.mosdepth.region.dist.txt`: Per-region coverage distribution  
  - `*.mosdepth.summary.txt`: Coverage summary statistics
  - `*.regions.bed.gz`: Per-region coverage depths
- `reports/cramino/`
  - `*_coverage_stats.txt`: Alignment and coverage metrics

</details>

Coverage analysis using [Mosdepth](https://github.com/brentp/mosdepth) and Cramino provides detailed depth and uniformity metrics:

**Coverage Metrics:**
- **Global Coverage**: Genome-wide depth distribution and statistics
- **Regional Coverage**: Coverage within specified genomic regions
- **Uniformity Assessment**: Coefficient of variation and coverage evenness
- **Quality Stratification**: Coverage by mapping quality and alignment flags

## Adaptive Sampling Analysis

<details markdown="1">
<summary>Output files</summary>

- `reports/`
  - `*_coverage_mapq.pdf`: Coverage visualization by mapping quality
  - `adaptive_sampling_summary.html`: Adaptive sampling efficiency report
- `alignments/`
  - `*_panel.bam`: On-target (panel) alignments
  - `*_bg.bam`: Off-target (background) alignments

</details>

Adaptive sampling analysis evaluates the efficiency of real-time target enrichment:

**On/Off-target Analysis:**
- **Panel Separation**: Splits alignments into on-target and off-target reads
- **Enrichment Metrics**: Calculates fold-enrichment and capture efficiency  
- **Coverage Distribution**: Analyzes depth across target regions
- **Quality Stratification**: Coverage by different mapping quality thresholds

**Visualization:**
- **Coverage Plots**: Target region coverage with background comparison
- **Enrichment Graphs**: Visual assessment of adaptive sampling performance
- **Quality Metrics**: Mapping quality distribution and filtering effects

## Time-series Analysis

<details markdown="1">
<summary>Output files</summary>

- `time_series/`
  - `*_t{timepoint}.bam`: BAM files at specific time points
  - `*_t{timepoint}_coverage.txt`: Coverage at each time point
  - `time_series_summary.html`: Progressive analysis report

</details>

Time-series analysis tracks sequencing progress over defined time points (adaptive mode only):

**Progressive Analysis:**
- **Time-point Splitting**: Creates BAM files for each specified time point
- **Coverage Tracking**: Monitors depth accumulation over time
- **Efficiency Assessment**: Evaluates adaptive sampling performance progression
- **Stopping Criteria**: Helps determine optimal sequencing duration

## MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: Comprehensive HTML report with interactive plots
  - `multiqc_data/`: Directory containing parsed statistics from all tools
  - `multiqc_plots/`: Static images from the report in various formats

</details>

[MultiQC](http://multiqc.info) generates a single HTML report summarizing all samples and analysis steps in your project. The report includes:

**Integrated Metrics:**
- **Basecalling Statistics**: Read counts, quality distributions, throughput
- **Alignment Metrics**: Mapping rates, coverage depth, quality scores  
- **Variant Calling Stats**: Variant counts, annotation summaries, quality metrics
- **Coverage Analysis**: Depth distributions, uniformity, target enrichment
- **Quality Control**: Pass/fail rates, filtering statistics, sample comparisons

**Interactive Features:**
- **Sample Filtering**: Filter plots by sample groups or quality metrics
- **Plot Customization**: Adjust scales, colors, and data representations
- **Data Export**: Download underlying data for custom analysis
- **Configuration**: Customizable report sections and plot parameters

Results from all pipeline tools are integrated including SeqKit, Cramino, Mosdepth, WhatsHap, and custom pipeline metrics.

## Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - `execution_report.html`: Nextflow execution report with resource usage
  - `execution_timeline.html`: Timeline of process execution
  - `execution_trace.txt`: Detailed trace of all process executions
  - `pipeline_dag.dot`/`pipeline_dag.svg`: Directed acyclic graph of pipeline processes
  - `pipeline_report.html`: Pipeline summary report (if `--email` used)
  - `software_versions.yml`: Versions of all software used in the pipeline
  - `samplesheet.valid.csv`: Validated and reformatted input samplesheet
  - `params.json`: All parameters used for the pipeline run

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides comprehensive execution reporting and monitoring:

**Execution Reports:**
- **Resource Usage**: CPU time, memory consumption, disk I/O per process
- **Timeline Visualization**: Process execution order and duration
- **Error Reporting**: Failed processes with error messages and context
- **Reproducibility**: Complete parameter sets and software versions

**Provenance Tracking:**
- **Software Versions**: All tools and their exact versions used
- **Parameter Documentation**: Complete record of pipeline configuration
- **Input Validation**: Processed samplesheet with validation results
- **Workflow Graph**: Visual representation of pipeline structure and dependencies

This information enables troubleshooting, performance optimization, and ensures full reproducibility of analysis results.

## Workflow-specific Outputs

### Adaptive Sampling Mode

When running with `--adaptive`, the pipeline generates additional outputs optimized for targeted sequencing analysis:

**Enhanced Coverage Analysis:**
- Separation of on-target vs off-target reads
- Target region enrichment metrics
- Coverage uniformity assessment across panel regions
- Background contamination analysis

**Time-series Support:**
- Progressive coverage accumulation tracking
- Adaptive sampling efficiency over time
- Optimal stopping point analysis

### Cell-free DNA Mode  

When running with `--cfdna`, the pipeline focuses on liquid biopsy analysis:

**CNV-focused Analysis:**
- Specialized copy number variation detection
- Fragment size distribution analysis  
- Low-input sample optimization
- cfDNA-specific quality metrics

### Whole Genome Sequencing Mode

When running with `--wgs`, the pipeline provides comprehensive genomic analysis:

**Complete Variant Analysis:**
- Both somatic and germline variant calling
- Structural variant detection across the genome
- Copy number variation analysis
- Comprehensive variant phasing

**Full Annotation:**
- Complete functional impact assessment
- Clinical significance annotation
- Population frequency data integration

## File Format Descriptions

### VCF Files (.vcf.gz)
Variant Call Format files containing detected variants with annotations:
- **Header**: Contains metadata, tool versions, and format descriptions  
- **Variants**: Genomic coordinates, reference/alternate alleles, quality scores
- **Annotations**: Functional impact, clinical significance, population frequencies
- **Genotype**: Sample-specific variant information including phase data

### BAM Files (.bam)
Binary Alignment Map files containing aligned sequencing reads:
- **Header**: Reference genome information and alignment parameters
- **Alignments**: Read sequences with mapping positions and quality scores
- **Tags**: Additional information like haplotype tags (HP), mapping quality (MQ)
- **Index**: Companion .bai files enable rapid region-specific access

### BED Files (.bed)
Browser Extensible Data format for genomic regions:
- **Coordinates**: Chromosome, start, end positions (0-based)
- **Annotations**: Region names, scores, strand information
- **Coverage**: Depth information for coverage analysis files

### Summary Statistics (.txt)
Tab-delimited text files with quantitative metrics:
- **Coverage**: Depth, uniformity, and quality metrics
- **Alignment**: Mapping rates, insert sizes, quality distributions  
- **Variants**: Counts, quality distributions, annotation summaries

## Data Interpretation Guidelines

### Quality Thresholds

**Read Quality:**
- Minimum Q10 for basic analysis
- Q20+ recommended for high-confidence variants
- Q30+ for clinical applications

**Coverage Depth:**
- 10x minimum for variant detection
- 30x+ recommended for reliable genotyping
- 50x+ for sensitive somatic variant detection

**Mapping Quality:**
- MAPQ ≥ 20 for reliable alignments
- MAPQ ≥ 60 for unique alignments
- Primary alignments preferred for variant calling

### Adaptive Sampling Metrics

**Enrichment Assessment:**
- On-target rate >50% indicates successful enrichment
- 5-10x fold enrichment typical for well-designed panels
- Monitor background contamination levels

**Time-series Analysis:**
- Compare coverage accumulation rates
- Identify plateau points for optimal stopping
- Assess real-time enrichment efficiency

### Clinical Interpretation

**Variant Classification:**
- Focus on ClinVar pathogenic/likely pathogenic variants
- Review variants of uncertain significance (VUS) carefully
- Consider population frequencies for interpretation

**CNV Analysis:**
- Validate large CNVs with orthogonal methods
- Consider copy number thresholds (>2.5 gains, <1.5 losses)
- Review segmentation quality in visualization plots
