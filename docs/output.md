## Pipeline output structure

<table>
<thead>
<tr>
<th>Type</th>
<th>File Path</th>
<th>Description</th>
<th>Condition</th>
</tr>
</thead>
<tbody>
<tr>
<td rowspan="2">Reads</td>
<td><code>reads/{sample}_pass.fq.gz</code></td>
<td>Merged raw reads passing QS filter ≥ <code>--minqs</code></td>
<td>Not generated with <code>--skip_basecalling</code></td>
</tr>
<tr>
<td><code>reads/{sample}_fail.fq.gz</code></td>
<td>Merged raw reads failing QS filter ≥ <code>--minqs</code></td>
<td>Not generated with <code>--skip_basecalling</code></td>
</tr>
<tr>
<td rowspan="2">Alignments</td>
<td>
<code>alignments/{sample}.bam</code><br>
<code>alignments/{sample}.bam.bai</code>
</td>
<td>Aligned and sorted BAM with index</td>
<td>If input is not aligned and sorted single bam</td>
</tr>
<tr>
</tr>
<tr>
<td rowspan="6">Variants</td>
<td>
<code>variants/{sample}_snp_phased.vcf.gz</code><br>
<td>Phased, VEP-annotated merged SNV/indels vcf file</td>
<td>Not generated with <code>--cfdna</code> or <code>--realtime</code></td>
</tr>
<tr>
<td>
<code>variants/SNPs/*vep.vcf.gz</code><br>
<code>variants/SNPs/*vep_filt.vcf.gz</code><br>
</td>
<td>VEP annotated and VEP filtered vcf files produced by SNP callers</a></td>
<td>For <code>--cfdna</code>, only if coverage ≥10×</td>
</tr>
<tr>
<td>
<code>variants/SVs/*_sv_vep.vcf.gz</code><br>
</td>
<td>VEP annotated and VEP filtered structural variants vcf files</td>
<td>Always</td>
</tr>
<tr>
<td>
<code>variants/CNVs/</code></td>
<td>Output files from CNV callers</td>
<td>Always</td>
</tr>
<tr>
<td><code>variants/pharmcat/</code></td>
<td>PharmCAT pharmacogenomic report</td>
<td>When SNPs are available</td>
</tr>
<tr>
</tr>
<tr>
<td rowspan="5">Report</td>
<td><code>reports/{sample}.html</code></td>
<td>Quarto report directory</td>
<td>Always</td>
</tr>
<tr>
<td><code>reports/tumor_classifiers/</code></td>
<td>Tumor classifier results (classifier JSON, figures and per-classifier <code>{sample}_{classifier}.tsv</code> tables)</td>
<td>If base modifications are present in input files or basecalling with modifications is turned on</td>
</tr>
<tr>
<td>
<code>reports/seqkit/{sample}_pass.tsv</code><br>
<code>reports/seqkit/{sample}_fail.tsv</code>
</td>
<td>Read QC stats (passed/failed)</td>
<td>Always</td>
</tr>
<tr>
<td><code>reports/mosdepth/</code></td>
<td>Coverage output from mosdepth (background and primary runs)</td>
<td>Always</td>
</tr>
<tr>
</tr>
<tr>
<td rowspan="2">Phasing</td>
<td><code>phasing/</code></td>
<td><code>Whatshap</code> happlotagged bam and haptoblocks</td>
<td>When <code>--cfdna</code> and <code>--realtime</code> is not used</td>
</tr>
<tr>
</tr>
<tr>
<td rowspan="2">Other</td>
<td><code>nasvar/</code></td>
<td>NASVAR pipeline output</td>
<td>If <code>--genome</code> is not hg19</td>
</tr>
<tr>
<td><code>program_outputs/{tool}/</code></td>
<td>Output of any process without its own <code>publishDir</code>, grouped by the first word of the process name (lowercased)</td>
<td>Always</td>
</tr>
</tbody>
</table>