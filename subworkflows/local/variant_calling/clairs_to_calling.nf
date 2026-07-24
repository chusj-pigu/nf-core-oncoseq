/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This section imports all necessary modules, subworkflows, and functions required
    for variant calling with ClairS-TO, annotation with SNPEff/SNPSift, and VCF
    processing with BCFtools.
*/
include { CLAIRS_TO_CALL               } from '../../../modules/local/clairsto/main.nf'
include { SNPEFF_ANNOTATE              } from '../../../modules/local/snpeff/main.nf'
include { SNPSIFT_ANNOTATE             } from '../../../modules/local/snpeff/main.nf'
include { BCFTOOLS_VIEW as BCFTOOLS_COUNT         } from '../../../modules/local/bcftools/main.nf'
include { BCFTOOLS_VIEW as BGZIP_VCF_FINAL        } from '../../../modules/local/bcftools/main.nf'
include { BCFTOOLS_VIEW as BGZIP_VCF_INTER        } from '../../../modules/local/bcftools/main.nf'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_HG38       } from '../../../modules/nf-core/ensemblvep/vep/main.nf'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_HG19       } from '../../../modules/nf-core/ensemblvep/vep/main.nf'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_HS1        } from '../../../modules/nf-core/ensemblvep/vep/main.nf'
include { ENSEMBLVEP_FILTERVEP            } from '../../../modules/nf-core/ensemblvep/filtervep/main.nf'
include { paramsSummaryMap             } from 'plugin/nf-schema'
include { paramsSummaryMultiqc         } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML       } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText       } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'
include { modifyMetaId                } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'
include { BCFTOOLS_CONCAT               } from '../../../modules/local/bcftools/main.nf'
include { BCFTOOLS_SORT                 } from '../../../modules/local/bcftools/main.nf'
include { BCFTOOLS_FILTER_REGION                 } from '../../../modules/local/bcftools/main.nf'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_SNV   } from '../../../modules/local/bcftools/main.nf'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_INDEL } from '../../../modules/local/bcftools/main.nf'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_RAW   } from '../../../modules/local/bcftools/main.nf'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_FINAL } from '../../../modules/local/bcftools/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Define named functions to replace anonymous closures for better readability
    and maintainability.
*/

// Extract reference FASTA and index from reference channel
def extractRefFiles(meta, _ref, ref_fasta, ref_fai) {
    return tuple(meta, ref_fasta, ref_fai)
}

// Extract reference genome type for database selection
def extractRefType(meta, refid, _ref_fasta, _ref_fai) {
    return tuple(meta, refid)
}

// Check if reference is hg38/GRCh38
def isHg38Reference(_meta, refid) {
    return refid.matches('hg38|GRCh38')
}

// Check if reference is hg19/GRCh37
def isHg19Reference(_meta, refid) {
    return refid.matches('hg19|GRCh37')
}

// Map hg38 reference to SNPEff database
def mapHg38Database(meta, _refid) {
    return tuple(meta, 'GRCh38.p14')
}

// Map hg19 reference to SNPEff database
def mapHg19Database(meta, _refid) {
    return tuple(meta, 'GRCh37.p13')
}

// Prepare VCF files for SNPSift annotation with clinical database
def prepareForSnpsift(meta, vcf, clin_db, clin_db_idx) {
    return tuple(meta, vcf, clin_db, clin_db_idx)
}

// Wrapper functions using the generalized modifyMetaId function
def restoreSnvMeta(meta, vcf, tbx) {
    def new_meta = modifyMetaId(meta, 'remove_suffix', '', '', '_snv')
    return tuple(new_meta, vcf, tbx)
}

def restoreIndelMeta(meta, vcf, tbx) {
    def new_meta = modifyMetaId(meta, 'remove_suffix', '', '', '_indel')
    return tuple(new_meta, vcf, tbx)
}

def prepareSomaticSnpAnnotation(meta, output, database) {
    def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_somatic_snp_snpeff')
    return tuple(new_meta, output, database)
}

def addSnvSuffix(meta, vcf) {
    def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_snv')
    return tuple(new_meta, vcf)
}

def addIndelSuffix(meta, vcf) {
    def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_indel')
    return tuple(new_meta, vcf)
}

def addClinvarSuffixWithTs(meta, vcf) {
    def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_clinvar')
    return tuple(new_meta, vcf)
}

workflow CLAIRS_TO_CALLING {

    //TODO Add reports for coverage stats figure and quality metrics visualization

    take:
    bam           // channel: from mapping workflow (tuple containing [meta, bam, bai])
    ref           // channel: reference genome information from input samplesheet (tuple containing [meta, refid, ref_fasta, ref_fai])
    model         // channel: basecalling model name used for platform-specific optimizations in variant calling
    bed
    vep_cache
    main:

    // Initialize empty channel for collecting software versions
    ch_versions = channel.empty()

    // Extract just the reference FASTA and its index from the reference channel
    ch_ref = ref.map { meta, reference, ref_fasta, ref_fai ->
        extractRefFiles(meta, reference, ref_fasta, ref_fai) }


    // Prepare input for ClairS-TO by joining BAM with reference and combining with chromosome list and model
    // Format: [meta, bam, bai, ref_fasta, ref_fai, chr_list, model]
    ch_input_clairs = bam
        .join(ch_ref)          // Join BAM with reference on meta
        .combine(model)        // Add model information for platform-specific optimizations

    // Run ClairS-TO variant caller on prepared input
    // Outputs SNV and indel VCFs separately
    CLAIRS_TO_CALL(ch_input_clairs)

    // Add indel and snv to meta to name index correctly while preserving ts
    ch_snv_to_index = CLAIRS_TO_CALL.out.snv
        .map { meta, vcf ->
            addSnvSuffix(meta, vcf)
        }

    ch_indel_to_index = CLAIRS_TO_CALL.out.indel
        .map { meta, vcf ->
            addIndelSuffix(meta, vcf)
        }

    BCFTOOLS_INDEX_SNV(ch_snv_to_index)
    BCFTOOLS_INDEX_INDEL(ch_indel_to_index)

    // Restore original meta id to join indel with snv together
    ch_snv_to_concat = BCFTOOLS_INDEX_SNV.out.vcf_tbi
        .map { meta, vcf, tbx ->
            restoreSnvMeta(meta, vcf, tbx)
        }

    ch_indel_to_concat = BCFTOOLS_INDEX_INDEL.out.vcf_tbi
        .map { meta, vcf, tbx ->
            restoreIndelMeta(meta, vcf, tbx)
        }

    // Merge SNV and INDEL together

    ch_to_concat =  ch_snv_to_concat
        .join(ch_indel_to_concat)

    BCFTOOLS_CONCAT(ch_to_concat)
    BCFTOOLS_SORT(BCFTOOLS_CONCAT.out.vcf)

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     VCF ANNOTATION
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

    BGZIP_VCF_INTER(BCFTOOLS_SORT.out.vcf)

    BCFTOOLS_INDEX_RAW(BGZIP_VCF_INTER.out.vcf)
    ch_in_filter_bcftools = BCFTOOLS_INDEX_RAW.out.vcf_tbi
        .join(bed)
    BCFTOOLS_FILTER_REGION(ch_in_filter_bcftools)

    ch_clairsto_out = BCFTOOLS_FILTER_REGION.out.filt_vcf

    BCFTOOLS_COUNT(ch_clairsto_out)

    ch_count_variant = BCFTOOLS_COUNT.out.vcf
        .branch { meta, vcf ->
            positive: vcf.size() > 0
                return meta
            negative: true
                return meta
        }

    ch_ref_type = ref
        .map { meta, refid, _ref_fasta, _ref_fai ->
            tuple(meta, refid) }


    ch_vep = ch_count_variant.positive
        .join(ch_clairsto_out)
        .join(ch_ref_type)
        .branch { meta, vcf, genome ->
            hg38: genome == "hg38"
            return tuple(modifyMetaId(meta, 'add_suffix', '', '', '_somatic_snp_vep'),vcf,[])
            hg19: genome == "hg19"
            return tuple(modifyMetaId(meta, 'add_suffix', '', '', '_somatic_snp_vep'),vcf,[])
            hs1: genome == "hs1"
                return tuple(modifyMetaId(meta, 'add_suffix', '', '', '_somatic_snp_vep'),vcf,[])
            }

    ch_fasta = ch_count_variant.positive
        .join(ch_clairsto_out)
        .join(ref)
        .branch { meta, _vcf, genome, ref_fasta, _ref_index ->
            hg38: genome == "hg38"
            return tuple(meta,ref_fasta)
            hg19: genome == "hg19"
            return tuple(meta,ref_fasta)
            hs1: genome == "hs1"
                return tuple(meta,ref_fasta)
            }

    def vep_cache_resolved = vep_cache ? vep_cache : []

    ENSEMBLVEP_HG38(
        ch_vep.hg38,
        "GRCh38",
        "homo_sapiens",
        params.vep_version,
        vep_cache_resolved,
        ch_fasta.hg38,
        []
    )

    ENSEMBLVEP_HG19(
        ch_vep.hg19,
        "GRCh37",
        "homo_sapiens",
        params.vep_version,
        vep_cache_resolved,
        ch_fasta.hg19,
        []
    )

    ENSEMBLVEP_HS1(
        ch_vep.hs1,
        "CHM13",
        "homo_sapiens",
        params.vep_version,
        vep_cache_resolved,
        ch_fasta.hs1,
        []
    )

    ch_vep_to_filter = ENSEMBLVEP_HG38.out.vcf
        .mix(ENSEMBLVEP_HG19.out.vcf)
        .mix(ENSEMBLVEP_HS1.out.vcf)

    ENSEMBLVEP_FILTERVEP(
        ch_vep_to_filter,
        []
    )

    ch_databases_hg38 = ch_ref_type
        .filter{meta, refid -> refid == "hg38"}
        .map { meta, _refid -> tuple(meta, 'GRCh38.p14') }
    ch_databases_hg19 = ch_ref_type
        .filter{meta, refid -> refid == "h19"}
        .map { meta, _refid -> tuple(meta, 'GRCh37.p13') }

    ch_databases_ref = ch_databases_hg38
        .mix(ch_databases_hg19)


    ch_snp_annotate = BCFTOOLS_SORT.out.vcf.join(ch_databases_ref)
        .map { meta, output, database ->
            // Prepare metadata for SNPEff annotation
            // Add reference genome type to metadata for database selection
            prepareSomaticSnpAnnotation(meta, output, database)
            }

    // -----------------------------------------------------------------------------

    // Run SNPEff annotation to add functional impact annotations to variants
    // Adds gene names, transcript IDs, effect predictions, etc.
    SNPEFF_ANNOTATE(ch_snp_annotate)

    // Combine both SNPEff-annotated and SNPSift-annotated VCFs for downstream processing
    // This preserves both annotation sets in separate files
    ch_vcf_final = SNPEFF_ANNOTATE.out.vcf

    // Compress and index VCF files for efficient storage and querying
    // Two-step process: compress with bgzip, then index with bcftools

    // Compress VCF files using bgzip (produces .vcf.gz files)
    // Bgzip creates block-compressed files that enable random access
    BGZIP_VCF_FINAL(ch_vcf_final)

    // Index compressed VCF files using BCFtools (produces .vcf.gz.tbi files)
    // Indexing allows fast retrieval of variants in specific genomic regions
    BCFTOOLS_INDEX_FINAL(BGZIP_VCF_FINAL.out.vcf)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COLLECT VERSIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Collect software version information from all processes for inclusion in reports
        and pipeline provenance tracking. Essential for reproducibility.
    */
    ch_versions = CLAIRS_TO_CALL.out.versions      // Version info from ClairS-TO
            .mix(SNPEFF_ANNOTATE.out.versions)     // Version info from SNPEff
            .mix(BGZIP_VCF_FINAL.out.versions)           // Version info from bgzip
            .mix(BCFTOOLS_INDEX_FINAL.out.versions) // Version info from BCFtools index
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        OUTPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Define the output channels that will be available to downstream processes.
        These include the annotated and indexed VCF files, version information,
        and any error channels for proper error handling.
    */

    emit:
    vcf_snpeff       = BCFTOOLS_INDEX_FINAL.out.vcf_tbi
    vcf_vep          = ENSEMBLVEP_FILTERVEP.out.output
    versions         = ch_versions                          // Output collected version information for MultiQC and reports

}
