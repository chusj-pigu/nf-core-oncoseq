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
include { BGZIP_VCF                    } from '../../../modules/local/bcftools/main.nf'
include { paramsSummaryMap             } from 'plugin/nf-schema'
include { paramsSummaryMultiqc         } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML       } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText       } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'
include { modifyMetaId                } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'
include { BCFTOOLS_CONCAT               } from '../../../modules/local/bcftools/main.nf'
include { BCFTOOLS_SORT                 } from '../../../modules/local/bcftools/main.nf'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_SNV   } from '../../../modules/local/bcftools/main.nf'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_INDEL } from '../../../modules/local/bcftools/main.nf'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_FINAL } from '../../../modules/local/bcftools/main.nf'
include { SUBCHROM_CALL_PANEL          } from '../../../modules/local/subchrom/main.nf'

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
    def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_somatic_snp')
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
    clinic_database // channel: path to clinical database (e.g., ClinVar) for variant annotation
    ch_panel_bin // subchrom panel bin file
    main:

    // Initialize empty channel for collecting software versions
    ch_versions = Channel.empty()

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

    ch_ref_type = ref
        .map { meta, refid, _ref_fasta, _ref_fai ->
            tuple(meta, refid) }

    // Branch ref channel to create database channel
    ch_databases = ch_ref_type.branch {
        hg38: { _meta, refid -> refid.matches('hg38|GRCh38') }
        hg19: { _meta, refid -> refid.matches('hg19|GRCh37') }
        other: true
            return 'Error'
    }

    // Generate error if reference is not hg38 nor hg19
    ch_error = ch_databases.other.map {
        throw new IllegalArgumentException("Unsupported reference genome: ${it.name}. Currently, only hg38/GRCh38 and hg19/GRCh37 are supported.")
    }

    ch_databases_hg38 = ch_databases.hg38
        .map { meta, _refid -> tuple(meta, 'GRCh38.p14') }
    ch_databases_hg19 = ch_databases.hg19
        .map { meta, _refid -> tuple(meta, 'GRCh37.p13') }

    ch_databases_ref = ch_databases_hg38
        .mix(ch_databases_hg19)

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

    // Convert clinical database channel to a sorted list for consistent processing
    // Ensures deterministic channel behavior when combining with multiple files
    ch_clin_db = clinic_database.toSortedList()

    // Prepare SNPEff-annotated VCFs for additional clinical annotation with SNPSift
    // Combines each annotated VCF with the clinical database (e.g., ClinVar)
    ch_snpsift_annotate = SNPEFF_ANNOTATE.out.vcf
        .combine(ch_clin_db)
        .map { meta, vcf, clin_db, clin_db_idx ->
            prepareForSnpsift(meta, vcf, clin_db, clin_db_idx) }

    // Run SNPSift annotation to add clinical significance information to variants
    // Adds ClinVar, OMIM, or other clinical database annotations
    SNPSIFT_ANNOTATE(ch_snpsift_annotate)

    // Add clinvar suffix to the sample ID in metadata for SNPSift output files
    // This helps distinguish between different annotation stages in output files
    ch_snipsift_out = SNPSIFT_ANNOTATE.out.vcf
        .map { meta, vcf ->
            addClinvarSuffixWithTs(meta, vcf)
        }

    // Combine both SNPEff-annotated and SNPSift-annotated VCFs for downstream processing
    // This preserves both annotation sets in separate files
    ch_vcf_final = SNPEFF_ANNOTATE.out.vcf
        .mix(ch_snipsift_out)

    // Compress and index VCF files for efficient storage and querying
    // Two-step process: compress with bgzip, then index with bcftools

    // Compress VCF files using bgzip (produces .vcf.gz files)
    // Bgzip creates block-compressed files that enable random access
    BGZIP_VCF(ch_vcf_final)

    // Index compressed VCF files using BCFtools (produces .vcf.gz.tbi files)
    // Indexing allows fast retrieval of variants in specific genomic regions
    BCFTOOLS_INDEX_FINAL(BGZIP_VCF.out.vcf_gz)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COLLECT VERSIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Collect software version information from all processes for inclusion in reports
        and pipeline provenance tracking. Essential for reproducibility.
    */
    ch_versions = CLAIRS_TO_CALL.out.versions      // Version info from ClairS-TO
            .mix(SNPEFF_ANNOTATE.out.versions)     // Version info from SNPEff
            .mix(BGZIP_VCF.out.versions)           // Version info from bgzip
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
    vcf              = BCFTOOLS_INDEX_FINAL.out.vcf_tbi  // Output compressed & indexed VCF files (tuple with [meta, vcf.gz, vcf.gz.tbi])
    versions         = ch_versions                          // Output collected version information for MultiQC and reports
    ch_error                                            // Output error channel for proper error handling

}
