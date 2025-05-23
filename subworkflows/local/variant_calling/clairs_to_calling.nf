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
include { BCFTOOLS_INDEX               } from '../../../modules/local/bcftools/main.nf'
include { paramsSummaryMap             } from 'plugin/nf-schema'
include { paramsSummaryMultiqc         } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML       } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText       } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CLAIRS_TO_CALLING {

    //TODO Add reports for coverage stats figure and quality metrics visualization

    take:
    bam           // channel: from mapping workflow (tuple containing [meta, bam, bai])
    ref           // channel: reference genome information from input samplesheet (tuple containing [meta, refid, ref_fasta, ref_fai])
    chr_list      // channel: comma-separated list of chromosomes to include for variant calling (from params.chr_list)
    model         // channel: basecalling model name used for platform-specific optimizations in variant calling
    clinic_database // channel: path to clinical database (e.g., ClinVar) for variant annotation
    main:

    // Initialize empty channel for collecting software versions
    ch_versions = Channel.empty()

    // Extract just the reference FASTA and its index from the reference channel
    ch_ref = ref
        .map { meta, _ref, ref_fasta, ref_fai ->
            tuple(meta, ref_fasta, ref_fai) }

    // Prepare input for ClairS-TO by joining BAM with reference and combining with chromosome list and model
    // Format: [meta, bam, bai, ref_fasta, ref_fai, chr_list, model]
    ch_input_clairs = bam
        .join(ch_ref)          // Join BAM with reference on meta
        .combine(chr_list)     // Add chromosome list
        .combine(model)        // Add model information for platform-specific optimizations

    // Run ClairS-TO variant caller on prepared input
    // Outputs SNV and indel VCFs separately
    CLAIRS_TO_CALL(ch_input_clairs)

    // Extract reference genome type (hg38/GRCh38 or hg19/GRCh37) for database selection
    ch_ref_type = ref
        .map { meta, refid, _ref_fasta, _ref_fai ->
            tuple(meta, refid) }

    // Branch reference channel based on genome build to select appropriate annotation database
    // SNPEff requires specific database names for different genome builds
    // This should be done early on to avoid unnecessary processing of unsupported references
    ch_databases = ch_ref_type.branch {
        hg38: { meta, refid -> refid.matches('hg38|GRCh38') }  // Match hg38/GRCh38 references
        hg19: { meta, refid -> refid.matches('hg19|GRCh37') }  // Match hg19/GRCh37 references
        other: true                                           // Catch any unsupported references
            return 'Error'
    }

    // Generate descriptive error if an unsupported reference genome is provided
    // ClairS-TO and annotation databases are optimized for human genomes
    ch_error = ch_databases.other.map {
        throw new IllegalArgumentException("Unsupported reference genome: ${it.name}. Currently, only hg38/GRCh38 and hg19/GRCh37 are supported.")
    }

    // Map reference genome builds to specific SNPEff database names
    // SNPEff requires exact database names that match installed databases
    ch_databases_hg38 = ch_databases.hg38
        .map { meta, _refid -> tuple(meta, 'GRCh38.p14') }  // Map hg38 to GRCh38.p14 SNPEff database
    ch_databases_hg19 = ch_databases.hg19
        .map { meta, _refid -> tuple(meta, 'GRCh37.p13') }  // Map hg19 to GRCh37.p13 SNPEff database

    // Combine both database channels for downstream processing
    ch_databases_ref = ch_databases_hg38
        .mix(ch_databases_hg19)

    // Prepare indel variant calls for annotation
    // Only pass downstream if the VCF file exists
    ch_clairs_indel = CLAIRS_TO_CALL.out.indel
        .filter { it[1] != null && file(it[1]).exists() } // Only keep if VCF exists
        .join(ch_databases_ref)                      // Join variant calls with database information
        .map { meta, output, database ->
            def meta_type = meta.id + '_indel'       // Add variant type to sample ID for output file naming
            def new_meta = meta.clone()
            new_meta.id = meta_type
            // preserve meta.ts if present
            if (meta.containsKey('ts')) new_meta.ts = meta.ts
                tuple(new_meta, output, database) }

    // Prepare SNV variant calls for annotation
    // Only pass downstream if the VCF file exists
    ch_clairs_snv = CLAIRS_TO_CALL.out.snv
        .filter { it[1] != null && file(it[1]).exists() } // Only keep if VCF exists
        .join(ch_databases_ref)                      // Join variant calls with database information
        .map { meta, output, database ->
            def meta_type = meta.id + '_snv'         // Add variant type to sample ID for output file naming
            def new_meta = meta.clone()
            new_meta.id = meta_type
            if (meta.containsKey('ts')) new_meta.ts = meta.ts
                tuple(new_meta, output, database) }

    // Combine indel and SNV channels for parallel annotation with SNPEff
    // Both variant types use the same annotation process
    ch_snp_annotate = ch_clairs_indel
        .mix(ch_clairs_snv)                 // Mix indel and SNV channels for parallel processing

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
        .map { vcf_tuple, clin_db ->
            def (meta, vcf) = vcf_tuple
            tuple(meta, vcf, clin_db)
        }

    // Run SNPSift annotation to add clinical significance information to variants
    // Adds ClinVar, OMIM, or other clinical database annotations
    SNPSIFT_ANNOTATE(ch_snpsift_annotate)

    // Add clinvar suffix to the sample ID in metadata for SNPSift output files
    // This helps distinguish between different annotation stages in output files
    ch_snipsift_out = SNPSIFT_ANNOTATE.out.vcf
        .map { meta, vcf ->
            def meta_type = meta.id + '_clinvar'     // Add clinvar suffix to sample ID for output file naming
                tuple(id:meta_type, vcf) }

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
    BCFTOOLS_INDEX(BGZIP_VCF.out.vcf_gz)

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
            .mix(BCFTOOLS_INDEX.out.versions)      // Version info from BCFtools

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        OUTPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Define the output channels that will be available to downstream processes.
        These include the annotated and indexed VCF files, version information,
        and any error channels for proper error handling.
    */

    emit:
    vcf              = BCFTOOLS_INDEX.out.vcf_tbi  // Output compressed & indexed VCF files (tuple with [meta, vcf.gz, vcf.gz.tbi])
    versions         = ch_versions                  // Output collected version information for MultiQC and reports
    ch_error                                        // Output error channel for proper error handling

}
