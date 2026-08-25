/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This section imports all necessary modules, subworkflows, and functions required
    for variant calling with ClairS-TO, annotation with SNPEff/SNPSift, and VCF
    processing with BCFtools.
*/
include { CLAIRS_TO_CALL               } from '../../../modules/local/clairsto/main.nf'
include { SNPEFF_SNPEFF as SNPEFF_HG38 } from '../../../modules/nf-core/snpeff/snpeff'
include { SNPEFF_SNPEFF as SNPEFF_HG19 } from '../../../modules/nf-core/snpeff/snpeff'
include { SNPEFF_DOWNLOAD              } from '../../../modules/nf-core/snpeff/download'
include { HTSLIB_BGZIPTABIX            } from '../../../modules/nf-core/htslib/bgziptabix'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_HG38       } from '../../../modules/nf-core/ensemblvep/vep/main.nf'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_HG19       } from '../../../modules/nf-core/ensemblvep/vep/main.nf'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_HS1        } from '../../../modules/nf-core/ensemblvep/vep/main.nf'
include { ENSEMBLVEP_FILTERVEP            } from '../../../modules/nf-core/ensemblvep/filtervep/main.nf'
include { paramsSummaryMap             } from 'plugin/nf-schema'
include { paramsSummaryMultiqc         } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML       } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText       } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'
include { modifyMetaId                 } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'
include { BCFTOOLS_CONCAT               } from '../../../modules/nf-core/bcftools/concat'  // TODO: Convert to nf-core module
include { BCFTOOLS_SORT                 } from '../../../modules/nf-core/bcftools/sort'
include { BCFTOOLS_VIEW                 } from '../../../modules/nf-core/bcftools/view'
include { BCFTOOLS_FILTER_REGION                 } from '../../../modules/local/bcftools/main.nf'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_SNV   } from '../../../modules/nf-core/bcftools/index'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_INDEL } from '../../../modules/nf-core/bcftools/index'

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

def prepareSomaticSnpAnnotation(meta, output, genome) {
    def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_somatic_snp_snpeff')
    return tuple(new_meta, output, genome)
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
    ch_ref = ref
        .map { _refid, ref_fasta, ref_fai ->
            tuple(ref_fasta, ref_fai)
        }


    // Prepare input for ClairS-TO by joining BAM with reference and combining with chromosome list and model
    // Format: [meta, bam, bai, ref_fasta, ref_fai, chr_list, model]
    ch_input_clairs = bam
        .combine(ch_ref)          // Join BAM with reference on meta
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
    ch_snv_to_concat = ch_snv_to_index
        .join(BCFTOOLS_INDEX_SNV.out.index)
        .map { meta, vcf, tbx ->
            restoreSnvMeta(meta, vcf, tbx)
        }

    ch_indel_to_concat = ch_indel_to_index
        .join(BCFTOOLS_INDEX_INDEL.out.index)
        .map { meta, vcf, tbx ->
            restoreIndelMeta(meta, vcf, tbx)
        }

    // Merge SNV and INDEL together

    ch_to_concat =  ch_snv_to_concat
        .mix(ch_indel_to_concat)
        .groupTuple()

    BCFTOOLS_CONCAT(ch_to_concat)
    BCFTOOLS_SORT(BCFTOOLS_CONCAT.out.vcf)

    ch_bcftools_sort_out = BCFTOOLS_SORT.out.vcf
        .join(BCFTOOLS_SORT.out.index)

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     VCF ANNOTATION
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

    ch_in_filter_bcftools = ch_bcftools_sort_out
        .join(bed)
    BCFTOOLS_FILTER_REGION(ch_in_filter_bcftools)

    ch_clairsto_out = BCFTOOLS_FILTER_REGION.out.filt_vcf
        .join(BCFTOOLS_FILTER_REGION.out.index)

    // Count variants remaining to avoid Ensembl annotation
    BCFTOOLS_VIEW(
        ch_clairsto_out,
        [],
        [],
        []
    )

    ch_bcftools_count = BCFTOOLS_VIEW.out.vcf
        .branch { meta, vcf ->
            pos: vcf.size() > 0
                return meta
            empty: true
                return meta
            }

    ch_ref_type = ref
        .map { meta, refid, _ref_fasta, _ref_fai ->
            tuple(meta, refid) }

    ch_vep = ch_bcftools_count.pos
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

    ch_fasta = ch_bcftools_count.pos
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

    // Only run snpeff annotation for available snpeff genome cache.
    ch_valid_genomes = ch_ref_type
        .filter { meta, genome -> genome == "hg38" || genome == "hg19" }

    ch_snp_annotate = BCFTOOLS_SORT.out.vcf
        .join(ch_valid_genomes)
        .map { meta, output, genome ->
            prepareSomaticSnpAnnotation(meta, output, genome)
        }
        .branch { meta, vcf, genome ->
            hg38: genome == "hg38"
                return tuple(meta,vcf)
            hg19: genome == "hg19"
                return tuple(meta,vcf) }

    // 1. Determine the required SnpEff database per sample, from the reference genome build
    ch_db_snpeff_in = BCFTOOLS_SORT.out.vcf
        .join(ch_ref_type)
        .map { meta, _vcf, genome ->
            def db = genome == "hg38" ? "GRCh38.p14" :
                    genome == "hg19" ? "GRCh37.p13" : null
            tuple(meta, db)
        }
        .filter { meta, db -> db != null }

    // 2. Resolve which directory to look in: snpeff_cache wins, else fall back to ref_cache/snpeff_cache
    def snpeff_cache_dir = null
    if (params.snpeff_cache) {
        snpeff_cache_dir = params.snpeff_cache.contains("snpeff_cache")
            ? file(params.snpeff_cache)
            : file("${params.snpeff_cache}/snpeff_cache").isDirectory()
                ? file("${params.snpeff_cache}/snpeff_cache")
                : file(params.snpeff_cache)
    } else if (params.ref_cache) {
        snpeff_cache_dir = file("${params.ref_cache}/snpeff_cache")
    }

    // 3. Dedup by db name so we don't repeat the disk check / download once per sample
    ch_db_unique = ch_db_snpeff_in
        .map { _meta, db -> db }
        .unique()
        .map { db -> tuple([id: db], db) }   // synthetic meta keyed by db name, for module I/O

    // 4. Check on disk whether that specific db is already present in the resolved cache dir
    ch_db_checked = ch_db_unique
        .map { meta, db ->
            def present = snpeff_cache_dir \
                && snpeff_cache_dir.isDirectory() \
                && file("${snpeff_cache_dir}/${db}").exists()
            tuple(meta, db, present)
        }
        .branch { meta, db, present ->
            to_download: !present
                return tuple(meta, db)
            cached: present
                return tuple(meta, db, snpeff_cache_dir)
        }

    // 5. Download only what's missing (once per distinct db)
    SNPEFF_DOWNLOAD(ch_db_checked.to_download)

    // 6. Unify already-cached dbs with newly-downloaded ones, keyed by db name
    ch_db_cache_by_name = ch_db_checked.cached
        .map { _meta, db, cache -> tuple(id:db, cache) }
        .mix(SNPEFF_DOWNLOAD.out.cache)
        .branch { meta, db ->
            hg38: meta.id == "GRCh38.p14"
            hg19: meta.id == "GRCh37.p13"
        }

    SNPEFF_HG38(
        ch_snp_annotate.hg38,
        "GRCh38.p14",
        ch_db_cache_by_name.hg38
    )

    SNPEFF_HG19(
        ch_snp_annotate.hg19,
        "GRCh37.p13",
        ch_db_cache_by_name.hg19
    )

    ch_snpeff_out = SNPEFF_HG38.out.vcf
        .mix(SNPEFF_HG19.out.vcf)

    // Combine both SNPEff-annotated and SNPSift-annotated VCFs for downstream processing
    // This preserves both annotation sets in separate files
    ch_snpeff_out = ch_snpeff_out
        .map { meta, vcf ->
            tuple(meta,vcf,[],[])}

    // Compress and index VCF files for efficient storage and querying
    // Two-step process: compress with bgzip, then index with bcftools

    // Compress VCF files using bgzip (produces .vcf.gz files)
    // Bgzip creates block-compressed files that enable random access
    HTSLIB_BGZIPTABIX (
        ch_snpeff_out,
        "compress",
        true,
        "vcf"
    )

    ch_vcf_final_tbi = HTSLIB_BGZIPTABIX.out.output
        .join(HTSLIB_BGZIPTABIX.out.index)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COLLECT VERSIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Collect software version information from all processes for inclusion in reports
        and pipeline provenance tracking. Essential for reproducibility.
    */
    ch_versions = CLAIRS_TO_CALL.out.versions      // Version info from ClairS-TO

    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        OUTPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Define the output channels that will be available to downstream processes.
        These include the annotated and indexed VCF files, version information,
        and any error channels for proper error handling.
    */

    emit:
    vcf_snpeff       = ch_vcf_final_tbi
    vcf_vep          = ENSEMBLVEP_FILTERVEP.out.output
    versions         = ch_collated_versions                       // Output collected version information for MultiQC and reports

}
