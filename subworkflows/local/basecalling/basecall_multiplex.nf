/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { DORADO_BASECALL                           } from '../../../modules/local/dorado/main.nf'
include { DORADO_DEMULTIPLEX                        } from '../../../modules/local/dorado/main.nf'
include { DORADO_TRIM                               } from '../../../modules/local/dorado/main.nf'
include { SAMTOOLS_QSFILTER                         } from '../../../modules/local/samtools/main.nf'
include { SAMTOOLS_TOFASTQ                          } from '../../../modules/local/samtools/main.nf'
include { HTSLIB_BGZIPTABIX                         } from '../../../modules/nf-core/htslib/bgziptabix/main'
include { SEQKIT_STATS                              } from '../../../modules/nf-core/seqkit/stats/main'

include { softwareVersionsToYAML                    } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { modifyMetaId                              } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BASECALL_MULTIPLEX {

    //TODO Add reports for read stats figure and tables

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_demux       // channel : demux samplesheet read in from --demux_samplesheet
    main:

    ch_versions = channel.empty()

    DORADO_BASECALL(ch_samplesheet)

    ch_to_dmux = ch_demux
        .map { meta, _sample, _barcode, kit ->
            tuple(meta,kit) }
        .join(DORADO_BASECALL.out.ubam)

    DORADO_DEMULTIPLEX(ch_to_dmux)

    split_bams_ch = DORADO_DEMULTIPLEX.out.demux
        .map { meta, dir ->
            def bams = []
            dir.toFile().eachFileRecurse { f ->
                if (f.name.endsWith('.bam')) {
                    bams << file(f)
                }
            }
            tuple(meta, bams)
        }
        .transpose()
        .filter { meta, bam ->
            bam.name.contains("barcode") }              // Exclude unclassified bam
        .map { meta, bam ->

            def barcode = bam.name.replaceFirst(/^.*?(barcode\d{2}).*$/, '$1')
            def new_meta = [project: meta.id, barcode: barcode ]

            tuple(new_meta, bam)
        }

    // Join back to samples and rename meta to sample
    ch_ubam_joined = ch_demux
        .map { meta, sample, barcode, kit ->
            def new_meta = [project: meta.id, barcode: barcode]
            tuple(new_meta, sample, kit)
        }
        .join(split_bams_ch)
        .map { _meta, sample, kit, bam ->
            tuple(id:sample, kit, bam)}     // Drop barcodes

    DORADO_TRIM(ch_ubam_joined)

    SAMTOOLS_QSFILTER(DORADO_TRIM.out.bam)

     // Get the sample_ids from the demux_samplesheet that specify each barcode to each sample_id
    ch_new_sample_ids_pass = SAMTOOLS_QSFILTER.out.ubam_pass
        .map { meta, ubam ->
            def meta_suffix = ubam.baseName.tokenize('_')[-1].replace('.bam', '')       // Add pass to meta in tuples for output naming
            def new_meta = [id : meta.id + '_' + meta_suffix]
            tuple(new_meta, ubam)
        }

    ch_new_sample_ids_fail = SAMTOOLS_QSFILTER.out.ubam_fail
        .map { meta, ubam ->
            def meta_suffix = ubam.baseName.tokenize('_')[-1].replace('.bam', '')       // Add fail to meta in tuples for output naming
            def new_meta = [id : meta.id + '_' + meta_suffix]
            tuple(new_meta, ubam)
        }

    SAMTOOLS_TOFASTQ(ch_new_sample_ids_pass.mix(ch_new_sample_ids_fail))

    SEQKIT_STATS(SAMTOOLS_TOFASTQ.out.fq)

    ch_fastq_pass = SAMTOOLS_TOFASTQ.out.fq
        .filter { meta, _reads -> meta.id.contains('pass')}
        .map { meta, fq ->
         def meta_restore = [id : meta.id.replace('_pass', '')]
         tuple(meta_restore, fq)}

    ch_stats_pass = SEQKIT_STATS.out.stats
        .filter { meta, _tbl -> meta.id.contains('pass')}
        .map { meta, table ->
        def meta_restore = modifyMetaId(meta, 'remove_suffix', '', '', '_pass')
        tuple(meta_restore, table)}

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COLLECT VERSIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    ch_versions = DORADO_BASECALL.out.versions
        .mix(DORADO_DEMULTIPLEX.out.versions)
        .mix(DORADO_TRIM.out.versions)
        .mix(SAMTOOLS_QSFILTER.out.versions)
        .mix(SAMTOOLS_TOFASTQ.out.versions)

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

    emit:
    fastq          = ch_fastq_pass
    stats_pass     = ch_stats_pass
    versions       = ch_collated_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
