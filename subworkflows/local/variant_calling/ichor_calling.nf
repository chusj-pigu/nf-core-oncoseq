/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { HMMCOPY_WIG       } from '../../../modules/local/ichorcna/main.nf'
include { ICHORCNA_DOWNLOAD } from '../../../modules/local/ichorcna/main.nf'
include { ICHORCNA          } from '../../../modules/local/ichorcna/main.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ICHORCNA_CALLING {

    //TODO Add reports for coverage stats figure ?

    take:
    bam                 // channel: from mapping workflow, includes index
    cfdna_samplesheet         // reference channel with index
    ichor_bin           // Channel from params.ichor_bin_size
    mapq_wig            // Channel from params.minmapq_wig
    main:

    ch_hmmwig_in = bam
        .combine(ichor_bin)
        .combine(mapq_wig)
    HMMCOPY_WIG(ch_hmmwig_in)

    ch_samplesheet_wig = cfdna_samplesheet
        .join(HMMCOPY_WIG.out.wig)
        .map { meta, purity, _filter, wig ->
            tuple(meta, wig, purity) }


    ICHORCNA_DOWNLOAD(ch_samplesheet_wig)

    ICHORCNA(ICHORCNA_DOWNLOAD.out.seq_info)

    ch_versions = HMMCOPY_WIG.out.versions
        .mix(ICHORCNA.out.versions)

    emit:
    ichorcna_plot       = ICHORCNA.out.ichor_dir   // TODO: Quarto report
    versions            = ch_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
