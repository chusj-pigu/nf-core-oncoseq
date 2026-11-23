/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { HMMCOPY_WIG       } from '../../../modules/local/ichorcna/main.nf'
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
    ref
    cfdna_samplesheet         // reference channel with index
    ichor_bin           // Channel from params.ichor_bin_size
    mapq_wig            // Channel from params.minmapq_wig
    main:

    ch_hmmwig_in = bam
        .combine(ichor_bin)
        .combine(mapq_wig)
    HMMCOPY_WIG(ch_hmmwig_in)

    ch_ichor_in = cfdna_samplesheet
        .combine(ref)
        .join(HMMCOPY_WIG.out.wig)
        .map { meta, purity, _filter, genome, _genome_fasta, wig ->
            tuple(meta, genome.id, wig, purity) }


    ICHORCNA(ch_ichor_in)

    ch_versions = HMMCOPY_WIG.out.versions
        .mix(ICHORCNA.out.versions)

    emit:
    ichorcna_plot       = ICHORCNA.out.plot   // TODO: Quarto report
    versions            = ch_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
