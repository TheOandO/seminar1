/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { STAR_ALIGN } from '../modules/nf-core/star/align/main'
include { TRIMGALORE } from '../modules/nf-core/trimgalore/main'
include { FASTQC } from '../modules/nf-core/fastqc/main'
include { MULTIQC } from '../modules/nf-core/multiqc/main'
include { SALMON_QUANT } from '../modules/nf-core/salmon/quant/main'
include { DUPRADAR } from '../modules/nf-core/dupradar/main'
include { QUALIMAP_RNASEQ } from '../modules/nf-core/qualimap/rnaseq/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SEMINAR1 {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()
    FASTQC(ch_samplesheet)

    // FastQC output files
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{ it[1] })

    // Trimgalore
    TRIMGALORE(ch_samplesheet)

    // Adapter-trimmed FASTQ files, FastQC files and logfiles
    ch_trimmed_reads = TRIMGALORE.out.reads
    ch_multiqc_files = ch_multiqc_files.mix(ch_trimmed_reads.collect{ it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.zip.collect{ it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.log.collect{ it[1] })


    // Run STAR
    ch_star_index = channel.value(file(params.star_index, checkIfExists: true))
    ch_gtf = channel.value(file(params.gtf, checkIfExists: true))
    STAR_ALIGN(
        ch_trimmed_reads,
        ch_star_index.map { [ [:], it ] },
        ch_gtf.map { [ [:], it ] },
        false
    
    )
    // STAR output files
    ch_star_bams = STAR_ALIGN.out.bam
    ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN.out.log_final.collect{ it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(ch_star_bams.collect{ it[1] })


    // Run SALMON quantification
    ch_salmon_index = channel.value(file(params.salmon_index, checkIfExists: true))
    ch_transcriptome = channel.value(file(params.transcriptome, checkIfExists: true))
    SALMON_QUANT(
        ch_trimmed_reads,
        ch_salmon_index,
        ch_gtf,
        ch_transcriptome,
        false,
        false
    )
    // SALMON output files
    ch_multiqc_files = ch_multiqc_files.mix(SALMON_QUANT.out.results.collect{ it[1] })

    // DupRadar
    DUPRADAR(
        ch_star_bams,
        ch_gtf
    )

    ch_multiqc_files = ch_multiqc_files.mix(DUPRADAR.out.multiqc.collect{ it[1] })

    // Qualimap RNA-seq
    QUALIMAP_RNASEQ(
        ch_star_bams,
        ch_gtf
    )

    ch_multiqc_files = ch_multiqc_files.mix(QUALIMAP_RNASEQ.out.results.collect{ it[1] })

    // Run MultiQC
    MULTIQC(
    ch_multiqc_files.collect(),
    [], [], [], [], []
    )
































    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
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

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'seminar1_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
