// ARTIC ncov workflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include {articStageScheme} from '../modules/artic.nf' 
include {articGuppyPlex} from '../modules/artic.nf' 
include {articMinIONNanopolish} from  '../modules/artic.nf' 
include {articMinIONMedaka} from  '../modules/artic.nf'
include {articRemoveUnmappedReads} from '../modules/artic.nf' 

include {makeQCCSV} from '../modules/qc.nf'
include {writeQCSummaryCSV} from '../modules/qc.nf'

include {collateSamples} from '../modules/upload.nf'


//CSB5
include {compressOutputs} from '../modules/send.nf'

// import subworkflows
include {CLIMBrsync} from './upload.nf'


// workflow component for artic pipeline
workflow sequenceAnalysisNanopolish {
    take:
      ch_runFastqDirs
      ch_fast5Pass
      ch_seqSummary
    
    main:
      articStageScheme( Channel.fromPath( "${params.schemeRepo}" ) )
      
      articGuppyPlex(ch_runFastqDirs.flatten())

      articMinIONNanopolish(articGuppyPlex.out.fastq
                                          .combine(articStageScheme.out.scheme)
                                          .combine(ch_fast5Pass)
                                          .combine(ch_seqSummary))

      articRemoveUnmappedReads(articMinIONNanopolish.out.mapped)

      makeQCCSV(articMinIONNanopolish.out.ptrim
                                     .join(articMinIONNanopolish.out.consensus_fasta, by: 0)
                                     .combine(articStageScheme.out.reffasta))

      makeQCCSV.out.csv.splitCsv()
                       .unique()
                       .branch {
                           header: it[-1] == 'qc_pass'
                           fail: it[-1] == 'FALSE'
                           pass: it[-1] == 'TRUE'
                       }
                       .set { qc }

     writeQCSummaryCSV(qc.header.concat(qc.pass).concat(qc.fail).toList())

     collateSamples(qc.pass.map{ it[0] }
                           .join(articMinIONNanopolish.out.consensus_fasta, by: 0)
                           .join(articRemoveUnmappedReads.out))


     compressOutputs(qc.pass.map{ it[0] }
                            .join(articMinIONNanopolish.out.consensus_fasta, by: 0)
			    .map {it[1]}
			    .collect()
			    )

    emit:
      qc_pass = collateSamples.out

}

workflow sequenceAnalysisMedaka {
    take:
      ch_runFastqDirs

    main:
      articStageScheme( Channel.fromPath( "${params.schemeRepo}" ) )
 
      articGuppyPlex(ch_runFastqDirs.flatten())

      articMinIONMedaka(articGuppyPlex.out.fastq
                                      .combine(articStageScheme.out.scheme))

      articRemoveUnmappedReads(articMinIONMedaka.out.mapped)

      makeQCCSV(articMinIONMedaka.out.ptrim.join(articMinIONMedaka.out.consensus_fasta, by: 0)
                           .combine(articStageScheme.out.reffasta))

      makeQCCSV.out.csv.splitCsv()
                       .unique()
                       .branch {
                           header: it[-1] == 'qc_pass'
                           fail: it[-1] == 'FALSE'
                           pass: it[-1] == 'TRUE'
                       }
                       .set { qc }

     writeQCSummaryCSV(qc.header.concat(qc.pass).concat(qc.fail).toList())

     collateSamples(qc.pass.map{ it[0] }
                           .join(articMinIONMedaka.out.consensus_fasta, by: 0)
                           .join(articRemoveUnmappedReads.out))

     compressOutputs(qc.pass.map{ it[0] }
                            .join(articMinIONMedaka.out.consensus_fasta, by: 0)
                            .map {it[1]}
                            .collect()
                            )

    emit:
      qc_pass = collateSamples.out

}


workflow articNcovNanopore {
    take:
      ch_fastqDirs
    
    main:
      if ( params.nanopolish ) {
          Channel.fromPath( "${params.fast5_pass}" )
                 .set{ ch_fast5Pass }

          Channel.fromPath( "${params.sequencing_summary}" )
                 .set{ ch_seqSummary }

          sequenceAnalysisNanopolish(ch_fastqDirs, ch_fast5Pass, ch_seqSummary)

      } else if ( params.medaka ) {
          sequenceAnalysisMedaka(ch_fastqDirs)
      }

      if ( params.upload ) {

        Channel.fromPath("${params.CLIMBkey}")
               .set{ ch_CLIMBkey }

        CLIMBrsync(sequenceAnalysis.out.qc_pass, ch_CLIMBkey )
      }
}

