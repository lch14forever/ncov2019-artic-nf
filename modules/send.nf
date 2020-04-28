process compressOutputs {

    publishDir "${params.outdir}/", mode: 'copy'

    input:
    file(fasta)

    output:
    file("*.tar.gz")

    script:
    """
    mkdir ${params.prefix}_consensus
    cp ${fasta} ${params.prefix}_consensus/
    tar -czf ${params.prefix}_consensus.tar.gz ${params.prefix}_consensus
    """
}