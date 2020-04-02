process illuminaDownloadScheme {
    tag params.scriptsRepoURL

    label 'internet'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "scripts", mode: "copy"

    output:
    path "scripts/${params.scriptsDir}/make_depth_mask.py" , emit: depthmask
    path "scripts/${params.scriptsDir}/vcftagprimersites.py" , emit: vcftagprimersites
    path "scripts/${params.scriptsDir}/*"

    script:
    """
    git clone ${params.scriptsRepoURL} scripts
    """
}

process makePythonPackage {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.lowcoverage.txt", mode: 'copy'

    input:
    tuple(sampleName, path(depthmask), path(vcftagprimersites))
    
    output:
    tuple path("artic/"), emit: articpackage
    
    script:
    """
    #!/usr/bin/env python3
    mkdir artic
    touch artic/__init__.py
    mv ${depthmask} artic/
    mv ${vcftagprimersites} artic/
    """
}

process readTrimming {
    /**
    * Trims paired fastq using trim_galore (https://github.com/FelixKrueger/TrimGalore)
    * @input tuple(sampleName, path(forward), path(reverse))
    * @output trimgalore_out tuple(sampleName, path("*_val_1.fq.gz"), path("*_val_2.fq.gz"))
    */

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: '*_val_{1,2}.fq.gz', mode: 'copy'

    cpus 2

    input:
    tuple(sampleName, path(forward), path(reverse))

    output:
    tuple(sampleName, path("*_val_1.fq.gz"), path("*_val_2.fq.gz")) optional true

    script:
    """
    if [[ \$(gunzip -c ${forward} | head -n4 | wc -l) -eq 0 ]]; then
      exit 0
    else
      trim_galore --paired $forward $reverse
    fi
    """
}

process readMapping {
    /**
    * Maps trimmed paired fastq using BWA (http://bio-bwa.sourceforge.net/)
    * Uses samtools to convert to BAM, sort and index sorted BAM (http://www.htslib.org/doc/samtools.html)
    * @input 
    * @output 
    */

    tag { sampleName }

    label 'largecpu'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.sorted.bam", mode: 'copy'

    input:
        tuple(path(schemeRepo), sampleName, path(forward), path(reverse))

    output:
        tuple(sampleName, path("ref.fa"), path("${sampleName}.sorted.bam"))

    script:
        """
        ln -s ${schemeRepo}/${params.schemeDir}/${params.scheme}/${params.schemeVersion}/*.reference.fasta ref.fa
        bwa index ref.fa
        bwa mem -t ${task.cpus} ref.fa ${forward} ${reverse} | samtools view -bS | \
        samtools sort -o ${sampleName}.sorted.bam
        """
}

process makeIvarBedfile {

    tag { schemeRepo }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "ivar.bed", mode: 'copy'

    input:
    file(schemeRepo)

    output:
    file("ivar.bed")

    script:
    """
    #!/usr/bin/env python3
  
    import csv
    bedrows = []
    with open("${schemeRepo}/${params.schemeDir}/${params.scheme}/${params.schemeVersion}/nCoV-2019.scheme.bed", newline='') as bedfile:
        reader = csv.reader(bedfile, delimiter='\t')
        for row in reader:
            row[4] = '60'
            if row[3].endswith('LEFT'):
                 row.append('+')
            else: 
                row.append('-')
            bedrows.append(row)
    with open('ivar.bed', 'w', newline='') as bedfile:
        writer = csv.writer(bedfile, delimiter='\t')
        for row in bedrows:
            writer.writerow(row)
    """
}

process trimPrimerSequences {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.mapped.bam", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.mapped.primertrimmed.sorted.bam", mode: 'copy'

    input:
    tuple(path(bedfile), sampleName, path(ref), path(bam))

    output:
    tuple sampleName, path("${sampleName}.mapped.bam"), emit: mapped
    tuple sampleName, path("${sampleName}.mapped.primertrimmed.sorted.bam" ), emit: ptrim

    script:
    if (params.allowNoprimer){
        ivarCmd = "ivar trim -e"
    } else {
        ivarCmd = "ivar trim"
    }
        """
        samtools view -F4 -o ${sampleName}.mapped.bam ${bam}
        samtools index ${sampleName}.mapped.bam
        ${ivarCmd} -i ${sampleName}.mapped.bam -b ${bedfile} -m ${params.illuminaKeepLen} -q ${params.illuminaQualThreshold} -p ivar.out
        samtools sort -o ${sampleName}.mapped.primertrimmed.sorted.bam ivar.out.bam
        """
}

process callVariants {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.variants.tsv", mode: 'copy'

    input:
    tuple(sampleName, path(bam), path(ref))

    output:
    tuple sampleName, path("${sampleName}.variants.tsv")

    script:
        """
        samtools mpileup -A -d 0 --reference ${ref} -B -Q 0 ${bam} |\
        ivar variants -r ${ref} -m ${params.ivarMinDepth} -p ${sampleName}.variants -q ${params.ivarMinVariantQuality} -t ${params.ivarMinFreqThreshold}
        """
}

process callVariantsLofreq {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.lofreq.vcf", mode: 'copy'

    input:
    tuple(sampleName, path(bam), path(ref))

    output:
    tuple sampleName, path("${sampleName}.lofreq.vcf")

    script:
        """
        lofreq indelqual --dindel ${bam} -f ${ref} |\
        lofreq call --call-indels --min-bq ${params.lofreqMinBaseQuality} --min-alt-bq ${params.lofreqMinBaseQuality} \
        --min-mq ${params.lofreqMinMapQuality} --no-default-filter --use-orphan --max-depth 1000000 \
        --min-cov ${params.lofreqMinCov} -f ${ref} -o ${sampleName}.lofreq.vcf -
        """
}

process findLowCoverageRegions {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.lowcoverage.txt", mode: 'copy'

    input:
    tuple(sampleName, path(bam), path(ref), path(depthmask), path(vcftagprimersites))

    output:
    tuple sampleName, path("${sampleName}.lowcoverage.txt")

    script:
        """
        sed 's/from .vcftagprimersites import read_bed_file/from vcftagprimersites import read_bed_file/g' ${depthmask} > edited_depth_mask.py
        
        #broken due to relative import
        python edited_depth_mask.py --depth ${params.minDepthThreshold} ${ref} \
        ${bam} ${sampleName}.lowcoverage.txt
        """
}

process makeConsensus {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.primertrimmed.consensus.fa", mode: 'copy'

    input:
        tuple(sampleName, path(bam))

    output:
        tuple(sampleName, path("${sampleName}.primertrimmed.consensus.fa"))

    script:
        """
        samtools mpileup -A -d ${params.mpileupDepth} -Q0 ${bam} | \
        ivar consensus -t ${params.ivarFreqThreshold} -m ${params.ivarMinDepth} \
        -n N -p ${sampleName}.primertrimmed.consensus
        """
}
