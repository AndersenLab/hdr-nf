
nextflow.enable.dsl=2

/*
    Parameters
*/

date = new Date().format('yyyyMMdd')
log.info("Source: ${params.source}.")

if (params.debug) {
    println """
    hdr-nf does not have a debug mode yet.....
    """
} else {
    println """
    Running hdr-nf in standard mode.
    """
}

def bam_dir = [
    c_elegans   : '/vast/eande106/data/c_elegans/WI/alignments',
    c_tropicalis: '/vast/eande106/data/c_tropicalis/WI/alignments',
    c_briggsae  : '/vast/eande106/data/c_briggsae/WI/alignments'
]

def ref_genome = [
    c_elegans   : '/vast/eande106/data/c_elegans/genomes/PRJNA13758/WS283/c_elegans.PRJNA13758.WS283.genome.fa',
    c_tropicalis: '/vast/eande106/data/c_tropicalis/genomes/NIC58_nanopore/June2021/c_tropicalis.NIC58_nanopore.June2021.genome.fa',
    c_briggsae  : '/vast/eande106/data/c_briggsae/genomes/QX1410_nanopore/Feb2020/c_briggsae.QX1410_nanopore.Feb2020.genome.fa'
]

def ref_vcf = [
    c_elegans   : '/vast/eande106/data/c_elegans/WI/variation/20250625/vcf/WI.20250625.hard-filter.isotype.vcf.gz',
    c_tropicalis: '/vast/eande106/data/c_tropicalis/WI/variation/20250627/vcf/WI.20250627.hard-filter.isotype.vcf.gz',
    c_briggsae  : '/vast/eande106/data/c_briggsae/WI/variation/20250626/vcf/WI.20250626.hard-filter.isotype.vcf.gz'
]

def ref_str = [
    c_elegans   : 'N2',
    c_tropicalis: 'NIC58',
    c_briggsae  : 'QX1410'
]

def invcf = params.vcf ?: ref_vcf[params.species]
def ingenome = params.reference ?: ref_genome[params.species]
def inbam = params.bam ?: bam_dir[params.species]
def refstrain = prams.ref ?: ref_str[params.species]

def log_summary() {
    // Corrected log summary function to print information instead of recursive call
    log.info("Workflow summary: \n" + 
             "Debug mode: ${params.debug}\n" + 
             "Output directory: ${params.output}\n")
    
    // Show help if requested
    if (params.help) {
        log.info("Help requested, exiting.")
        exit 1
    }
}

workflow {

    ch_vcf    = channel.fromPath(invcf, checkIfExists: true)
    ch_genome = channel.fromPath(ingenome, checkIfExists: true)
    ch_bam_dir = channel.fromPath(inbam, checkIfExists: true)

    GENERATE_SAMPLE_LIST_AND_WINDOWS(ch_vcf, ch_genome)

    ch_samples = GENERATE_SAMPLE_LIST_AND_WINDOWS.out.sample_list
    ch_windows = GENERATE_SAMPLE_LIST_AND_WINDOWS.out.windows_bed

    // Fan out: one strain per emission, trimming the newline
    ch_strains = ch_samples
        .splitText()
        .map { line -> line.trim() }

    // One task per strain
    COUNT_VARIANTS_PER_WINDOW(
        ch_strains,
        ch_vcf.first(), //first allows for channel re-usage across ch_strain lines
        ch_windows.first()
    )

    ch_all_counts = COUNT_VARIANTS_PER_WINDOW.out.variant_counts
        .map { strain, tsv -> tsv }
        .collect()

    MERGE_VARIANT_COUNTS(ch_all_counts)

    MOSDEPTH_COVERAGE(
        ch_strains,
        ch_bam_dir.first(),
        ch_windows.first()
    )

    ch_all_thresholds = MOSDEPTH_COVERAGE.out.thresholds_bed
        .map { strain, bed -> bed }
        .collect()

    MERGE_THRESHOLDS(ch_all_thresholds)

    coverage_ch = MERGE_THRESHOLDS.out.merged_thresholds
    varct_ch    = MERGE_VARIANT_COUNTS.out.merged_counts

    CALL_HDRS(
        MERGE_THRESHOLDS.out.merged_thresholds,
        MERGE_VARIANT_COUNTS.out.merged_counts,
        ch_windows,
        params.refstrain,
        params.covthresh,
        params.vcthresh
    )

}


process GENERATE_SAMPLE_LIST_AND_WINDOWS {
    tag "${invcf.simpleName}"
    label 'process_low'
    container 'docker://quay.io/biocontainers/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127e1b5a6d81d97e2b7e4d53b3b6b0fe1e6023e-0'

    input:
    path invcf
    path ingenome

    output:
    path "samples.txt",  emit: sample_list
    path "windows.bed",  emit: windows_bed

    script:
    """
    bcftools query -l ${invcf} > samples.txt

    samtools faidx ${ingenome}
    cut -f1,2 ${ingenome}.fai > genome.txt
    bedtools makewindows -g genome.txt -w 1000 > windows.bed
    """
}

process COUNT_VARIANTS_PER_WINDOW {
    tag "${strain}_var"
    label 'process_med'
    container 'docker://quay.io/biocontainers/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127e1b5a6d81d97e2b7e4d53b3b6b0fe1e6023e-0'

    input:
    val strain
    path invcf
    path windows_bed

    output:
    tuple val(strain), path("${strain}.variant_counts.tsv"), emit: variant_counts

    script:
    """
    bcftools view -s ${strain} ${invcf} | \
        bcftools filter -i 'GT="alt"' -Oz -o ${strain}.vcf.gz

    bedtools coverage -a ${windows_bed} -b ${strain}.vcf.gz -counts > ${strain}.variant_counts.tsv
    """
}

process MERGE_VARIANT_COUNTS {
    tag "merge_var"
    label 'process_low'
    publishDir "${params.output}", mode: 'copy'
    container 'docker://quay.io/biocontainers/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127e1b5a6d81d97e2b7e4d53b3b6b0fe1e6023e-0'

    input:
    path variant_count_files

    output:
    path "all_variant_counts.tsv", emit: merged_counts

    script:
    """
    for f in ${variant_count_files}; do
        strain=\$(basename \$f .variant_counts.tsv)
        awk -v s="\$strain" 'BEGIN{OFS="\\t"} {print s, \$0}' \$f
    done > all_variant_counts.tsv
    """
}

process MOSDEPTH_COVERAGE {
    tag "${strain}_cov"
    label 'process_med'
    container 'docker://quay.io/biocontainers/mosdepth:0.3.8--hd299d5a_0'

    input:
    val strain
    path bam_dir
    path windows_bed

    output:
    tuple val(strain), path("${strain}.thresholds.bed"), emit: thresholds_bed

    script:
    """
    mosdepth ${strain} ${bam_dir}/${strain}.bam -b ${windows_bed} -t 4 -T 1,2,5 -n
    gunzip ${strain}.thresholds.bed.gz
    """
}

process MERGE_THRESHOLDS {
    tag "merge_thresh"
    label 'process_low'
    publishDir "${params.output}", mode: 'copy'
    container 'docker://quay.io/biocontainers/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127e1b5a6d81d97e2b7e4d53b3b6b0fe1e6023e-0'

    input:
    path threshold_files

    output:
    path "all_thresholds.tsv", emit: merged_thresholds

    script:
    """
    for f in ${threshold_files}; do
        strain=\$(basename \$f .thresholds.bed)
        awk -v s="\$strain" 'BEGIN{OFS="\\t"} !/^#/ {print s, \$0}' \$f
    done > all_thresholds.tsv
    """
}

process CALL_HDRS {
    tag "call_hdrs"
    label 'process_med'
    publishDir "${params.output}", mode: 'copy'
    container 'docker://docker.io/nicmoya/hdr_r_image:2026_07_24'

    input:
    path coverage_df
    path varct_df
    path bins_1kb_stripped
    val refstrain
    val covthresh
    val vcthresh

    output:
    path "hdrs.tsv", emit: hdrs

    script:
    """
    export OMP_NUM_THREADS=${task.cpus}
    Rscript call_hdr.R \
        ${coverage_df} \
        ${varct_df} \
        ${bins_1kb_stripped} \
        ${refstrain} \
        ${covthresh} \
        ${vcthresh}
    """
}