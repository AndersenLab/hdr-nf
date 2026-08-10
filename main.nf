
nextflow.enable.dsl=2

/*
    Parameters
*/

date = new Date().format('yyyyMMdd')
log.info("Species: ${params.species}.")

if (params.debug) {
    println """
    hdr-nf does not have a debug mode yet.....
    """
} else {
    println """
    Running hdr-nf in standard mode.
    """
}

if (params.species == 'c_briggsae' && !params.samplesheet) {
    error """
    A samplesheet is required for c_briggsae.

    Example:
      --samplesheet samples.csv
    """.stripIndent()
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

def cov_thresh = [
    c_elegans   : '0.8',
    c_tropicalis: '0.9',
    c_briggsae  : '0.9'
]

def var_thresh = [
    c_elegans   : '11',
    c_tropicalis: '9',
    c_briggsae  : '11'
]

def invcf = params.vcf ?: ref_vcf[params.species]
def ingenome = params.refgen ?: ref_genome[params.species]
def inbam = params.bam ?: bam_dir[params.species]
def refstrain = params.str ?: ref_str[params.species]
def covthresh = params.pbt ?: cov_thresh[params.species]
def vcthresh = params.vct ?: var_thresh[params.species]

def paramSummary = [
    'Species'       : params.species,
    'VCF'           : invcf,
    'Reference genome'     : ingenome,
    'BAM directory' : inbam,
    'REF strain'    : refstrain,
    'Coverage threshold' : covthresh,
    'Variant threshold' : vcthresh,
    'Output directory'    : params.output
]

if (covthresh == null || vcthresh == null) {
    error "Thresholds not defined for species: ${params.species}"
}

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

def maxLen = paramSummary.keySet().collect { k -> k.size() }.max()
    def summary = paramSummary.collect { k, v ->
        "${k.padRight(maxLen)} : ${v}"
    }.join('\n    ')

    log.info """
    =========================================
      Pipeline Parameters
    =========================================
    ${summary}
    =========================================
    """.stripIndent()

workflow {
    
    if (params.species == 'c_briggsae') {

        ch_samples = channel
            .fromPath(params.samplesheet, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                tuple(
                    row.group,
                    row.strain,
                    file(row.vcf),
                    file(row.ref),
                    file(row.bam_path),
                    row.refstrain
                )
            }

        ch_group_refs = ch_samples
            .map { group, strain, vcf, ref, bam, group_refstrain ->
                tuple(group, ref, group_refstrain)
            }
            .unique()

        GENERATE_WINDOWS_GROUP(ch_group_refs)

        ch_group_metadata = GENERATE_WINDOWS_GROUP.out.windows_bed

        ch_samples_for_join = ch_samples
            .map { group, strain, vcf, ref, bam_path, group_refstrain ->
                tuple(group, strain, vcf, bam_path)
            }

        ch_samples_with_windows = ch_samples_for_join
            .join(GENERATE_WINDOWS_GROUP.out.windows_bed)

        ch_variant_inputs = ch_samples_with_windows
            .map { group, strain, vcf, bam_path, group_refstrain, windows ->
                tuple(group, strain, vcf, windows)
            }

        ch_coverage_inputs = ch_samples_with_windows
            .map { group, strain, vcf, bam_path, group_refstrain, windows ->
                tuple(group, strain, bam_path, windows)
            }

    } else {

        ch_vcf     = channel.fromPath(invcf, checkIfExists: true)
        ch_genome  = channel.fromPath(ingenome, checkIfExists: true)
        ch_bam_dir = channel.fromPath(inbam, checkIfExists: true)

        GENERATE_SAMPLE_LIST_AND_WINDOWS(
            ch_vcf,
            ch_genome
        )

        ch_shared_windows = GENERATE_SAMPLE_LIST_AND_WINDOWS.out.windows_bed

        ch_group_metadata = ch_shared_windows.map { windows ->
            tuple("GLOBAL", refstrain, windows)
        }

        ch_strains = GENERATE_SAMPLE_LIST_AND_WINDOWS.out.sample_list
            .splitText()
            .map { line ->
                tuple("GLOBAL", line.trim())
        }

        ch_variant_inputs = ch_strains
            .combine(ch_vcf)
            .combine(ch_shared_windows)

        ch_coverage_inputs = ch_strains
            .combine(ch_bam_dir)
            .combine(ch_shared_windows)
 
    }

    COUNT_VARIANTS_PER_WINDOW(ch_variant_inputs)
    MOSDEPTH_COVERAGE(ch_coverage_inputs)

    ch_all_counts = COUNT_VARIANTS_PER_WINDOW.out.variant_counts
        .map { group, strain, tsv -> tsv }
        .collect()

    MERGE_VARIANT_COUNTS(ch_all_counts)

    ch_all_thresholds = MOSDEPTH_COVERAGE.out.thresholds_bed
        .map { group, strain, bed -> bed }
        .collect()

    MERGE_THRESHOLDS(ch_all_thresholds)

    ch_coverage_merged = MERGE_THRESHOLDS.out.merged_thresholds
    ch_varct_merged = MERGE_VARIANT_COUNTS.out.merged_counts

    ch_hdr_inputs = ch_group_metadata
        .combine(ch_coverage_merged)
        .combine(ch_varct_merged)

    CALL_HDRS(
        ch_hdr_inputs,
        covthresh,
        vcthresh
    )

}

process GENERATE_SAMPLE_LIST_AND_WINDOWS {
    tag "${vcf_file.simpleName}"
    label 'process_low'
    container 'docker://docker.io/nicmoya/bedvcf_hdr_image:2026_07_24'
    beforeScript   = 'module load singularity'

    input:
    path vcf_file
    path genome_file

    output:
    path "samples.txt",  emit: sample_list
    path "windows.bed",  emit: windows_bed

    script:
    """
    bcftools query -l ${vcf_file} > samples.txt

    samtools faidx ${genome_file}
    cut -f1,2 ${genome_file}.fai > genome.txt
    bedtools makewindows -g genome.txt -w 1000 > windows.bed
    """
}

process GENERATE_WINDOWS_GROUP {
    tag "${group}"

    label 'process_low'
    container 'docker://docker.io/nicmoya/bedvcf_hdr_image:2026_07_24'
    beforeScript = 'module load singularity'

    input:
    tuple val(group), path(genome_file), val(refstrain)

    output:
    tuple val(group), val(refstrain), path("windows.bed"),
        emit: windows_bed

    script:
    """
    samtools faidx ${genome_file}
    cut -f1,2 ${genome_file}.fai > genome.txt
    bedtools makewindows -g genome.txt -w 1000 > windows.bed
    """
}

process COUNT_VARIANTS_PER_WINDOW {
    tag "${strain}_var"
    label 'process_med'
    container 'docker://docker.io/nicmoya/bedvcf_hdr_image:2026_07_24'
    beforeScript = 'module load singularity'

    input:
    tuple val(group), val(strain), path(vcf_file), path(windows_bed)

    output:
    tuple val(group), val(strain), path("${strain}.variant_counts.tsv"),
        emit: variant_counts

    script:
    """
    bcftools view -s ${strain} ${vcf_file} | \
        bcftools filter -i 'GT="alt"' -Oz -o ${strain}.vcf.gz

    bedtools coverage \
        -a ${windows_bed} \
        -b ${strain}.vcf.gz \
        -counts | \
        awk -v g="${group}" -v s="${strain}" \
            'BEGIN{OFS="\t"} {print g, s, \$0}' \
        > ${strain}.variant_counts.tsv
    """
}

process MOSDEPTH_COVERAGE {
    tag "${strain}_cov"
    label 'process_med'
    container 'docker://docker.io/nicmoya/mosdepth_hdr_image:2026_07_24'
    beforeScript = 'module load singularity'

    input:
    tuple val(group), val(strain), path(bamdir), path(windows_bed)

    output:
    tuple val(group), val(strain), path("${strain}.thresholds.bed"),
        emit: thresholds_bed

    script:
    """
    mosdepth \
        -b ${windows_bed} \
        -t 4 \
        -T 1,2,5 \
        -n \
        ${strain} \
        ${bamdir}/${strain}.bam

    gunzip ${strain}.thresholds.bed.gz

    awk -v g="${group}" -v s="${strain}" \
        'BEGIN{OFS="\t"} !/^#/ {print g, s, \$0}' \
        ${strain}.thresholds.bed \
        > ${strain}.thresholds.tmp

    mv ${strain}.thresholds.tmp ${strain}.thresholds.bed
    """
}

process MERGE_VARIANT_COUNTS {
    tag "merge_var"
    label 'process_low'
    publishDir "${params.output}/variants", mode: 'copy'
    container 'docker://docker.io/nicmoya/bedvcf_hdr_image:2026_07_24'
    beforeScript = 'module load singularity'

    input:
    path variant_count_files

    output:
    path "all_variant_counts.tsv", emit: merged_counts

    script:
    """
    cat ${variant_count_files} > all_variant_counts.tsv
    """
}

process MERGE_THRESHOLDS {
    tag "merge_thresh"
    label 'process_low'
    publishDir "${params.output}/coverage", mode: 'copy'
    container 'docker://docker.io/nicmoya/bedvcf_hdr_image:2026_07_24'
    beforeScript = 'module load singularity'

    input:
    path threshold_files

    output:
    path "all_thresholds.tsv", emit: merged_thresholds

    script:
    """
    cat ${threshold_files} > all_thresholds.tsv
    """
}

process CALL_HDRS {
    tag "call_hdrs_${group}"
    label 'process_med'
    publishDir "${params.output}/hdrs", mode: 'copy'
    container 'docker://docker.io/nicmoya/hdr_r_image:2026_07_24'
    beforeScript = 'module load singularity'

    input:
    tuple val(group),
          val(refstrain),
          path(windows_file),
          path(coverage_df),
          path(varct_df)

    val covthresh
    val vcthresh

    output:
    tuple val(group), path("${group}.hdrs.tsv"),
        emit: hdrs

    script:
    """
    export OMP_NUM_THREADS=${task.cpus}

    call_hdr.R \
        ${coverage_df} \
        ${varct_df} \
        ${windows_file} \
        ${refstrain} \
        ${group} \
        ${covthresh} \
        ${vcthresh}

    mv hdrs.tsv ${group}.hdrs.tsv
    """
}