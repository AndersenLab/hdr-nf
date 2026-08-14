
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
} else if (params.species == 'c_briggsae') {

    if (params.samplesheet) {
        println """
        Running hdr-nf in relatedness group mode.
        Using the provided sample sheet:
          ${params.samplesheet}

        Coordinate transformation will be performed to the Tropical reference.
        """
    } else {
        println """
        Running hdr-nf in relatedness group (RG) mode.
        No RG sample sheet was provided.
        A RG sample sheet template will be along with relatedness group lists necessary
        to run alignment-nf and wi-gatk (required to fill the sample sheet).
        Coordinate transformation will be performed to the Tropical reference.
        """
    }
} else {
    println """
    Running hdr-nf in standard mode (without relatedness groups).
    Using the standard VCF, reference genome, and BAM inputs.
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

def rpo_thresh = [
    c_elegans   : '25',
    c_tropicalis: '25',
    c_briggsae  : '25'
]

def ref_gtcheck = [
    c_elegans   : '/vast/eande106/data/c_elegans/WI/concordance/20250625/gtcheck.txt',
    c_tropicalis: '/vast/eande106/data/c_tropicalis/WI/concordance/20250627/gtcheck.txt',
    c_briggsae  : '/vast/eande106/data/c_briggsae/WI/concordance/20250626/gtcheck.txt'
]

def ref_groups = [
    c_elegans   : '/vast/eande106/data/c_elegans/WI/concordance/20250625/isotype_groups.tsv',
    c_tropicalis: '/vast/eande106/data/c_tropicalis/WI/concordance/20250627/isotype_groups.tsv',
    c_briggsae  : '/vast/eande106/data/c_briggsae/WI/concordance/20250626/isotype_groups.tsv'
]

def invcf = params.vcf ?: ref_vcf[params.species]
def ingenome = params.refgen ?: ref_genome[params.species]
def inbam = params.bam ?: bam_dir[params.species]
def refstrain = params.str ?: ref_str[params.species]
def covthresh = params.pbt ?: cov_thresh[params.species]
def vcthresh = params.vct ?: var_thresh[params.species]
def rpothresh = params.rpo ?: rpo_thresh[params.species]
def ingtcheck = params.gtcheck ?: ref_gtcheck[params.species]
def ingroups = params.groups ?:ref_groups[params.species]

def paramSummary = [
    'Species'            : params.species,
    'VCF'                : invcf,
    'Reference genome'   : ingenome,
    'BAM directory'      : inbam,
    'REF strain'         : refstrain,
    'Coverage threshold' : covthresh,
    'Variant threshold'  : vcthresh,
    'Output directory'   : params.output
]

if (params.species == 'c_briggsae') {

    if (params.samplesheet) {
        paramSummary['Sample Sheet'] = params.samplesheet
    } else {
        paramSummary['Sample Sheet']    = params.samplesheet
        paramSummary['GTCHECK']         = params.gtcheck
        paramSummary['Isotype groups']  = params.groups
    }
}

if (covthresh == null || vcthresh == null) {
    error "Thresholds not defined for species: ${params.species}"
}

def maxLen = paramSummary.keySet().collect { k -> k.size() }.max()

def summary = paramSummary.collect { k, v ->
    "${k.padRight(maxLen)} : ${v}"
}.join('\n')

log.info """
=========================================
  Pipeline Parameters
=========================================
${summary}
=========================================
""".stripIndent()

workflow {

    if (params.species == 'c_briggsae' && !params.samplesheet) {
        if (!params.gtcheck || !params.groups) {

            def default_inputs = []

            if (!params.gtcheck) {
                default_inputs << [
                    "  GTCHECK:",
                    "    ${ingtcheck}"
                ].join('\n')
            }

            if (!params.groups) {
                default_inputs << [
                    "  Isotype groups:",
                    "    ${ingroups}"
                ].join('\n')
            }

            def warning_msg = [
                "Sample sheet is missing.",
                "One or more isotype input files were not provided with --gtcheck/--groups.",
                "The following default file(s) for ${params.species} will be used:",
                "",
                default_inputs.join('\n\n'),
                "",
                "Please ensure these files are appropriate for your dataset.",
                "Mismatches between the isotype groups, GTCHECK, and VCF samples may",
                "generate an incomplete/truncated sample sheet template."
            ].join('\n')

            log.warn warning_msg
}
        
        ch_gtcheck = channel.fromPath(
            ingtcheck,
            checkIfExists: true
        )

        ch_isogroups = channel.fromPath(
            ingroups,
            checkIfExists: true
        )

        GENERATE_TEMPLATE(
            ch_gtcheck,
            ch_isogroups
        )

    } else {

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
                .combine(GENERATE_WINDOWS_GROUP.out.windows_bed, by: 0)

            ch_variant_inputs = ch_samples_with_windows
                .map { group, strain, vcf, bam_path, group_refstrain, windows ->
                    tuple(group, strain, vcf, windows)
                }

            ch_coverage_inputs = ch_samples_with_windows
                .map { group, strain, vcf, bam_path, group_refstrain, windows ->
                    tuple(group, strain, bam_path, windows)
                }

            ch_tropical_ref = ch_group_refs
                .filter { group, ref, group_refstrain ->
                    group == 'Tropical'
                }
                .map { group, ref, group_refstrain ->
                    ref
                }
                .first()

            ch_nucmer_inputs = ch_group_refs
                .filter { group, ref, group_refstrain ->
                    group != 'Tropical'
                }
                .map { group, ref, group_refstrain ->
                    tuple(group, ref)
                }
                .combine(ch_tropical_ref)

            ALIGN_TO_TROPICAL(ch_nucmer_inputs)

            ch_all_group_alignments = ALIGN_TO_TROPICAL.out.coords
                .map { group, coords ->
                    coords
                }
                .collect()

            MERGE_GROUP_ALIGNMENTS(ch_all_group_alignments)

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

        if (params.species == 'c_briggsae') {

            ch_hdrs_to_transform = CALL_HDRS.out.hdrs
                .filter { group, hdrs ->
                    group != 'Tropical'
                }

            TRANSFORM_HDRS(
                ch_hdrs_to_transform,
                MERGE_GROUP_ALIGNMENTS.out.merged_coords
            )

            ch_transformed_hdrs = TRANSFORM_HDRS.out.transformed_hdrs
                .map { group, hdr -> hdr }
                .collect()

            ch_tropical_hdr = CALL_HDRS.out.hdrs
                .filter { group, hdr ->
                    group == 'Tropical'
                }
                .map { group, hdr ->
                    hdr
                }
                .first()

            MERGE_TRANSFORMED_HDRS(
                ch_transformed_hdrs,
                ch_tropical_hdr
            )
        }
    }
}

process GENERATE_SAMPLE_LIST_AND_WINDOWS {
    tag "${vcf_file.simpleName}"
    label 'process_low'
    container 'library://docker.io/nicmoya/bedvcf_hdr_image:2026_07_24'
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
    container 'library://docker.io/nicmoya/bedvcf_hdr_image:2026_07_24'
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
    container 'library://docker.io/nicmoya/bedvcf_hdr_image:2026_07_24'
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
    container 'library://docker.io/nicmoya/mosdepth_hdr_image:2026_07_24'
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
    container 'library://docker.io/nicmoya/bedvcf_hdr_image:2026_07_24'
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
    container 'library://docker.io/nicmoya/bedvcf_hdr_image:2026_07_24'
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
    label 'process_hi'
    publishDir(
        params.species == 'c_briggsae'
            ? "${params.output}/hdr_untransformed"
            : "${params.output}/hdrs",
        mode: 'copy'
    )
    container 'library://docker.io/nicmoya/hdr_r_image:2026_07_24'
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
        ${vcthresh} \
        ${rpothresh}

    mv hdrs.tsv ${group}.hdrs.tsv
    """
}

process ALIGN_TO_TROPICAL {
    tag "${group}_vs_Tropical"
    label 'process_med'
    container 'library://docker.io/nicmoya/nucmer_hdr_image:2026_08_12'
    //publishDir "${params.output}/alignments_raw", mode: 'copy'
    beforeScript = 'module load singularity'

    input:
    tuple val(group),
          path(query_genome),
          path(tropical_genome)

    output:
    tuple val(group),
          path("${group}_transformed.tsv"),
          emit: coords

    script:
    """
    nucmer \
        --maxgap=500 \
        --mincluster=100 \
        --prefix=${group} \
        --coords \
        ${tropical_genome} \
        ${query_genome}

    show-coords \
        -r \
        -l \
        -T \
        ${group}.delta | \
        awk -v g="${group}" \
            'BEGIN{OFS="\\t"} \$5 > 1000 {print g, \$0}' \
        > ${group}_transformed.tsv
    """
}

process MERGE_GROUP_ALIGNMENTS {
    tag "merge_group_alignments"
    label 'process_low'
    publishDir "${params.output}/genome_alignments", mode: 'copy'

    input:
    path alignment_files

    output:
    path "all_groups.vs_Tropical.coords.tsv",
        emit: merged_coords

    script:
    """
    cat ${alignment_files} | grep -v "\\[S1\\]" > all_groups.vs_Tropical.coords.tsv
    """
}

process TRANSFORM_HDRS {
    tag "transform_coords_${group}"
    label 'process_med'
    //publishDir "${params.output}/transformed_coords", mode: 'copy'
    container 'library://docker.io/nicmoya/hdr_r_image:2026_07_24'
    beforeScript = 'module load singularity'

    input:
    tuple val(group), path(hdr_file)
    path merged_coords

    output:
    tuple val(group), path("${group}.transformed_hdrs.tsv"),
        emit: transformed_hdrs

    script:
    """
    transform_hdr.R \
        ${hdr_file} \
        ${merged_coords} \
        ${group} \
        ${group}.transformed_hdrs.tsv
    """
}

process MERGE_TRANSFORMED_HDRS {
    tag "merge_transformed_hdrs"
    label 'process_low'
    publishDir "${params.output}/hdrs_transformed", mode: 'copy'

    input:
    path transformed_hdrs
    path tropical_hdr

    output:
    path "all_groups.transformed_hdrs.tsv",
        emit: merged_hdrs

    script:
    """
    awk 'BEGIN{FS=OFS="\\t"}
        NR==1 {print "CHROM","start","end","group","STRAIN","size","source"; next}
        {print \$1,\$2,\$3,\$4,\$5,\$6,"${refstrain}"}' \
        ${tropical_hdr} > tropical.normalized.tsv

    head -n 1 tropical.normalized.tsv > all_groups.unsorted.tsv

    tail -n +2 tropical.normalized.tsv >> all_groups.unsorted.tsv

    for f in ${transformed_hdrs}; do
        tail -n +2 "\$f" >> all_groups.unsorted.tsv
    done

    head -n 1 all_groups.unsorted.tsv > all_groups.transformed_hdrs.tsv

    tail -n +2 all_groups.unsorted.tsv | \
    awk 'BEGIN{FS=OFS="\\t"}
        \$1=="I"   {chr=1}
        \$1=="II"  {chr=2}
        \$1=="III" {chr=3}
        \$1=="IV"  {chr=4}
        \$1=="V"   {chr=5}
        \$1=="X"   {chr=6}
        {print \$4,\$5,chr,\$2,\$0}' | \
    sort -t \$'\\t' -k1,1 -k2,2 -k3,3n -k4,4n | \
    cut -f5- >> all_groups.transformed_hdrs.tsv
    """
}

process GENERATE_TEMPLATE {
    tag "generate_template"
    label 'process_low'
    publishDir "${params.output}/template_sample_sheet", mode: 'copy'
    container 'library://docker.io/nicmoya/hdr_r_image:2026_07_24'
    beforeScript = 'module load singularity'

    input:
    path gtcheck
    path isogroups

    output:
    path "c_briggsae.samplesheet.template.csv",
        emit: samplesheet_template

    path "*_alignment_sample_sheet.txt",
        emit: alignment_sample_sheets

    script:
    """
    generate_template.R \
        ${gtcheck} \
        ${isogroups}
    """
}
