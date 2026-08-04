```text
╔═══════════════════════════════════════════════╗
║                                               ║
║  ██╗  ██╗██████╗ ██████╗      ███╗   ██╗███████╗
║  ██║  ██║██╔══██╗██╔══██╗     ████╗  ██║██╔════╝
║  ███████║██║  ██║██████╔╝     ██╔██╗ ██║█████╗
║  ██╔══██║██║  ██║██╔══██╗     ██║╚██╗██║██╔══╝
║  ██║  ██║██████╔╝██║  ██║     ██║ ╚████║██║
║  ╚═╝  ╚═╝╚═════╝ ╚═╝  ╚═╝     ╚═╝  ╚═══╝╚═╝
║                                               ║
║    Hyper Divergent Region Nextflow Pipeline   ║
║                                               ║
╚═══════════════════════════════════════════════╝

A reproducible Nextflow pipeline to call HDRs in selfing nematodes.
```

## Pipeline Parameters

Most input parameters are optional. If they are not provided, the pipeline automatically selects the appropriate defaults based on the value of `--species`.

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--species <string>` | Target species. | NULL. Accepted values: `c_elegans`, `c_briggsae`, `c_tropicalis`. |
| `--vcf <full_path>` | Reference VCF file used for extracting variant information. | Inferred from `--species`. |
| `--refgen <full_path>` | Reference genome (FASTA). | Inferred from `--species`. |
| `--bam <full_path>` | Directory containing input BAM files. | Inferred from `--species`. |
| `--str <string>` | Reference strain identifier. | Inferred from `--species`. |
| `--pbt <float>` | Percent bases covered threshold for calling HDRs. | Inferred from `--species`. |
| `--vct <int>` | Variant count threshold for calling HDRs. | Inferred from `--species`. |
| `--output <string>` | Name of the output directory. | `results_MM-DD-YY` |
| `-profile <string>` | Nextflow configuration profile. | NULL. Only `-profile rockfish` is accepted currently. |

> **Note:** Unless explicitly specified, the values for `--vcf`, `--reference`, `--bam`, `--ref`, `--pbt`, and `--vct` are automatically determined from the selected `--species`. These parameters are determined from the most recent releases as of 07/30/26. Manually define each custom parameter (except for HDR calling thresholds) when using a newer release.

## Usage

> **Disclaimer**
>
> This pipeline is designed to run within the `nf24_env` environment. Ensure this environment is activated before executing any of the commands below.

### Standard run

Run the pipeline using the default reference files for the selected species.

```bash
nextflow run main.nf --species c_elegans -profile rockfish
```

This is the recommended command for most users. The pipeline automatically uses the most recent released VCF and BAM files available as of **07/30/26** for the specified species.

### Using an updated VCF

Run the pipeline with a custom VCF while using the default reference genome and BAM directory.

```bash
nextflow run main.nf --species c_elegans --vcf <path_to_vcf_file> -profile rockfish
```

When a custom VCF is provided, the pipeline assumes that any newly introduced samples have corresponding BAM files located in:

```text
/vast/<allocation>/data/<species>/WI/alignments
```

### Custom HDR thresholds

Run the pipeline with user-defined HDR calling thresholds.

```bash
nextflow run main.nf \
    --species c_elegans \
    --vct 15 \
    --pbt 0.7 \
    --output CE_VCT15_PBT70 \
    -profile rockfish
```

This command specifies custom variant count (`--vct`) and percent bases covered (`--pbt`) thresholds, as well as a custom output directory. Modifying these thresholds is intended only for experimental analyses and is **strongly discouraged unless there is a well-justified reason**.

### Development workflows

The following parameters are reserved for alternative workflows that are currently under development for **_C. briggsae_**:

- `--bam`
- `--refgen`
- `--str`

These parameters are not intended for routine use.

---

## Current Release Files (Revised: 08/04/26)

The pipeline automatically selects the appropriate reference files based on the value of `--species`.

| Species | Reference Genome | Reference VCF |
|---------|------------------|---------------|
| `c_elegans` | `c_elegans.PRJNA13758.WS283.genome.fa` | `WI.20250625.hard-filter.isotype.vcf.gz` |
| `c_briggsae` | `c_briggsae.QX1410_nanopore.Feb2020.genome.fa` | `WI.20250626.hard-filter.isotype.vcf.gz` |
| `c_tropicalis` | `c_tropicalis.NIC58_nanopore.June2021.genome.fa` | `WI.20250627.hard-filter.isotype.vcf.gz` |
