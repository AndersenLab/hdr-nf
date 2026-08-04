# hdr-nf
A reproducible nextflow to call HDRs in selfing nematodes.

## Pipeline Parameters

Most input parameters are optional. If they are not provided, the pipeline automatically selects the appropriate defaults based on the value of `--species`.

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--species` | Target species selected by user | Default is null. Options are "c_elegans", "c_briggsae", or "c_tropicalis" |
| `--vcf` | Reference VCF file used for extracting variant information. | Default is null.Inferred from `--species`. |
| `--refgen` | Reference genome (FASTA). | Default is null. Inferred from `--species`. |
| `--bam` | Directory containing input BAM files. | Default is null. Inferred from `--species`. |
| `--str` | Reference strain identifier. | Default is null. Inferred from `--species`. |
| `--pbt` | Percent bases covered threshold for calling HDRs. | Default is null. Inferred from `--species`. |
| `--vct` | Variant count threshold for calling HDRs. | Default is null. Inferred from `--species`. |
| `--output` | Output directory for all pipeline results. | Default is "results". Uses the value provided by the user. |
| `-profile` | Select the nextflow configuration profile. | Only "rockfish" is configured currently. |

> **Note:** Unless explicitly specified, the values for `--vcf`, `--reference`, `--bam`, `--ref`, `--pbt`, and `--vct` are automatically determined from the selected `--species`. These parameters are determined from the most recent releases as of 07/30/26. Manually define each custom parameter (except for HDR calling thresholds) when using a newer release.

