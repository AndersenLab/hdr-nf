# hdr-nf
A reproducible nextflow to call HDRs in selfing nematodes.

## Pipeline Parameters

Most input parameters are optional. If they are not provided, the pipeline automatically selects the appropriate defaults based on the value of `--species`.

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--species` | Target species selected by user | Options are "c_elegans", "c_briggsae", or "c_tropicalis" |
| `--vcf` | Reference VCF file used for extracting variant information. | Inferred from `--species`. |
| `--reference` | Reference genome (FASTA). | Inferred from `--species`. |
| `--bam` | Directory containing input BAM files. | Inferred from `--species`. |
| `--ref` | Reference strain identifier. | Inferred from `--species`. |
| `--pbt` | Percent bases covered threshold for calling HDRs. | Inferred from `--species`. |
| `--vct` | Variant count threshold for calling HDRs. | Inferred from `--species`. |
| `--output` | Output directory for all pipeline results. | Uses the value provided by the user. Defaults to  "results". |
| `-profile` | Select the nextflow configuration profile. | Only "rockfish" is configured currently. |

> **Note:** Unless explicitly specified, the values for `--vcf`, `--reference`, `--bam`, `--ref`, `--pbt`, and `--vct` are automatically determined from the selected `--species`. These parameters are determined from the most recent releases as of 07/30/26. Manually define each custom parameter when using a newer release.

