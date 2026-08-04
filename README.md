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
║      Hyper Divergent Region Nextflow Pipeline ║
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

