#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(tibble)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: generate_template.R <gtcheck.tsv> <isogroups.tsv>")
}

gtcheck_file <- args[1]
isogroups_file <- args[2]

gt <- readr::read_tsv(
  gtcheck_file,
  show_col_types = FALSE
)

isogroups <- readr::read_tsv(
  isogroups_file,
  show_col_types = FALSE
)

required_isogroup_cols <- c(
  "group",
  "strain",
  "isotype",
  "isotype_ref_strain"
)

if (!all(required_isogroup_cols %in% colnames(isogroups))) {
  stop(
    "Isogroups file is missing required column(s): ",
    paste(
      setdiff(required_isogroup_cols, colnames(isogroups)),
      collapse = ", "
    )
  )
}

strains_of_interest <- unique(
  stats::na.omit(isogroups$isotype_ref_strain)
)

gt <- gt %>%
  dplyr::filter(
    i %in% strains_of_interest,
    j %in% strains_of_interest
  )

if (nrow(gt) == 0) {
  stop(
    "No GTCHECK comparisons remain after filtering i and j ",
    "to strains in isotype_ref_strain"
  )
}

ref_groups <- c(
  "BRC20492" = "TD2",
  "BRC20530" = "TD1",
  "ECA2670"  = "AD",
  "JU1348"   = "KD",
  "JU2536"   = "Temperate",
  "NIC1660"  = "TH",
  "QX1410"   = "Tropical"
)

refs <- names(ref_groups)

missing_refs <- setdiff(
  refs,
  strains_of_interest
)

if (length(missing_refs) > 0) {
  stop(
    "Required relatedness-group reference strain(s) are absent from ",
    "isotype_ref_strain: ",
    paste(missing_refs, collapse = ", ")
  )
}

matches <- gt %>%
  dplyr::filter(
    (i %in% refs & !j %in% refs) |
      (j %in% refs & !i %in% refs)
  ) %>%
  dplyr::mutate(
    strain = dplyr::if_else(
      i %in% refs,
      j,
      i
    ),
    refstrain = dplyr::if_else(
      i %in% refs,
      i,
      j
    ),
    concordance = (sites - discordance) / sites,
    group = unname(ref_groups[refstrain])
  ) %>%
  dplyr::filter(
    concordance > 0.95
  ) %>%
  dplyr::group_by(
    strain
  ) %>%
  dplyr::slice_max(
    order_by = concordance,
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    strain,
    refstrain,
    group
  )

ref_matches <- tibble::tibble(
  strain = names(ref_groups),
  refstrain = names(ref_groups),
  group = unname(ref_groups)
)

matches_bound <- dplyr::bind_rows(
  matches,
  ref_matches
) %>%
  dplyr::arrange(
    group,
    strain
  )

ref_paths <- c(
  "BRC20492" = "/vast/eande106/data/c_briggsae/WI/relatedness_group/20260324/references/BRC20492.ragtag.fa",
  "BRC20530" = "/vast/eande106/data/c_briggsae/WI/relatedness_group/20260324/references/BRC20530.ragtag.fa",
  "ECA2670"  = "/vast/eande106/data/c_briggsae/WI/relatedness_group/20260324/references/ECA2670.ragtag.fa",
  "JU1348"   = "/vast/eande106/data/c_briggsae/WI/relatedness_group/20260324/references/JU1348.ragtag.fa",
  "JU2536"   = "/vast/eande106/data/c_briggsae/WI/relatedness_group/20260324/references/JU2536.ragtag.fa",
  "NIC1660"  = "/vast/eande106/data/c_briggsae/WI/relatedness_group/20260324/references/NIC1660.ragtag.fa",
  "QX1410"   = "/vast/eande106/data/c_briggsae/genomes/QX1410_nanopore/Feb2020/c_briggsae.QX1410_nanopore.Feb2020.genome"
)

template <- matches_bound %>%
  dplyr::mutate(
    vcf = "[REQUIRES_WI_GATK_NF]",
    ref = unname(
      ref_paths[refstrain]
    ),
    bam_path = "[REQUIRES_ALIGNMENT_NF]"
  ) %>%
  dplyr::select(
    group,
    strain,
    vcf,
    ref,
    bam_path,
    refstrain
  )

if (any(is.na(template$ref))) {
  stop(
    "Missing reference path for refstrain(s): ",
    paste(
      unique(
        template$refstrain[
          is.na(template$ref)
        ]
      ),
      collapse = ", "
    )
  )
}

readr::write_csv(
  template,
  "c_briggsae.samplesheet.template.csv"
)

template %>%
  dplyr::group_by(
    group,
    refstrain
  ) %>%
  dplyr::group_walk(
    ~ {
      strains <- sort(
        unique(.x$strain)
      )

      readr::write_lines(
        strains,
        paste0(
          .y$group,
          "_",
          .y$refstrain,
          "_alignment_sample_sheet.txt"
        )
      )
    }
  )