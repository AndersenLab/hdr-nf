#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: generate_template.R <gtcheck.tsv>")
}

gtcheck_file <- args[1]

gt <- readr::read_tsv(
  gtcheck_file,
  show_col_types = FALSE
)

excluded_strains <- c(
  "MY681",
  "ECA1146",
  "JU356",
  "ECA1503"
)

gt <- gt %>%
  dplyr::filter(
    !i %in% excluded_strains,
    !j %in% excluded_strains
  )

readr::write_tsv(
  gt,
  "gtcheck.filtered.tsv"
)

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

matches <- gt %>%
  dplyr::filter(
    (i %in% refs & !j %in% refs) |
      (j %in% refs & !i %in% refs)
  ) %>%
  dplyr::mutate(
    strain = dplyr::if_else(i %in% refs, j, i),
    refstrain = dplyr::if_else(i %in% refs, i, j),
    concordance = (sites - discordance) / sites,
    group = unname(ref_groups[refstrain])
  ) %>%
  dplyr::filter(concordance > 0.95) %>%
  dplyr::group_by(strain) %>%
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
    ref = unname(ref_paths[refstrain]),
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
      unique(template$refstrain[is.na(template$ref)]),
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
  dplyr::group_walk(~ {
    strains <- sort(unique(.x$strain))
    
    readr::write_lines(
      strains,
      paste0(
        .y$group,
        "_",
        .y$refstrain,
        "_alignment_sample_sheet.txt"
      )
    )
  })