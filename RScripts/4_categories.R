## ============================================================
## 4_categories.R — UpSet + categories (+ gate + shifted-baseline)
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexUpset)
})

## ---------------- Parameters ----------------
analysis_level <- "gene"  # "gene" or "isoform"

## ---------------- Helper functions ----------------
strip_version <- function(x) sub("\\..*$", "", x)

load_id_table <- function(file, analysis_level) {
  if (!file.exists(file)) return(data.frame())
  df <- read.delim(file, check.names = FALSE)
  
  if (analysis_level == "gene") {
    if (!"gene_id" %in% names(df)) return(data.frame())
    df <- df %>% dplyr::select(gene_id) %>% distinct()
    df$gene_id <- strip_version(df$gene_id)
  } else if (analysis_level == "isoform") {
    if (!all(c("isoform_id", "gene_id") %in% names(df))) return(data.frame())
    df <- df %>% dplyr::select(isoform_id, gene_id) %>% distinct()
    df$isoform_id <- strip_version(df$isoform_id)
    df$gene_id    <- strip_version(df$gene_id)
  }
  df
}

## ---------------- Directories ----------------
base_root <- getwd()

de_tbl   <- file.path(base_root, "results/1_DEAnalysis/tables/shrunk")
int_tbl  <- file.path(base_root, "results/2_interaction/tables")
out_tbl  <- file.path(base_root, "results/4_categories/tables")
out_fig  <- file.path(base_root, "results/4_categories/figures")
dir.create(out_tbl, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)

## ---------------- Annotation map: gene_id -> SYMBOL ----------------
annotation <- read.delim(file.path(base_root, "data", "annotation.txt"), check.names = FALSE)
sym_col <- if ("symbol" %in% names(annotation)) "symbol" else if ("gene_name" %in% names(annotation)) "gene_name" else NA
if (is.na(sym_col)) stop("annotation.txt must have symbol or gene_name column.")

gene2sym <- annotation %>%
  dplyr::select(gene_id, SYMBOL = all_of(sym_col)) %>%
  dplyr::mutate(gene_id = strip_version(gene_id)) %>%
  dplyr::distinct(gene_id, .keep_all = TRUE)

## ---------------- Drugs ----------------
drugs <- c("11j", "KVS", "11j_PlaB", "KVS_PlaB", "PlaB")

## ============================================================
## PART 1 — UpSet plots
## ============================================================

for (drug in drugs) {
  message(">>> [UpSet] Processing: ", drug)
  
  up_drug       <- load_id_table(file.path(de_tbl,  paste0("up_",   drug, ".tsv")), analysis_level)
  up_drug_OHT   <- load_id_table(file.path(de_tbl,  paste0("up_",   drug, "_OHT.tsv")), analysis_level)
  down_drug     <- load_id_table(file.path(de_tbl,  paste0("down_", drug, ".tsv")), analysis_level)
  down_drug_OHT <- load_id_table(file.path(de_tbl,  paste0("down_", drug, "_OHT.tsv")), analysis_level)
  up_int        <- load_id_table(file.path(int_tbl, paste0("up_int_",   drug, ".tsv")), analysis_level)
  down_int      <- load_id_table(file.path(int_tbl, paste0("down_int_", drug, ".tsv")), analysis_level)
  
  ## gate sets (drug_OHT vs drug) in DE (shrunk)
  gate_up       <- load_id_table(file.path(de_tbl, paste0("up_",   drug, "_OHT__", drug, ".tsv")), analysis_level)
  gate_down     <- load_id_table(file.path(de_tbl, paste0("down_", drug, "_OHT__", drug, ".tsv")), analysis_level)
  
  if (analysis_level == "gene") {
    all_ids <- unique(na.omit(c(
      up_drug$gene_id, up_drug_OHT$gene_id, down_drug$gene_id, down_drug_OHT$gene_id,
      up_int$gene_id, down_int$gene_id,
      gate_up$gene_id, gate_down$gene_id
    )))
    if (length(all_ids) == 0) {
      message("  [skip] No genes found for ", drug)
      next
    }
    
    upset_df <- tibble(gene_id = all_ids)
    ref_col <- upset_df$gene_id
    
    upset_df[[paste0("up_", drug)]]                 <- ref_col %in% up_drug$gene_id
    upset_df[[paste0("up_", drug, "_OHT")]]         <- ref_col %in% up_drug_OHT$gene_id
    upset_df[[paste0("down_", drug)]]               <- ref_col %in% down_drug$gene_id
    upset_df[[paste0("down_", drug, "_OHT")]]       <- ref_col %in% down_drug_OHT$gene_id
    upset_df[[paste0("up_int_", drug)]]             <- ref_col %in% up_int$gene_id
    upset_df[[paste0("down_int_", drug)]]           <- ref_col %in% down_int$gene_id
    
    ## explicit gate membership columns (kept for Part 2)
    upset_df[[paste0("up_", drug, "_OHT__", drug)]]   <- ref_col %in% gate_up$gene_id
    upset_df[[paste0("down_", drug, "_OHT__", drug)]] <- ref_col %in% gate_down$gene_id
    
  } else {
    all_iso <- unique(na.omit(c(
      up_drug$isoform_id, up_drug_OHT$isoform_id, down_drug$isoform_id, down_drug_OHT$isoform_id,
      up_int$isoform_id, down_int$isoform_id,
      gate_up$isoform_id, gate_down$isoform_id
    )))
    if (length(all_iso) == 0) {
      message("  [skip] No isoforms found for ", drug)
      next
    }
    
    map_df <- bind_rows(up_drug, up_drug_OHT, down_drug, down_drug_OHT, up_int, down_int, gate_up, gate_down) %>%
      distinct(isoform_id, gene_id)
    
    upset_df <- tibble(isoform_id = all_iso) %>%
      left_join(map_df, by = "isoform_id")
    
    ref_col <- upset_df$isoform_id
    
    upset_df[[paste0("up_", drug)]]                 <- ref_col %in% up_drug$isoform_id
    upset_df[[paste0("up_", drug, "_OHT")]]         <- ref_col %in% up_drug_OHT$isoform_id
    upset_df[[paste0("down_", drug)]]               <- ref_col %in% down_drug$isoform_id
    upset_df[[paste0("down_", drug, "_OHT")]]       <- ref_col %in% down_drug_OHT$isoform_id
    upset_df[[paste0("up_int_", drug)]]             <- ref_col %in% up_int$isoform_id
    upset_df[[paste0("down_int_", drug)]]           <- ref_col %in% down_int$isoform_id
    
    ## explicit gate membership columns (kept for Part 2)
    upset_df[[paste0("up_", drug, "_OHT__", drug)]]   <- ref_col %in% gate_up$isoform_id
    upset_df[[paste0("down_", drug, "_OHT__", drug)]] <- ref_col %in% gate_down$isoform_id
  }
  
  write.table(
    upset_df,
    file.path(out_tbl, paste0("upset_data_", drug, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  sets <- setdiff(names(upset_df), c("gene_id", "isoform_id"))
  
  ## hide gate columns from UpSet plot (still present in upset_df for Part 2)
  gate_cols <- c(
    paste0("up_",   drug, "_OHT__", drug),
    paste0("down_", drug, "_OHT__", drug)
  )
  sets <- setdiff(sets, gate_cols)
  
  p <- upset(
    upset_df,
    intersect = sets,
    base_annotations = list(
      "Intersection size" = intersection_size(counts = TRUE, text = list(size = 3.5))
    ),
    set_sizes = upset_set_size(),
    width_ratio = 0.3,
    height_ratio = 0.3
  ) + theme_bw(base_size = 9)
  
  ggsave(file.path(out_fig, paste0("UpSet_", drug, ".png")), p, width = 9, height = 6, dpi = 300)
}

message("\n>>> Part 1 complete: UpSet plots and tables saved.\n")

## ============================================================
## PART 2 — Category detection + save IDs + SYMBOL
## ============================================================

for (drug in drugs) {
  message(">>> [Categories] Processing: ", drug)
  
  fp <- file.path(out_tbl, paste0("upset_data_", drug, ".tsv"))
  if (!file.exists(fp)) {
    message("  [skip] No UpSet data for ", drug)
    next
  }
  df <- read.delim(fp, check.names = FALSE)
  
  up_noOHT   <- paste0("up_", drug)
  up_OHT     <- paste0("up_", drug, "_OHT")
  down_noOHT <- paste0("down_", drug)
  down_OHT   <- paste0("down_", drug, "_OHT")
  up_int     <- paste0("up_int_", drug)
  down_int   <- paste0("down_int_", drug)
  
  ## gate columns
  up_gate    <- paste0("up_",   drug, "_OHT__", drug)
  down_gate  <- paste0("down_", drug, "_OHT__", drug)
  
  has_col <- function(x) if (x %in% names(df)) df[[x]] else rep(FALSE, nrow(df))
  U1 <- has_col(up_noOHT); U2 <- has_col(up_OHT)
  D1 <- has_col(down_noOHT); D2 <- has_col(down_OHT)
  UI <- has_col(up_int); DI <- has_col(down_int)
  UG <- has_col(up_gate)
  DG <- has_col(down_gate)
  
  category <- rep(NA_character_, nrow(df))
  
  ## enhanced
  category[(U1 & U2 & UI) | (!U1 & U2 & UI)] <- "enhanced_up"
  category[(D1 & D2 & DI) | (!D1 & D2 & DI)] <- "enhanced_down"
  
  ## suppressed
  category[(U1 & U2 & DI) | (U1 & !U2 & DI)] <- "suppressed_up"
  category[(D1 & D2 & UI) | (D1 & !D2 & UI)] <- "suppressed_down"
  
  ## switched
  category[(D1 & U2 & UI)] <- "switched_positive"
  category[(U1 & D2 & DI)] <- "switched_negative"
  
  ## independent (true independent = no baseline shift)
  
  is_unassigned <- is.na(category)
  
  ind_up  <- is_unassigned & (U1 & U2)
  ind_down <- is_unassigned & (D1 & D2)
  
  ## true independent: no significant change in drug+OHT/drug
  category[ind_up  & !UG] <- "independent_up"
  category[ind_down & !DG] <- "independent_down"
  
  ## shifted baseline independent
  category[ind_up  & UG] <- "shifted_baseline_independent_up"
  category[ind_down & DG] <- "shifted_baseline_independent_down"
  
  ## ---------------- Split by gate (FIXED) ----------------
  shifted_sub <- rep(NA_character_, nrow(df))
  
  shifted_sub[category == "enhanced_up"       & !UG] <- "shifted_baseline_enhanced_up"
  shifted_sub[category == "enhanced_down"     & !DG] <- "shifted_baseline_enhanced_down"
  shifted_sub[category == "suppressed_up"     & !UG] <- "shifted_baseline_suppressed_up"
  shifted_sub[category == "suppressed_down"   & !DG] <- "shifted_baseline_suppressed_down"
  
  category[!is.na(shifted_sub)] <- shifted_sub[!is.na(shifted_sub)]
  df$category <- category
  
  ## Drop unclassified
  df2 <- df[!is.na(df$category), , drop = FALSE]
  if (nrow(df2) == 0) next

  cat_list <- split(df2, df2$category)
  
  for (nm in names(cat_list)) {
    subdf <- cat_list[[nm]]
    if (is.null(subdf) || nrow(subdf) == 0) {
      message("  (no gene detected in ", nm, " category)")
      next
    }
    
    if (analysis_level == "gene") {
      out_df <- subdf %>%
        dplyr::select(gene_id) %>%
        dplyr::mutate(gene_id = strip_version(gene_id)) %>%
        dplyr::distinct() %>%
        dplyr::left_join(gene2sym, by = "gene_id")
      
      if (nrow(out_df) == 0) {
        message("  (no gene detected in ", nm, " category)")
        next
      }
      
    } else {
      out_df <- subdf %>%
        dplyr::select(isoform_id, gene_id) %>%
        dplyr::mutate(
          isoform_id = strip_version(isoform_id),
          gene_id    = strip_version(gene_id)
        ) %>%
        dplyr::distinct() %>%
        dplyr::left_join(gene2sym, by = "gene_id")
      
      if (nrow(out_df) == 0) {
        message("  (no gene detected in ", nm, " category)")
        next
      }
    }
    
    write.table(
      out_df,
      file.path(out_tbl, paste0(nm, "_", drug, ".tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
  
  ## Also report categories that were not created at all for this drug
  expected_cats <- c(
    "enhanced_up","enhanced_down","suppressed_up","suppressed_down",
    "switched_positive","switched_negative","independent_up","independent_down",
    "shifted_baseline_enhanced_up","shifted_baseline_enhanced_down",
    "shifted_baseline_suppressed_up","shifted_baseline_suppressed_down"
  )
  missing_cats <- setdiff(expected_cats, names(cat_list))
  for (mc in missing_cats) {
    message("  (no gene detected in ", mc, " category)")
  }
}

cat("\n>>> Part 2 complete: categories saved under results/4_categories/tables\n")