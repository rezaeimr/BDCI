## ============================================================
## 4_categories.R — UpSet + categories + ADDITIVE (simple + correct)
## - Categories: enhanced / suppressed / switched / independent
## - Shifted-baseline: enhanced / suppressed / independent
## - Additive: subset of TRUE independent (up/down), using shrunk LFC difference
## - Saves: tables (.tsv) + figures (.pdf + .tiff)  [NO png]
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(ComplexUpset)
  library(ggplot2)
})

## ---------------- Parameters ----------------
analysis_level  <- "gene"   # "gene" or "isoform"

## ---------------- Helpers ----------------
strip_version <- function(x) sub("\\..*$", "", x)

save_both <- function(path_no_ext, plot, w = 9, h = 6) {
  ggsave(paste0(path_no_ext, ".pdf"), plot, width = w, height = h)
  ggsave(paste0(path_no_ext, ".tiff"), plot, width = w, height = h,
         dpi = 600, device = "tiff", compression = "lzw")
}

load_id_table <- function(file, analysis_level) {
  if (!file.exists(file)) return(data.frame())
  df <- read.delim(file, check.names = FALSE)
  
  if (analysis_level == "gene") {
    if (!"gene_id" %in% names(df)) return(data.frame())
    df %>%
      transmute(gene_id = strip_version(gene_id)) %>%
      distinct()
  } else {
    if (!all(c("isoform_id", "gene_id") %in% names(df))) return(data.frame())
    df %>%
      transmute(
        isoform_id = strip_version(isoform_id),
        gene_id    = strip_version(gene_id)
      ) %>%
      distinct()
  }
}

load_lfc_table <- function(file, analysis_level) {
  if (!file.exists(file)) return(data.frame())
  df <- read.delim(file, check.names = FALSE)
  if (!"log2FoldChange" %in% names(df)) return(data.frame())
  
  if (analysis_level == "gene") {
    if (!"gene_id" %in% names(df)) return(data.frame())
    df %>% transmute(gene_id = strip_version(gene_id), lfc = log2FoldChange) %>% distinct()
  } else {
    if (!all(c("isoform_id","gene_id") %in% names(df))) return(data.frame())
    df %>% transmute(
      isoform_id = strip_version(isoform_id),
      gene_id    = strip_version(gene_id),
      lfc        = log2FoldChange
    ) %>% distinct()
  }
}

pick_id_col <- function(df) {
  if ("isoform_id" %in% names(df)) return("isoform_id")
  if ("gene_id" %in% names(df)) return("gene_id")
  stop("No ID column found in UpSet table.")
}

## ---------------- Paths ----------------
base_root <- getwd()

de_tbl  <- file.path(base_root, "results/1_DEAnalysis/tables/shrunk")
int_tbl <- file.path(base_root, "results/2_interaction/tables")
out_tbl <- file.path(base_root, "results/4_categories/tables")
out_fig <- file.path(base_root, "results/4_categories/figures")
dir.create(out_tbl, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)

## ---------------- Annotation: gene_id -> SYMBOL ----------------
annotation <- read.delim(file.path(base_root, "data", "annotation.txt"), check.names = FALSE)
sym_col <- if ("symbol" %in% names(annotation)) "symbol" else if ("gene_name" %in% names(annotation)) "gene_name" else NA
if (is.na(sym_col)) stop("annotation.txt must have symbol or gene_name column.")

gene2sym <- annotation %>%
  dplyr::select(gene_id, SYMBOL = all_of(sym_col)) %>%
  mutate(gene_id = strip_version(gene_id)) %>%
  distinct(gene_id, .keep_all = TRUE)

## ---------------- Drugs ----------------
drugs <- c("11j", "KVS", "11j_PlaB", "KVS_PlaB", "PlaB")

## ============================================================
## PART 1 — Build upset_data_<drug>.tsv + UpSet figure
## ============================================================

for (drug in drugs) {
  message(">>> [UpSet] ", drug)
  
  up_no   <- load_id_table(file.path(de_tbl,  paste0("up_",   drug, ".tsv")), analysis_level)
  up_oht  <- load_id_table(file.path(de_tbl,  paste0("up_",   drug, "_OHT.tsv")), analysis_level)
  dn_no   <- load_id_table(file.path(de_tbl,  paste0("down_", drug, ".tsv")), analysis_level)
  dn_oht  <- load_id_table(file.path(de_tbl,  paste0("down_", drug, "_OHT.tsv")), analysis_level)
  
  up_int  <- load_id_table(file.path(int_tbl, paste0("up_int_",   drug, ".tsv")), analysis_level)
  dn_int  <- load_id_table(file.path(int_tbl, paste0("down_int_", drug, ".tsv")), analysis_level)
  
  gate_up <- load_id_table(file.path(de_tbl,  paste0("up_",   drug, "_OHT__", drug, ".tsv")), analysis_level)
  gate_dn <- load_id_table(file.path(de_tbl,  paste0("down_", drug, "_OHT__", drug, ".tsv")), analysis_level)
  
  if (analysis_level == "gene") {
    
    ids <- unique(na.omit(c(up_no$gene_id, up_oht$gene_id, dn_no$gene_id, dn_oht$gene_id,
                            up_int$gene_id, dn_int$gene_id,
                            gate_up$gene_id, gate_dn$gene_id)))
    if (length(ids) == 0) next
    
    upset_df <- tibble(gene_id = ids)
    ref <- upset_df$gene_id
    
    upset_df[[paste0("up_", drug)]]                    <- ref %in% up_no$gene_id
    upset_df[[paste0("up_", drug, "_OHT")]]            <- ref %in% up_oht$gene_id
    upset_df[[paste0("down_", drug)]]                  <- ref %in% dn_no$gene_id
    upset_df[[paste0("down_", drug, "_OHT")]]          <- ref %in% dn_oht$gene_id
    upset_df[[paste0("up_int_", drug)]]                <- ref %in% up_int$gene_id
    upset_df[[paste0("down_int_", drug)]]              <- ref %in% dn_int$gene_id
    upset_df[[paste0("up_", drug, "_OHT__", drug)]]   <- ref %in% gate_up$gene_id
    upset_df[[paste0("down_", drug, "_OHT__", drug)]] <- ref %in% gate_dn$gene_id
    
    id_cols <- c("gene_id")
    
  } else {
    
    ids <- unique(na.omit(c(up_no$isoform_id, up_oht$isoform_id, dn_no$isoform_id, dn_oht$isoform_id,
                            up_int$isoform_id, dn_int$isoform_id,
                            gate_up$isoform_id, gate_dn$isoform_id)))
    if (length(ids) == 0) next
    
    map_df <- bind_rows(up_no, up_oht, dn_no, dn_oht, up_int, dn_int, gate_up, gate_dn) %>%
      distinct(isoform_id, gene_id)
    
    upset_df <- tibble(isoform_id = ids) %>% left_join(map_df, by = "isoform_id")
    ref <- upset_df$isoform_id
    
    upset_df[[paste0("up_", drug)]]                    <- ref %in% up_no$isoform_id
    upset_df[[paste0("up_", drug, "_OHT")]]            <- ref %in% up_oht$isoform_id
    upset_df[[paste0("down_", drug)]]                  <- ref %in% dn_no$isoform_id
    upset_df[[paste0("down_", drug, "_OHT")]]          <- ref %in% dn_oht$isoform_id
    upset_df[[paste0("up_int_", drug)]]                <- ref %in% up_int$isoform_id
    upset_df[[paste0("down_int_", drug)]]              <- ref %in% dn_int$isoform_id
    upset_df[[paste0("up_", drug, "_OHT__", drug)]]   <- ref %in% gate_up$isoform_id
    upset_df[[paste0("down_", drug, "_OHT__", drug)]] <- ref %in% gate_dn$isoform_id
    
    id_cols <- c("isoform_id","gene_id")
  }
  
  write.table(upset_df, file.path(out_tbl, paste0("upset_data_", drug, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  sets_all <- setdiff(names(upset_df), id_cols)
  gate_cols <- c(paste0("up_", drug, "_OHT__", drug),
                 paste0("down_", drug, "_OHT__", drug))
  sets_plot <- setdiff(sets_all, gate_cols)
  
  p <- upset(
    upset_df,
    intersect = sets_plot,
    base_annotations = list(
      "Intersection size" = intersection_size(counts = TRUE, text = list(size = 3.5))
    ),
    set_sizes = upset_set_size(),
    width_ratio = 0.3,
    height_ratio = 0.3
  ) + theme_bw(base_size = 9)
  
  save_both(file.path(out_fig, paste0("UpSet_", drug)), p, w = 9, h = 6)
}

message(">>> Part 1 done.")

## ============================================================
## PART 2 — Categories + save tables (+ ADDITIVE)
## ============================================================

for (drug in drugs) {
  message(">>> [Categories] ", drug)
  
  fp <- file.path(out_tbl, paste0("upset_data_", drug, ".tsv"))
  if (!file.exists(fp)) next
  df <- read.delim(fp, check.names = FALSE)
  
  ## column names in upset table
  up_no   <- paste0("up_", drug)
  up_oht  <- paste0("up_", drug, "_OHT")
  dn_no   <- paste0("down_", drug)
  dn_oht  <- paste0("down_", drug, "_OHT")
  up_int  <- paste0("up_int_", drug)
  dn_int  <- paste0("down_int_", drug)
  
  up_gate <- paste0("up_", drug, "_OHT__", drug)
  dn_gate <- paste0("down_", drug, "_OHT__", drug)
  
  has <- function(nm) if (nm %in% names(df)) df[[nm]] else rep(FALSE, nrow(df))
  
  U1 <- has(up_no);  U2 <- has(up_oht)
  D1 <- has(dn_no);  D2 <- has(dn_oht)
  UI <- has(up_int); DI <- has(dn_int)
  UG <- has(up_gate); DG <- has(dn_gate)
  
  category <- rep(NA_character_, nrow(df))
  
  ## Enhanced: drug and drug+OHT significant in same direction; interaction same direction
  category[(U1 & U2 & UI) | (!U1 & U2 & UI)] <- "enhanced_up"
  category[(D1 & D2 & DI) | (!D1 & D2 & DI)] <- "enhanced_down"
  
  ## Suppressed: drug effect reduced under OHT (interaction opposite) or lost significance
  category[(U1 & U2 & DI) | (U1 & !U2 & DI)] <- "suppressed_up"
  category[(D1 & D2 & UI) | (D1 & !D2 & UI)] <- "suppressed_down"
  
  ## Switched: sign flips with OHT + significant interaction consistent with flip
  category[(D1 & U2 & UI)] <- "switched_positive"
  category[(U1 & D2 & DI)] <- "switched_negative"
  
  ## Independent: significant in both conditions, not assigned above
  un <- is.na(category)
  ind_up   <- un & (U1 & U2)
  ind_down <- un & (D1 & D2)
  
  category[ind_up]   <- "independent_up"
  category[ind_down] <- "independent_down"
  
  ## Shifted-baseline reclassification for enhanced/suppressed
  ## enhanced_up:   OHT lowers baseline → drug+OHT < drug → DG
  ## enhanced_down: OHT raises baseline → drug+OHT > drug → UG
  ## suppressed_up:   OHT raises baseline → drug+OHT > drug → UG
  ## suppressed_down: OHT lowers baseline → drug+OHT < drug → DG
  category[category == "enhanced_up"    & DG] <- "shifted_baseline_enhanced_up"
  category[category == "enhanced_down"  & UG] <- "shifted_baseline_enhanced_down"
  category[category == "suppressed_up"  & UG] <- "shifted_baseline_suppressed_up"
  category[category == "suppressed_down"& DG] <- "shifted_baseline_suppressed_down"
  
  ## Assign to df
  df$category <- category
  df2 <- df[!is.na(df$category), , drop = FALSE]
  if (nrow(df2) == 0) next
  
  ## -------- save category tables --------
  id_col <- pick_id_col(df2)
  
  for (nm in sort(unique(df2$category))) {
    subdf <- df2[df2$category == nm, , drop = FALSE]
    if (nrow(subdf) == 0) next
    
    if (analysis_level == "gene") {
      out_df <- subdf %>%
        transmute(gene_id = strip_version(.data[[id_col]])) %>%
        distinct() %>%
        left_join(gene2sym, by = "gene_id")
    } else {
      out_df <- subdf %>%
        transmute(
          isoform_id = strip_version(isoform_id),
          gene_id    = strip_version(gene_id)
        ) %>%
        distinct() %>%
        left_join(gene2sym, by = "gene_id")
    }
    
    write.table(out_df, file.path(out_tbl, paste0(nm, "_", drug, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  
  ## ============================================================
  ## ADDITIVE (final step): among TRUE independent only
  ## Definition:
  ## Additive_up   = independent_up   & up in OHT/DMSO
  ## Additive_down = independent_down & down in OHT/DMSO
  ## Interaction filter is kept as safety check
  ## ============================================================
  
  ## load OHT/DMSO DE genes
  up_oht_only <- load_id_table(
    file.path(de_tbl, "up_DMSO_OHT_vs_DMSO.tsv"),
    analysis_level
  )
  
  dn_oht_only <- load_id_table(
    file.path(de_tbl, "down_DMSO_OHT_vs_DMSO.tsv"),
    analysis_level
  )
  
  if (analysis_level == "gene") {
    
    indep_df <- df2 %>%
      filter(category %in% c("independent_up","independent_down")) %>%
      transmute(
        gene_id  = strip_version(gene_id),
        category = category
      ) %>%
      distinct()
    
    if (nrow(indep_df) == 0) next
    
    ## interaction genes
    up_int_ids <- load_id_table(
      file.path(int_tbl, paste0("up_int_", drug, ".tsv")),
      analysis_level
    )
    
    dn_int_ids <- load_id_table(
      file.path(int_tbl, paste0("down_int_", drug, ".tsv")),
      analysis_level
    )
    
    int_ids <- unique(c(up_int_ids$gene_id, dn_int_ids$gene_id))
    
    add_df <- indep_df %>%
      filter(!gene_id %in% int_ids)
    
    add_up <- add_df %>%
      filter(category == "independent_up",
             gene_id %in% up_oht_only$gene_id) %>%
      dplyr::select(gene_id) %>%
      distinct() %>%
      left_join(gene2sym, by = "gene_id")
    
    add_dn <- add_df %>%
      filter(category == "independent_down",
             gene_id %in% dn_oht_only$gene_id) %>%
      dplyr::select(gene_id) %>%
      distinct() %>%
      left_join(gene2sym, by = "gene_id")
    
  } else {
    
    indep_df <- df2 %>%
      filter(category %in% c("independent_up","independent_down")) %>%
      transmute(
        isoform_id = strip_version(isoform_id),
        gene_id    = strip_version(gene_id),
        category   = category
      ) %>%
      distinct()
    
    if (nrow(indep_df) == 0) next
    
    up_int_ids <- load_id_table(
      file.path(int_tbl, paste0("up_int_", drug, ".tsv")),
      analysis_level
    )
    
    dn_int_ids <- load_id_table(
      file.path(int_tbl, paste0("down_int_", drug, ".tsv")),
      analysis_level
    )
    
    int_ids <- unique(c(up_int_ids$isoform_id, dn_int_ids$isoform_id))
    
    add_df <- indep_df %>%
      filter(!isoform_id %in% int_ids)
    
    add_up <- add_df %>%
      filter(category == "independent_up",
             isoform_id %in% up_oht_only$isoform_id) %>%
      dplyr::select(isoform_id, gene_id) %>%
      distinct() %>%
      left_join(gene2sym, by = "gene_id")
    
    add_dn <- add_df %>%
      filter(category == "independent_down",
             isoform_id %in% dn_oht_only$isoform_id) %>%
      dplyr::select(isoform_id, gene_id) %>%
      distinct() %>%
      left_join(gene2sym, by = "gene_id")
  }
  
  if (nrow(add_up) > 0)
    write.table(add_up,
                file.path(out_tbl, paste0("additive_up_", drug, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  
  if (nrow(add_dn) > 0)
    write.table(add_dn,
                file.path(out_tbl, paste0("additive_down_", drug, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  
}
message(">>> Part 2 done. Outputs in results/4_categories/{tables,figures}")