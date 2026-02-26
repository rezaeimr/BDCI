## ============================================================
## 6_scatterplots.R — Two scatterplots per drug (shrunk log2FC)
## 1) All genes gray + enhanced/suppressed/switched colored (legend only, no labels)
## 2) All genes gray + independent_up/down colored (legend only, no labels)
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

## -------------------- Parameters --------------------
analysis_level <- "gene"   # or "isoform" (must match your 4_categories output)
padj_cutoff <- 0.05

## -------------------- Paths --------------------
base_root <- getwd()
de_tbl   <- file.path(base_root, "results/1_DEAnalysis/tables/shrunk")
cat_tbl  <- file.path(base_root, "results/4_categories/tables")
fig_out  <- file.path(base_root, "results/6_scatterplots/figures")
dir.create(fig_out, recursive = TRUE, showWarnings = FALSE)

## -------------------- Drug list --------------------
drug_list <- c("11j", "KVS", "11j_PlaB", "KVS_PlaB", "PlaB")

## -------------------- Colors --------------------
col_gray <- "gray88"

## Main categories (match 4_categories naming)
cols_main <- c(
  enhanced_up       = "#D62728",  # red
  enhanced_down     = "#1F77B4",  # blue
  suppressed_up     = "#FF7F0E",  # orange
  suppressed_down   = "#9467BD",  # purple
  switched_positive = "#2CA02C",  # green
  switched_negative = "#17BECF"   # cyan
)

## Independent: near red / near blue
cols_ind <- c(
  independent_up   = "#E31A1C",   # red-ish
  independent_down = "#1F78B4"    # blue-ish
)

## -------------------- Dot sizes (panel-friendly) --------------------
pt_other <- 0.85
pt_cat   <- 1.05

## -------------------- Helper functions --------------------
safe_num <- function(x) x[is.finite(x) & !is.na(x)]

pick_id_col <- function(df) {
  cand <- intersect(c("isoform_id", "gene_id", "feature_id"), names(df))
  if (length(cand) == 0) stop("No ID column found in DE table.")
  cand[1]
}

load_tt <- function(name) {
  fp <- file.path(de_tbl, paste0("tT_", name, ".tsv"))
  if (!file.exists(fp)) stop("Missing DE table: ", fp)
  read.delim(fp, check.names = FALSE)
}

load_ids <- function(tag, drug, analysis_level) {
  fn <- file.path(cat_tbl, paste0(tag, "_", drug, ".tsv"))
  if (!file.exists(fn)) return(character(0))
  df <- read.delim(fn, check.names = FALSE)
  
  if (analysis_level == "gene" && "gene_id" %in% names(df)) {
    unique(na.omit(df$gene_id))
  } else if (analysis_level == "isoform" && "isoform_id" %in% names(df)) {
    unique(na.omit(df$isoform_id))
  } else {
    character(0)
  }
}

make_merged_fc <- function(drug) {
  t_no  <- load_tt(drug)
  t_yes <- load_tt(paste0(drug, "_OHT"))
  
  id_no  <- pick_id_col(t_no)
  id_yes <- pick_id_col(t_yes)
  t_no  <- dplyr::rename(t_no,  feature_id = !!id_no)
  t_yes <- dplyr::rename(t_yes, feature_id = !!id_yes)
  
  full_join(
    t_no  %>% dplyr::select(feature_id, log2FoldChange, padj) %>%
      dplyr::rename(log2FC_noOHT = log2FoldChange, padj_noOHT = padj),
    t_yes %>% dplyr::select(feature_id, log2FoldChange, padj) %>%
      dplyr::rename(log2FC_withOHT = log2FoldChange, padj_withOHT = padj),
    by = "feature_id"
  )
}

axis_x <- function(drug) bquote(atop(.(drug) ~ "/ DMSO", "(" * log[2] * "FC" * ")"))
axis_y <- function(drug) bquote(atop(.(drug) ~ "+ OHT / DMSO + OHT", "(" * log[2] * "FC" * ")"))

theme_scatter <- function() {
  theme_bw(base_size = 10) +
    theme(
      panel.grid   = element_blank(),
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 11),
      axis.title   = element_text(size = 10),
      axis.text    = element_text(size = 9),
      legend.title = element_text(size = 9),
      legend.text  = element_text(size = 8),
      legend.position = "right"
    )
}

## ============================================================
## Pre-pass: compute GLOBAL x/y limits across all drugs
## ============================================================
scatter_list <- lapply(drug_list, make_merged_fc)
names(scatter_list) <- drug_list

all_x <- safe_num(unlist(lapply(scatter_list, function(d) d$log2FC_noOHT),  use.names = FALSE))
all_y <- safe_num(unlist(lapply(scatter_list, function(d) d$log2FC_withOHT), use.names = FALSE))

lim_max <- max(abs(c(all_x, all_y)), na.rm = TRUE)
xy_lim <- c(-lim_max, lim_max)

## ============================================================
## Main loop: TWO plots per drug
## ============================================================
for (drug in drug_list) {
  message("Plotting: ", drug)
  
  df <- scatter_list[[drug]]
  df$category_main <- "Other"
  df$category_ind  <- "Other"
  
  ## ---- Load IDs ----
  enh_up   <- load_ids("enhanced_up",       drug, analysis_level)
  enh_down <- load_ids("enhanced_down",     drug, analysis_level)
  sup_up   <- load_ids("suppressed_up",     drug, analysis_level)
  sup_down <- load_ids("suppressed_down",   drug, analysis_level)
  swi_pos  <- load_ids("switched_positive", drug, analysis_level)
  swi_neg  <- load_ids("switched_negative", drug, analysis_level)
  
  ind_up   <- load_ids("independent_up",    drug, analysis_level)
  ind_down <- load_ids("independent_down",  drug, analysis_level)
  
  ## ---- Assign for plot 1 (main categories only) ----
  df$category_main[df$feature_id %in% enh_up]   <- "enhanced_up"
  df$category_main[df$feature_id %in% enh_down] <- "enhanced_down"
  df$category_main[df$feature_id %in% sup_up]   <- "suppressed_up"
  df$category_main[df$feature_id %in% sup_down] <- "suppressed_down"
  df$category_main[df$feature_id %in% swi_pos]  <- "switched_positive"
  df$category_main[df$feature_id %in% swi_neg]  <- "switched_negative"
  df$category_main <- factor(df$category_main, levels = c("Other", names(cols_main)))
  
  ## ---- Assign for plot 2 (independent only) ----
  df$category_ind[df$feature_id %in% ind_up]   <- "independent_up"
  df$category_ind[df$feature_id %in% ind_down] <- "independent_down"
  df$category_ind <- factor(df$category_ind, levels = c("Other", names(cols_ind)))
  
  ## ========================================================
  ## Plot 1: gray all + enhanced/suppressed/switched colored
  ## ========================================================
  p1 <- ggplot(df, aes(x = log2FC_noOHT, y = log2FC_withOHT)) +
    ## --- reference lines ---
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
    ## --- points ---
    geom_point(
      data  = df[df$category_main == "Other" | is.na(df$category_main), , drop = FALSE],
      color = col_gray, alpha = 0.55, size = pt_other
    ) +
    geom_point(
      data = df[df$category_main != "Other" & !is.na(df$category_main), , drop = FALSE],
      aes(color = category_main),
      alpha = 0.9, size = pt_cat
    ) +
    ## --- y = x diagonal (on top of points) ---
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "black", linewidth = 0.35) +
    scale_color_manual(values = cols_main, drop = FALSE, name = "Category") +
    labs(title = drug, x = axis_x(drug), y = axis_y(drug)) +
    coord_fixed(ratio = 1) +
    coord_cartesian(xlim = xy_lim, ylim = xy_lim) +
    theme_scatter()
  
  ggsave(file.path(fig_out, paste0(drug, ".png")),
         p1, width = 6, height = 5, dpi = 300)
  
  ## ========================================================
  ## Plot 2: gray all + independent_up/down colored
  ## ========================================================
  p2 <- ggplot(df, aes(x = log2FC_noOHT, y = log2FC_withOHT)) +
    ## --- reference lines ---
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
    ## --- points ---
    geom_point(
      data  = df[df$category_ind == "Other" | is.na(df$category_ind), , drop = FALSE],
      color = col_gray, alpha = 0.55, size = pt_other
    ) +
    geom_point(
      data = df[df$category_ind != "Other" & !is.na(df$category_ind), , drop = FALSE],
      aes(color = category_ind),
      alpha = 0.9, size = pt_cat
    ) +
    ## --- y = x diagonal (on top of points) ---
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "black", linewidth = 0.35) +
    scale_color_manual(values = cols_ind, drop = FALSE, name = "Independent") +
    labs(title = drug, x = axis_x(drug), y = axis_y(drug)) +
    coord_fixed(ratio = 1) +
    coord_cartesian(xlim = xy_lim, ylim = xy_lim) +
    theme_scatter()
  
  ggsave(file.path(fig_out, paste0(drug, "_independent.png")),
         p2, width = 6, height = 5, dpi = 300)
}

message("Scatterplots complete: ", fig_out)