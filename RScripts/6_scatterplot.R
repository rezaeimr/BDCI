## 6_scatterplots.R — scatter + correlation (PDF + TIFF)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
})

## -------------------- Settings --------------------
analysis_level <- "gene"  # "gene" or "isoform"

drug_list <- c("11j", "KVS", "11j_PlaB", "KVS_PlaB", "PlaB")

## Required DE tables (shrunk):
## - tT_<drug>.tsv                   : drug / DMSO
## - tT_<drug>_OHT.tsv               : (drug+OHT) / (DMSO+OHT)
## - tT_DMSO_OHT_vs_DMSO.tsv         : OHT / DMSO
## - tT_<drug>_OHT_vs_DMSO.tsv       : (drug+OHT) / DMSO   <-- MUST exist for CORR_OHT y-axis

## Required category tables:
## - results/4_categories/tables/<category>_<drug>.tsv
##   includes additive_up / additive_down

## -------------------- Paths --------------------
base_root <- getwd()
de_tbl   <- file.path(base_root, "results/1_DEAnalysis/tables/shrunk")
cat_tbl  <- file.path(base_root, "results/4_categories/tables")
fig_out  <- file.path(base_root, "results/6_scatterplots/figures")
dir.create(fig_out, recursive = TRUE, showWarnings = FALSE)

stopifnot(dir.exists(de_tbl), dir.exists(cat_tbl))

## -------------------- Colors --------------------
col_gray <- "gray88"

cols_main <- c(
  enhanced_up       = "#D62728",
  enhanced_down     = "#1F77B4",
  suppressed_up     = "#FF7F0E",
  suppressed_down   = "#9467BD",
  switched_positive = "#2CA02C",
  switched_negative = "#17BECF"
)

cols_ind <- c(
  independent_up   = "#E31A1C",
  independent_down = "#1F78B4"
)

cols_shift <- c(
  shifted_baseline_enhanced_up     = "#D62728",
  shifted_baseline_enhanced_down   = "#1F77B4",
  shifted_baseline_suppressed_up   = "#FF7F0E",
  shifted_baseline_suppressed_down = "#9467BD"
)

## Additive: red/blue identical to independent
cols_add <- c(
  additive_up   = cols_ind[["independent_up"]],
  additive_down = cols_ind[["independent_down"]]
)

## -------------------- Plot sizing --------------------
pt_other <- 1.2
pt_cat   <- 1.6

## -------------------- Helpers --------------------
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

load_ids <- function(tag, drug) {
  fp <- file.path(cat_tbl, paste0(tag, "_", drug, ".tsv"))
  if (!file.exists(fp)) return(character(0))
  df <- read.delim(fp, check.names = FALSE)
  
  if (analysis_level == "gene" && "gene_id" %in% names(df)) {
    unique(na.omit(df$gene_id))
  } else if (analysis_level == "isoform" && "isoform_id" %in% names(df)) {
    unique(na.omit(df$isoform_id))
  } else {
    character(0)
  }
}

theme_scatter <- function() {
  theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 15),
      axis.text  = element_text(face = "bold", size = 13),
      legend.title = element_text(face = "bold", size = 14),
      legend.text  = element_text(face = "bold", size = 12),
      legend.position = "right"
    )
}

save_both <- function(path_no_ext, plot, w, h) {
  ggsave(paste0(path_no_ext, ".pdf"), plot, width = w, height = h)
  ggsave(
    paste0(path_no_ext, ".tiff"),
    plot, width = w, height = h, dpi = 600,
    device = "tiff", compression = "lzw"
  )
}

get_legend <- function(p) {
  g <- ggplotGrob(p)
  idx <- which(vapply(g$grobs, \(x) x$name, character(1)) == "guide-box")
  if (length(idx) == 0) return(NULL)
  g$grobs[[idx[1]]]
}

make_vertical_legend <- function(p, title) {
  p +
    guides(color = guide_legend(ncol = 1, byrow = TRUE)) +
    theme(
      legend.position  = "right",
      legend.direction = "vertical",
      legend.box       = "vertical"
    ) +
    labs(color = title)
}

plot_xy <- function(df, xcol, ycol, cat_col, cols, title, xlab, ylab, lim, legend_title) {
  ggplot(df, aes(x = .data[[xcol]], y = .data[[ycol]])) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
    geom_point(data = df[df[[cat_col]] == "Other",], color = col_gray, size = pt_other, alpha = 0.6) +
    geom_point(data = df[df[[cat_col]] != "Other",], aes(color = .data[[cat_col]]), size = pt_cat, alpha = 0.9) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.8) +
    scale_color_manual(values = cols, drop = FALSE, name = legend_title) +
    labs(title = title, x = xlab, y = ylab) +
    coord_fixed() +
    coord_cartesian(xlim = lim, ylim = lim) +
    theme_scatter()
}

## -------------------- Axis labels --------------------
axis_x <- function(drug)
  bquote(bold(atop(.(drug) ~ "/ DMSO", "(" * log[2] * "FC" * ")")))

axis_y <- function(drug)
  bquote(bold(atop(.(drug) ~ "+ OHT / DMSO + OHT", "(" * log[2] * "FC" * ")")))

corr_axis_x <- function()
  bquote(bold(atop("OHT / DMSO", "(" * log[2] * "FC" * ")")))

corr_axis_y <- function(drug)
  bquote(bold(atop(.(drug) ~ "+ OHT / DMSO", "(" * log[2] * "FC" * ")")))

## ============================================================
## PART 1) Scatter: x = drug/DMSO, y = (drug+OHT)/(DMSO+OHT)
## ============================================================

make_scatter_df <- function(drug) {
  t_no  <- load_tt(drug)
  t_yes <- load_tt(paste0(drug, "_OHT"))
  
  id_no  <- pick_id_col(t_no)
  id_yes <- pick_id_col(t_yes)
  
  t_no  <- dplyr::rename(t_no,  feature_id = !!id_no)
  t_yes <- dplyr::rename(t_yes, feature_id = !!id_yes)
  
  full_join(
    t_no  %>% select(feature_id, log2FoldChange) %>% rename(x = log2FoldChange),
    t_yes %>% select(feature_id, log2FoldChange) %>% rename(y = log2FoldChange),
    by = "feature_id"
  )
}

scatter_list <- lapply(drug_list, make_scatter_df)
names(scatter_list) <- drug_list

sx <- safe_num(unlist(lapply(scatter_list, \(d) d$x), use.names = FALSE))
sy <- safe_num(unlist(lapply(scatter_list, \(d) d$y), use.names = FALSE))
lim_scatter <- max(abs(c(sx, sy)), na.rm = TRUE)
scatter_lim <- c(-lim_scatter, lim_scatter)

p_main_list  <- list()
p_ind_list   <- list()
p_shift_list <- list()

for (drug in drug_list) {
  
  df <- scatter_list[[drug]]
  df$main  <- "Other"
  df$ind   <- "Other"
  df$shift <- "Other"
  
  df$main[df$feature_id %in% load_ids("enhanced_up", drug)]       <- "enhanced_up"
  df$main[df$feature_id %in% load_ids("enhanced_down", drug)]     <- "enhanced_down"
  df$main[df$feature_id %in% load_ids("suppressed_up", drug)]     <- "suppressed_up"
  df$main[df$feature_id %in% load_ids("suppressed_down", drug)]   <- "suppressed_down"
  df$main[df$feature_id %in% load_ids("switched_positive", drug)] <- "switched_positive"
  df$main[df$feature_id %in% load_ids("switched_negative", drug)] <- "switched_negative"
  
  df$ind[df$feature_id %in% load_ids("independent_up", drug)]   <- "independent_up"
  df$ind[df$feature_id %in% load_ids("independent_down", drug)] <- "independent_down"
  
  df$shift[df$feature_id %in% load_ids("shifted_baseline_enhanced_up", drug)]       <- "shifted_baseline_enhanced_up"
  df$shift[df$feature_id %in% load_ids("shifted_baseline_enhanced_down", drug)]     <- "shifted_baseline_enhanced_down"
  df$shift[df$feature_id %in% load_ids("shifted_baseline_suppressed_up", drug)]     <- "shifted_baseline_suppressed_up"
  df$shift[df$feature_id %in% load_ids("shifted_baseline_suppressed_down", drug)]   <- "shifted_baseline_suppressed_down"
  
  p1 <- plot_xy(df, "x", "y", "main",  cols_main,  drug, axis_x(drug), axis_y(drug), scatter_lim, "Category")
  p2 <- plot_xy(df, "x", "y", "ind",   cols_ind,   drug, axis_x(drug), axis_y(drug), scatter_lim, "Independent")
  p3 <- plot_xy(df, "x", "y", "shift", cols_shift, drug, axis_x(drug), axis_y(drug), scatter_lim, "Shifted")
  
  save_both(file.path(fig_out, drug),                        p1 + theme(legend.position = "none"), 6, 6)
  save_both(file.path(fig_out, paste0(drug, "_independent")), p2 + theme(legend.position = "none"), 6, 6)
  save_both(file.path(fig_out, paste0(drug, "_shifted")),     p3 + theme(legend.position = "none"), 6, 6)
  
  p_main_list[[drug]]  <- p1
  p_ind_list[[drug]]   <- p2
  p_shift_list[[drug]] <- p3
}

save_both(file.path(fig_out, "ALL_main"),
          wrap_plots(p_main_list, ncol = 3) & theme(legend.position = "none"), 18, 12)
save_both(file.path(fig_out, "ALL_independent"),
          wrap_plots(p_ind_list, ncol = 3) & theme(legend.position = "none"), 18, 12)
save_both(file.path(fig_out, "ALL_shifted"),
          wrap_plots(p_shift_list, ncol = 3) & theme(legend.position = "none"), 18, 12)

legend_types <- list(
  list(p_main_list[[1]],  "Category",    "LEGEND_main"),
  list(p_ind_list[[1]],   "Independent", "LEGEND_independent"),
  list(p_shift_list[[1]], "Shifted",     "LEGEND_shifted")
)

for (lg in legend_types) {
  p_leg <- make_vertical_legend(lg[[1]], lg[[2]])
  leg <- get_legend(p_leg)
  if (!is.null(leg)) {
    save_both(file.path(fig_out, lg[[3]]), patchwork::wrap_elements(full = leg), 4, 8)
  }
}

## ============================================================
## PART 2) CORR_OHT: x = OHT/DMSO, y = (drug+OHT)/DMSO
## + Additive overlay (additive_up/down only)
## ============================================================

corr_make_df <- function(drug) {
  t_x <- load_tt("DMSO_OHT_vs_DMSO")
  t_y <- load_tt(paste0(drug, "_OHT_vs_DMSO"))
  
  id_x <- pick_id_col(t_x)
  id_y <- pick_id_col(t_y)
  
  t_x <- dplyr::rename(t_x, feature_id = !!id_x)
  t_y <- dplyr::rename(t_y, feature_id = !!id_y)
  
  full_join(
    t_x %>% select(feature_id, log2FoldChange) %>% rename(x = log2FoldChange),
    t_y %>% select(feature_id, log2FoldChange) %>% rename(y = log2FoldChange),
    by = "feature_id"
  )
}

corr_list <- lapply(drug_list, corr_make_df)
names(corr_list) <- drug_list

cx <- safe_num(unlist(lapply(corr_list, \(d) d$x), use.names = FALSE))
cy <- safe_num(unlist(lapply(corr_list, \(d) d$y), use.names = FALSE))
lim_corr <- max(abs(c(cx, cy)), na.rm = TRUE)
corr_lim <- c(-lim_corr, lim_corr)

corr_main_list  <- list()
corr_ind_list   <- list()
corr_shift_list <- list()
corr_add_list   <- list()

for (drug in drug_list) {
  
  df <- corr_list[[drug]]
  
  df$main  <- "Other"
  df$ind   <- "Other"
  df$shift <- "Other"
  df$add   <- "Other"
  
  df$main[df$feature_id %in% load_ids("enhanced_up", drug)]       <- "enhanced_up"
  df$main[df$feature_id %in% load_ids("enhanced_down", drug)]     <- "enhanced_down"
  df$main[df$feature_id %in% load_ids("suppressed_up", drug)]     <- "suppressed_up"
  df$main[df$feature_id %in% load_ids("suppressed_down", drug)]   <- "suppressed_down"
  df$main[df$feature_id %in% load_ids("switched_positive", drug)] <- "switched_positive"
  df$main[df$feature_id %in% load_ids("switched_negative", drug)] <- "switched_negative"
  
  df$ind[df$feature_id %in% load_ids("independent_up", drug)]   <- "independent_up"
  df$ind[df$feature_id %in% load_ids("independent_down", drug)] <- "independent_down"
  
  df$shift[df$feature_id %in% load_ids("shifted_baseline_enhanced_up", drug)]       <- "shifted_baseline_enhanced_up"
  df$shift[df$feature_id %in% load_ids("shifted_baseline_enhanced_down", drug)]     <- "shifted_baseline_enhanced_down"
  df$shift[df$feature_id %in% load_ids("shifted_baseline_suppressed_up", drug)]     <- "shifted_baseline_suppressed_up"
  df$shift[df$feature_id %in% load_ids("shifted_baseline_suppressed_down", drug)]   <- "shifted_baseline_suppressed_down"
  
  df$add[df$feature_id %in% load_ids("additive_up", drug)]   <- "additive_up"
  df$add[df$feature_id %in% load_ids("additive_down", drug)] <- "additive_down"
  
  p1 <- plot_xy(df, "x", "y", "main",  cols_main,  drug, corr_axis_x(), corr_axis_y(drug), corr_lim, "Category")
  p2 <- plot_xy(df, "x", "y", "ind",   cols_ind,   drug, corr_axis_x(), corr_axis_y(drug), corr_lim, "Independent")
  p3 <- plot_xy(df, "x", "y", "shift", cols_shift, drug, corr_axis_x(), corr_axis_y(drug), corr_lim, "Shifted")
  p5 <- plot_xy(df, "x", "y", "add",   cols_add,   drug, corr_axis_x(), corr_axis_y(drug), corr_lim, "Additive")
  
  save_both(file.path(fig_out, paste0("CORR_OHT_", drug)),
            p1 + theme(legend.position = "none"), 6, 6)
  save_both(file.path(fig_out, paste0("CORR_OHT_", drug, "_independent")),
            p2 + theme(legend.position = "none"), 6, 6)
  save_both(file.path(fig_out, paste0("CORR_OHT_", drug, "_shifted")),
            p3 + theme(legend.position = "none"), 6, 6)
  save_both(file.path(fig_out, paste0("CORR_OHT_", drug, "_additive")),
            p5 + theme(legend.position = "none"), 6, 6)
  
  corr_main_list[[drug]]  <- p1
  corr_ind_list[[drug]]   <- p2
  corr_shift_list[[drug]] <- p3
  corr_add_list[[drug]]   <- p5
}

save_both(file.path(fig_out, "CORR_OHT_ALL_main"),
          wrap_plots(corr_main_list, ncol = 3) & theme(legend.position = "none"), 18, 12)
save_both(file.path(fig_out, "CORR_OHT_ALL_independent"),
          wrap_plots(corr_ind_list, ncol = 3) & theme(legend.position = "none"), 18, 12)
save_both(file.path(fig_out, "CORR_OHT_ALL_shifted"),
          wrap_plots(corr_shift_list, ncol = 3) & theme(legend.position = "none"), 18, 12)
save_both(file.path(fig_out, "CORR_OHT_ALL_additive"),
          wrap_plots(corr_add_list, ncol = 3) & theme(legend.position = "none"), 18, 12)

corr_legend_types <- list(
  list(corr_main_list[[1]],  "Category",    "CORR_OHT_LEGEND_main"),
  list(corr_ind_list[[1]],   "Independent", "CORR_OHT_LEGEND_independent"),
  list(corr_shift_list[[1]], "Shifted",     "CORR_OHT_LEGEND_shifted"),
  list(corr_add_list[[1]],   "Additive",    "CORR_OHT_LEGEND_additive")
)

for (lg in corr_legend_types) {
  p_leg <- make_vertical_legend(lg[[1]], lg[[2]])
  leg <- get_legend(p_leg)
  if (!is.null(leg)) {
    save_both(file.path(fig_out, lg[[3]]), patchwork::wrap_elements(full = leg), 4, 8)
  }
}

message("Done. Figures saved in: ", fig_out)