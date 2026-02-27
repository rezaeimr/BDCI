## =========================
## 6_scatterplots.R
## =========================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
})

## -------------------- Parameters --------------------
analysis_level <- "gene"

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

cols_shift_ind <- c(
  shifted_baseline_independent_up   = "#E31A1C",
  shifted_baseline_independent_down = "#1F78B4"
)

## -------------------- Sizes --------------------
pt_other <- 1.2
pt_cat   <- 1.6

## -------------------- Helpers --------------------
safe_num <- function(x) x[is.finite(x) & !is.na(x)]

pick_id_col <- function(df) {
  cand <- intersect(c("isoform_id", "gene_id", "feature_id"), names(df))
  if (length(cand) == 0) stop("No ID column found.")
  cand[1]
}

load_tt <- function(name) {
  read.delim(file.path(de_tbl, paste0("tT_", name, ".tsv")), check.names = FALSE)
}

load_ids <- function(tag, drug) {
  fn <- file.path(cat_tbl, paste0(tag, "_", drug, ".tsv"))
  if (!file.exists(fn)) return(character(0))
  df <- read.delim(fn, check.names = FALSE)
  if ("gene_id" %in% names(df)) unique(na.omit(df$gene_id)) else character(0)
}

make_merged_fc <- function(drug) {
  t_no  <- load_tt(drug)
  t_yes <- load_tt(paste0(drug, "_OHT"))
  
  id_no  <- pick_id_col(t_no)
  id_yes <- pick_id_col(t_yes)
  
  t_no  <- dplyr::rename(t_no,  feature_id = !!id_no)
  t_yes <- dplyr::rename(t_yes, feature_id = !!id_yes)
  
  full_join(
    t_no  %>% dplyr::select(feature_id, log2FoldChange) %>%
      rename(log2FC_noOHT = log2FoldChange),
    t_yes %>% dplyr::select(feature_id, log2FoldChange) %>%
      rename(log2FC_withOHT = log2FoldChange),
    by = "feature_id"
  )
}

axis_x <- function(drug)
  bquote(bold(atop(.(drug) ~ "/ DMSO", "(" * log[2] * "FC" * ")")))

axis_y <- function(drug)
  bquote(bold(atop(.(drug) ~ "+ OHT / DMSO + OHT", "(" * log[2] * "FC" * ")")))

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

save_both <- function(name, plot, w, h) {
  ggsave(paste0(name, ".pdf"), plot, width = w, height = h)
  ggsave(paste0(name, ".tiff"), plot,
         width = w, height = h, dpi = 600,
         device = "tiff", compression = "lzw")
}

get_legend <- function(p) {
  g <- ggplotGrob(p)
  idx <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
  if (length(idx) == 0) return(NULL)
  g$grobs[[idx[1]]]
}

make_vertical_legend <- function(p, title) {
  p +
    guides(color = guide_legend(ncol = 1, byrow = TRUE)) +
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      legend.box = "vertical"
    ) +
    labs(color = title)
}

## -------------------- Global limits --------------------
scatter_list <- lapply(drug_list, make_merged_fc)
names(scatter_list) <- drug_list

all_x <- safe_num(unlist(lapply(scatter_list, function(d) d$log2FC_noOHT)))
all_y <- safe_num(unlist(lapply(scatter_list, function(d) d$log2FC_withOHT)))

lim_max <- max(abs(c(all_x, all_y)), na.rm = TRUE)
xy_lim <- c(-lim_max, lim_max)

## -------------------- Lists --------------------
p_main_list      <- list()
p_ind_list       <- list()
p_shift_list     <- list()
p_shift_ind_list <- list()

## ==================== LOOP ====================
for (drug in drug_list) {
  
  df <- scatter_list[[drug]]
  
  df$main  <- "Other"
  df$ind   <- "Other"
  df$shift <- "Other"
  df$shift_ind <- "Other"
  
  ## Load IDs
  df$main[df$feature_id %in% load_ids("enhanced_up", drug)]       <- "enhanced_up"
  df$main[df$feature_id %in% load_ids("enhanced_down", drug)]     <- "enhanced_down"
  df$main[df$feature_id %in% load_ids("suppressed_up", drug)]     <- "suppressed_up"
  df$main[df$feature_id %in% load_ids("suppressed_down", drug)]   <- "suppressed_down"
  df$main[df$feature_id %in% load_ids("switched_positive", drug)] <- "switched_positive"
  df$main[df$feature_id %in% load_ids("switched_negative", drug)] <- "switched_negative"
  
  df$ind[df$feature_id %in% load_ids("independent_up", drug)]   <- "independent_up"
  df$ind[df$feature_id %in% load_ids("independent_down", drug)] <- "independent_down"
  
  df$shift[df$feature_id %in% load_ids("shifted_baseline_enhanced_up", drug)]     <- "shifted_baseline_enhanced_up"
  df$shift[df$feature_id %in% load_ids("shifted_baseline_enhanced_down", drug)]   <- "shifted_baseline_enhanced_down"
  df$shift[df$feature_id %in% load_ids("shifted_baseline_suppressed_up", drug)]   <- "shifted_baseline_suppressed_up"
  df$shift[df$feature_id %in% load_ids("shifted_baseline_suppressed_down", drug)] <- "shifted_baseline_suppressed_down"
  
  df$shift_ind[df$feature_id %in% load_ids("shifted_baseline_independent_up", drug)]   <- "shifted_baseline_independent_up"
  df$shift_ind[df$feature_id %in% load_ids("shifted_baseline_independent_down", drug)] <- "shifted_baseline_independent_down"
  
  plot_fun <- function(cat_col, cols, legend_title) {
    ggplot(df, aes(x = log2FC_noOHT, y = log2FC_withOHT)) +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
      geom_point(data = df[df[[cat_col]] == "Other",], color = col_gray, size = pt_other) +
      geom_point(data = df[df[[cat_col]] != "Other",],
                 aes(color = .data[[cat_col]]), size = pt_cat) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.8) +
      scale_color_manual(values = cols, drop = FALSE, name = legend_title) +
      labs(title = drug, x = axis_x(drug), y = axis_y(drug)) +
      coord_fixed() +
      coord_cartesian(xlim = xy_lim, ylim = xy_lim) +
      theme_scatter()
  }
  
  p1 <- plot_fun("main",      cols_main,      "Category")
  p2 <- plot_fun("ind",       cols_ind,       "Independent")
  p3 <- plot_fun("shift",     cols_shift,     "Shifted")
  p4 <- plot_fun("shift_ind", cols_shift_ind, "Shifted Independent")
  
  save_both(file.path(fig_out, drug), p1 + theme(legend.position = "none"), 6, 6)
  save_both(file.path(fig_out, paste0(drug,"_independent")), p2 + theme(legend.position = "none"), 6, 6)
  save_both(file.path(fig_out, paste0(drug,"_shifted")), p3 + theme(legend.position = "none"), 6, 6)
  save_both(file.path(fig_out, paste0(drug,"_shifted_independent")), p4 + theme(legend.position = "none"), 6, 6)
  
  p_main_list[[drug]]      <- p1
  p_ind_list[[drug]]       <- p2
  p_shift_list[[drug]]     <- p3
  p_shift_ind_list[[drug]] <- p4
}

## ---------------- Combined Panels (no legend) ----------------
save_both(file.path(fig_out,"ALL_main"), wrap_plots(p_main_list,ncol=3)&theme(legend.position="none"),18,12)
save_both(file.path(fig_out,"ALL_independent"), wrap_plots(p_ind_list,ncol=3)&theme(legend.position="none"),18,12)
save_both(file.path(fig_out,"ALL_shifted"), wrap_plots(p_shift_list,ncol=3)&theme(legend.position="none"),18,12)
save_both(file.path(fig_out,"ALL_shifted_independent"), wrap_plots(p_shift_ind_list,ncol=3)&theme(legend.position="none"),18,12)

## ---------------- Vertical Legends ----------------
legend_types <- list(
  list(p_main_list[[1]], "Category", "LEGEND_main"),
  list(p_ind_list[[1]], "Independent", "LEGEND_independent"),
  list(p_shift_list[[1]], "Shifted", "LEGEND_shifted"),
  list(p_shift_ind_list[[1]], "Shifted Independent", "LEGEND_shifted_independent")
)

for (lg in legend_types) {
  p_leg <- make_vertical_legend(lg[[1]], lg[[2]])
  leg <- get_legend(p_leg)
  if (!is.null(leg)) {
    save_both(file.path(fig_out, lg[[3]]),
              patchwork::wrap_elements(full = leg),
              4, 8)
  }
}

message("Done.")