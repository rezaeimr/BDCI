## ============================================================
## 2_interaction.R — Interaction analysis (drug-OHT)
## Reference levels: drug = DMSO | OHT = OFF
## Outputs:
##   results/2_interaction/tables/: tT_int_<drug>.tsv, up_int_<drug>.tsv, down_int_<drug>.tsv
##   results/2_interaction/figures/: MA_int_<drug>.png, Volcano_int_<drug>.png
## Notes:
##   - Interaction term tests whether the drug effect changes when OHT is ON vs OFF
##   - Uses user-defined n_cores (no auto-detection)
##   - Global volcano + MA scales are shared across all drugs (interaction only)
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(ggplot2)
  library(BiocParallel)
})

set.seed(1)
message(">>> 2_interaction: baseline-OHT interaction mode")

## -------------------- USER SETTINGS --------------------
n_cores <- 4                 
analysis_level <- "gene"     # "gene" or "isoform"
padj_cutoff <- 0.01
lfc_cutoff  <- 0

## -------------------- Parallel backend (RStudio-safe) -----
BiocParallel::register(SnowParam(workers = n_cores, type = "SOCK"))
options(mc.cores = n_cores)

## -------------------- Directories --------------------
base_dir    <- getwd()
data_dir    <- file.path(base_dir, "data")
results_dir <- file.path(base_dir, "results", "2_interaction")
tbl_dir     <- file.path(results_dir, "tables")
fig_dir     <- file.path(results_dir, "figures")
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

## ============================================================
## PART 1 — Interaction DESeq2 Calculation
## ============================================================

## ---- Load metadata ----
metadata <- read.delim(file.path(data_dir, "metadata.txt"), check.names = FALSE)
stopifnot(all(c("SampleID", "drug", "OHT") %in% colnames(metadata)))

metadata$OHT  <- factor(metadata$OHT, levels = c("OFF", "ON"))
metadata$drug <- factor(metadata$drug)
metadata$drug <- relevel(metadata$drug, ref = "DMSO")

## ---- Load counts ----
if (analysis_level == "gene") {
  count_file <- file.path(data_dir, "gene_counts.tsv")
  id_col <- "gene_id"
} else if (analysis_level == "isoform") {
  count_file <- file.path(data_dir, "isoform_counts.tsv")
  id_col <- "isoform_id"
} else {
  stop("analysis_level must be 'gene' or 'isoform'")
}

counts <- readr::read_tsv(count_file, show_col_types = FALSE) %>% as.data.frame()
rownames(counts) <- counts[[1]]

count_mat <- counts[, intersect(colnames(counts), metadata$SampleID), drop = FALSE]
metadata  <- metadata[match(colnames(count_mat), metadata$SampleID), ]
stopifnot(all(colnames(count_mat) == metadata$SampleID))

## ---- Load annotation (adds SYMBOL; keeps gene_id for enrichment compatibility) ----
annotation <- read.delim(file.path(data_dir, "annotation.txt"), check.names = FALSE)
sym_col <- if ("symbol" %in% names(annotation)) "symbol" else if ("gene_name" %in% names(annotation)) "gene_name" else NA

ann_keep_cols <- unique(c("gene_id", id_col, sym_col))
ann_keep_cols <- ann_keep_cols[ann_keep_cols %in% names(annotation)]
annotation2 <- annotation[, ann_keep_cols, drop = FALSE]
if (!is.na(sym_col)) names(annotation2)[names(annotation2) == sym_col] <- "SYMBOL"

if (analysis_level == "gene") {
  annotation2 <- annotation2 %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(SYMBOL = dplyr::first(na.omit(SYMBOL)), .groups = "drop")
}

if (!"gene_id" %in% names(annotation2)) stop("Annotation must include gene_id column.")

## ---- Build DESeq2 dataset ----
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(count_mat)),
  colData   = metadata,
  design    = ~ drug + OHT + drug:OHT
)

dds <- DESeq(dds, parallel = TRUE)

## ---- Identify interaction terms robustly ----
res_names <- resultsNames(dds)
interaction_terms <- res_names[grepl("^drug", res_names) & grepl("OHT", res_names)]
if (length(interaction_terms) == 0) stop("No interaction terms detected in resultsNames(dds).")

message("Detected interaction terms: ", paste(interaction_terms, collapse = ", "))

term_to_drug <- function(term) {
  x <- term
  x <- sub("^drug_", "", x)
  x <- sub("^drug",  "", x)
  x <- sub("_vs_.*$", "", x)
  x <- sub("\\.OHT.*$", "", x)
  x
}

## ---- Run + save tables ----
for (term in interaction_terms) {
  drug <- term_to_drug(term)
  message("\n>>> Interaction term: ", term, " | drug: ", drug)

  res <- results(dds, name = term)
  res_df <- as.data.frame(res)

  # ID column (strip version)
  res_df[[id_col]] <- sub("\\..*$", "", rownames(res_df))

  # Merge annotation
  res_df <- merge(res_df, annotation2, by = id_col, all.x = TRUE, sort = FALSE)

  # Significance
  res_df$sig <- ifelse(!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange >  lfc_cutoff, "up",
                       ifelse(!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange < -lfc_cutoff, "down", "ns"))

  write.table(res_df, file.path(tbl_dir, paste0("tT_int_", drug, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(subset(res_df, sig == "up"), file.path(tbl_dir, paste0("up_int_", drug, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(subset(res_df, sig == "down"), file.path(tbl_dir, paste0("down_int_", drug, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

message("\n>>> Part 1 complete — Interaction tables saved to: ", tbl_dir)

## ============================================================
## PART 2 — Plotting (MA & Volcano) with shared auto-scales
## ============================================================

col_scale <- c(up = "#D62728", down = "#1F77B4", ns = "gray80")

files <- list.files(tbl_dir, pattern = "^tT_int_.*\\.tsv$", full.names = TRUE)
if (length(files) == 0) stop("No tT_int_*.tsv files found in: ", tbl_dir)

safe_num <- function(x) x[is.finite(x) & !is.na(x)]

## ---- Global limits across all interaction tables ----
all_tt <- lapply(files, function(fp) {
  df <- tryCatch(read.delim(fp, check.names = FALSE), error = function(e) NULL)
  if (is.null(df)) return(NULL)
  if (!all(c("baseMean","log2FoldChange","padj","sig") %in% names(df))) return(NULL)
  df <- df[!is.na(df$log2FoldChange) & !is.na(df$padj) & !is.na(df$baseMean), , drop = FALSE]
  if (nrow(df) == 0) return(NULL)
  df
})
all_tt <- Filter(Negate(is.null), all_tt)

if (length(all_tt) > 0) {
  all_lfc <- safe_num(unlist(lapply(all_tt, `[[`, "log2FoldChange"), use.names = FALSE))
  all_maX <- safe_num(unlist(lapply(all_tt, function(d) log10(d$baseMean + 1)), use.names = FALSE))
  all_maY <- safe_num(unlist(lapply(all_tt, `[[`, "log2FoldChange"), use.names = FALSE))
  all_vY  <- safe_num(unlist(lapply(all_tt, function(d) -log10(pmax(d$padj, .Machine$double.xmin))), use.names = FALSE))

  # Volcano x/y
  x_max_v <- max(abs(all_lfc), na.rm = TRUE)
  x_lim_v <- c(-x_max_v, x_max_v)
  y_max_v <- max(all_vY, na.rm = TRUE)
  y_lim_v <- c(0, y_max_v)

  # MA x/y
  x_max_ma <- max(all_maX, na.rm = TRUE)
  x_lim_ma <- c(0, x_max_ma)
  y_max_ma <- max(abs(all_maY), na.rm = TRUE)
  y_lim_ma <- c(-y_max_ma, y_max_ma)
} else {
  x_lim_v <- y_lim_v <- x_lim_ma <- y_lim_ma <- NULL
}

for (fp in files) {
  nm <- sub("^tT_int_", "", tools::file_path_sans_ext(basename(fp)))
  res_df <- tryCatch(read.delim(fp, check.names = FALSE), error = function(e) NULL)
  if (is.null(res_df) || nrow(res_df) == 0) {
    message(" [skip] ", nm, " — empty or unreadable.")
    next
  }

  needed <- c("baseMean", "log2FoldChange", "padj", "sig")
  if (!all(needed %in% names(res_df))) {
    message(" [skip] ", nm, " — missing required columns: ", paste(setdiff(needed, names(res_df)), collapse = ", "))
    next
  }

  res_df$cat <- factor(res_df$sig, levels = c("ns", "down", "up"))
  plot_df <- res_df %>% arrange(cat)

  ## --- MA plot ---
  p_ma <- ggplot(plot_df, aes(x = log10(baseMean + 1), y = log2FoldChange, color = cat)) +
    geom_point(size = 0.6, alpha = 0.8) +
    scale_color_manual(values = col_scale,
                       breaks = c("up","down","ns"),
                       labels = c("Up","Down","NS"),
                       name = "DE category") +
    geom_hline(yintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    labs(title = paste("MA (interaction):", nm),
         x = "log10(baseMean + 1)",
         y = "log2(Fold Change)") +
    theme_bw(base_size = 10)

  if (!is.null(x_lim_ma) && !is.null(y_lim_ma)) {
    p_ma <- p_ma + coord_cartesian(xlim = x_lim_ma, ylim = y_lim_ma)
  }

  ## --- Volcano plot ---
  p_vol <- ggplot(plot_df, aes(x = log2FoldChange,
                               y = -log10(pmax(padj, .Machine$double.xmin)),
                               color = cat)) +
    geom_point(size = 0.6, alpha = 0.8) +
    scale_color_manual(values = col_scale,
                       breaks = c("up","down","ns"),
                       labels = c("Up","Down","NS"),
                       name = "DE category") +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
    labs(title = paste("Volcano (interaction):", nm),
         x = "log2(Fold Change)",
         y = "-log10(adjusted p-value)") +
    theme_bw(base_size = 10)

  if (!is.null(x_lim_v) && !is.null(y_lim_v)) {
    p_vol <- p_vol + coord_cartesian(xlim = x_lim_v, ylim = y_lim_v)
  }

  ggsave(file.path(fig_dir, paste0("MA_int_", nm, ".png")), p_ma, width = 6, height = 5, dpi = 300)
  ggsave(file.path(fig_dir, paste0("Volcano_int_", nm, ".png")), p_vol, width = 6.5, height = 5, dpi = 300)
}

message("\n>>> Part 2 complete — MA & Volcano plots saved to: ", fig_dir)