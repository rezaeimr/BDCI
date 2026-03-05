## ============================================================
## 5_ORA_Categories.R — ORA per category
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(msigdbr)
  library(org.Mm.eg.db)
})

set.seed(1)
message(">>> 5_ORA_Categories running")

## -------------------- Parameters --------------------
analysis_level <- "gene"
padj_cutoff    <- 0.05
species_name   <- "Mus musculus"

## Colors (KEEP CONSISTENT with 3_enrichment)
col_up    <- "#D62728"
col_down  <- "#1F77B4"
col_other <- "#7F7F7F"

strip_version <- function(x) sub("\\..*$", "", x)

## -------------------- Paths --------------------
project_root <- getwd()
data_dir     <- file.path(project_root, "data")

cat_tbl  <- file.path(project_root, "results/4_categories/tables")
out_root <- file.path(project_root, "results/5_ORA_categories")
tbl_dir  <- file.path(out_root, "tables")
fig_dir  <- file.path(out_root, "figures")
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

## -------------------- Annotation (optional symbols) --------------------
annotation_fp <- file.path(data_dir, "annotation.txt")
id2sym <- NULL
if (file.exists(annotation_fp)) {
  annotation <- read.delim(annotation_fp, check.names = FALSE)
  sym_col <- if ("symbol" %in% names(annotation)) "symbol" else if ("gene_name" %in% names(annotation)) "gene_name" else NA
  if (!is.na(sym_col) && "gene_id" %in% names(annotation)) {
    ann2 <- annotation %>%
      dplyr::transmute(gene_id = strip_version(gene_id), SYMBOL = .data[[sym_col]]) %>%
      dplyr::distinct(gene_id, .keep_all = TRUE)
    id2sym <- setNames(ann2$SYMBOL, ann2$gene_id)
  }
}

## -------------------- Universe (REQUIRED) --------------------
universe_fp <- file.path(project_root, "results", "1_DEAnalysis", "tables",
                         paste0("expressed_universe_", analysis_level, ".tsv"))
if (!file.exists(universe_fp)) stop("Missing universe: ", universe_fp)

u <- read.delim(universe_fp, check.names = FALSE)
if (ncol(u) < 1) stop("Universe file has no columns: ", universe_fp)

universe_ens <- unique(na.omit(strip_version(u[[1]])))
if (length(universe_ens) < 50) stop("Universe too small (<50): ", universe_fp)

## Universe (ENTREZ) for KEGG ORA
universe_entrez <- mapIds(
  org.Mm.eg.db,
  keys      = universe_ens,
  column    = "ENTREZID",
  keytype   = "ENSEMBL",
  multiVals = "first"
) %>% unname() %>% na.omit() %>% unique()

## -------------------- MSigDB gene sets --------------------
message("Loading MSigDB...")
msigdb <- list(
  HALLMARK = msigdbr::msigdbr(species = species_name, category = "H"),
  GOBP     = msigdbr::msigdbr(species = species_name, category = "C5", subcategory = "BP"),
  C2CP     = msigdbr::msigdbr(species = species_name, category = "C2", subcategory = "CP")
)

term2gene_list <- lapply(msigdb, function(df) {
  df %>%
    dplyr::transmute(gs_name, ensembl_gene = strip_version(ensembl_gene)) %>%
    dplyr::distinct()
})

term_sizes_u <- lapply(term2gene_list, function(t2g) {
  t2g %>%
    dplyr::filter(ensembl_gene %in% universe_ens) %>%
    dplyr::count(gs_name, name = "setSize")
})

## -------------------- Theme --------------------
theme_enrich <- function() {
  theme_bw(9) +
    theme(
      panel.grid  = element_blank(),
      plot.title  = element_text(face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 8)
    )
}

## -------------------- Plot --------------------
plot_ora_bar <- function(df_sig, out_base, fill_col) {
  if (nrow(df_sig) == 0) return(invisible(NULL))
  
  dfp <- df_sig %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = 20) %>%
    dplyr::mutate(Description = forcats::fct_reorder(Description, GeneRatio))
  
  p <- ggplot(dfp, aes(x = GeneRatio, y = Description)) +
    geom_col(fill = fill_col) +
    geom_text(aes(x = GeneRatio / 2, label = RatioLabel),
              size = 3.2, color = "black", fontface = "bold") +
    labs(x = "GeneRatio", y = NULL) +
    theme_enrich()
  
  ggsave(paste0(out_base, ".pdf"),  p, width = 7, height = 5)
  ggsave(paste0(out_base, ".tiff"), p, width = 7, height = 5, dpi = 300,
         device = "tiff", compression = "lzw")
}

## -------------------- Fill color from tag --------------------
pick_fill <- function(tag) {
  if (grepl("positive", tag, ignore.case = TRUE)) return(col_up)
  if (grepl("negative", tag, ignore.case = TRUE)) return(col_down)
  if (grepl("up",       tag, ignore.case = TRUE)) return(col_up)
  if (grepl("down",     tag, ignore.case = TRUE)) return(col_down)
  col_other
}

## ============================================================
## ORA (MSigDB + KEGG) for each category file
## (includes shifted_baseline_* automatically)
## upset_data_* files are excluded
## ============================================================
cat_files <- list.files(cat_tbl, pattern = "\\.tsv$", full.names = TRUE)
cat_files <- cat_files[!grepl("upset_data_", basename(cat_files))]
if (length(cat_files) == 0) stop("No category files found in ", cat_tbl)

for (fp in cat_files) {
  tag <- tools::file_path_sans_ext(basename(fp))
  message("ORA: ", tag)
  
  df <- read.delim(fp, check.names = FALSE)
  
  id_col <- if ("gene_id" %in% names(df)) "gene_id" else
    if ("feature_id" %in% names(df)) "feature_id" else
      if ("isoform_id" %in% names(df)) "isoform_id" else NULL
  if (is.null(id_col)) { message("  [skip] No ID column: ", tag); next }
  
  genes_ens <- unique(na.omit(strip_version(df[[id_col]])))
  if (length(genes_ens) < 1) { message("  [skip] 0 genes: ", tag); next }
  
  fill_col <- pick_fill(tag)
  inputN   <- length(genes_ens)
  
  ## ---- MSigDB ORA ----
  for (set_name in names(term2gene_list)) {
    t2g <- term2gene_list[[set_name]]
    
    res <- tryCatch(
      clusterProfiler::enricher(
        gene          = genes_ens,
        universe      = universe_ens,
        TERM2GENE     = t2g,
        pAdjustMethod = "BH",
        qvalueCutoff  = padj_cutoff
      ),
      error = function(e) NULL
    )
    if (is.null(res)) next
    
    df_out <- as.data.frame(res)
    if (nrow(df_out) == 0) next
    
    df_out <- df_out %>%
      dplyr::left_join(term_sizes_u[[set_name]], by = c("ID" = "gs_name")) %>%
      dplyr::mutate(
        setSize    = tidyr::replace_na(setSize, 0L),
        GeneRatio  = Count / pmax(inputN, 1L),
        RatioLabel = paste0(Count, "/", setSize)
      )
    
    if (!is.null(id2sym)) {
      df_out$SYMBOLS <- vapply(
        strsplit(df_out$geneID, "/"),
        function(ids) paste(na.omit(id2sym[strip_version(ids)]), collapse = ";"),
        character(1)
      )
    }
    
    out_tsv <- file.path(tbl_dir, paste0("ORA_", set_name, "_", tag, ".tsv"))
    write.table(df_out, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    
    df_sig <- df_out %>% dplyr::filter(p.adjust < padj_cutoff)
    if (nrow(df_sig) == 0) next
    
    out_base <- file.path(fig_dir, paste0("ORA_", set_name, "_", tag))
    plot_ora_bar(df_sig, out_base, fill_col = fill_col)
  }
  
  ## ---- KEGG ORA ----
  message("  KEGG ORA: ", tag)
  
  genes_entrez <- mapIds(
    org.Mm.eg.db,
    keys      = genes_ens,
    column    = "ENTREZID",
    keytype   = "ENSEMBL",
    multiVals = "first"
  ) %>% unname() %>% na.omit() %>% unique()
  
  if (length(genes_entrez) < 1 || length(universe_entrez) < 50) next
  
  kegg <- tryCatch(
    clusterProfiler::enrichKEGG(
      gene          = genes_entrez,
      universe      = universe_entrez,
      organism      = "mmu",
      pAdjustMethod = "BH",
      qvalueCutoff  = padj_cutoff
    ),
    error = function(e) NULL
  )
  if (is.null(kegg)) next
  
  df_out <- as.data.frame(kegg)
  if (nrow(df_out) == 0) next
  
  df_out <- df_out %>%
    dplyr::mutate(
      GeneRatio  = Count / pmax(length(genes_entrez), 1L),
      setSize    = as.integer(sub("/.*$", "", BgRatio)),
      RatioLabel = paste0(Count, "/", setSize)
    )
  
  out_tsv <- file.path(tbl_dir, paste0("ORA_KEGG_", tag, ".tsv"))
  write.table(df_out, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  
  df_sig <- df_out %>% dplyr::filter(p.adjust < padj_cutoff)
  if (nrow(df_sig) == 0) next
  
  out_base <- file.path(fig_dir, paste0("ORA_KEGG_", tag))
  plot_ora_bar(df_sig, out_base, fill_col = fill_col)
}

## -------------------- Reproducibility --------------------
writeLines(capture.output(sessionInfo()), file.path(out_root, "sessionInfo.txt"))

message(">>> done")
message("Tables: ", tbl_dir)
message("Figures: ", fig_dir)