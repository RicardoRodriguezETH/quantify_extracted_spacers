suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(patchwork)
  library(edgeR)
  library(logger)
})

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
# Usage:
#   Rscript process_alignments.R <aligned_tsv> <counts_tsv> <output_dir>
#             [genome_seqname]  [cpm_max_thresh]  [cpm_range_thresh]
#
#   genome_seqname   : seqnames value for the chromosome (default: auto-detect
#                      as the seqname with the most unique alignments)
#   cpm_max_thresh   : minimum max log2-CPM for a gene to be a "hit" (default: 7)
#   cpm_range_thresh : minimum log2-CPM range across samples for a "hit" (default: 0.5)
#                      set to 0 when only one sample is present

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript process_alignments.R <aligned_tsv> <counts_tsv> <output_dir> [genome_seqname] [cpm_max_thresh] [cpm_range_thresh]")
}

aligned_path     <- args[1]
counts_path      <- args[2]
output_dir       <- args[3]
genome_seqname   <- if (length(args) >= 4 && nchar(args[4]) > 0) args[4] else NULL
cpm_max_thresh   <- if (length(args) >= 5) as.numeric(args[5]) else 7
cpm_range_thresh <- if (length(args) >= 6) as.numeric(args[6]) else 0.5

# Known non-sample columns that may appear in the counts file
COUNTS_META_COLS <- c("target_seq")

# Required columns in the aligned file
REQUIRED_ALIGNED_COLS <- c("match_type", "seqnames", "gene_name")

# ---------------------------------------------------------------------------
# Directories & logging
# ---------------------------------------------------------------------------
spacer_dir <- file.path(output_dir, "spacer_level")
gene_dir   <- file.path(output_dir, "gene_level")
plots_dir  <- file.path(output_dir, "plots")
logs_dir   <- file.path(output_dir, "logs")

for (d in c(spacer_dir, gene_dir, plots_dir, logs_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_appender(appender_file(file.path(logs_dir, "process_alignments_log.txt")))
log_layout(layout_glue_generator(format = "{time} {level}: {msg}"))

# ---------------------------------------------------------------------------
# Helper: save a ggplot with consistent sizing
# ---------------------------------------------------------------------------
save_plot <- function(p, filename, width = 8, height = 5, dpi = 150) {
  ggsave(file.path(plots_dir, filename), p,
         width = width, height = height, dpi = dpi)
  log_info("  Plot saved: plots/{filename}")
}

# ---------------------------------------------------------------------------
# 1. Load & validate data
# ---------------------------------------------------------------------------
log_info("Loading aligned spacers from: {aligned_path}")
aligned <- fread(aligned_path)
log_info("  {nrow(aligned)} rows, {ncol(aligned)} columns")

missing_cols <- setdiff(REQUIRED_ALIGNED_COLS, colnames(aligned))
if (length(missing_cols) > 0) {
  stop("aligned_tsv is missing required columns: ", paste(missing_cols, collapse = ", "))
}

log_info("Loading raw counts from: {counts_path}")
counts <- fread(counts_path)
log_info("  {nrow(counts)} spacers, {ncol(counts)} columns")

if (!"target_seq" %in% colnames(counts)) {
  stop("counts_tsv must contain a 'target_seq' column")
}

sample_cols <- setdiff(colnames(counts), COUNTS_META_COLS)
if (length(sample_cols) == 0) {
  stop("No sample columns detected in counts_tsv")
}
log_info("  {length(sample_cols)} sample(s): {paste(sample_cols, collapse = ', ')}")

n_samples <- length(sample_cols)

# ---------------------------------------------------------------------------
# 2. Merge counts with alignment results
# ---------------------------------------------------------------------------
log_info("Merging on target_seq ...")
merged <- merge(counts, aligned, by = "target_seq", all.x = TRUE)
log_info("  Merged: {nrow(merged)} rows, {ncol(merged)} columns")

# ---------------------------------------------------------------------------
# 3. Per-spacer max count; filter and sort
# ---------------------------------------------------------------------------
merged[, max_count := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = sample_cols]
n_before <- nrow(merged)
merged   <- merged[max_count > 0]
setorder(merged, -max_count, target_seq)
merged[, rank := .I]
log_info("  {nrow(merged)} / {n_before} spacers retained after max_count > 0 filter")

# ---------------------------------------------------------------------------
# 4. Auto-detect genome seqname if not supplied
# ---------------------------------------------------------------------------
if (is.null(genome_seqname)) {
  seqname_counts <- merged[match_type == "unique" & !is.na(seqnames),
                            .N, by = seqnames][order(-N)]
  genome_seqname <- seqname_counts[1, seqnames]
  log_info("  Auto-detected genome seqname: '{genome_seqname}' ({seqname_counts[1,N]} unique hits)")
} else {
  log_info("  Using supplied genome seqname: '{genome_seqname}'")
}
other_seqnames <- merged[match_type == "unique" & !is.na(seqnames) &
                           seqnames != genome_seqname, unique(seqnames)]
log_info("  Other seqname(s): {paste(other_seqnames, collapse = ', ')}")

# Reduce range threshold to 0 for single-sample runs
if (n_samples == 1 && cpm_range_thresh > 0) {
  cpm_range_thresh <- 0
  log_info("  Single sample detected — cpm_range_thresh set to 0")
}

# ---------------------------------------------------------------------------
# 5. Match-type cumulative stats & summary table
# ---------------------------------------------------------------------------
log_info("Computing match-type statistics ...")
merged[, unique_cum   := cumsum(match_type == "unique")]
merged[, multi_cum    := cumsum(match_type == "multi")]
merged[, unmapped_cum := cumsum(match_type == "unmapped")]

mt_table <- table(merged$match_type)
log_info("  {paste(names(mt_table), mt_table, sep = '=', collapse = ', ')}")

all_totals      <- colSums(merged[, ..sample_cols], na.rm = TRUE)
unmapped_totals <- colSums(merged[match_type == "unmapped", ..sample_cols], na.rm = TRUE)
multi_totals    <- colSums(merged[match_type == "multi",    ..sample_cols], na.rm = TRUE)
unique_totals   <- colSums(merged[match_type == "unique",   ..sample_cols], na.rm = TRUE)

match_type_summary <- data.table(
  sample        = sample_cols,
  total_counts  = all_totals,
  unique_counts = unique_totals,
  multi_counts  = multi_totals,
  unmapped      = unmapped_totals,
  pct_unique    = round(unique_totals   / all_totals * 100, 2),
  pct_multi     = round(multi_totals    / all_totals * 100, 2),
  pct_unmapped  = round(unmapped_totals / all_totals * 100, 2)
)
fwrite(match_type_summary, file.path(output_dir, "match_type_summary.tsv"), sep = "\t")
log_info("  Written match_type_summary.tsv")

# ---------------------------------------------------------------------------
# 6. Unique genome-mapped spacers
# ---------------------------------------------------------------------------
log_info("Filtering to unique '{genome_seqname}' spacers ...")
unique_genome <- merged[match_type == "unique" & seqnames == genome_seqname]
setorder(unique_genome, -max_count, target_seq)
log_info("  {nrow(unique_genome)} unique genome spacers")

# Label intergenic spacers
unique_genome[, gene_label := fifelse(is.na(gene_name), "intergenic", gene_name)]

fwrite(unique_genome, file.path(spacer_dir, "spacers_unique_genome.tsv"), sep = "\t")
log_info("  Written spacer_level/spacers_unique_genome.tsv")

n_intergenic <- unique_genome[is.na(gene_name), .N]
n_genic      <- unique_genome[!is.na(gene_name), .N]
log_info("  Genic: {n_genic} ({round(n_genic/nrow(unique_genome)*100,1)}%)  Intergenic: {n_intergenic}")

# Other seqnames (plasmid etc.)
for (sn in other_seqnames) {
  dt_other <- merged[match_type == "unique" & seqnames == sn]
  pct      <- round(colSums(dt_other[, ..sample_cols], na.rm = TRUE) / all_totals * 100, 2)
  log_info("  {sn}: {nrow(dt_other)} unique spacers | % of total: {paste(names(pct), pct, sep='=', collapse=', ')}")
}

# ---------------------------------------------------------------------------
# 7. Gene-level detection count (spacers with count > 0 per gene per sample)
# ---------------------------------------------------------------------------
log_info("Building gene-level detection count table ...")
gene_counts <-
  unique_genome[!is.na(gene_name),
    lapply(.SD, function(x) sum(x > 0, na.rm = TRUE)),
    .SDcols = sample_cols,
    by = .(seqnames, gene_name)]
setorderv(gene_counts, sample_cols[1], order = -1L)
fwrite(gene_counts, file.path(gene_dir, "gene_counts.tsv"), sep = "\t")
log_info("  {nrow(gene_counts)} genes written to gene_level/gene_counts.tsv")

# ---------------------------------------------------------------------------
# 8. Gene-level raw count sums per sample
# ---------------------------------------------------------------------------
log_info("Building gene-level count sum table ...")
gene_count_sums <-
  unique_genome[!is.na(gene_name),
    lapply(.SD, sum, na.rm = TRUE),
    .SDcols = sample_cols,
    by = .(seqnames, gene_name)]
setorderv(gene_count_sums, sample_cols[1], order = -1L)
fwrite(gene_count_sums, file.path(gene_dir, "gene_count_sums.tsv"), sep = "\t")
log_info("  Written gene_level/gene_count_sums.tsv")

# ---------------------------------------------------------------------------
# 9. CPM normalization (on raw count sums — correct basis for abundance)
# ---------------------------------------------------------------------------
log_info("Computing log2 CPM normalization on raw count sums ...")
gene_matrix           <- as.matrix(gene_count_sums[, ..sample_cols])
rownames(gene_matrix) <- gene_count_sums$gene_name

cpm_matrix <- cpm(gene_matrix, log = TRUE, prior.count = 1)
cpm_dt     <- as.data.table(cpm_matrix, keep.rownames = "gene_name")

stat_cols <- if (n_samples > 1) sample_cols else character(0)
cpm_dt[, min_cpm   := do.call(pmin, c(.SD, na.rm = TRUE)), .SDcols = sample_cols]
cpm_dt[, max_cpm   := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = sample_cols]
cpm_dt[, mean_cpm  := rowMeans(.SD, na.rm = TRUE),         .SDcols = sample_cols]
cpm_dt[, range_cpm := max_cpm - min_cpm]

fwrite(cpm_dt, file.path(gene_dir, "gene_counts_cpm.tsv"), sep = "\t")
log_info("  Written gene_level/gene_counts_cpm.tsv")

cpm_hits <- cpm_dt[max_cpm > cpm_max_thresh & range_cpm > cpm_range_thresh]
fwrite(cpm_hits, file.path(gene_dir, "gene_counts_cpm_hits.tsv"), sep = "\t")
log_info("  CPM hits (max>{cpm_max_thresh}, range>{cpm_range_thresh}): {nrow(cpm_hits)} genes")

# ---------------------------------------------------------------------------
# 10. Plots — spacer level
# ---------------------------------------------------------------------------
log_info("Generating spacer-level plots ...")

MATCH_COLORS <- c(unique = "#2196F3", multi = "#FF9800", unmapped = "#9E9E9E")

# S1 — Spacer abundance distribution by match type
p <- ggplot(merged, aes(x = max_count, fill = match_type)) +
  geom_histogram(bins = 60, alpha = 0.8, position = "identity") +
  scale_x_log10(labels = scales::label_comma()) +
  scale_fill_manual(values = MATCH_COLORS) +
  facet_wrap(~match_type, ncol = 1, scales = "free_y") +
  labs(x = "Max count (log10)", y = "Number of spacers",
       title = "Spacer abundance distribution by match type") +
  theme_bw() + theme(legend.position = "none")
save_plot(p, "s01_spacer_abundance_distribution.png", width = 7, height = 7)

# S2 — Match type composition per sample (stacked %)
comp_long <- match_type_summary[,
  .(sample, unique = pct_unique, multi = pct_multi, unmapped = pct_unmapped)] |>
  melt(id.vars = "sample", variable.name = "match_type", value.name = "pct")
comp_long[, match_type := factor(match_type,
  levels = c("unique", "multi", "unmapped"),
  labels = c("Unique", "Multi-mapped", "Unmapped"))]
comp_long[, sample := factor(sample, levels = rev(sample_cols))]

p <- ggplot(comp_long, aes(x = sample, y = pct, fill = match_type)) +
  geom_col() +
  scale_fill_manual(values = c(Unique = "#2196F3", `Multi-mapped` = "#FF9800", Unmapped = "#9E9E9E")) +
  coord_flip() +
  labs(x = NULL, y = "Percentage (%)", fill = NULL,
       title = "Match type composition per sample") +
  theme_bw() + theme(legend.position = "bottom")
save_plot(p, "s02_match_type_composition.png",
          width = 8, height = max(4, n_samples * 0.35 + 2))

# S3 — Cumulative match types by rank
p <- merged |>
  select(rank, unique_cum, multi_cum, unmapped_cum) |>
  pivot_longer(-rank, names_to = "type", values_to = "cumulative") |>
  mutate(type = recode(type,
    unique_cum = "unique", multi_cum = "multi", unmapped_cum = "unmapped")) |>
  ggplot(aes(rank, cumulative, color = type)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = MATCH_COLORS) +
  scale_x_continuous(labels = scales::label_comma()) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(x = "Rank (by max count, descending)", y = "Cumulative spacers",
       title = "Cumulative match types by rank", color = NULL) +
  theme_bw()
save_plot(p, "s03_cumulative_match_types.png", width = 8, height = 4)

# S4 — Spacer length distribution (if len column present)
if ("len" %in% colnames(merged)) {
  p <- ggplot(merged, aes(x = len, fill = match_type)) +
    geom_histogram(binwidth = 1, alpha = 0.8, position = "identity") +
    scale_fill_manual(values = MATCH_COLORS) +
    facet_wrap(~match_type, ncol = 1, scales = "free_y") +
    labs(x = "Spacer length (bp)", y = "Number of spacers",
         title = "Spacer length distribution by match type") +
    theme_bw() + theme(legend.position = "none")
  save_plot(p, "s04_spacer_length_distribution.png", width = 7, height = 6)
}

# S5 — Strand distribution for unique genome spacers
if ("strand" %in% colnames(unique_genome) && nrow(unique_genome) > 0) {
  strand_tally <- unique_genome[, .N, by = strand]
  strand_tally[, pct := round(N / sum(N) * 100, 1)]
  p <- ggplot(strand_tally, aes(x = strand, y = N, fill = strand)) +
    geom_col(width = 0.5) +
    geom_text(aes(label = paste0(pct, "%")), vjust = -0.4, size = 4) +
    scale_fill_brewer(palette = "Set2") +
    labs(x = "Strand", y = "Number of spacers",
         title = paste("Strand distribution —", genome_seqname, "unique spacers")) +
    theme_bw() + theme(legend.position = "none")
  save_plot(p, "s05_strand_distribution.png", width = 4, height = 4)
}

# S6 — Distance to nearest gene for unique genome spacers
if ("distance" %in% colnames(unique_genome) && nrow(unique_genome) > 0) {
  dist_dt <- unique_genome[!is.na(distance) & distance > 0]
  if (nrow(dist_dt) > 0) {
    p <- ggplot(dist_dt, aes(x = distance + 1)) +
      geom_histogram(bins = 60, fill = "#2196F3", alpha = 0.8) +
      scale_x_log10(labels = scales::label_comma()) +
      labs(x = "Distance to nearest gene + 1 (log10 bp)", y = "Number of spacers",
           title = "Distance to nearest gene (unique genome spacers)") +
      theme_bw()
    save_plot(p, "s06_distance_to_gene.png", width = 7, height = 4)
  }
}

# S7 — Genome coverage: spacer positions along chromosome
if (all(c("start", "strand") %in% colnames(unique_genome)) && nrow(unique_genome) > 0) {
  p <- ggplot(unique_genome, aes(x = start, fill = strand)) +
    geom_histogram(bins = 200, alpha = 0.8, position = "identity") +
    scale_fill_brewer(palette = "Set1") +
    scale_x_continuous(labels = scales::label_comma()) +
    labs(x = paste(genome_seqname, "position (bp)"), y = "Spacer count",
         fill = "Strand", title = "Spacer coverage along chromosome") +
    theme_bw()
  save_plot(p, "s07_genome_coverage.png", width = 10, height = 4)
}

# S8 — Genic vs intergenic fraction by abundance bin
overlap_tally <- unique_genome[,
  .(has_gene = !is.na(gene_name)),
  by = max_count][,
  .(genic = sum(has_gene), total = .N), by = max_count][,
  prop_genic := genic / total]

p <- ggplot(overlap_tally, aes(max_count, prop_genic)) +
  geom_point(alpha = 0.5, size = 1.5, color = "#2196F3") +
  geom_smooth(method = "loess", se = FALSE, color = "firebrick", linewidth = 0.8) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(x = "Max spacer count (log10)", y = "Fraction genic",
       title = "Genic fraction vs. spacer abundance") +
  theme_bw()
save_plot(p, "s08_genic_fraction_vs_abundance.png", width = 7, height = 4)

# ---------------------------------------------------------------------------
# 11. Plots — gene level
# ---------------------------------------------------------------------------
log_info("Generating gene-level plots ...")

# Spacers per gene: sum across samples for ranking
gene_counts[, total_spacers := rowSums(.SD, na.rm = TRUE), .SDcols = sample_cols]

# G1 — Top genes by total spacer detection
top_n <- min(30, nrow(gene_counts))
top_genes <- gene_counts[order(-total_spacers)][seq_len(top_n)]
top_genes[, gene_name := factor(gene_name, levels = rev(gene_name))]

p <- ggplot(top_genes, aes(x = gene_name, y = total_spacers)) +
  geom_col(fill = "#2196F3", alpha = 0.85) +
  coord_flip() +
  labs(x = NULL, y = "Total spacers detected (across samples)",
       title = paste("Top", top_n, "genes by spacer detection")) +
  theme_bw()
save_plot(p, "g01_top_genes_spacer_count.png",
          width = 7, height = max(5, top_n * 0.25 + 2))

# G2 — Distribution of spacer count per gene
p <- ggplot(gene_counts, aes(x = total_spacers + 1)) +
  geom_histogram(bins = 60, fill = "#2196F3", alpha = 0.8) +
  scale_x_log10(labels = scales::label_comma()) +
  labs(x = "Total spacers per gene + 1 (log10)", y = "Number of genes",
       title = "Spacers per gene distribution") +
  theme_bw()
save_plot(p, "g02_spacers_per_gene_distribution.png", width = 6, height = 4)

# G3 — Sample pairwise scatter (all pairs if ≤ 8 samples, else first vs others)
if (n_samples >= 2) {
  scatter_data <- gene_count_sums[, c("gene_name", sample_cols), with = FALSE]

  if (n_samples <= 8) {
    pairs_grid <- as.data.table(t(combn(sample_cols, 2)))
    colnames(pairs_grid) <- c("s1", "s2")
  } else {
    pairs_grid <- data.table(s1 = sample_cols[1], s2 = sample_cols[-1])
    log_info("  >8 samples: showing first sample vs all others in pair scatter")
  }

  plot_list <- lapply(seq_len(nrow(pairs_grid)), function(i) {
    s1 <- pairs_grid[i, s1]; s2 <- pairs_grid[i, s2]
    ggplot(scatter_data, aes(x = .data[[s1]] + 1, y = .data[[s2]] + 1)) +
      geom_point(alpha = 0.3, size = 0.8, color = "#2196F3") +
      scale_x_log10() + scale_y_log10() +
      labs(x = s1, y = s2) +
      theme_bw(base_size = 8)
  })

  ncols  <- min(3, length(plot_list))
  nrows  <- ceiling(length(plot_list) / ncols)
  p <- wrap_plots(plot_list, ncol = ncols) +
    plot_annotation(title = "Gene count sums: pairwise sample scatters (log10)")
  save_plot(p, "g03_sample_scatter_pairs.png",
            width = ncols * 3.5, height = nrows * 3.5)
}

# G4 — Sample correlation heatmap
if (n_samples >= 2) {
  cor_mat <- cor(cpm_matrix, method = "pearson", use = "pairwise.complete.obs")

  cor_long <- as.data.table(as.table(cor_mat))
  colnames(cor_long) <- c("s1", "s2", "r")

  p <- ggplot(cor_long, aes(x = s1, y = s2, fill = r)) +
    geom_tile() +
    geom_text(aes(label = round(r, 2)), size = 2.5) +
    scale_fill_gradient2(low = "#D32F2F", mid = "white", high = "#1976D2",
                         midpoint = 0.9, limits = c(max(0, min(cor_long$r) - 0.05), 1)) +
    labs(x = NULL, y = NULL, fill = "Pearson r",
         title = "Sample correlation (log2 CPM)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p, "g04_sample_correlation_heatmap.png",
            width = max(5, n_samples * 0.7 + 2),
            height = max(4, n_samples * 0.7 + 1.5))
}

# G5 — CPM variability: max_cpm vs range_cpm, hits highlighted
if (nrow(cpm_dt) > 0) {
  cpm_plot_dt <- copy(cpm_dt)
  cpm_plot_dt[, is_hit := max_cpm > cpm_max_thresh & range_cpm > cpm_range_thresh]

  p <- ggplot(cpm_plot_dt, aes(x = mean_cpm, y = range_cpm, color = is_hit)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "#E53935"),
                       labels = c(`FALSE` = "background", `TRUE` = "hit")) +
    geom_hline(yintercept = cpm_range_thresh, linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = cpm_max_thresh,   linetype = "dashed", color = "grey40") +
    labs(x = "Mean log2 CPM", y = "CPM range (max - min)",
         color = NULL, title = "Gene CPM variability",
         subtitle = paste0("Hits: max > ", cpm_max_thresh, ", range > ", cpm_range_thresh)) +
    theme_bw()
  save_plot(p, "g05_cpm_variability.png", width = 6, height = 5)
}

# G6 — CPM heatmap of top variable genes (pheatmap if available, else base)
if (nrow(cpm_hits) > 0 && n_samples >= 2) {
  n_show     <- min(50, nrow(cpm_hits))
  top_hits   <- cpm_hits[order(-range_cpm)][seq_len(n_show)]
  plot_mat   <- as.matrix(top_hits[, ..sample_cols])
  rownames(plot_mat) <- top_hits$gene_name

  has_pheatmap <- requireNamespace("pheatmap", quietly = TRUE)
  if (has_pheatmap) {
    png(file.path(plots_dir, "g06_cpm_heatmap.png"),
        width  = max(800, n_samples * 55),
        height = max(600, n_show * 14 + 120),
        res    = 150)
    pheatmap::pheatmap(plot_mat,
      cluster_rows = TRUE, cluster_cols = TRUE,
      scale = "row",
      color = colorRampPalette(c("#1976D2", "white", "#E53935"))(100),
      fontsize_row = 7, fontsize_col = 8,
      main = paste0("Top ", n_show, " variable genes (log2 CPM, row-scaled)"))
    dev.off()
  } else {
    png(file.path(plots_dir, "g06_cpm_heatmap.png"),
        width = max(800, n_samples * 55), height = max(600, n_show * 14 + 120), res = 150)
    heatmap(plot_mat, scale = "row", margins = c(10, 8),
            main = paste0("Top ", n_show, " variable genes (log2 CPM, row-scaled)"))
    dev.off()
  }
  log_info("  Plot saved: plots/g06_cpm_heatmap.png")
}

# G7 — Sense vs antisense spacers for top genes
if (all(c("strand", "gene_strand") %in% colnames(unique_genome))) {
  sense_dt <- unique_genome[!is.na(gene_name) & !is.na(gene_strand),
    .(n_sense    = sum(strand == gene_strand),
      n_antisense = sum(strand != gene_strand)),
    by = gene_name]
  sense_dt[, total := n_sense + n_antisense]
  sense_dt[, pct_sense := n_sense / total]
  setorder(sense_dt, -total)

  top_s  <- min(30, nrow(sense_dt))
  top_sa <- sense_dt[seq_len(top_s)]
  top_sa[, gene_name := factor(gene_name, levels = rev(gene_name))]

  sa_long <- melt(top_sa[, .(gene_name, Sense = n_sense, Antisense = n_antisense)],
                  id.vars = "gene_name", variable.name = "orientation", value.name = "n")

  p <- ggplot(sa_long, aes(x = gene_name, y = n, fill = orientation)) +
    geom_col(position = "fill") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey30") +
    coord_flip() +
    scale_y_continuous(labels = scales::label_percent()) +
    scale_fill_brewer(palette = "Set2") +
    labs(x = NULL, y = "Proportion of spacers", fill = NULL,
         title = paste("Sense/antisense orientation — top", top_s, "genes")) +
    theme_bw()
  save_plot(p, "g07_sense_antisense.png",
            width = 7, height = max(5, top_s * 0.25 + 2))
}

# ---------------------------------------------------------------------------
# 12. Summary stats to log
# ---------------------------------------------------------------------------
log_info("Summary statistics:")
log_info("  Total spacers (max_count > 0): {nrow(merged)}")
log_info("  Unique {genome_seqname} spacers: {nrow(unique_genome)}")
log_info("  Genes with ≥1 spacer: {nrow(gene_counts)}")
log_info("  Median spacers per gene: {median(gene_counts$total_spacers)}")
log_info("  CPM hits: {nrow(cpm_hits)}")
log_info("  Intergenic spacers: {n_intergenic} ({round(n_intergenic/nrow(unique_genome)*100,1)}%)")
log_info("Pipeline completed successfully.")
