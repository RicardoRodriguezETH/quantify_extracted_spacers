suppressPackageStartupMessages({
  library(tidyverse)
  library(Biostrings)
  library(BiocParallel)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(rtracklayer)
  library(stringdist)
  library(progressr)
  library(logger)
  library(lubridate)
})



convert_coordinates <- function(out) {
  genome_len <- length(genome_seq)
  plasmid_len <- length(plasmid_seq)
  out <-
    out %>%
    mutate(
      rev_start = if_else(seqnames == "MG1655",
        genome_len - end + 1, plasmid_len - end + 1
      ),
      rev_end = if_else(seqnames == "MG1655",
        genome_len - start + 1,
        plasmid_len - start + 1
      )
    ) %>%
    mutate(
      start = if_else(strand == "-", rev_start, start),
      end = if_else(strand == "-", rev_end, end)
    ) %>%
    select(-rev_start, -rev_end)
  return(out)
}

match_sequence <- function(sequence, genome_sequence, min_mismatch = 0, max_mismatch = 2) {
  out <- matchPattern(sequence,
    genome_sequence,
    with.indels = TRUE,
    min.mismatch = min_mismatch,
    max.mismatch = max_mismatch
  )
  out <- as.data.frame(out)
  out
}

match_to_genome <- function(sequences, reference_seq, seqname, strand) {
  # Perform parallel alignment

  start_time <- now()

  matches <-
    bplapply(sequences,
      match_sequence,
      genome_sequence = reference_seq
    )

  matches_dt <-
    enframe(matches) %>%
    unnest(value) %>%
    mutate(
      seqnames = seqname,
      strand = strand
    )

  log_info("Matched spacers to {seqname} {strand} strand in  {time_length(now() - start_time, 'seconds')} seconds")

  return(matches_dt)
}

add_edit_distance <- function(out, sequences) {
  num_mismatch <-
    stringdist::stringdist(as.character(sequences[out$name]),
      out$seq,
      method = "lv"
    )

  out$num_mismatch <- num_mismatch
  return(out)
}

add_id <- function(out) {
  out <-
    out %>%
    group_by(name, seqnames, strand, num_mismatch) %>%
    mutate(n_match = row_number()) %>%
    ungroup() %>%
    mutate(id = paste0(
      "s", name, "_",
      n_match, "_",
      seqnames, "_",
      strand, "_",
      num_mismatch
    ))
  return(out)
}



find_nearest_gene <- function(out, genome_genes) {
  granges <-
    GRanges(
      seqnames = out$seqnames,
      ranges = IRanges(
        start = out$start,
        end = out$end,
        width = out$width,
        names = out$id
      ),
      strand = out$strand
    )
  
  dist_to_nearest <-
    distanceToNearest(granges, genome_genes,
                      ignore.strand = TRUE,
                      select="arbitrary")  %>%
    as_tibble()
  
  
  subject_df <-
    as.data.frame(genome_genes) %>%
    select(start, end, width, strand, Name) %>%
    set_names(c("gene_start",
                "gene_end",
                "gene_width",
                "gene_strand",
                "gene_name")) %>%
    mutate(subjectHits= row_number())
  
  dist_to_nearest <-
    left_join(dist_to_nearest, subject_df, by="subjectHits") %>%
    select(-subjectHits) %>%
    relocate(distance, .after="gene_name")
  
  query_df <- 
    as.data.frame(granges) %>%
    mutate(queryHits= row_number())
  query_df$id <- rownames(query_df)
  
  query_df <-
    left_join(query_df, dist_to_nearest, by="queryHits") %>%
    select(id:distance)
  
  
  out <-
    left_join(out, query_df, by = "id")
  
  return(out)
}

add_alignment_metadata_columns <- function(spacers_dt, out, seq_name, plasmid_name) {
  # Every batch must emit the same 12 category columns (2 seqnames x 2 strands x
  # {0,1,2} mismatches), even when a batch has zero matches in some category —
  # otherwise batches get written with different column counts and appending
  # them to one TSV silently corrupts or truncates the output.
  all_categs <-
    apply(
      expand.grid(seqnames = c(seq_name, plasmid_name), strand = c("+", "-"), num_mismatch = 0:2),
      1, function(x) paste0(x[1], "_", x[2], "_", x[3])
    )

  out_wide <-
    out %>%
    mutate(categ = paste0(seqnames, "_", strand, "_", num_mismatch)) %>%
    group_by(name, categ) %>%
    summarize(num_align = n()) %>%
    pivot_wider(names_from = "categ",
                values_from = "num_align",
                values_fill = 0) %>%
    ungroup()

  missing_categs <- setdiff(all_categs, colnames(out_wide))
  out_wide[missing_categs] <- 0L

  col_names <-
    c(colnames(out_wide)[1], sort(all_categs))

  out_wide <-
    out_wide %>%
    select(!!!col_names)

  log_info("out_wide")


  spacers_dt <-
    left_join(spacers_dt, out_wide, by = "name")

  spacers_dt[is.na(spacers_dt)] <- 0


  spacers_dt <-
    spacers_dt %>%
    mutate(across(!!col_names[-1], ~ if_else(. > 1, 2, .),
                  .names = "categ_{.col}")) %>%
    unite(col = "group", starts_with("categ_"), sep = "") %>%
    mutate(
      num_singles = str_count(group, "1"),
      num_multi = str_count(group, "2")
    ) %>%
    mutate(match_type = case_when(
      num_singles > 1 | num_multi > 0 ~ "multi",
      num_singles == 1 & num_multi == 0 ~ "unique",
      num_singles == 0 & num_multi == 0 ~ "unmapped"
    )) %>%
    select(-num_singles, -num_multi) %>%
    mutate(group= paste0("s", group))


  unique_matched_spacers <-
    spacers_dt %>%
    filter(match_type == "unique") %>%
    pull(name)

  spacers_dt <-
    left_join(spacers_dt,
      out %>%
        filter(name %in% unique_matched_spacers) %>%
        select(-n_match, -id, -num_mismatch),
      by = "name"
    ) %>%
    select(-name)

  return(spacers_dt)
}

match_spacers_in_genome <- function(idx_vector,
                                    genome_seq,
                                    plasmid_seq,
                                    seq_name,
                                    plasmid_name) {
  genome_seq_rc <- reverseComplement(genome_seq)
  plasmid_seq_rc <- reverseComplement(plasmid_seq)

  genome_p <- list(idx_vector, genome_seq, seq_name, "+")
  genome_n <- list(idx_vector, genome_seq_rc, seq_name, "-")
  plasmid_p <- list(idx_vector, plasmid_seq, plasmid_name, "+")
  plasmid_n <- list(idx_vector, plasmid_seq_rc, plasmid_name, "-")

  args_list <- list(genome_p, genome_n, plasmid_p, plasmid_n)

  out <- with_progress({
    p <- progressor(along = args_list)
    lapply(args_list, function(args) {
      p()
      do.call(match_to_genome, args)
    })
  })

  out <- do.call(rbind, out)

  return(out)
}






process_batch <- function(spacers_dt, batch_index,
                           genome_seq, plasmid_seq,
                           genome_genes,
                           seq_name, plasmid_name) {
  sequences <-
    spacers_dt %>%
    pull(target_seq) %>%
    DNAStringSet(.)

  start_time <- now()
  log_info("Batch {batch_index}: aligning {nrow(spacers_dt)} spacers")

  out <- match_spacers_in_genome(
    sequences,
    genome_seq,
    plasmid_seq,
    seq_name,
    plasmid_name
  )

  out <- add_edit_distance(out, sequences)
  out <- add_id(out)
  out <- convert_coordinates(out)
  log_info("Batch {batch_index}: aligned in {time_length(now() - start_time, 'seconds')} seconds")

  out <- find_nearest_gene(out, genome_genes)
  spacers_dt <- add_alignment_metadata_columns(spacers_dt, out, seq_name, plasmid_name)

  return(spacers_dt)
}


args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("No file provided. Usage: Rscript align_spacers.R <input> <genome> <plasmid> <genome_gff> <plasmid_gff> <cores> [<chunk_size>]")
}

input_file              <- args[1]
genome_filepath         <- args[2]
plasmid_genome_filepath <- args[3]
genome_gff_filepath     <- args[4]
plasmid_gff_filepath    <- args[5]
K                       <- as.numeric(args[6])
chunk_size              <- if (length(args) >= 7 && !is.na(args[7])) as.numeric(args[7]) else 100000
output_dir_arg          <- if (length(args) >= 8 && nchar(args[8]) > 0) args[8] else NULL

if (!is.null(output_dir_arg)) {
  spacer_dir   <- file.path(output_dir_arg, "spacer_level")
  logs_dir     <- file.path(output_dir_arg, "logs")
  for (d in c(spacer_dir, logs_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  output_file  <- file.path(spacer_dir, "spacers_aligned.tsv")
  log_filepath <- file.path(logs_dir,   "align_spacers_log.txt")
} else {
  output_file  <- sub("\\.tsv$", "_aligned.tsv", input_file)
  log_filepath <- sub("\\.tsv$", "_log.txt",     input_file)
}

set.seed(100)
register(MulticoreParam(workers = K))

log_appender(appender_file(log_filepath))
log_layout(layout_glue_generator(
  format = "{node}/{pid}/{namespace}/{fn} {time} {level}: {msg}"
))
handlers(global = TRUE)

start_time_0 <- now()
log_info("Pipeline started. Workers: {K}, chunk_size: {chunk_size}")

# Load reference data once
log_info("Reading genome sequences and annotations")
genome_seq  <- readDNAStringSet(genome_filepath)[[1]]
plasmid_seq <- readDNAStringSet(plasmid_genome_filepath)[[1]]

genome_gff  <- import(genome_gff_filepath)
genome_genes <- genome_gff[genome_gff$type == "gene"]
genome_genes <- renameSeqlevels(genome_genes,
  c("MG1655%20(U00096)%20DSMZ" = "MG1655x",
    "MG1655_(U00096)_DSMZ"     = "MG1655"))
genome_genes <- dropSeqlevels(genome_genes, "MG1655x")
genome_genes$Name <- str_remove(genome_genes$Name, " gene")

plasmid_gff   <- import(plasmid_gff_filepath)

seq_name     <- "MG1655"
plasmid_name <- "plasmid"

# Stream input in batches — never load the full file into memory
con <- file(input_file, open = "r")
header <- readLines(con, n = 1)
col_names <- strsplit(header, "\t")[[1]]

header_written <- FALSE
batch_index    <- 1
total_rows     <- 0

repeat {
  raw_lines <- readLines(con, n = chunk_size)
  if (length(raw_lines) == 0) break

  batch_text   <- paste(c(header, raw_lines), collapse = "\n")
  spacers_dt   <- read_tsv(I(batch_text), show_col_types = FALSE) %>%
                    mutate(name = row_number())

  result <- process_batch(spacers_dt, batch_index,
                          genome_seq, plasmid_seq, genome_genes,
                          seq_name, plasmid_name)

  write_tsv(result, output_file,
            append = header_written,
            col_names = !header_written)
  header_written <- TRUE
  total_rows     <- total_rows + nrow(result)

  log_info("Batch {batch_index}: wrote {nrow(result)} rows (total so far: {total_rows})")
  batch_index <- batch_index + 1
}

close(con)
log_info("Pipeline completed. {total_rows} rows written to {output_file}. Runtime: {time_length(now() - start_time_0, 'seconds')} seconds.")
