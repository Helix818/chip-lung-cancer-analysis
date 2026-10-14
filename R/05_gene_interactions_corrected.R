# Rebuild gene-specific interaction inputs from each gene's own variant VAF.
# Temporary participant-level inputs are removed when this script exits.
# ==================== 变量定义（需要自定义） ====================
data_dir <- Sys.getenv("CHIP_LC_DATA_DIR")
output_root <- Sys.getenv("CHIP_LC_OUTPUT_DIR")
cohort_file <- Sys.getenv("CHIP_LC_COHORT_FILE",
                          unset = file.path(data_dir, "CHIP_ALL.csv"))
variant_file <- Sys.getenv("CHIP_LC_VARIANT_FILE",
                           unset = file.path(data_dir, "CH_matched.csv"))
snp_file <- Sys.getenv("CHIP_LC_SNP_FILE",
                       unset = file.path(data_dir, "SNP645人加个人信息.csv"))
genes <- c("DNMT3A", "TET2", "PPM1D")
chip_vaf_threshold <- 0.02
expected_cohort_n <- 765L
expected_snp_n <- 645L
expected_cohort_gene_n <- c(DNMT3A = 96L, TET2 = 35L, PPM1D = 21L)
expected_snp_gene_n <- c(DNMT3A = 92L, TET2 = 34L, PPM1D = 19L)
# ================================================================

main <- function() {
  suppressWarnings(try(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"),
                       silent = TRUE))
  if (!nzchar(data_dir) || !nzchar(output_root)) {
    stop("Set CHIP_LC_DATA_DIR and CHIP_LC_OUTPUT_DIR.")
  }
  paths <- c(cohort_file, variant_file, snp_file)
  if (!all(file.exists(paths))) stop("Missing restricted analysis input: ",
                                    paste(paths[!file.exists(paths)], collapse = ", "))
  script_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(script_arg) != 1L) stop("Run with Rscript.")
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg),
                                      winslash = "/"))
  cohort <- read.csv(cohort_file, stringsAsFactors = FALSE)
  variants <- read.csv(variant_file, stringsAsFactors = FALSE)
  snp <- read.csv(snp_file, stringsAsFactors = FALSE)
  check_cols <- function(dat, cols, label) {
    missing <- setdiff(cols, names(dat))
    if (length(missing)) stop(label, " lacks: ", paste(missing, collapse = ", "))
  }
  check_cols(cohort, c("Filename", "TYPE", "CHIPstatus"), "cohort")
  check_cols(variants, c("Filename", "GENE", "VAF"), "variants")
  check_cols(snp, c("Filename", "TYPE", "CHIPstatus"), "SNP cohort")
  if (nrow(cohort) != expected_cohort_n || nrow(snp) != expected_snp_n ||
      anyDuplicated(cohort$Filename) || anyDuplicated(snp$Filename) ||
      !all(snp$Filename %in% cohort$Filename)) {
    stop("Cohort size or Filename linkage differs from validated inputs.")
  }
  variants$VAF <- suppressWarnings(as.numeric(variants$VAF))
  if (anyNA(variants$VAF) || any(variants$VAF < 0 | variants$VAF > 1)) {
    stop("Variant VAF must be numeric, complete, and between 0 and 1.")
  }
  qualified <- variants[variants$VAF >= chip_vaf_threshold, , drop = FALSE]
  chip_ids <- unique(qualified$Filename)
  expected_cohort_chip <- ifelse(cohort$Filename %in% chip_ids, "CHIP", "Non_CHIP")
  expected_snp_chip <- ifelse(snp$Filename %in% chip_ids, "CHIP", "Non_CHIP")
  if (!identical(cohort$CHIPstatus, expected_cohort_chip) ||
      !identical(snp$CHIPstatus, expected_snp_chip) ||
      !identical(snp$TYPE, cohort$TYPE[match(snp$Filename, cohort$Filename)])) {
    stop("Existing CHIPstatus or case/control labels disagree with source records.")
  }
  cohort_gene_n <- vapply(genes, function(g)
    length(unique(qualified$Filename[qualified$GENE == g])), integer(1))
  snp_gene_n <- vapply(genes, function(g)
    sum(snp$Filename %in% qualified$Filename[qualified$GENE == g]), integer(1))
  if (!identical(cohort_gene_n, expected_cohort_gene_n) ||
      !identical(snp_gene_n, expected_snp_gene_n) ||
      sum(expected_cohort_chip == "Non_CHIP") != 582L ||
      sum(expected_snp_chip == "Non_CHIP") != 473L) {
    stop("Corrected gene or No CHIP counts differ from validated analysis.")
  }
  audit <- data.frame(
    gene = genes,
    cohort_gene_carriers = unname(cohort_gene_n),
    snp_gene_carriers = unname(snp_gene_n),
    cohort_no_chip = rep(sum(expected_cohort_chip == "Non_CHIP"), length(genes)),
    snp_no_chip = rep(sum(expected_snp_chip == "Non_CHIP"), length(genes)))

  temporary_data_dir <- tempfile("chip_gene_inputs_")
  dir.create(temporary_data_dir)
  on.exit(unlink(temporary_data_dir, recursive = TRUE, force = TRUE), add = TRUE)
  for (g in genes) {
    ids <- unique(qualified$Filename[qualified$GENE == g])
    d <- snp
    d[[g]] <- ifelse(d$Filename %in% ids, "Yes",
                     ifelse(d$Filename %in% chip_ids, "No", "NCHIP"))
    d <- d[d[[g]] %in% c("NCHIP", "Yes"), , drop = FALSE]
    if (sum(d[[g]] == "Yes") != expected_snp_gene_n[[g]] ||
        sum(d[[g]] == "NCHIP") != 473L) stop("Unexpected model input for ", g)
    write.csv(d, file.path(temporary_data_dir,
                           paste0("CHIP_", g, "_SNP.csv")), row.names = FALSE)
  }
  old_data_dir <- Sys.getenv("CHIP_LC_DATA_DIR", unset = NA_character_)
  on.exit(if (is.na(old_data_dir)) Sys.unsetenv("CHIP_LC_DATA_DIR") else
    Sys.setenv(CHIP_LC_DATA_DIR = old_data_dir), add = TRUE)
  Sys.setenv(CHIP_LC_DATA_DIR = temporary_data_dir)
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  for (g in genes) {
    script <- file.path(script_dir, paste0("04_interaction_", g, ".R"))
    if (!file.exists(script)) stop("Missing original interaction script: ", script)
    message("Running corrected ", g, " model input with original analysis method.")
    sys.source(script, envir = new.env(parent = globalenv()))
  }
  write.csv(audit, file.path(output_root, "gene_input_counts_corrected.csv"),
            row.names = FALSE)
  writeLines(capture.output(sessionInfo()),
             file.path(output_root, "R_session_info_gene_interactions.txt"))
  message("Gene-specific interaction analyses completed: ", output_root)
}
main()
