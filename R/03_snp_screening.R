# Audits the original seven-SNP selection against an omnibus categorical test.
# Requires only base R; no participant-level output is written.
# ==================== 变量定义（需要自定义） ====================
data_dir <- Sys.getenv("CHIP_LC_DATA_DIR")
output_root <- Sys.getenv("CHIP_LC_OUTPUT_DIR")
# ================================================================
suppressWarnings(try(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"),
                     silent = TRUE))
if (!nzchar(data_dir) || !nzchar(output_root)) {
  stop("Set CHIP_LC_DATA_DIR and CHIP_LC_OUTPUT_DIR before running.")
}
path <- file.path(data_dir, "SNP645人加个人信息.csv")
if (!file.exists(path)) stop("Missing restricted input: SNP645人加个人信息.csv")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dat <- read.csv(path, stringsAsFactors = FALSE)
required <- c("TYPE", "CHIPstatus", "sex", "age", "Smoke", "Family_tumor",
              "BingLang", "Drink")
if (length(setdiff(required, names(dat)))) {
  stop("Missing required columns: ", paste(setdiff(required, names(dat)), collapse = ", "))
}
dat$TYPE <- factor(dat$TYPE, levels = c("Non_LC", "LC"))
dat$CHIPstatus <- factor(dat$CHIPstatus, levels = c("Non_CHIP", "CHIP"))
dat$sex <- factor(dat$sex, levels = c("Male", "Female"))
for (v in c("Smoke", "Family_tumor", "BingLang", "Drink")) {
  dat[[v]] <- factor(dat[[v]], levels = c("No", "Yes"))
}
dat$age <- as.numeric(dat$age)
snps <- grep("^rs", names(dat), value = TRUE)
rows <- list()
for (snp in snps) {
  x <- dat[[snp]]
  genotype_n <- sum(!is.na(x))
  levels_present <- sort(unique(stats::na.omit(as.character(x))))
  status <- "analyzed"
  first_p <- omnibus_p <- NA_real_
  n_model <- NA_integer_
  if (genotype_n < nrow(dat) * 0.8) {
    status <- "excluded_gt_20pct_genotype_missing"
  } else if (length(levels_present) < 2L) {
    status <- "excluded_fewer_than_two_genotypes"
  } else {
    dat[[snp]] <- factor(dat[[snp]], levels = levels_present)
    full_formula <- stats::reformulate(c(snp, "CHIPstatus", "sex", "age", "Smoke",
                                         "Family_tumor", "BingLang", "Drink"),
                                       response = "TYPE")
    null_formula <- stats::reformulate(c("CHIPstatus", "sex", "age", "Smoke",
                                         "Family_tumor", "BingLang", "Drink"),
                                       response = "TYPE")
    vars <- all.vars(full_formula)
    fit_data <- dat[stats::complete.cases(dat[, vars, drop = FALSE]), , drop = FALSE]
    n_model <- nrow(fit_data)
    fit <- tryCatch(stats::glm(full_formula, data = fit_data, family = binomial()),
                    error = function(e) NULL)
    null <- tryCatch(stats::glm(null_formula, data = fit_data, family = binomial()),
                     error = function(e) NULL)
    if (is.null(fit) || is.null(null)) {
      status <- "model_failed"
    } else {
      coef_table <- summary(fit)$coefficients
      snp_terms <- grep(paste0("^", snp), rownames(coef_table), value = TRUE)
      if (!length(snp_terms)) {
        status <- "genotype_term_not_estimable"
      } else {
        # Matches original SNPnew.R: first nonreference genotype contrast only.
        first_p <- coef_table[snp_terms[1], "Pr(>|z|)"]
        lrt <- stats::anova(null, fit, test = "Chisq")
        omnibus_p <- lrt[2, "Pr(>Chi)"]
      }
    }
  }
  rows[[length(rows) + 1L]] <- data.frame(
    SNP = snp, genotype_nonmissing_n = genotype_n, model_n = n_model,
    genotype_levels = length(levels_present), status = status,
    first_contrast_P = first_p, omnibus_LRT_P = omnibus_p,
    selected_original_rule = !is.na(first_p) && first_p < 0.1,
    selected_omnibus_rule = !is.na(omnibus_p) && omnibus_p < 0.1)
}
out <- do.call(rbind, rows)
write.csv(out, file.path(output_root, "snp_screening_audit.csv"), row.names = FALSE)
message("Wrote SNP screen audit; original selection uses one genotype coefficient.")
