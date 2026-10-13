# Reconstructs the main logistic models from restricted participant-level inputs.
# The original exploratory scripts remain in the private project; this file is
# organized as a single, non-interactive analysis script.
data_dir <- Sys.getenv("CHIP_LC_DATA_DIR")
output_root <- Sys.getenv("CHIP_LC_OUTPUT_DIR")
if (!nzchar(data_dir) || !nzchar(output_root)) {
  stop("Set CHIP_LC_DATA_DIR and CHIP_LC_OUTPUT_DIR before running.")
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
input <- function(name) {
  path <- file.path(data_dir, name)
  if (!file.exists(path)) stop("Missing restricted input: ", name)
  path
}
required <- function(dat, cols, label) {
  missing <- setdiff(cols, names(dat))
  if (length(missing)) stop(label, " lacks: ", paste(missing, collapse = ", "))
}
estimate <- function(fit, term, analysis) {
  coefs <- summary(fit)$coefficients
  if (!term %in% rownames(coefs)) stop("Missing model term: ", term, " in ", analysis)
  b <- coefs[term, "Estimate"]
  se <- coefs[term, "Std. Error"]
  profile <- tryCatch(suppressMessages(stats::confint(fit, parm = term)),
                      error = function(e) c(NA_real_, NA_real_))
  data.frame(analysis = analysis, term = term, n = stats::nobs(fit),
             OR = exp(b), CI_low_wald = exp(b - 1.96 * se),
             CI_high_wald = exp(b + 1.96 * se),
             CI_low_profile = exp(profile[1]), CI_high_profile = exp(profile[2]),
             P_wald = coefs[term, "Pr(>|z|)"], row.names = NULL)
}

all <- read.csv(input("CHIP_ALL.csv"), stringsAsFactors = FALSE)
variants <- read.csv(input("CH_matched.csv"), stringsAsFactors = FALSE)
required(all, c("Filename", "TYPE", "CHIPstatus", "sex", "age", "Smoke",
                "Family_tumor", "BingLang", "Drink"), "CHIP_ALL.csv")
required(variants, c("Filename", "GENE", "max_vaf_value"), "CH_matched.csv")
if (anyDuplicated(all$Filename)) stop("CHIP_ALL.csv must have one row per participant.")
if (!base::all(all$TYPE %in% c("Non_LC", "LC"))) stop("Unexpected TYPE values.")
if (!base::all(all$CHIPstatus %in% c("Non_CHIP", "CHIP"))) stop("Unexpected CHIPstatus values.")
all$TYPE <- factor(all$TYPE, levels = c("Non_LC", "LC"))
all$CHIPstatus <- factor(all$CHIPstatus, levels = c("Non_CHIP", "CHIP"))
all$sex <- factor(all$sex, levels = c("Male", "Female"))
all$Smoke <- factor(all$Smoke, levels = c("No", "Yes"))
all$Drink <- factor(all$Drink, levels = c("No", "Yes"))
all$BingLang <- factor(all$BingLang, levels = c("No", "Yes"))
all$Family_tumor <- factor(all$Family_tumor, levels = c("No", "Yes"))
all$age <- as.numeric(all$age)
variants$max_vaf_value <- as.numeric(variants$max_vaf_value)

results <- list()
fit_main <- glm(TYPE ~ CHIPstatus + sex + age + Smoke + Family_tumor +
                  BingLang + Drink, data = all, family = binomial())
results[[length(results) + 1L]] <- estimate(fit_main, "CHIPstatusCHIP", "main_CHIP")

# Source Fig. S2A: CHIP carriage as the outcome.
fit_chip <- glm(CHIPstatus ~ TYPE + sex + age + Smoke + Family_tumor +
                  BingLang + Drink, data = all, family = binomial())
results[[length(results) + 1L]] <- estimate(fit_chip, "TYPELC", "CHIP_determinants_LC")
results[[length(results) + 1L]] <- estimate(fit_chip, "age", "CHIP_determinants_age")

# Source stage/histology models. The "No" reference indicator contains
# control participants and cases in non-target/unknown categories; this is
# retained to reproduce the reported estimates, not endorsed as final wording.
stage <- read.csv(input("ALL因素(分期).csv"), stringsAsFactors = FALSE)
required(stage, c("Filename", "cancer_type", "stage"), "ALL因素(分期).csv")
stage <- stage[, c("Filename", "cancer_type", "stage"), drop = FALSE]
if (anyDuplicated(stage$Filename)) stop("Stage file has duplicate Filename values.")
all <- merge(all, stage, by = "Filename", all.x = TRUE, sort = FALSE)
if (anyNA(all$cancer_type) || anyNA(all$stage)) {
  stop("Stage/histology records do not cover all cohort participants.")
}
for (stage_level in c("I", "II", "III", "IV")) {
  all[[paste0("stage_", stage_level)]] <- factor(
    ifelse(all$stage == stage_level, "Yes", "No"), levels = c("No", "Yes"))
}
all$non_small_cell <- factor(
  ifelse(all$cancer_type == "Non-small cell carcinoma", "Yes", "No"),
  levels = c("No", "Yes"))
all$small_cell <- factor(
  ifelse(all$cancer_type == "Small cell carcinoma", "Yes", "No"),
  levels = c("No", "Yes"))
fit_stage <- glm(CHIPstatus ~ stage_I + stage_II + stage_III + stage_IV +
                   sex + age + Smoke, data = all, family = binomial())
results[[length(results) + 1L]] <- estimate(fit_stage, "stage_IVYes", "stage_IV_source_coding")
fit_histology <- glm(CHIPstatus ~ non_small_cell + small_cell + sex + age + Smoke,
                     data = all, family = binomial())
results[[length(results) + 1L]] <- estimate(fit_histology, "non_small_cellYes", "NSCLC_source_coding")
results[[length(results) + 1L]] <- estimate(fit_histology, "small_cellYes", "SCLC_source_coding")

# Match the source script's three-level gene exposure: no called CH, other CH,
# and carrier of the specified CHIP gene. This is NOT identical to contrasting
# carriers with every participant labelled Non_CHIP; see AUDIT_FINDINGS.md.
all_ch_ids <- unique(variants$Filename)
chip_variants <- variants[!is.na(variants$max_vaf_value) &
                            variants$max_vaf_value >= 0.02, , drop = FALSE]
for (gene in c("DNMT3A", "TET2", "PPM1D")) {
  ids <- unique(chip_variants$Filename[chip_variants$GENE == gene])
  all$gene_group <- ifelse(all$Filename %in% ids, "Yes",
                           ifelse(all$Filename %in% all_ch_ids, "Other_CH", "No_CH"))
  all$gene_group <- factor(all$gene_group,
                           levels = c("No_CH", "Other_CH", "Yes"))
  fit <- glm(TYPE ~ gene_group + sex + age + Smoke,
             data = all, family = binomial())
  results[[length(results) + 1L]] <- estimate(fit, "gene_groupYes", paste0(gene, "_CHIP"))
}

# These subgroup models reproduce the covariate sets actually used in the
# exploratory source script, which are narrower than the current figure legend.
all$age_group <- cut(all$age, breaks = c(40, 50, 60, 70, 80),
                     right = FALSE, labels = c("40-49", "50-59", "60-69", "70-79"))
for (group in levels(all$age_group)) {
  sub <- all[!is.na(all$age_group) & all$age_group == group, , drop = FALSE]
  fit <- glm(TYPE ~ CHIPstatus + sex + Smoke, data = sub, family = binomial())
  results[[length(results) + 1L]] <- estimate(fit, "CHIPstatusCHIP", paste0("age_", group))
}
for (group in c("Female", "Male")) {
  sub <- all[!is.na(all$sex) & all$sex == group, , drop = FALSE]
  fit <- glm(TYPE ~ CHIPstatus + age + Smoke, data = sub, family = binomial())
  results[[length(results) + 1L]] <- estimate(fit, "CHIPstatusCHIP", paste0("sex_", group))
}
for (group in c("Yes", "No")) {
  sub <- all[!is.na(all$Smoke) & all$Smoke == group, , drop = FALSE]
  fit <- glm(TYPE ~ CHIPstatus + age + sex, data = sub, family = binomial())
  results[[length(results) + 1L]] <- estimate(fit, "CHIPstatusCHIP", paste0("smoking_", group))
}
out <- do.call(rbind, results)
write.csv(out, file.path(output_root, "core_model_estimates.csv"), row.names = FALSE)
message("Wrote aggregate model estimates. No participant-level rows were exported.")

# Figure 3B-E. Each row of CHIP_matched.csv is one mutation with VAF >= 2%.
# The clone-size percentages in D/E use mutations, not unique participants.
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Install ggplot2 to generate Figure 3B-E and 3G-J.")
}
fig3_dir <- file.path(output_root, "figure3")
dir.create(fig3_dir, recursive = TRUE, showWarnings = FALSE)
chip_mut <- read.csv(input("CHIP_matched.csv"), stringsAsFactors = FALSE)
required(chip_mut, c("GENE", "TYPE", "VAF", "HGVS_P", "HGVS_C"),
         "CHIP_matched.csv")
chip_mut <- chip_mut[chip_mut$GENE %in% c("PPM1D", "TET2"), , drop = FALSE]
chip_mut$VAF <- as.numeric(chip_mut$VAF)
if (!nrow(chip_mut) || anyNA(chip_mut$VAF) ||
    !base::all(chip_mut$TYPE %in% c("LC", "Non_LC")) ||
    !base::all(chip_mut$VAF >= 0.02 & chip_mut$VAF <= 1)) {
  stop("Check the gene, TYPE, and VAF values in CHIP_matched.csv.")
}
chip_mut$clone_size <- ifelse(chip_mut$VAF >= 0.10, "large", "small")
chip_mut$site <- ifelse(!is.na(chip_mut$HGVS_P) & nzchar(chip_mut$HGVS_P),
                        chip_mut$HGVS_P, chip_mut$HGVS_C)
if (anyNA(chip_mut$site) || base::any(!nzchar(chip_mut$site))) {
  stop("PPM1D/TET2 mutations need an HGVS_P or HGVS_C site label.")
}

case_colors <- c("LC" = "#ea5506", "Non_LC" = "#274a78")
clone_colors <- c("large" = "#5688ac", "small" = "#c5e2e7")
clone_counts <- list()
for (gene in c("PPM1D", "TET2")) {
  dat <- chip_mut[chip_mut$GENE == gene, , drop = FALSE]
  if (!nrow(dat)) stop("No CHIP mutations found for ", gene)

  # Protein positions determine the x-axis order; c. annotations follow them.
  protein_pos <- suppressWarnings(as.integer(
    sub("^[^0-9]*([0-9]+).*$", "\\1", dat$site)))
  protein_pos[!grepl("^p\\.", dat$site)] <- NA_integer_
  site_order <- unique(dat$site[order(ifelse(is.na(protein_pos), Inf,
                                          protein_pos), dat$site)])
  dat$site <- factor(dat$site, levels = site_order)
  dat$clone_size <- factor(dat$clone_size, levels = c("small", "large"))
  dat$TYPE <- factor(dat$TYPE, levels = c("LC", "Non_LC"))
  dat$VAF_percent <- 100 * dat$VAF

  panel <- if (gene == "PPM1D") "B" else "C"
  lollipop <- ggplot2::ggplot(dat, ggplot2::aes(x = site, y = VAF_percent,
                                                color = TYPE)) +
    ggplot2::geom_segment(ggplot2::aes(xend = site, yend = 0), linewidth = 0.5) +
    ggplot2::geom_point(ggplot2::aes(shape = clone_size), size = 2.2) +
    ggplot2::scale_color_manual(values = case_colors,
                                labels = c("Cases", "Controls")) +
    ggplot2::scale_shape_manual(values = c("small" = 16, "large" = 17)) +
    ggplot2::labs(title = gene, x = "Mutation site", y = "VAF (%)",
                  color = NULL, shape = "Clone size") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1))
  ggplot2::ggsave(file.path(fig3_dir, paste0("Fig3", panel, "_", gene,
                                             "_mutation_sites.pdf")),
                  lollipop, width = if (gene == "PPM1D") 7 else 14,
                  height = 4.5)

  tab <- as.data.frame(table(TYPE = dat$TYPE,
                             clone_size = factor(dat$clone_size,
                                                 levels = c("large", "small"))))
  names(tab)[3] <- "n"
  tab$percent <- 100 * tab$n / ave(tab$n, tab$TYPE, FUN = sum)
  tab$gene <- gene
  clone_counts[[gene]] <- tab
  panel <- if (gene == "PPM1D") "D" else "E"
  bars <- ggplot2::ggplot(tab, ggplot2::aes(x = TYPE, y = percent,
                                            fill = clone_size)) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::scale_x_discrete(labels = c("Cases", "Controls")) +
    ggplot2::scale_fill_manual(values = clone_colors,
                               breaks = c("small", "large")) +
    ggplot2::scale_y_continuous(limits = c(0, 100),
                                breaks = seq(0, 100, 20)) +
    ggplot2::labs(title = gene, x = NULL, y = "Frequency (%)",
                  fill = "Clone size") +
    ggplot2::theme_classic(base_size = 11)
  ggplot2::ggsave(file.path(fig3_dir, paste0("Fig3", panel, "_", gene,
                                             "_clone_size.pdf")),
                  bars, width = 4, height = 4.5)
}
write.csv(do.call(rbind, clone_counts),
          file.path(fig3_dir, "Fig3DE_clone_size_counts.csv"), row.names = FALSE)

# Figure 3G-J. Include lung cancer cases with stage I-IV only. Gene groups
# are CHIP-positive people with the named mutation in CH_matched.csv; the
# mutation itself need not have VAF >= 2% when another clone establishes
# CHIP status.
stage_levels <- c("I", "II", "III", "IV")
cases <- all[all$TYPE == "LC" & all$stage %in% stage_levels,
             c("Filename", "CHIPstatus", "stage")]
stage_colors <- c("I" = "#e99cad", "II" = "#bd5276",
                  "III" = "#9dc5e2", "IV" = "#5788ac")
chip_cases <- cases[cases$CHIPstatus == "CHIP", , drop = FALSE]
stage_groups <- list(
  G = cases[cases$CHIPstatus == "Non_CHIP", , drop = FALSE],
  H = chip_cases,
  I = chip_cases[chip_cases$Filename %in%
                   variants$Filename[variants$GENE == "PPM1D"], , drop = FALSE],
  J = chip_cases[chip_cases$Filename %in%
                   variants$Filename[variants$GENE == "TET2"], , drop = FALSE]
)
stage_titles <- c(G = "Non-CHIP", H = "CHIP",
                  I = "PPM1D-CHIP", J = "TET2-CHIP")
stage_counts <- list()
for (panel in names(stage_groups)) {
  group <- stage_groups[[panel]]
  if (!nrow(group)) stop("No cases in Figure 3", panel)
  n <- as.integer(table(factor(group$stage, levels = stage_levels)))
  tab <- data.frame(panel = panel, group = stage_titles[[panel]],
                    stage = factor(stage_levels, levels = stage_levels),
                    n = n, percent = 100 * n / sum(n))
  tab$label <- ifelse(tab$n > 0,
                      sprintf("%d (%.1f%%)", tab$n, tab$percent), "")
  stage_counts[[panel]] <- tab
  pie <- ggplot2::ggplot(tab, ggplot2::aes(x = 1, y = n, fill = stage)) +
    ggplot2::geom_col(width = 1, color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = label),
                       position = ggplot2::position_stack(vjust = 0.5),
                       size = 3) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::scale_fill_manual(values = stage_colors, drop = FALSE) +
    ggplot2::labs(title = stage_titles[[panel]], fill = "Stage") +
    ggplot2::theme_void(base_size = 11)
  ggplot2::ggsave(file.path(fig3_dir,
                            paste0("Fig3", panel, "_stage_distribution.pdf")),
                  pie, width = 4, height = 4)
}
stage_counts <- do.call(rbind, stage_counts)
stage_counts$label <- NULL
write.csv(stage_counts, file.path(fig3_dir, "Fig3GJ_stage_counts.csv"),
          row.names = FALSE)
message("Wrote Figure 3B-E and 3G-J panels and aggregate counts.")
