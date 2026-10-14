
library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(forestploter)
library(grid)
library(stringr)
library(readr)
library(flextable)
library(officer)

base_dir <- Sys.getenv("CHIP_LC_DATA_DIR")
if (!nzchar(base_dir)) stop("Set CHIP_LC_DATA_DIR.")
input_file <- file.path(base_dir, "CHIP_PPM1D_SNP.csv")
output_root <- Sys.getenv("CHIP_LC_OUTPUT_DIR")
if (!nzchar(output_root)) stop("Set CHIP_LC_OUTPUT_DIR.")
out_dir <- file.path(output_root, "interaction_PPM1D")

analysis_label <- "PPM1D-CHIP"
chip_var <- "PPM1D"

chip_no  <- "No"
chip_yes <- "Yes"


snp_list <- c("rs3743078","rs7178270","rs1292037",
              "rs2222722","rs2228001","rs2239611","rs4911154")

snp_levels <- list(
  rs3743078 = c("CC","GC","GG"),
  rs7178270 = c("CC","CG","GG"),
  rs1292037 = c("TT","TC","CC"),
  rs2222722 = c("GG","GA","AA"),
  rs2228001 = c("GG","GT","TT"),
  rs2239611 = c("GG","AG","AA"),
  rs4911154 = c("GG","GA","AA")
)

ref_genotypes <- list(
  rs3743078 = "CC",
  rs7178270 = "CC",
  rs1292037 = "TT",
  rs2222722 = "GG",
  rs2228001 = "GG",
  rs2239611 = "GG",
  rs4911154 = "GG"
)

outcome_var <- "TYPE"
outcome_levels <- c("Non_LC", "LC")

covar_candidates <- c("sex", "age", "Smoke", "Family_tumor", "BingLang", "Drink")

sig_cutoff <- 0.05

min_case_or_control_n   <- 5
extreme_or_cutoff       <- 100
extreme_ci_high_cutoff  <- 1000
extreme_abs_reri_cutoff <- 100

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_file)) stop("Input file not found: ", input_file)

dat <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

required_base_vars <- c(outcome_var, chip_var)
missing_base_vars <- setdiff(required_base_vars, names(dat))
if (length(missing_base_vars) > 0) {
  stop("缺少以下核心变量：", paste(missing_base_vars, collapse = ", "))
}

available_snps <- intersect(snp_list, names(dat))
missing_snps <- setdiff(snp_list, names(dat))
if (length(missing_snps) > 0) {
  warning("以下 SNP 在数据中不存在，已自动跳过：", paste(missing_snps, collapse = ", "))
}
snp_list <- available_snps
if (length(snp_list) == 0) stop("没有可分析的 SNP。")

snp_levels    <- snp_levels[names(snp_levels) %in% snp_list]
ref_genotypes <- ref_genotypes[names(ref_genotypes) %in% snp_list]

standardize_binary_var <- function(x, no_candidates, yes_candidates, var_name) {
  x_chr <- trimws(as.character(x))
  x_chr[x_chr %in% c("", "NA", "NaN")] <- NA_character_
  
  out <- rep(NA_character_, length(x_chr))
  out[tolower(x_chr) %in% tolower(no_candidates)]  <- chip_no
  out[tolower(x_chr) %in% tolower(yes_candidates)] <- chip_yes
  
  if (all(is.na(out))) {
    stop("变量 ", var_name, " 无法识别为二分类暴露。请检查其取值并修改 chip_no_candidates / chip_yes_candidates。")
  }
  
  unmapped_vals <- unique(na.omit(x_chr[is.na(out)]))
  if (length(unmapped_vals) > 0) {
    warning(
      "变量 ", var_name, " 存在未映射值：",
      paste(unmapped_vals, collapse = ", "),
      "；这些值将被设为 NA。"
    )
  }
  
  factor(out, levels = c(chip_no, chip_yes))
}

dat[[chip_var]] <- standardize_binary_var(dat[[chip_var]], chip_no_candidates, chip_yes_candidates, chip_var)

dat[[outcome_var]] <- factor(dat[[outcome_var]], levels = outcome_levels)
if (!all(outcome_levels %in% levels(dat[[outcome_var]]))) {
  stop("结局变量 ", outcome_var, " 的 levels 需要包含：", paste(outcome_levels, collapse = ", "))
}

if ("sex" %in% names(dat))          dat$sex <- factor(dat$sex, levels = c("Male", "Female"))
if ("Smoke" %in% names(dat))        dat$Smoke <- factor(dat$Smoke, levels = c("No", "Yes"))
if ("BingLang" %in% names(dat))     dat$BingLang <- factor(dat$BingLang, levels = c("No", "Yes"))
if ("Family_tumor" %in% names(dat)) dat$Family_tumor <- factor(dat$Family_tumor, levels = c("No", "Yes"))
if ("Drink" %in% names(dat))        dat$Drink <- factor(dat$Drink, levels = c("No", "Yes"))
if ("age" %in% names(dat))          dat$age <- as.numeric(dat$age)

covars <- intersect(covar_candidates, names(dat))

fmt_p <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

fmt_or_ci <- function(est, low, high, ref = FALSE) {
  if (ref) return("1.000 (ref.)")
  if (any(is.na(c(est, low, high)))) return("NE")
  
  if (is.finite(est) && is.finite(low) && is.finite(high)) {
    return(sprintf("%.3f (%.3f–%.3f)", est, low, high))
  }
  if (is.finite(est) && is.finite(low) && !is.finite(high)) {
    return(sprintf("%.3f (%.3f–Inf)", est, low))
  }
  "NE"
}

fmt_num_ci <- function(est, low, high) {
  if (any(is.na(c(est, low, high)))) return("NE")
  if (is.finite(est) && is.finite(low) && is.finite(high)) {
    return(sprintf("%.3f (%.3f–%.3f)", est, low, high))
  }
  if (is.finite(est) && is.finite(low) && !is.finite(high)) {
    return(sprintf("%.3f (%.3f–Inf)", est, low))
  }
  "NE"
}

n_fit <- function(fit) {
  if (is.null(fit)) return(NA_integer_)
  if (!is.null(fit$model)) return(nrow(fit$model))
  if (!is.null(fit$fitted.values)) return(length(fit$fitted.values))
  NA_integer_
}

build_formula <- function(snp, with_interaction = TRUE, covars = NULL) {
  rhs <- if (with_interaction) {
    paste0(snp, " * ", chip_var)
  } else {
    paste0(snp, " + ", chip_var)
  }
  if (!is.null(covars) && length(covars) > 0) {
    rhs <- paste(rhs, paste(covars, collapse = " + "), sep = " + ")
  }
  as.formula(paste(outcome_var, "~", rhs))
}

safe_glm <- function(formula, data) {
  tryCatch(
    suppressWarnings(glm(formula, data = data, family = binomial(link = "logit"))),
    error = function(e) NULL
  )
}

overall_lrt_p <- function(fit0, fit1) {
  if (is.null(fit0) || is.null(fit1)) return(NA_real_)
  tryCatch({
    anova(fit0, fit1, test = "LRT")[2, "Pr(>Chi)"]
  }, error = function(e) NA_real_)
}

make_profile_row <- function(data, vars) {
  nd <- data[1, vars, drop = FALSE]
  for (v in vars) {
    x <- data[[v]]
    if (is.factor(x)) {
      nd[[v]] <- factor(levels(x)[1], levels = levels(x))
    } else if (is.numeric(x)) {
      nd[[v]] <- mean(x, na.rm = TRUE)
    } else {
      nd[[v]] <- x[which(!is.na(x))[1]]
    }
  }
  nd
}

make_nd <- function(data, snp, genotype, chip_status) {
  vars <- unique(c(snp, chip_var, covars))
  nd <- make_profile_row(data, vars)
  nd[[snp]] <- factor(genotype, levels = levels(data[[snp]]))
  nd[[chip_var]] <- factor(chip_status, levels = levels(data[[chip_var]]))
  nd
}

get_mm <- function(fit, nd) {
  tt <- delete.response(terms(fit))
  X <- model.matrix(tt, nd)
  X <- X[, names(coef(fit)), drop = FALSE]
  X
}

make_contrast <- function(fit, nd1, nd0) {
  x1 <- get_mm(fit, nd1)
  x0 <- get_mm(fit, nd0)
  L <- as.numeric(x1 - x0)
  names(L) <- colnames(x1)
  L
}

estimate_from_L <- function(fit, L) {
  b <- coef(fit)
  V <- vcov(fit)
  Lm <- matrix(L, nrow = 1)
  est <- as.numeric(Lm %*% b)
  se  <- suppressWarnings(sqrt(as.numeric(Lm %*% V %*% t(Lm))))
  z   <- ifelse(is.na(se) || se == 0, NA_real_, est / se)
  p   <- ifelse(is.na(z), NA_real_, 2 * pnorm(-abs(z)))
  tibble(
    log_est = est,
    se_log = se,
    z = z,
    p = p,
    OR = exp(est),
    CI_low = exp(est - 1.96 * se),
    CI_high = exp(est + 1.96 * se)
  )
}

collapse_reasons <- function(reasons) {
  reasons <- reasons[!is.na(reasons) & nzchar(reasons)]
  if (length(reasons) == 0) return("")
  paste(unique(reasons), collapse = "; ")
}

flag_from_reasons <- function(reasons) {
  reasons <- reasons[!is.na(reasons) & nzchar(reasons)]
  if (length(reasons) == 0) {
    list(flag = "Stable", reason = "")
  } else {
    list(flag = "Unstable", reason = paste(unique(reasons), collapse = "; "))
  }
}

get_group_counts <- function(data_snp, snp) {
  tab <- table(data_snp[[snp]], data_snp[[chip_var]], data_snp[[outcome_var]], useNA = "no")
  levs <- dimnames(tab)[[1]]
  chips <- dimnames(tab)[[2]]
  
  expand.grid(genotype = levs, CHIP = chips, stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    rowwise() %>%
    mutate(
      case_n    = as.integer(tab[genotype, CHIP, outcome_levels[2]]),
      control_n = as.integer(tab[genotype, CHIP, outcome_levels[1]]),
      total_n   = case_n + control_n
    ) %>%
    ungroup()
}

get_one_count <- function(count_df, genotype, chip_status) {
  out <- count_df %>% filter(genotype == !!genotype, CHIP == !!chip_status)
  if (nrow(out) == 0) {
    tibble(case_n = 0L, control_n = 0L, total_n = 0L)
  } else {
    out %>% select(case_n, control_n, total_n)
  }
}

joint_effects_from_model <- function(fit, data_snp, snp, count_df) {
  if (is.null(fit)) return(tibble())
  
  levs <- levels(data_snp[[snp]])
  ref  <- ref_genotypes[[snp]]
  nd00 <- make_nd(data_snp, snp, ref, chip_no)
  
  grid_df <- expand.grid(
    genotype = levs,
    CHIP = c(chip_no, chip_yes),
    stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    mutate(
      row_order = (match(genotype, levs) - 1L) * 2L + ifelse(CHIP == chip_no, 1L, 2L),
      group_type = case_when(
        genotype == ref & CHIP == chip_no ~ "reference",
        genotype == ref & CHIP == chip_yes ~ "chip_only",
        genotype != ref & CHIP == chip_no ~ "genotype_only",
        TRUE ~ "joint_exposure"
      )
    ) %>%
    left_join(count_df, by = c("genotype", "CHIP"))
  
  map_dfr(seq_len(nrow(grid_df)), function(i) {
    g  <- grid_df$genotype[i]
    ch <- grid_df$CHIP[i]
    
    if (g == ref && ch == chip_no) {
      tibble(
        SNP = snp,
        n = n_fit(fit),
        ref_genotype = ref,
        genotype = g,
        CHIP = ch,
        group_type = grid_df$group_type[i],
        row_order = grid_df$row_order[i],
        case_n = grid_df$case_n[i],
        control_n = grid_df$control_n[i],
        total_n = grid_df$total_n[i],
        log_est = 0,
        se_log = NA_real_,
        z = NA_real_,
        p = NA_real_,
        OR = 1,
        CI_low = 1,
        CI_high = 1,
        OR_95CI = "1.000 (ref.)",
        P_display = "",
        unstable_flag = "Stable",
        unstable_reason = ""
      )
    } else {
      nd <- make_nd(data_snp, snp, g, ch)
      L  <- make_contrast(fit, nd, nd00)
      est <- estimate_from_L(fit, L)
      
      reasons <- c(
        ifelse(grid_df$case_n[i] < min_case_or_control_n,
               paste0("case_n<", min_case_or_control_n), ""),
        ifelse(grid_df$control_n[i] < min_case_or_control_n,
               paste0("control_n<", min_case_or_control_n), ""),
        ifelse(any(!is.finite(c(est$OR, est$CI_low, est$CI_high))),
               "non-finite OR/CI", ""),
        ifelse(is.finite(est$OR) && est$OR > extreme_or_cutoff,
               paste0("OR>", extreme_or_cutoff), ""),
        ifelse(is.finite(est$CI_high) && est$CI_high > extreme_ci_high_cutoff,
               paste0("CI_high>", extreme_ci_high_cutoff), "")
      )
      fr <- flag_from_reasons(reasons)
      
      tibble(
        SNP = snp,
        n = n_fit(fit),
        ref_genotype = ref,
        genotype = g,
        CHIP = ch,
        group_type = grid_df$group_type[i],
        row_order = grid_df$row_order[i],
        case_n = grid_df$case_n[i],
        control_n = grid_df$control_n[i],
        total_n = grid_df$total_n[i],
        log_est = est$log_est,
        se_log = est$se_log,
        z = est$z,
        p = est$p,
        OR = est$OR,
        CI_low = est$CI_low,
        CI_high = est$CI_high,
        OR_95CI = fmt_or_ci(est$OR, est$CI_low, est$CI_high),
        P_display = fmt_p(est$p),
        unstable_flag = fr$flag,
        unstable_reason = fr$reason
      )
    }
  })
}

interaction_metrics_from_model <- function(fit, data_snp, snp, overall_p, count_df) {
  if (is.null(fit)) return(tibble())
  
  levs <- levels(data_snp[[snp]])
  ref  <- ref_genotypes[[snp]]
  nd00 <- make_nd(data_snp, snp, ref, chip_no)
  nd01 <- make_nd(data_snp, snp, ref, chip_yes)
  
  V <- vcov(fit)
  
  map_dfr(levs[levs != ref], function(g) {
    nd10 <- make_nd(data_snp, snp, g, chip_no)
    nd11 <- make_nd(data_snp, snp, g, chip_yes)
    
    L10  <- make_contrast(fit, nd10, nd00)
    L01  <- make_contrast(fit, nd01, nd00)
    L11  <- make_contrast(fit, nd11, nd00)
    Lint <- L11 - L10 - L01
    
    est10  <- estimate_from_L(fit, L10)
    est01  <- estimate_from_L(fit, L01)
    est11  <- estimate_from_L(fit, L11)
    estInt <- estimate_from_L(fit, Lint)
    
    Lmat  <- rbind(L10, L01, L11)
    Sigma <- Lmat %*% V %*% t(Lmat)
    eta   <- c(est10$log_est, est01$log_est, est11$log_est)
    
    OR10 <- exp(eta[1])
    OR01 <- exp(eta[2])
    OR11 <- exp(eta[3])
    
    RERI <- OR11 - OR10 - OR01 + 1
    g_reri <- c(-OR10, -OR01, OR11)
    var_reri <- as.numeric(t(g_reri) %*% Sigma %*% g_reri)
    se_reri <- ifelse(is.na(var_reri) || var_reri < 0, NA_real_, sqrt(var_reri))
    RERI_low  <- ifelse(is.na(se_reri), NA_real_, RERI - 1.96 * se_reri)
    RERI_high <- ifelse(is.na(se_reri), NA_real_, RERI + 1.96 * se_reri)
    p_reri    <- ifelse(is.na(se_reri) || se_reri == 0, NA_real_, 2 * pnorm(-abs(RERI / se_reri)))
    
    AP <- RERI / OR11
    g_ap <- c(
      -OR10 / OR11,
      -OR01 / OR11,
      (OR10 + OR01 - 1) / OR11
    )
    var_ap <- as.numeric(t(g_ap) %*% Sigma %*% g_ap)
    se_ap <- ifelse(is.na(var_ap) || var_ap < 0, NA_real_, sqrt(var_ap))
    AP_low  <- ifelse(is.na(se_ap), NA_real_, AP - 1.96 * se_ap)
    AP_high <- ifelse(is.na(se_ap), NA_real_, AP + 1.96 * se_ap)
    p_ap    <- ifelse(is.na(se_ap) || se_ap == 0, NA_real_, 2 * pnorm(-abs(AP / se_ap)))
    
    c00 <- get_one_count(count_df, ref, chip_no)
    c01 <- get_one_count(count_df, ref, chip_yes)
    c10 <- get_one_count(count_df, g, chip_no)
    c11 <- get_one_count(count_df, g, chip_yes)
    
    min_case_ctrl_cell <- min(c(
      c00$case_n, c00$control_n,
      c01$case_n, c01$control_n,
      c10$case_n, c10$control_n,
      c11$case_n, c11$control_n
    ), na.rm = TRUE)
    
    reasons <- c(
      ifelse(min_case_ctrl_cell < min_case_or_control_n,
             paste0("min case/control cell n<", min_case_or_control_n), ""),
      ifelse(any(!is.finite(c(estInt$OR, estInt$CI_low, estInt$CI_high))),
             "non-finite interaction OR/CI", ""),
      ifelse(is.finite(estInt$OR) && estInt$OR > extreme_or_cutoff,
             paste0("interaction OR>", extreme_or_cutoff), ""),
      ifelse(is.finite(estInt$CI_high) && estInt$CI_high > extreme_ci_high_cutoff,
             paste0("interaction CI_high>", extreme_ci_high_cutoff), ""),
      ifelse(any(!is.finite(c(RERI, RERI_low, RERI_high))),
             "non-finite RERI/CI", ""),
      ifelse(is.finite(RERI) && abs(RERI) > extreme_abs_reri_cutoff,
             paste0("|RERI|>", extreme_abs_reri_cutoff), "")
    )
    fr <- flag_from_reasons(reasons)
    
    tibble(
      SNP = snp,
      n = n_fit(fit),
      ref_genotype = ref,
      genotype = g,
      genotype_contrast = paste0(g, " vs ", ref),
      
      n_case_00 = c00$case_n,
      n_ctrl_00 = c00$control_n,
      n_case_01 = c01$case_n,
      n_ctrl_01 = c01$control_n,
      n_case_10 = c10$case_n,
      n_ctrl_10 = c10$control_n,
      n_case_11 = c11$case_n,
      n_ctrl_11 = c11$control_n,
      min_case_ctrl_cell = min_case_ctrl_cell,
      
      OR10 = OR10,
      OR10_low = est10$CI_low,
      OR10_high = est10$CI_high,
      OR10_95CI = fmt_or_ci(OR10, est10$CI_low, est10$CI_high),
      P_OR10 = est10$p,
      P_OR10_display = fmt_p(est10$p),
      
      OR01 = OR01,
      OR01_low = est01$CI_low,
      OR01_high = est01$CI_high,
      OR01_95CI = fmt_or_ci(OR01, est01$CI_low, est01$CI_high),
      P_OR01 = est01$p,
      P_OR01_display = fmt_p(est01$p),
      
      OR11 = OR11,
      OR11_low = est11$CI_low,
      OR11_high = est11$CI_high,
      OR11_95CI = fmt_or_ci(OR11, est11$CI_low, est11$CI_high),
      P_OR11 = est11$p,
      P_OR11_display = fmt_p(est11$p),
      
      interaction_OR = estInt$OR,
      interaction_low = estInt$CI_low,
      interaction_high = estInt$CI_high,
      interaction_OR_95CI = fmt_or_ci(estInt$OR, estInt$CI_low, estInt$CI_high),
      P_interaction = estInt$p,
      P_interaction_display = fmt_p(estInt$p),
      
      overall_multiplicative_p = overall_p,
      overall_multiplicative_p_display = fmt_p(overall_p),
      
      RERI = RERI,
      RERI_low = RERI_low,
      RERI_high = RERI_high,
      RERI_95CI = fmt_num_ci(RERI, RERI_low, RERI_high),
      P_RERI = p_reri,
      P_RERI_display = fmt_p(p_reri),
      
      AP = AP,
      AP_low = AP_low,
      AP_high = AP_high,
      AP_95CI = fmt_num_ci(AP, AP_low, AP_high),
      P_AP = p_ap,
      P_AP_display = fmt_p(p_ap),
      
      unstable_flag = fr$flag,
      unstable_reason = fr$reason
    )
  })
}

analyze_snp <- function(snp) {
  message("Analyzing ", snp, " ...")
  
  analysis_vars <- unique(c(outcome_var, chip_var, snp, covars))
  data_snp <- dat[complete.cases(dat[, analysis_vars, drop = FALSE]), , drop = FALSE]
  
  if (nrow(data_snp) == 0) {
    return(list(
      coef_table = tibble(),
      counts_table = tibble(),
      joint_table = tibble(),
      interaction_table = tibble()
    ))
  }
  
  data_snp[[snp]] <- factor(data_snp[[snp]], levels = snp_levels[[snp]])
  count_df <- get_group_counts(data_snp, snp) %>% mutate(SNP = snp, .before = 1)
  
  f_red  <- build_formula(snp, with_interaction = FALSE, covars = covars)
  f_full <- build_formula(snp, with_interaction = TRUE,  covars = covars)
  
  fit_red  <- safe_glm(f_red, data_snp)
  fit_full <- safe_glm(f_full, data_snp)
  
  overall_p <- overall_lrt_p(fit_red, fit_full)
  
  coef_table <- if (!is.null(fit_full)) {
    sm <- summary(fit_full)$coefficients
    tibble(
      SNP = snp,
      n = n_fit(fit_full),
      term = rownames(sm),
      beta = sm[, "Estimate"],
      SE = sm[, "Std. Error"],
      z = sm[, "z value"],
      p = sm[, "Pr(>|z|)"],
      OR = exp(beta),
      CI_low = exp(beta - 1.96 * SE),
      CI_high = exp(beta + 1.96 * SE),
      OR_95CI = mapply(fmt_or_ci, OR, CI_low, CI_high),
      P_display = fmt_p(p),
      overall_multiplicative_p = overall_p,
      overall_multiplicative_p_display = fmt_p(overall_p)
    )
  } else {
    tibble()
  }
  
  joint_table <- joint_effects_from_model(fit_full, data_snp, snp, count_df %>% select(-SNP))
  interaction_table <- interaction_metrics_from_model(fit_full, data_snp, snp, overall_p, count_df %>% select(-SNP))
  
  list(
    coef_table = coef_table,
    counts_table = count_df,
    joint_table = joint_table,
    interaction_table = interaction_table
  )
}

all_res <- map(snp_list, analyze_snp)

coef_table_all  <- bind_rows(map(all_res, "coef_table"))
counts_all      <- bind_rows(map(all_res, "counts_table"))
joint_all       <- bind_rows(map(all_res, "joint_table"))
interaction_all <- bind_rows(map(all_res, "interaction_table"))

joint_main_snps <- joint_all %>%
  filter(group_type == "joint_exposure", CHIP == chip_yes, !is.na(p), p < sig_cutoff) %>%
  distinct(SNP) %>%
  pull(SNP)

joint_supp_snps <- setdiff(snp_list, joint_main_snps)

interaction_main_snps <- interaction_all %>%
  filter(!is.na(overall_multiplicative_p), overall_multiplicative_p < sig_cutoff) %>%
  distinct(SNP) %>%
  pull(SNP)

interaction_supp_snps <- setdiff(snp_list, interaction_main_snps)

message("Joint-effect main-figure SNPs: ",
        ifelse(length(joint_main_snps) == 0, "None", paste(joint_main_snps, collapse = ", ")))
message("Interaction main-table SNPs: ",
        ifelse(length(interaction_main_snps) == 0, "None", paste(interaction_main_snps, collapse = ", ")))

write.csv(coef_table_all, file.path(out_dir, "00_model_coefficients_adjusted.csv"), row.names = FALSE, na = "")
write.csv(counts_all,     file.path(out_dir, "01_joint_group_counts_by_case_control.csv"), row.names = FALSE, na = "")
write.csv(joint_all,      file.path(out_dir, "02_joint_effects_adjusted_long.csv"), row.names = FALSE, na = "")
write.csv(interaction_all,file.path(out_dir, "03_interaction_metrics_adjusted_long.csv"), row.names = FALSE, na = "")

joint_table_export <- joint_all %>%
  transmute(
    SNP, ref_genotype, genotype, CHIP, group_type, row_order,
    case_n, control_n, total_n,
    OR, CI_low, CI_high, OR_95CI, p, P_display,
    unstable_flag, unstable_reason
  ) %>%
  arrange(SNP, row_order)
write.csv(joint_table_export, file.path(out_dir, "04_joint_effects_adjusted_table.csv"), row.names = FALSE, na = "")

interaction_all_info <- interaction_all %>%
  arrange(SNP, genotype)
write.csv(interaction_all_info, file.path(out_dir, "05_interaction_all_information.csv"), row.names = FALSE, na = "")

interaction_compact_all <- interaction_all %>%
  transmute(
    SNP,
    genotype_contrast,
    interaction_OR_95CI = ifelse(unstable_flag == "Unstable",
                                 paste0(interaction_OR_95CI, " †"),
                                 interaction_OR_95CI),
    P_interaction_display,
    overall_multiplicative_p_display,
    RERI_95CI = ifelse(unstable_flag == "Unstable",
                       paste0(RERI_95CI, " †"),
                       RERI_95CI),
    P_RERI_display,
    AP_95CI = ifelse(unstable_flag == "Unstable",
                     paste0(AP_95CI, " †"),
                     AP_95CI),
    P_AP_display,
    unstable_flag,
    unstable_reason
  ) %>%
  arrange(SNP, genotype_contrast)

interaction_compact_main <- interaction_compact_all %>%
  filter(SNP %in% interaction_main_snps)

interaction_compact_supp <- interaction_compact_all %>%
  filter(SNP %in% interaction_supp_snps)

write.csv(interaction_compact_all,  file.path(out_dir, "06_interaction_table_compact_all_SNPs.csv"), row.names = FALSE, na = "")
write.csv(interaction_compact_main, file.path(out_dir, "07_interaction_table_compact_main_overall_mult_sig_SNPs.csv"), row.names = FALSE, na = "")
write.csv(interaction_compact_supp, file.path(out_dir, "08_interaction_table_compact_supplementary_SNPs.csv"), row.names = FALSE, na = "")

format_joint_for_forest <- function(joint_df, snps_to_plot) {
  if (length(snps_to_plot) == 0) return(tibble())
  
  df <- joint_df %>%
    filter(SNP %in% snps_to_plot) %>%
    arrange(factor(SNP, levels = snps_to_plot), row_order) %>%
    mutate(
      SNPs = ifelse(duplicated(SNP), "", SNP),
      Genotype = genotype,
      CHIP = ifelse(CHIP == chip_yes, "YES", "NO"),
      OR_95CI_plot = OR_95CI,
      P = P_display,
      is_reference = group_type == "reference"
    ) %>%
    select(SNPs, Genotype, CHIP, OR_95CI_plot, P,
           OR, CI_low, CI_high, is_reference, SNP, row_order)
  
  df$` ` <- paste(rep(" ", 40), collapse = "")
  
  df %>%
    select(SNPs, Genotype, CHIP, OR_95CI_plot, ` `, P,
           OR, CI_low, CI_high, is_reference, SNP, row_order)
}

get_forest_axis <- function(plot_df) {
  finite_upper <- plot_df$CI_high[is.finite(plot_df$CI_high) & !is.na(plot_df$CI_high)]
  upper_raw <- if (length(finite_upper) == 0) 4 else max(finite_upper)
  
  upper_use <- 4
  for (x in c(4, 8, 16, 32, 64)) {
    if (upper_raw <= x) {
      upper_use <- x
      break
    }
    upper_use <- x
  }
  ticks <- c(0.25, 0.5, 1, 2, 4, 8, 16, 32, 64)
  ticks <- ticks[ticks <= upper_use]
  list(xlim = c(0.25, upper_use), ticks = ticks)
}

make_forest_theme_custom <- function() {
  forest_theme(
    base_size = 11,
    refline_col = "#7F7F7F",
    refline_lty = "dashed",
    ci_pch = 15,
    ci_col = "#274a78",
    ci_fill = "#274a78",
    ci_alpha = 1,
    ci_lty = 1,
    ci_lwd = 1.4,
    ci_Theight = 0,
    footnote_cex = 0.8,
    core = list(bg_params = list(fill = c("#eaf4fc", "#b9d7ea"), col = NA)),
    colhead = list(bg_params = list(fill = "#FFFFFF", col = NA))
  )
}

save_forest_plot <- function(plot_obj, pdf_file, png_file, width = 11, height = 8) {
  pdf(pdf_file, width = width, height = height, onefile = FALSE)
  grid.newpage()
  grid.draw(plot_obj)
  dev.off()
  
  png(png_file, width = width, height = height, units = "in", res = 300)
  grid.newpage()
  grid.draw(plot_obj)
  dev.off()
}

make_joint_forest <- function(joint_df, snps_to_plot, title_text, file_prefix) {
  plot_df <- format_joint_for_forest(joint_df, snps_to_plot)
  if (nrow(plot_df) == 0) {
    message("No rows available for ", file_prefix, ". Skipped.")
    return(invisible(NULL))
  }
  
  plot_cap <- 100
  
  plot_df <- plot_df %>%
    mutate(
      OR_plot = OR,
      CI_low_plot = CI_low,
      CI_high_plot = CI_high
    )
  
  plot_df$OR_plot[!is.finite(plot_df$OR_plot) | plot_df$OR_plot <= 0] <- 0.001
  plot_df$CI_low_plot[!is.finite(plot_df$CI_low_plot) | plot_df$CI_low_plot <= 0] <- 0.001
  plot_df$OR_plot[!is.finite(plot_df$OR_plot) | plot_df$OR_plot > plot_cap] <- plot_cap
  plot_df$CI_high_plot[!is.finite(plot_df$CI_high_plot) |
                         plot_df$CI_high_plot <= 0 |
                         plot_df$CI_high_plot > plot_cap] <- plot_cap
  
  axis_info <- get_forest_axis(plot_df %>%
                                 mutate(
                                   OR = OR_plot,
                                   CI_low = CI_low_plot,
                                   CI_high = CI_high_plot
                                 ))
  
  table_text <- plot_df %>%
    select(SNPs, Genotype, CHIP, OR_95CI_plot, ` `, P)
  
  p <- forest(
    table_text,
    est = plot_df$OR_plot,
    lower = plot_df$CI_low_plot,
    upper = plot_df$CI_high_plot,
    ci_column = 5,
    ref_line = 1,
    xlim = axis_info$xlim,
    ticks_at = axis_info$ticks,
    x_trans = "log",
    size = 0.5,
    theme = make_forest_theme_custom()
  )
  
  p <- edit_plot(p, col = 1:6, row = 1, which = "text", gp = gpar(fontface = "bold"))
  p <- add_text(p, text = title_text, row = 0, col = 1:6,
                gp = gpar(fontsize = 12, fontface = "bold"), just = "left")
  
  pdf_file <- file.path(out_dir, paste0(file_prefix, ".pdf"))
  png_file <- file.path(out_dir, paste0(file_prefix, ".png"))
  height_use <- max(4.5, 0.33 * nrow(plot_df) + 1.8)
  save_forest_plot(p, pdf_file, png_file, width = 11, height = height_use)
  invisible(p)
}

clean_p <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- str_replace_all(x, regex("^p\\s*<\\s*", ignore_case = TRUE), "<")
  x
}

clean_ci <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x
}

is_sig_p <- function(x, cutoff = 0.05) {
  x <- clean_p(x)
  out <- rep(FALSE, length(x))
  
  is_lt <- str_detect(x, "^<")
  out[is_lt] <- TRUE
  
  is_num <- !is_lt & !is.na(suppressWarnings(as.numeric(x)))
  out[is_num] <- suppressWarnings(as.numeric(x[is_num]) < cutoff)
  
  out[is.na(out)] <- FALSE
  out
}

make_empty_docx <- function(path, title_text, body_text) {
  doc <- read_docx()
  doc <- body_add_par(doc, title_text, style = "heading 1")
  doc <- body_add_par(doc, body_text, style = "Normal")
  print(doc, target = path)
}

make_interaction_table_docx <- function(compact_df, caption_text, out_path, empty_text) {
  if (nrow(compact_df) == 0) {
    make_empty_docx(
      path = out_path,
      title_text = caption_text,
      body_text = empty_text
    )
    return(invisible(NULL))
  }
  
  tbl_long <- bind_rows(
    compact_df %>%
      mutate(row_id = row_number()) %>%
      transmute(
        row_id,
        SNP_full = SNP,
        Contrast_full = genotype_contrast,
        Metric = "Multiplicative OR",
        `Estimate (95% CI)` = clean_ci(interaction_OR_95CI),
        P = clean_p(P_interaction_display),
        `Overall P(mult.)` = clean_p(overall_multiplicative_p_display)
      ),
    compact_df %>%
      mutate(row_id = row_number()) %>%
      transmute(
        row_id,
        SNP_full = SNP,
        Contrast_full = genotype_contrast,
        Metric = "RERI",
        `Estimate (95% CI)` = clean_ci(RERI_95CI),
        P = clean_p(P_RERI_display),
        `Overall P(mult.)` = ""
      ),
    compact_df %>%
      mutate(row_id = row_number()) %>%
      transmute(
        row_id,
        SNP_full = SNP,
        Contrast_full = genotype_contrast,
        Metric = "AP",
        `Estimate (95% CI)` = clean_ci(AP_95CI),
        P = clean_p(P_AP_display),
        `Overall P(mult.)` = ""
      )
  ) %>%
    mutate(
      Metric = factor(Metric, levels = c("Multiplicative OR", "RERI", "AP"))
    ) %>%
    arrange(factor(SNP_full, levels = unique(compact_df$SNP)), row_id, Metric) %>%
    group_by(SNP_full) %>%
    mutate(
      `Overall P(mult.)` = if_else(
        row_number() == 1,
        first(`Overall P(mult.)`[`Overall P(mult.)` != ""]),
        ""
      )
    ) %>%
    ungroup() %>%
    group_by(SNP_full) %>%
    mutate(SNP = if_else(row_number() == 1, SNP_full, "")) %>%
    ungroup() %>%
    group_by(SNP_full, Contrast_full) %>%
    mutate(Contrast = if_else(row_number() == 1, Contrast_full, "")) %>%
    ungroup() %>%
    select(SNP, Contrast, Metric, `Estimate (95% CI)`, P, `Overall P(mult.)`)
  
  ft <- flextable(tbl_long)
  
  ft <- set_header_labels(
    ft,
    SNP = "SNP",
    Contrast = "Contrast",
    Metric = "Metric",
    `Estimate (95% CI)` = "Estimate (95% CI)",
    P = "P",
    `Overall P(mult.)` = "Overall P(mult.)"
  )
  
  ft <- align(ft, j = c("SNP", "Contrast", "Metric"), align = "left", part = "all")
  ft <- align(ft, j = c("Estimate (95% CI)", "P", "Overall P(mult.)"), align = "center", part = "all")
  
  ft <- font(ft, fontname = "Times New Roman", part = "all")
  ft <- fontsize(ft, size = 10.5, part = "all")
  ft <- padding(ft, padding = 4, part = "all")
  ft <- bold(ft, part = "header")
  
  big_border   <- fp_border(color = "black", width = 1.2)
  small_border <- fp_border(color = "black", width = 0.6)
  
  ft <- border_remove(ft)
  ft <- hline_top(ft, border = big_border, part = "header")
  ft <- hline(ft, border = small_border, part = "header")
  
  block_end_rows <- seq(3, nrow(tbl_long), by = 3)
  if (length(block_end_rows) > 1) {
    ft <- hline(ft, i = block_end_rows[-length(block_end_rows)], border = small_border, part = "body")
  }
  ft <- hline_bottom(ft, border = big_border, part = "body")
  
  sig_p_body <- is_sig_p(tbl_long$P)
  sig_p_overall <- is_sig_p(tbl_long$`Overall P(mult.)`)
  
  if (any(sig_p_body)) {
    ft <- bold(ft, i = which(sig_p_body), j = "P", bold = TRUE, part = "body")
  }
  if (any(sig_p_overall)) {
    ft <- bold(ft, i = which(sig_p_overall), j = "Overall P(mult.)", bold = TRUE, part = "body")
  }
  
  ft <- width(ft, j = "SNP", width = 0.95)
  ft <- width(ft, j = "Contrast", width = 1.00)
  ft <- width(ft, j = "Metric", width = 1.25)
  ft <- width(ft, j = "Estimate (95% CI)", width = 2.80)
  ft <- width(ft, j = "P", width = 0.80)
  ft <- width(ft, j = "Overall P(mult.)", width = 1.10)
  
  ft <- set_caption(ft, caption = caption_text)
  
  ft <- add_footer_lines(
    ft,
    values = c(
      paste0("Multiplicative OR represents the genotype-specific SNP × ", analysis_label, " interaction term from logistic regression."),
      "Overall P(mult.) was obtained by likelihood ratio test comparing models with and without the interaction term.",
      "RERI and AP reflect additive interaction on the OR scale approximation. P < 0.05 is shown in bold.",
      "† Unstable estimate (sparse data and/or extreme/non-finite interval); interpret cautiously."
    )
  )
  
  ft <- fontsize(ft, size = 9, part = "footer")
  ft <- autofit(ft)
  
  save_as_docx(`Interaction table` = ft, path = out_path)
  invisible(ft)
}

make_joint_forest(
  joint_df = joint_all,
  snps_to_plot = snp_list,
  title_text = paste0("Figure 1. Adjusted joint effects of SNP genotype and ", analysis_label, " on lung cancer risk (all SNPs)"),
  file_prefix = "Figure1_joint_effects_all_SNPs"
)

make_joint_forest(
  joint_df = joint_all,
  snps_to_plot = joint_main_snps,
  title_text = paste0("Figure 2. Adjusted joint effects of SNP genotype and ", analysis_label, " on lung cancer risk (main-figure SNPs)"),
  file_prefix = "Figure2_joint_effects_main_figure_SNPs"
)

make_joint_forest(
  joint_df = joint_all,
  snps_to_plot = joint_supp_snps,
  title_text = paste0("Figure 3. Adjusted joint effects of SNP genotype and ", analysis_label, " on lung cancer risk (supplementary SNPs)"),
  file_prefix = "Figure3_joint_effects_supplementary_SNPs"
)

write.csv(
  tibble(joint_main_figure_snps = joint_main_snps),
  file.path(out_dir, "09_joint_main_figure_SNPs.csv"),
  row.names = FALSE, na = ""
)

write.csv(
  tibble(joint_supplementary_snps = joint_supp_snps),
  file.path(out_dir, "10_joint_supplementary_figure_SNPs.csv"),
  row.names = FALSE, na = ""
)

write.csv(
  tibble(interaction_main_table_snps = interaction_main_snps),
  file.path(out_dir, "11_interaction_main_table_SNPs.csv"),
  row.names = FALSE, na = ""
)

write.csv(
  tibble(interaction_supplementary_table_snps = interaction_supp_snps),
  file.path(out_dir, "12_interaction_supplementary_table_SNPs.csv"),
  row.names = FALSE, na = ""
)

make_interaction_table_docx(
  compact_df = interaction_compact_main,
  caption_text = paste0("Table. Multiplicative and additive interactions between SNPs and ", analysis_label, " in relation to lung cancer risk"),
  out_path = file.path(out_dir, "13_interaction_table_main_text.docx"),
  empty_text = "No SNP met the prespecified criterion of overall multiplicative interaction P < 0.05."
)

make_interaction_table_docx(
  compact_df = interaction_compact_supp,
  caption_text = paste0("Supplementary Table. Multiplicative and additive interactions between SNPs and ", analysis_label, " in relation to lung cancer risk"),
  out_path = file.path(out_dir, "14_interaction_table_supplementary.docx"),
  empty_text = "All analyzed SNPs met the prespecified criterion of overall multiplicative interaction P < 0.05; therefore no supplementary interaction table was generated."
)


