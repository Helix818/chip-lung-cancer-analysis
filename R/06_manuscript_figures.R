result_dir <- Sys.getenv("CHIP_LC_OUTPUT_DIR")
out_dir <- file.path(result_dir, "figures_matching_submission_style")
if (!nzchar(result_dir)) stop("Set CHIP_LC_OUTPUT_DIR.")
suppressWarnings(try(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"), silent=TRUE))
library(forestploter)
library(grid)
if (dir.exists(out_dir)) stop("Refusing to overwrite existing figure directory: ", out_dir)
dir.create(out_dir)

matching_theme <- function(footnote=FALSE) {
  forest_theme(
    base_size=11,
    ci_pch=15, ci_col="#274a78", ci_fill="#274a78", ci_lwd=1.4,
    ci_Theight=0,
    xaxis_gp=gpar(fontsize=9),
    refline_gp=gpar(lwd=1, lty="dashed", col="#7F7F7F"),
    core=list(bg_params=list(fill=c("#eaf4fc", "#b9d7ea"), col=NA)),
    colhead=list(bg_params=list(fill="#FFFFFF", col=NA)),
    footnote_gp=gpar(fontsize=8, col="#444444"),
    footnote_parse=FALSE
  )
}

save_forest <- function(plot_obj, stem, width, height) {
  pdf_file <- file.path(out_dir, paste0(stem, ".pdf"))
  png_file <- file.path(out_dir, paste0(stem, ".png"))
  pdf(pdf_file, width=width, height=height, onefile=FALSE, useDingbats=FALSE)
  grid.newpage(); grid.draw(plot_obj); dev.off()
  png(png_file, width=width, height=height, units="in", res=300, bg="white")
  grid.newpage(); grid.draw(plot_obj); dev.off()
}

gene <- read.csv(file.path(result_dir, "Fig3A_gene_OR_corrected.csv"),
                 stringsAsFactors=FALSE)
gene <- gene[match(c("DNMT3A", "TET2", "PPM1D"), gene$Gene), , drop=FALSE]
stopifnot(nrow(gene)==3L, !anyNA(gene$Gene))
gene_text <- data.frame(
  "CHIP Genes"=gene$Gene,
  "OR (95% CI)"=sprintf("%.2f (%.2f–%.2f)", gene$OR,
                         gene$CI_lower, gene$CI_upper),
  " "=rep(strrep(" ", 35), nrow(gene)),
  "p value"=paste0("  p = ", sprintf("%.3f", gene$P)),
  check.names=FALSE
)
fig3a <- forest(gene_text, est=gene$OR, lower=gene$CI_lower,
  upper=gene$CI_upper, sizes=0.5, ci_column=3, ref_line=1,
  xlim=c(0,7), ticks_at=0:7, x_trans="none",
  xlab="Odds Ratio (95% CI)", theme=matching_theme())
fig3a <- edit_plot(fig3a, col=1, row=1:3, which="text",
                   gp=gpar(fontface="italic"))
save_forest(fig3a, "Fig3A_matching_submission_style", 10, 2.7)

snps <- c("rs3743078", "rs7178270", "rs1292037", "rs2222722",
          "rs2228001", "rs2239611", "rs4911154")
read_joint <- function(gene_name) {
  p <- file.path(result_dir, paste0("interaction_", gene_name),
                 "04_joint_effects_adjusted_table.csv")
  read.csv(p, stringsAsFactors=FALSE, check.names=FALSE)
}

format_p <- function(x) {
  ifelse(is.na(x) | !nzchar(x), "",
         ifelse(grepl("^<", x), paste0("  p ", x), paste0("  p = ", x)))
}

draw_joint <- function(gene_name, selected_snps, stem) {
  d <- read_joint(gene_name)
  d <- d[d$SNP %in% selected_snps, , drop=FALSE]
  d <- d[order(match(d$SNP, selected_snps), d$row_order), , drop=FALSE]
  stopifnot(nrow(d)==6L*length(selected_snps),
            all(selected_snps %in% d$SNP))
  invalid <- is.na(d$OR) | !is.finite(d$OR) | d$OR<=0 | d$OR>100 |
             is.na(d$CI_low) | !is.finite(d$CI_low) | d$CI_low<=0 |
             is.na(d$CI_high) | !is.finite(d$CI_high) | d$CI_high>1000
  sparse <- d$unstable_flag=="Unstable"
  or_text <- d$OR_95CI
  or_text[invalid] <- "NE"
  or_text[sparse] <- paste0(or_text[sparse], " †")
  p_text <- format_p(d$P_display)
  p_text[invalid] <- "NE"
  displayed_snp <- ifelse(duplicated(d$SNP), "", d$SNP)
  table_text <- data.frame(
    "SNPs"=displayed_snp,
    "Genotype"=d$genotype,
    "CHIP"=ifelse(d$CHIP=="Yes", "YES", "NO"),
    "OR (95% CI)"=or_text,
    " "=rep(strrep(" ", if(length(selected_snps)==1L) 40 else 65), nrow(d)),
    "p value"=p_text,
    check.names=FALSE
  )
  names(table_text)[3] <- paste0(gene_name, "-CHIP")
  d$OR[invalid] <- NA_real_
  d$CI_low[invalid] <- NA_real_
  d$CI_high[invalid] <- NA_real_
  finite_upper <- d$CI_high[is.finite(d$CI_high)]
  upper_raw <- if(length(finite_upper)) max(finite_upper) else 4
  axis_upper <- 64
  for(candidate in c(4,8,16,32,64)) {
    if(upper_raw<=candidate) { axis_upper <- candidate; break }
  }
  ticks <- c(.25,.5,1,2,4,8,16,32,64)
  ticks <- ticks[ticks<=axis_upper]
  footnote <- "† Sparse or unstable estimate. NE: extreme/non-finite estimate; see full model CSV."
  p <- forest(table_text, est=d$OR, lower=d$CI_low, upper=d$CI_high,
    sizes=0.5, ci_column=5, ref_line=1,
    xlim=c(.25,axis_upper), ticks_at=ticks, x_trans="log",
    xlab="Odds Ratio (95% CI)", footnote=footnote,
    theme=matching_theme(footnote=TRUE))
  height <- if(length(selected_snps)==1L) 3.1 else 0.32*nrow(d)+1.7
  save_forest(p, stem, 11, height)
}


