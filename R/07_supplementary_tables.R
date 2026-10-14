res_dir <- Sys.getenv("CHIP_LC_OUTPUT_DIR")
out_xlsx <- file.path(res_dir,"Supplementary_Tables_S3-S5_corrected.xlsx")
out_table2 <- file.path(res_dir,"Table2_DNMT3A_replacement_rows.csv")
if (!nzchar(res_dir)) stop("Set CHIP_LC_OUTPUT_DIR.")
suppressWarnings(try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE))
if(file.exists(out_xlsx) || file.exists(out_table2)) stop("Refusing to overwrite a corrected output.")
spec <- data.frame(Gene=c("DNMT3A","TET2","PPM1D"),
                   Sheet=c("Table S3","Table S4","Table S5"),
                   stringsAsFactors=FALSE)
fmt <- function(g, supplementary=TRUE) {
  base <- file.path(res_dir,paste0("interaction_",g))
  fname <- if(supplementary) "08_interaction_table_compact_supplementary_SNPs.csv" else "06_interaction_table_compact_all_SNPs.csv"
  x <- read.csv(file.path(base,fname),stringsAsFactors=FALSE,check.names=FALSE)
  if(!supplementary) x <- x[x$SNP=="rs7178270",,drop=FALSE]
  ans <- vector("list",nrow(x)*3L)
  k <- 0L
  for(i in seq_len(nrow(x))) {
    snp_label <- if(i==1L || x$SNP[i]!=x$SNP[i-1L]) x$SNP[i] else ""
    vals <- c(x$interaction_OR_95CI[i],x$RERI_95CI[i],x$AP_95CI[i])
    ps <- c(x$P_interaction_display[i],x$P_RERI_display[i],x$P_AP_display[i])
    for(j in seq_len(3L)) {
      k <- k+1L
      ans[[k]] <- data.frame(SNP=if(j==1L) snp_label else "",
        Contrast=if(j==1L) x$genotype_contrast[i] else "",
        Metric=c("Multiplicative OR","RERI","AP")[j],
        Estimate=vals[j],P=ps[j],
        Overall=if(j==1L && nzchar(snp_label)) x$overall_multiplicative_p_display[i] else "",
        stringsAsFactors=FALSE)
    }
  }
  y <- do.call(rbind,ans)
  names(y) <- c("SNP","Contrast","Metric","Estimate (95% CI)","P","Overall P(mult.)")
  y
}
sheets <- setNames(lapply(spec$Gene,fmt,supplementary=TRUE),spec$Sheet)
sheets[["Notes"]] <- data.frame(Note=c(
  "Dagger marks an unstable estimate (sparse case/control cells, extreme estimate, or non-finite interval); see the full model CSV for its reason.",
  "OR-scale RERI/AP use the original delta-method calculation.",
  "Models adjust for sex, age, smoking, family history, areca use, and alcohol consumption.",
  "Table S3 excludes rs7178270, which remains in main Table 2 under the original overall interaction P < 0.05 selection rule."))
writexl::write_xlsx(sheets,out_xlsx)
write.csv(fmt("DNMT3A",FALSE),out_table2,row.names=FALSE)
