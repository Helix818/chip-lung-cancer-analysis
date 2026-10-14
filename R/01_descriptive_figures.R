# Descriptive analyses from CHIPanalysis/CHIP_LC.R, with portable I/O.
# Panel filenames reflect the historical analysis script.
# Input tables remain outside this repository.
# ==================== 变量定义（需要自定义） ====================
data_dir <- Sys.getenv("CHIP_LC_DATA_DIR")
output_root <- Sys.getenv("CHIP_LC_OUTPUT_DIR")
# ================================================================
suppressWarnings(try(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"),
                     silent = TRUE))

library(dplyr)
library(tidyr)
library(ggplot2)
if (!nzchar(data_dir) || !nzchar(output_root)) stop("Set CHIP_LC_DATA_DIR and CHIP_LC_OUTPUT_DIR.")
out_dir <- file.path(output_root, "descriptive_figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(out_dir)
CHIP_matched <- read.csv(file.path(data_dir, "CHIP_matched.csv"), header = TRUE)
CHIP_ALL <- read.csv(file.path(data_dir, "CHIP_ALL.csv"), header = T)

CHIPALL11table <- table(CHIP_ALL$TYPE, CHIP_ALL$CHIPstatus)
CHIPALL11table <- as.matrix(CHIPALL11table)

# 运行卡方检验
chi_test <- chisq.test(CHIPALL11table)
print(chi_test)

library(tidyverse)

# 按 TYPE 和 CHIPstatus 分组计数并计算百分比
CHIPALL11_prop <- CHIP_ALL %>%
  group_by(TYPE, CHIPstatus) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(TYPE) %>%
  mutate(percent = n / sum(n) * 100)

# 新增组合填色变量：CHIPstatus.TYPE
CHIPALL11_prop <- CHIPALL11_prop %>%
  mutate(fill_group = paste(CHIPstatus, TYPE, sep = "."))

# 绘制堆叠柱状图并添加 P 值
p_value <- chi_test$p.value
p_label <- paste0("p = ", formatC(p_value, format = "e", digits = 2))

ggplot(
  CHIPALL11_prop,
  aes(
    x = TYPE,
    y = percent,
    fill = factor(fill_group, levels = c("Non_CHIP.LC", "CHIP.LC", "Non_CHIP.Non_LC", "CHIP.Non_LC"))
  )
) +
  geom_bar(
    stat = "identity", position = "stack", width = 0.6, linewidth = 0.3
  ) +
  scale_fill_manual(
    values = c(
      "CHIP.LC"         = "#ea5506",   # LC组  CHIP
      "Non_CHIP.LC"     = "#fff1cf",   # LC组  Non_CHIP
      "CHIP.Non_LC"     = "#274a78",   # Non_LC组 CHIP
      "Non_CHIP.Non_LC" = "#eaf4fc"    # Non_LC组 Non_CHIP
    ),
    labels = c(
      "CHIP.LC"         = "CHIP",
      "Non_CHIP.LC"     = "Non_CHIP",
      "CHIP.Non_LC"     = "CHIP",
      "Non_CHIP.Non_LC" = "Non_CHIP"
    ),
    name = "Category"
  ) +
  scale_y_continuous(
    limits = c(0, 110),
    breaks = seq(0, 100, 25),
    expand = c(0, 0)
  ) +
  labs(
    title = "Frequency of CHIP and Non-CHIP",
    x = "Type",
    y = "Frequency (%)",
    fill = "Category"
  ) +
  annotate(
    "text",
    x = 1.5,
    y = 105,
    label = p_label,
    size = 5,
    fontface = "bold",
    color = "black"
  ) +
  theme_classic() +
  theme(
    text = element_text(color = "black"),
    plot.title = element_text(
      size = 16, face = "bold", hjust = 0.5, color = "black"
    ),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    legend.title = element_text(size = 13, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank()
  )

ggsave("CHIP人数占比带P值.pdf",width = 5,height = 5)


#### ===================================================================== ####
####            14 张图（CHIP, VAF >= 2%）  Fig1 / FigS1 / FigS2             ####
####  分子来自 CHIP_matched（携带者/突变人数），分母来自 CHIP_ALL（队列总人数）   ####
####  配色：LC=Case=#ea5506 ，Non_LC=Control=#274a78                          ####
#### ===================================================================== ####
library(tidyr)
library(scales)

## ---- 通用约定 ----
grp_col    <- c("LC" = "#ea5506", "Non_LC" = "#274a78")
age_breaks <- c(39, 49, 59, 69, 79)
age_labs   <- c("40-49", "50-59", "60-69", "70-79")

base_theme <- theme_classic() +
  theme(
    text            = element_text(color = "black"),
    plot.title      = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title      = element_text(size = 13, color = "black"),
    axis.text       = element_text(size = 11, color = "black"),
    axis.line       = element_line(color = "black", linewidth = 0.7),
    axis.ticks      = element_line(color = "black", linewidth = 0.7),
    legend.title    = element_blank(),
    legend.text     = element_text(size = 11, color = "black"),
    legend.position = "top",
    panel.grid      = element_blank()
  )

## ---- 数据预处理 ----
CHIP_matched$TYPE    <- factor(CHIP_matched$TYPE, levels = c("LC", "Non_LC"))
CHIP_ALL$TYPE        <- factor(CHIP_ALL$TYPE,     levels = c("LC", "Non_LC"))
CHIP_ALL$age_grp     <- cut(CHIP_ALL$age,     breaks = age_breaks, labels = age_labs)
CHIP_matched$age_grp <- cut(CHIP_matched$age, breaks = age_breaks, labels = age_labs)

# 每个 CHIP 携带者一行（个体级聚合）
indiv <- CHIP_matched %>%
  group_by(Filename) %>%
  summarise(TYPE    = first(TYPE),
            age     = first(age),
            age_grp = first(age_grp),
            mut_count = n(),
            max_vaf   = max(VAF),
            .groups = "drop")

# 各 (年龄段 × 组) 的队列总人数（频率分母）
denom <- CHIP_ALL %>% count(age_grp, TYPE, name = "total")

# 基因按携带者人数排序（多图复用）
gene_indiv <- CHIP_matched %>% distinct(Filename, GENE, TYPE)
gene_rank  <- gene_indiv %>% count(GENE, sort = TRUE)
gene_order <- gene_rank$GENE
nLC  <- sum(CHIP_ALL$TYPE == "LC")
nNon <- sum(CHIP_ALL$TYPE == "Non_LC")


#### ---- Task 1  Fig1B：CHIP(VAF>=2%) 患病率随年龄段 ----
carr2 <- CHIP_ALL %>% filter(CHIPstatus == "CHIP") %>% count(age_grp, TYPE, name = "n_chip")
fig1b <- denom %>%
  left_join(carr2, by = c("age_grp", "TYPE")) %>%
  mutate(n_chip = ifelse(is.na(n_chip), 0, n_chip),
         freq   = n_chip / total * 100) %>%
  filter(!is.na(age_grp))

p1b <- ggplot(fig1b, aes(age_grp, freq, color = TYPE, group = TYPE)) +
  geom_line(linewidth = 1) + geom_point(size = 2.2) +
  scale_color_manual(values = grp_col) +
  scale_y_continuous(limits = c(0, 60), breaks = seq(0, 60, 20)) +
  labs(x = "Age Group", y = "Frequency (%)", title = "CHIP (VAF ≥ 2%)") +
  base_theme +
  theme(axis.line   = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8))
ggsave("Fig1B.pdf", p1b, width = 5, height = 4)


#### ---- Task 2  Fig1C：大克隆(VAF>=10%) 患病率随年龄段 ----
lc10 <- CHIP_matched %>% filter(VAF >= 0.10) %>%
  distinct(Filename, age_grp, TYPE) %>%
  count(age_grp, TYPE, name = "n10")
fig1c <- denom %>%
  left_join(lc10, by = c("age_grp", "TYPE")) %>%
  mutate(n10  = ifelse(is.na(n10), 0, n10),
         freq = n10 / total * 100) %>%
  filter(!is.na(age_grp))

p1c <- ggplot(fig1c, aes(age_grp, freq, color = TYPE, group = TYPE)) +
  geom_line(linewidth = 1) + geom_point(size = 2.2) +
  scale_color_manual(values = grp_col) +
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, 5)) +
  labs(x = "Age Group", y = "Frequency (%)", title = "Large clone (VAF ≥ 10%)") +
  base_theme +
  theme(axis.line   = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8))
ggsave("Fig1C.pdf", p1c, width = 5, height = 4)


#### ---- Task 3  Fig1D：single- vs multi-mutated 组别构成（堆叠柱）----
indiv <- indiv %>%
  mutate(mut_grp = factor(ifelse(mut_count == 1, "single-mutated", "multi-mutated"),
                          levels = c("single-mutated", "multi-mutated")))
tab_d <- table(indiv$mut_grp, indiv$TYPE)
pval_d <- chisq.test(tab_d)$p.value
d_df <- indiv %>% count(mut_grp, TYPE, name = "n") %>%
  group_by(mut_grp) %>% mutate(pct = n / sum(n) * 100) %>% ungroup()

p1d <- ggplot(d_df, aes(mut_grp, pct, fill = TYPE)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = grp_col) +
  scale_y_continuous(limits = c(0, 107), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  annotate("text", x = 1.5, y = 103, label = paste0("p = ", signif(pval_d, 3)),
           fontface = "italic", size = 4.5) +
  labs(x = NULL, y = "Frequency (%)") + base_theme
ggsave("Fig1D.pdf", p1d, width = 4, height = 4.5)


#### ---- Task 4  Fig1E：max VAF(%) ~ Age（个体级；分面散点 + lm 回归 + 95%CI + Spearman R/p） ----
# 统计单位为个体（每位患者一点：最大克隆 VAF ~ 年龄），避免同一人多突变的伪重复。
# LC 最大克隆随年龄弱-显著增大（Spearman ρ≈0.21, p≈0.015），Non_LC 几乎无趋势
# （ρ≈0.06, p≈0.67）。注：去除最高 2–3 位患者后 LC 显著性减弱，属弱效应，报告时宜注明。
e_df <- indiv %>% mutate(VAFpct = max_vaf * 100)

cc_e <- e_df %>% group_by(TYPE) %>%
  summarise(R = cor(age, VAFpct, method = "spearman"),
            p = suppressWarnings(cor.test(age, VAFpct, method = "spearman")$p.value),
            .groups = "drop") %>%
  mutate(lab = sprintf("Spearman R = %.2f\np = %s", R, signif(p, 2)),
         x   = min(e_df$age),
         y   = max(e_df$VAFpct) * 0.97)

p1e <- ggplot(e_df, aes(age, VAFpct, color = TYPE, fill = TYPE)) +
  geom_jitter(width = 0.6, height = 0, alpha = 0.35, size = 1.3) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1, alpha = 0.18) +
  geom_text(data = cc_e, aes(x = x, y = y, label = lab, color = TYPE),
            hjust = 0, vjust = 1, size = 4, fontface = "italic",
            inherit.aes = FALSE, show.legend = FALSE) +
  facet_wrap(~ TYPE, nrow = 1) +
  scale_color_manual(values = grp_col) +
  scale_fill_manual(values = grp_col) +
  guides(color = "none", fill = "none") +
  labs(x = "Age", y = "Max VAF(%)") + base_theme +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold"))
ggsave("Fig1E.pdf", p1e, width = 7, height = 4)


#### ---- Task 4b  Fig1E_gam：max VAF(%) ~ Age（GAM 样条平滑，两组同面板）----
# 方法可行：用惩罚回归样条（mgcv，REML）替代 loess —— 惩罚项自动控制曲线粗糙度，
# 不会像 loess 那样出现多拐点/多交点。小样本(Non_LC n≈61)+弱信号下，固定低基维 k=4
# 防过拟合；penalty 通常会把 edf 收缩到 ~1（近线性），属诚实结果而非缺陷。
library(mgcv)

eg_df <- indiv %>% mutate(VAFpct = max_vaf * 100)

# 曲线用 GAM 样条展示，统计推断用 Spearman（与 Task 4 一致），标注 Spearman R/p
sp_e <- eg_df %>% group_by(TYPE) %>%
  summarise(R = cor(age, VAFpct, method = "spearman"),
            p = suppressWarnings(cor.test(age, VAFpct, method = "spearman")$p.value),
            .groups = "drop") %>%
  mutate(lab = sprintf("%s: Spearman R = %.2f, p = %s", TYPE, R, signif(p, 2)))

p1e_gam <- ggplot(eg_df, aes(age, VAFpct, color = TYPE, fill = TYPE)) +
  geom_jitter(width = 0.6, height = 0, alpha = 0.3, size = 1.3) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 4),
              method.args = list(method = "REML"),
              se = TRUE, linewidth = 1, alpha = 0.15) +
  scale_color_manual(values = grp_col) +
  scale_fill_manual(values = grp_col) +
  annotate("text", x = min(eg_df$age), y = max(eg_df$VAFpct) * 0.98, hjust = 0, size = 3.8,
           color = unname(grp_col["LC"]),
           label = sp_e$lab[sp_e$TYPE == "LC"]) +
  annotate("text", x = min(eg_df$age), y = max(eg_df$VAFpct) * 0.90, hjust = 0, size = 3.8,
           color = unname(grp_col["Non_LC"]),
           label = sp_e$lab[sp_e$TYPE == "Non_LC"]) +
  labs(x = "Age", y = "Max VAF(%)") + base_theme
ggsave("Fig1E_gam.pdf", p1e_gam, width = 5.5, height = 4.2)


#### ---- Task 4c  Fig1E_q75：max VAF 第75百分位 ~ 年龄段（上尾趋势版配套）----
# 信号在 VAF 上尾，故只追踪各年龄段 maxVAF 的第75百分位；bootstrap 95%CI。
# 仅保留 n>=5 的年龄段（40-49 太小，不可解释）。不替换 GAM 版，仅新增。
set.seed(1)
boot_q <- function(x, prob = 0.75, B = 2000) {
  qs <- replicate(B, quantile(sample(x, replace = TRUE), prob, names = FALSE))
  c(q  = quantile(x, prob, names = FALSE),
    lo = quantile(qs, 0.025, names = FALSE),
    hi = quantile(qs, 0.975, names = FALSE))
}
q75 <- indiv %>% mutate(VAFpct = max_vaf * 100) %>%
  filter(!is.na(age_grp)) %>%
  group_by(TYPE, age_grp) %>%
  group_modify(~{
    b <- boot_q(.x$VAFpct)
    tibble::tibble(n = nrow(.x), q = b["q"], lo = b["lo"], hi = b["hi"])
  }) %>% ungroup() %>% filter(n >= 5)

p1e_q75 <- ggplot(q75, aes(age_grp, q, color = TYPE, group = TYPE)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.7,
                position = position_dodge(0.25)) +
  geom_line(linewidth = 1, position = position_dodge(0.25)) +
  geom_point(size = 2.4, position = position_dodge(0.25)) +
  scale_color_manual(values = grp_col) +
  labs(x = "Age group", y = "Max VAF 75th pct (%)") + base_theme
ggsave("Fig1E_q75.pdf", p1e_q75, width = 5, height = 4)


#### ---- Task 5  Fig1F：基因突变频率（分组柱）+ 每基因 Fisher p ----
top_genes <- head(gene_order, 10)
f_df <- gene_indiv %>% filter(GENE %in% top_genes) %>%
  count(GENE, TYPE, name = "n") %>%
  complete(GENE, TYPE, fill = list(n = 0)) %>%
  mutate(total = ifelse(TYPE == "LC", nLC, nNon),
         freq  = n / total * 100)
fish_p <- sapply(top_genes, function(g) {
  a <- f_df$n[f_df$GENE == g & f_df$TYPE == "LC"]
  b <- f_df$n[f_df$GENE == g & f_df$TYPE == "Non_LC"]
  fisher.test(matrix(c(a, nLC - a, b, nNon - b), nrow = 2))$p.value
})
f_df$GENE <- factor(f_df$GENE, levels = top_genes)
plab_f <- data.frame(
  GENE  = factor(top_genes, levels = top_genes),
  p     = paste0("p=", signif(fish_p, 3)),
  y     = sapply(top_genes, function(g) max(f_df$freq[f_df$GENE == g]))
)

p1f <- ggplot(f_df, aes(GENE, freq, fill = TYPE)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  scale_fill_manual(values = grp_col) +
  geom_text(data = plab_f, aes(GENE, y + 1, label = p), inherit.aes = FALSE,
            angle = 45, hjust = 0, vjust = 0, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(x = NULL, y = "Frequency(%)") + base_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Fig1F.pdf", p1f, width = 7, height = 4.5)


#### ---- Task 6  FigS1A：Mutation counts ~ Age（Poisson GAM 样条平滑）+ Spearman R/p ----
# 图注：GAM curves with spline smoothing；计数响应 → 用 Poisson GAM（family=poisson）更严谨；
# 推断仍用 Spearman（与图注一致）。低基维 k=4 防过拟合。
cc_a <- indiv %>% group_by(TYPE) %>%
  summarise(R = cor(age, mut_count, method = "spearman"),
            p = suppressWarnings(cor.test(age, mut_count, method = "spearman")$p.value),
            .groups = "drop")
ymax_a <- max(indiv$mut_count)

pS1a <- ggplot(indiv, aes(age, mut_count, color = TYPE)) +
  geom_jitter(height = 0.06, width = 0, alpha = 0.4, size = 1.2) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 4),
              method.args = list(family = poisson(), method = "REML"),
              se = FALSE, linewidth = 1) +
  scale_color_manual(values = grp_col) +
  annotate("text", x = min(indiv$age), y = ymax_a * 0.98, hjust = 0, size = 4,
           color = unname(grp_col["LC"]),
           label = sprintf("LC: R = %.2f  p = %s",
                           cc_a$R[cc_a$TYPE == "LC"], signif(cc_a$p[cc_a$TYPE == "LC"], 2))) +
  annotate("text", x = min(indiv$age), y = ymax_a * 0.88, hjust = 0, size = 4,
           color = unname(grp_col["Non_LC"]),
           label = sprintf("Non_LC: R = %.2f  p = %s",
                           cc_a$R[cc_a$TYPE == "Non_LC"], signif(cc_a$p[cc_a$TYPE == "Non_LC"], 2))) +
  labs(x = "Age", y = "Mutation Counts") + base_theme
ggsave("FigS1A.pdf", pS1a, width = 5, height = 4)


#### ---- Task 6b  FigS1A_prop：multi-mutated 比例 ~ 年龄段（比例版配套）----
# 把离散计数(mut_count)改成"携带≥2个突变的人群比例"，趋势更清晰；binom 精确95%CI。
prop_a <- indiv %>% filter(!is.na(age_grp)) %>%
  group_by(TYPE, age_grp) %>%
  summarise(n = n(), k = sum(mut_count >= 2), .groups = "drop") %>%
  mutate(prop = k / n * 100,
         lo = mapply(function(k, n) binom.test(k, n)$conf.int[1], k, n) * 100,
         hi = mapply(function(k, n) binom.test(k, n)$conf.int[2], k, n) * 100)

pS1a_prop <- ggplot(prop_a, aes(age_grp, prop, color = TYPE, group = TYPE)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.7,
                position = position_dodge(0.25)) +
  geom_line(linewidth = 1, position = position_dodge(0.25)) +
  geom_point(size = 2.4, position = position_dodge(0.25)) +
  scale_color_manual(values = grp_col) +
  labs(x = "Age group", y = "Multi-mutated (%)") + base_theme
ggsave("FigS1A_prop.pdf", pS1a_prop, width = 5, height = 4)


#### ---- Task 7  FigS1B：VAF(%) 箱线图 + 抖动点 + Wilcoxon p ----
b_df  <- CHIP_matched %>% mutate(VAFpct = VAF * 100)
pv_b  <- wilcox.test(VAFpct ~ TYPE, data = b_df)$p.value
# 计算两组 VAF 中位数
median_vaf <- b_df %>%
  group_by(TYPE) %>%
  summarise(median_VAFpct = median(VAFpct, na.rm = TRUE), .groups = "drop")

print(median_vaf)
pS1b <- ggplot(b_df, aes(TYPE, VAFpct, color = TYPE)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, fill = NA, linewidth = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  scale_color_manual(values = grp_col) +
  annotate("text", x = 1.5, y = max(b_df$VAFpct) * 0.97,
           label = paste0("p = ", signif(pv_b, 2)), fontface = "italic", size = 4.5) +
  labs(x = NULL, y = "VAF(%)") + base_theme + theme(legend.position = "none")
ggsave("FigS1B.pdf", pS1b, width = 4, height = 4.5)


#### ---- Task 8  FigS1C：VAF 密度图 + 放大插图(0.1-0.2) ----
dmax <- max(sapply(split(CHIP_matched$VAF, CHIP_matched$TYPE),
                   function(x) max(density(x)$y)))
main8 <- ggplot(CHIP_matched, aes(VAF, fill = TYPE, color = TYPE)) +
  geom_density(alpha = 0.4, linewidth = 0.7) +
  scale_fill_manual(values = grp_col) + scale_color_manual(values = grp_col) +
  labs(x = "VAF", y = "Density") + base_theme
inset8 <- ggplot(CHIP_matched, aes(VAF, fill = TYPE, color = TYPE)) +
  geom_density(alpha = 0.4, linewidth = 0.5) +
  scale_fill_manual(values = grp_col) + scale_color_manual(values = grp_col) +
  coord_cartesian(xlim = c(0.1, 0.2), ylim = c(0, 1.5)) +
  labs(x = NULL, y = NULL) + theme_classic() +
  theme(legend.position = "none", axis.text = element_text(size = 7),
        panel.grid = element_blank())
pS1c <- main8 +
  annotation_custom(ggplotGrob(inset8),
                    xmin = 0.15, xmax = 0.42,
                    ymin = dmax * 0.30, ymax = dmax * 0.95)
ggsave("FigS1C.pdf", pS1c, width = 5.5, height = 4)


#### ---- Task 9  FigS1D：Mutation count ~ log10(maxVAF)（Poisson GAM 样条平滑）+ Spearman R/p ----
# 图注：smoothed curves（与 S1A 同属一图 → GAM 样条）；计数响应 → Poisson GAM；推断用 Spearman。
d9 <- indiv %>% mutate(logmax = log10(max_vaf * 100))
cc_d <- d9 %>% group_by(TYPE) %>%
  summarise(R = cor(logmax, mut_count, method = "spearman"),
            p = suppressWarnings(cor.test(logmax, mut_count, method = "spearman")$p.value),
            .groups = "drop")
ymax_d <- max(d9$mut_count)

pS1d <- ggplot(d9, aes(logmax, mut_count, color = TYPE)) +
  geom_jitter(height = 0.08, width = 0, alpha = 0.5, size = 1.2) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 4),
              method.args = list(family = poisson(), method = "REML"),
              se = FALSE, linewidth = 1) +
  scale_color_manual(values = grp_col) +
  annotate("text", x = 0, y = ymax_d * 0.98, hjust = 0, size = 4,
           color = unname(grp_col["LC"]),
           label = sprintf("LC: R = %.2f  p = %s",
                           cc_d$R[cc_d$TYPE == "LC"], signif(cc_d$p[cc_d$TYPE == "LC"], 2))) +
  annotate("text", x = 0, y = ymax_d * 0.88, hjust = 0, size = 4,
           color = unname(grp_col["Non_LC"]),
           label = sprintf("Non_LC: R = %.2f  p = %s",
                           cc_d$R[cc_d$TYPE == "Non_LC"], signif(cc_d$p[cc_d$TYPE == "Non_LC"], 2))) +
  labs(x = "log10(maxVAF)", y = "Mutation count") + base_theme
ggsave("FigS1D.pdf", pS1d, width = 5, height = 4)


#### ---- Task 9b  FigS1D_prop：multi-mutated 比例 ~ 克隆大小分箱（比例版配套）----
# 把"计数 ~ maxVAF"散点改成"按克隆大小分箱的 multi-mutated 比例"，趋势更清晰；binom 95%CI。
vaf_breaks <- c(0, 2, 5, 10, 100)
vaf_labs   <- c("<2", "2-5", "5-10", ">=10")
prop_d <- indiv %>%
  mutate(vaf_bin = cut(max_vaf * 100, breaks = vaf_breaks, labels = vaf_labs,
                       right = TRUE, include.lowest = TRUE)) %>%
  filter(!is.na(vaf_bin)) %>%
  group_by(TYPE, vaf_bin) %>%
  summarise(n = n(), k = sum(mut_count >= 2), .groups = "drop") %>%
  mutate(prop = k / n * 100,
         lo = mapply(function(k, n) binom.test(k, n)$conf.int[1], k, n) * 100,
         hi = mapply(function(k, n) binom.test(k, n)$conf.int[2], k, n) * 100)

pS1d_prop <- ggplot(prop_d, aes(vaf_bin, prop, color = TYPE, group = TYPE)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.7,
                position = position_dodge(0.25)) +
  geom_line(linewidth = 1, position = position_dodge(0.25)) +
  geom_point(size = 2.4, position = position_dodge(0.25)) +
  scale_color_manual(values = grp_col) +
  labs(x = "Max VAF (%) bin", y = "Multi-mutated (%)") + base_theme
ggsave("FigS1D_prop.pdf", pS1d_prop, width = 5, height = 4)


#### ---- Task 10  FigS2A_bar：基因突变计数柱状图（整体，不分组）----
# 整体观察每个基因的出现情况，统一淡绿色填充
gcount <- CHIP_matched %>% count(GENE, name = "n")
gcount$GENE <- factor(gcount$GENE, levels = rev(gene_order))
pS2a_bar <- ggplot(gcount, aes(n, GENE)) +
  geom_bar(stat = "identity", fill = "#a8d5a2", width = 0.75) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Mutation count", y = NULL) + base_theme +
  theme(legend.position = "none")
ggsave("FigS2A_bar.pdf", pS2a_bar, width = 4.5, height = 8)


#### ---- Task 11  FigS2A_waterfall：简化版瀑布图(oncoprint) ----
# onco <- CHIP_matched %>% distinct(Filename, GENE, TYPE) %>%
#   mutate(gi = as.integer(factor(GENE, levels = gene_order)))
# samp_order <- onco %>% group_by(Filename) %>%
#   summarise(top = min(gi), TYPE = first(TYPE), .groups = "drop") %>%
#   arrange(TYPE, top)
# onco$Filename <- factor(onco$Filename, levels = samp_order$Filename)
# onco$GENE     <- factor(onco$GENE, levels = rev(gene_order))
# 
# pS2a_w <- ggplot(onco, aes(Filename, GENE, fill = TYPE)) +
#   geom_tile(width = 0.9, height = 0.8) +
#   scale_fill_manual(values = grp_col) +
#   labs(x = NULL, y = NULL) + base_theme +
#   theme(axis.text.x  = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.line    = element_blank(),
#         panel.background = element_rect(fill = "grey92", color = NA))
# ggsave("FigS2A_waterfall.pdf", pS2a_w, width = 10, height = 7)
library(dplyr)
library(ggplot2)
library(forcats)

gcount <- CHIP_matched %>% 
  count(GENE, name = "n") %>%
  arrange(n)

gcount$GENE <- factor(gcount$GENE, levels = gcount$GENE)

pS2a_bar <- ggplot(gcount, aes(x = n, y = GENE)) +
  geom_bar(
    stat = "identity",
    fill = "#a8d5a2",
    width = 0.75
  ) +
  geom_text(
    aes(label = n),
    hjust = -0.2,
    size = 3.5,
    color = "black"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(
    x = "Mutation count",
    y = NULL
  ) +
  base_theme +
  theme(
    legend.position = "none"
  )

ggsave("FigS2A_bar.pdf", pS2a_bar, width = 4.5, height = 8)

#### ---- Task 12  FigS2B：各基因 VAF(%) 散点（两组着色）----
b12 <- CHIP_matched %>% mutate(VAFpct = VAF * 100)
b12$GENE <- factor(b12$GENE, levels = gene_order)
pS2b <- ggplot(b12, aes(GENE, VAFpct, color = TYPE)) +
  geom_vline(xintercept = seq_along(gene_order), color = "grey85", linewidth = 0.3) +
  geom_jitter(width = 0.22, height = 0, size = 1, alpha = 0.7) +
  scale_color_manual(values = grp_col) +
  labs(x = NULL, y = "VAF(%)") + base_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("FigS2B.pdf", pS2b, width = 11, height = 4.5)


#### ---- Task 13  FigS2C：突变类型频率（分组柱）----
mut_cats <- c("missense_variant", "stopgain", "frameshift_deletion",
              "frameshift_insertion", "splice")
tot13 <- CHIP_matched %>% count(TYPE, name = "tot")
c13 <- CHIP_matched %>% filter(mutTYPE %in% mut_cats) %>%
  count(TYPE, mutTYPE, name = "n") %>%
  complete(TYPE, mutTYPE = mut_cats, fill = list(n = 0)) %>%
  left_join(tot13, by = "TYPE") %>%
  mutate(freq = n / tot * 100,
         mutTYPE = factor(mutTYPE, levels = mut_cats))
pS2c <- ggplot(c13, aes(mutTYPE, freq, fill = TYPE)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  scale_fill_manual(values = grp_col) +
  labs(x = NULL, y = "Frequency(%)") + base_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("FigS2C.pdf", pS2c, width = 6, height = 4.5)


#### ---- Task 14  FigS2D：碱基替换类型频率（分组柱，嘧啶折叠）----
comp <- c(A = "T", T = "A", C = "G", G = "C")
sub_levels <- c("C>T", "T>C", "C>A", "T>A", "C>G", "T>G")
snv <- CHIP_matched %>%
  filter(REF %in% names(comp), ALT %in% names(comp)) %>%
  mutate(r   = ifelse(REF %in% c("A", "G"), comp[REF], REF),
         a   = ifelse(REF %in% c("A", "G"), comp[ALT], ALT),
         sub = paste0(r, ">", a))
tot14 <- snv %>% count(TYPE, name = "tot")
d14 <- snv %>% filter(sub %in% sub_levels) %>%
  count(TYPE, sub, name = "n") %>%
  complete(TYPE, sub = sub_levels, fill = list(n = 0)) %>%
  left_join(tot14, by = "TYPE") %>%
  mutate(freq = n / tot * 100,
         sub  = factor(sub, levels = sub_levels))
pS2d <- ggplot(d14, aes(sub, freq, fill = TYPE)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  scale_fill_manual(values = grp_col) +
  labs(x = NULL, y = "Frequency (%)") + base_theme
ggsave("FigS2D.pdf", pS2d, width = 6, height = 4.5)

message("已生成 18 张图：Fig1B/C/D/E/F, Fig1E_gam, Fig1E_q75, ",
        "FigS1A/B/C/D, FigS1A_prop, FigS1D_prop, FigS2A_bar/A_waterfall/B/C/D")


















