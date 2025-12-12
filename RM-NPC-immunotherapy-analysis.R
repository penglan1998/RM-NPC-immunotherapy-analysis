# This script performs Inverse Probability of Treatment Weighting and survival analyses
# for patients with recurrent/metastatic nasopharyngeal carcinoma (RM-NPC).
# Input data should include treatment group, baseline covariates,
# and survival time and event indicators (PFS/OS).
# The analyses include time-dependent ROC analysis, IPTW,
# and weighted Cox proportional hazards models.

set.seed(12345)

library(timeROC)
library(readxl)
library(survival)
library(dplyr)
library(ggplot2)
library (ggsurvfit)
library(nnet)
library(cobalt)
library(survey)
library(emmeans)

data <- read_excel("xxx.xlsx")

time_roc <- timeROC(
  T = data$PFS,
  delta = data$PFS_event,
  marker = data$NLR,
  cause = 1,
  weighting = "marginal",
  times = c(12),
  ROC = TRUE,
  iid = TRUE
)

# Determine optimal cut-off using Youden index
cutoff_NLR <- time_roc$cut.values[
  which.max(time_roc$TP - time_roc$FP)
]

cutoff_NLR


# Multinomial logistic regression was used to calculate propensity scores

data2 <- read_excel("xxx.xlsx")

fit_multinom <- multinom(
  Group ~ Sex + Age + Body_Mass_Index + Smoking + Drinking + Family_history + EBV_DNA + Prior_anti_PD1 + Prior_antiangiogenic + Stage + Number_of_metastatic_lesions + Prior_lines + Lactate_dehydrogenase + Karnofsky_performance_status_score + Albumin + Neutrophil_to_Lymphocyte_Ratio + Prior_local_regional_radiotherapy,
  data = data2,
  trace = FALSE 
)

summary(fit_multinom)

ps_mat <- fitted(fit_multinom)

dim(ps_mat)
head(ps_mat)
colnames(ps_mat)

p_treat <- prop.table(table(data2$Group))
p_treat

w_iptw <- rep(NA, nrow(data2))

treat_levels <- levels(data2$Group)

for (k in treat_levels) {
  p_k <- p_treat[k]
  idx_k <- which(data2$Group == k)
  ps_k <- ps_mat[idx_k, k]
  w_iptw[idx_k] <- p_k / ps_k
}

data2$w_iptw <- w_iptw

summary(data2$w_iptw)
quantile(data2$w_iptw, probs = c(0.01, 0.1, 0.5, 0.9, 0.99))

hist(data2$w_iptw, breaks = 50, main = "Distribution of IPTW weights")

# weights were truncated at the 1st and 99th percentiles
wt_lo <- quantile(data2$w_iptw, 0.01)
wt_hi <- quantile(data2$w_iptw, 0.99)

data2$w_iptw_trim <- data2$w_iptw
data2$w_iptw_trim[data2$w_iptw < wt_lo] <- wt_lo
data2$w_iptw_trim[data2$w_iptw > wt_hi] <- wt_hi

summary(data2$w_iptw_trim)

bal_unw <- bal.tab(
  Group ~ Sex + Age + Body_Mass_Index + Smoking + Drinking + Family_history + EBV_DNA + Prior_anti_PD1 + Prior_antiangiogenic + Stage + Number_of_metastatic_lesions + Prior_lines + Lactate_dehydrogenase + Karnofsky_performance_status_score + Albumin + Neutrophil_to_Lymphocyte_Ratio + Prior_local_regional_radiotherapy,
  data = data2,
  un = TRUE, 
  weights = NULL,
  method = "weighting"
)
bal_unw

bal_w <- bal.tab(
  Group ~ Sex + Age + Body_Mass_Index + Smoking + Drinking + Family_history + EBV_DNA + Prior_anti_PD1 + Prior_antiangiogenic + Stage + Number_of_metastatic_lesions + Prior_lines + Lactate_dehydrogenase + Karnofsky_performance_status_score + Albumin + Neutrophil_to_Lymphocyte_Ratio + Prior_local_regional_radiotherapy,
  data = data2,
  weights = data2$w_iptw_trim,
  method = "weighting",
  estimand = "ATE",
  un = TRUE 
)
bal_w

var_labels <- c(
  "Sex" = "Sex",
  "Age" = "Age",
  "Body_Mass_Index" = "Body Mass Index",
  "Smoking" = "Smoking",
  "Drinking" = "Drinking",
  "Family_history" = "Family history",
  "EBV_DNA" = "EBV-DNA",
  "Prior_anti_PD1" = "Prior anti-PD-1",
  "Prior_antiangiogenic" = "Prior antiangiogenic",
  "Stage" = "Stage",
  "Number_of_metastatic_lesions" = "Number of metastatic lesions",
  "Prior_lines" = "Prior lines",
  "Lactate_dehydrogenase" = "Lactate dehydrogenase",
  "Karnofsky_performance_status_score" = "Karnofsky performance status score",
  "Albumin" = "Albumin",
  "Neutrophil_to_Lymphocyte_Ratio" = "Neutrophil-to-Lymphocyte Ratio",
  "Prior_local_regional_radiotherapy" = "Prior local-regional radiotherapy"
)

p <- love.plot(
  bal_w,
  stats = "mean.diffs",
  abs = TRUE,
  thresholds = c(m = 0.1),
  var.order = "unadjusted",
  var.names = var_labels 
) +
  theme_minimal(base_family = "Arial") +
  theme(
    text = element_text(family = "Arial", size = 40),
    axis.text.y = element_text(family = "Arial", size = 40),
    axis.text.x = element_text(family = "Arial", size = 40),
    legend.text = element_text(family = "Arial", size = 40),
    legend.title = element_text(family = "Arial", size = 40),
    strip.text = element_text(family = "Arial", size = 40)
  )

p

# Weighted descriptive statistics and weighted Rao-Scott adjusted chi-square tests

design_w <- svydesign(ids = ~1,
                      data = data2,
                      weights = ~w_iptw_trim)

make_weighted_table <- function(var){
  tab <- svytable(as.formula(paste("~ Group +", var)), design = design_w)
  as.data.frame(tab)
}
make_weighted_table("Sex")
make_weighted_table("Age")
make_weighted_table("Body_Mass_Index")
make_weighted_table("Smoking")
make_weighted_table("Drinking")
make_weighted_table("Family_history")
make_weighted_table("EBV_DNA")
make_weighted_table("Prior_anti_PD1")
make_weighted_table("Prior_antiangiogenic")
make_weighted_table("Stage")
make_weighted_table("Number_of_metastatic_lesions")
make_weighted_table("Prior_lines")
make_weighted_table("Lactate_dehydrogenase")
make_weighted_table("Karnofsky_performance_status_score")
make_weighted_table("Albumin")
make_weighted_table("Neutrophil_to_Lymphocyte_Ratio")
make_weighted_table("Prior_local_regional_radiotherapy")

svychisq(~ Sex + Group, design = design_w)
svychisq(~ Age + Group, design = design_w)
svychisq(~ Body_Mass_Index + Group, design = design_w)
svychisq(~ Smoking + Group, design = design_w)
svychisq(~ Drinking + Group, design = design_w)
svychisq(~ Family_history + Group, design = design_w)
svychisq(~ EBV_DNA + Group, design = design_w)
svychisq(~ Prior_anti_PD1 + Group, design = design_w)
svychisq(~ Prior_antiangiogenic + Group, design = design_w)
svychisq(~ Stage + Group, design = design_w)
svychisq(~ Number_of_metastatic_lesions + Group, design = design_w)
svychisq(~ Prior_lines + Group, design = design_w)
svychisq(~ Lactate_dehydrogenase + Group, design = design_w)
svychisq(~ Karnofsky_performance_status_score + Group, design = design_w)
svychisq(~ Albumin + Group, design = design_w)
svychisq(~ Neutrophil_to_Lymphocyte_Ratio + Group, design = design_w)
svychisq(~ Prior_local_regional_radiotherapy + Group, design = design_w)

# Reverse K-M method to calculate the median follow-up time
surv_follow <- Surv(time = data2$`OS`,
                    event = (data2$`OS_event` == 0))

fit_follow_ipw <- survfit(
  surv_follow ~ 1,
  data    = data2,
  weights = data2$w_iptw_trim
)
summary(fit_follow_ipw)$table

# PFS analysis

cox_iptw <- coxph(
  Surv(`PFS`, PFS_event) ~ Group,
  data = data2,
  weights = data2$w_iptw_trim,
  robust = TRUE
)
summary(cox_iptw)

p.global <- summary(cox_iptw)$logtest["pvalue"]
p.global

if (p.global < 0.001) {
  p.text <- "Weighted Cox p < 0.001"
} else {
  p.text <- paste0("Weighted Cox p = ", signif(p.global, 3))
}
p.text

pair_results <- emmeans(cox_iptw, pairwise ~ Group, type = "response")
pair_results

fit.ipw <- survfit(Surv(`PFS`, PFS_event) ~ Group, data = data2,weights = data2$w_iptw_trim)
print(fit.ipw)
summary(fit.ipw)

sf <- fit.ipw
sf$n.risk <- round(sf$n.risk)

ggsurvplot_obj <- sf %>% 
  ggsurvfit(linewidth = 0.8) +
  add_risktable(
    risktable_height = 0.23, 
    risktable_stats = c("{n.risk}"),
    stats_label = list(n.risk = "Number at risk"),
    size = 12,
    theme = list(
      theme_risktable_default(
        axis.text.y.size = 36, 
        plot.title.size  = 36
      ),
      theme(
        plot.title = element_text(face = "bold", family = "Arial"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(family = "Arial"),
        text        = element_text(family = "Arial")
      )
    )
  ) +
  add_risktable_strata_symbol(symbol = "\U25CF", size = 40) +
  labs(
    title = "", 
    x = "Months", 
    y = "Progression-free survival (%)"
  ) +
  scale_x_continuous(
    breaks = seq(0, 40, 4),
    expand = c(0.05, 0),
    limits = c(0, 40)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    labels = seq(0, 100, 25),
    limits = c(0, 1)
  ) + 
  scale_color_manual(
    values = c("#FB6A4A","#2171B5","red","purple","darkgreen"),
    labels = c("CAP","C2P","C1P","AP","P")
  ) +
  scale_fill_manual(
    values = c("#FB6A4A","#2171B5","red","purple","darkgreen"),
    labels = c("CAP","C2P","C1P","AP","P")
  ) +
  theme_classic(base_family = "Arial")
ggsurvplot_obj <- ggsurvplot_obj +
  annotate(
    "text",
    x      = 28,
    y      = 0.98,
    label  = p.text,
    size   = 15,
    family = "Arial"
  ) +
  theme(
    text = element_text(family = "Arial"),
    axis.text = element_text(size = 36, color = "black"),
    axis.title.x = element_text(
      size   = 36,
      color  = "black",
      margin = margin(t = -15)
    ),
    axis.title.y = element_text(
      size   = 36,
      color  = "black",
      margin = margin(r = -200)
    ),
    legend.position    = "none",
    legend.title       = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

print(ggsurvplot_obj)

# OS analysis

cox_iptw <- coxph(
  Surv(`OS`, OS_event) ~ Group,
  data = data2,
  weights = data2$w_iptw_trim,
  robust = TRUE
)
summary(cox_iptw)

p.global <- summary(cox_iptw)$logtest["pvalue"]
p.global

if (p.global < 0.001) {
  p.text <- "Weighted Cox p < 0.001"
} else {
  p.text <- paste0("Weighted Cox p = ", signif(p.global, 3))
}
p.text

pair_results <- emmeans(cox_iptw, pairwise ~ Group, type = "response")
pair_results

fit.ipw <- survfit(Surv(`OS`, OS_event) ~ Group, data = data2,weights = data2$w_iptw_trim)
print(fit.ipw)
summary(fit.ipw)

sf <- fit.ipw
sf$n.risk <- round(sf$n.risk)

ggsurvplot_obj <- sf %>% 
  ggsurvfit(linewidth = 0.8) +
  add_risktable(
    risktable_height = 0.23, 
    risktable_stats = c("{n.risk}"),
    stats_label = list(n.risk = "Number at risk"),
    size = 12,
    theme = list(
      theme_risktable_default(
        axis.text.y.size = 36, 
        plot.title.size  = 36
      ),
      theme(
        plot.title = element_text(face = "bold", family = "Arial"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(family = "Arial"),
        text        = element_text(family = "Arial")
      )
    )
  ) +
  add_risktable_strata_symbol(symbol = "\U25CF", size = 40) +
  labs(
    title = "", 
    x = "Months", 
    y = "Overall survival (%)"
  ) +
  scale_x_continuous(
    breaks = seq(0, 40, 4),
    expand = c(0.05, 0),
    limits = c(0, 40)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    labels = seq(0, 100, 25),
    limits = c(0, 1)
  ) + 
  scale_color_manual(
    values = c("#FB6A4A","#2171B5","red","purple","darkgreen"),
    labels = c("CAP","C2P","C1P","AP","P")
  ) +
  scale_fill_manual(
    values = c("#FB6A4A","#2171B5","red","purple","darkgreen"),
    labels = c("CAP","C2P","C1P","AP","P")
  ) +
  theme_classic(base_family = "Arial")
ggsurvplot_obj <- ggsurvplot_obj +
  annotate(
    "text",
    x      = 28,
    y      = 0.98,
    label  = p.text,
    size   = 15,
    family = "Arial"
  ) +
  theme(
    text = element_text(family = "Arial"),
    axis.text = element_text(size = 36, color = "black"),
    axis.title.x = element_text(
      size   = 36,
      color  = "black",
      margin = margin(t = -15)
    ),
    axis.title.y = element_text(
      size   = 36,
      color  = "black",
      margin = margin(r = -200)
    ),
    legend.position    = "none",
    legend.title       = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

print(ggsurvplot_obj)

sessionInfo()
