###############################################################################
# Visual Imagery Measures Relate Differently to Episodic Detail Production
# and Subjective Vividness: analysis code (OSF)
#
# Sections follow the order of the manuscript and the Supplementary Material:
#   0. Setup
#   1. Data and exclusions (careless responding / longstring)
#   2. Variable preparation
#   3. Preliminary analyses: Table 1, Appendix Table 1, manipulation checks
#   4. Episodic details:  baseline (App. Table 2), Model 1 (App. Table 2a),
#                         Model 2 (Table 2), hurdle model (App. Table 3x), Figure 1
#   5. Subjective vividness: baseline (App. Table 3), Model 3 (Table 3),
#                         Model 4 (App. Table 3a), Model 5 (App. Table 4), Figure 2
#   6. Sensitivity analysis without zero-episodic-detail memories
#                         (App. Tables 1x, 2x, 4x-10x, Figures 1x-2x)
#
# Input:  data_osf.xlsx (Sheet1), one row per memory.
###############################################################################


###############################################################################
# 0. SETUP
###############################################################################

# Run once if needed:
# install.packages(c("readxl", "dplyr", "tidyr", "ggplot2", "psych", "coin",
#                    "glmmTMB", "ordinal", "emmeans", "performance"))

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(psych)        # corr.test
library(coin)         # wilcoxsign_test (z values)
library(glmmTMB)      # negative binomial (GLMM) and hurdle models
library(ordinal)      # clmm for vividness ratings
library(emmeans)      # pairwise comparisons
library(performance)  # compare_performance (AIC weights)

# Helper: as.numeric for columns stored with a decimal comma ("3,25")
to_num <- function(x) as.numeric(gsub(",", ".", as.character(x), fixed = TRUE))

# Helper: M, SD, min, max
describe_var <- function(x) {
  data.frame(M = round(mean(x, na.rm = TRUE), 2), SD = round(sd(x, na.rm = TRUE), 2),
             Min = min(x, na.rm = TRUE), Max = max(x, na.rm = TRUE))
}

# Helper: fixed-effect table (B, SE, z, p) for glmmTMB or clmm models
coef_table <- function(model, component = "cond") {
  s <- summary(model)$coefficients
  if (is.list(s)) s <- s[[component]]                      # glmmTMB
  s <- s[!grepl("\\|", rownames(s)), , drop = FALSE]       # drop clmm thresholds
  out <- data.frame(B = s[, 1], SE = s[, 2], z = s[, 3], p = s[, 4])
  round(out, 3)
}

# Helper: Akaike weights for a set of models
akaike_weights <- function(...) {
  a <- sapply(list(...), AIC)
  w <- exp(-0.5 * (a - min(a)))
  round(data.frame(AIC = a, weight = w / sum(w)), 3)
}


###############################################################################
# 1. DATA AND EXCLUSIONS
###############################################################################

data <- read_excel("data_osf.xlsx", sheet = "Sheet1")
names(data) <- trimws(names(data))     # several column names contain spaces

# Data file already excludes the 4 participants described in the Method
# (outliers on age/imagery scales; memories for 3 or fewer cue words).
nrow(data)              # 670 memories
n_distinct(data$ID)     # 115 participants

# --- 1.1 Careless responding: longstring index (Johnson, 2005) -------------
# Longest run of identical consecutive answers on the 15 memory-characteristics
# items (1-7 scale), per memory, averaged per participant.
# Month_ago (between Talk and Observer) is not a 1-7 rating and breaks the run,
# so the longest possible run is the first 14 items.
mc_items <- c("Vividness", "Reliving", "Visual", "Auditory", "Odortaste", "Tactile",
              "Verbal", "Access", "MTT", "Pos_Emo", "Neg_Emo", "Intensity",
              "Importance", "Talk", "Observer")

longstring <- function(x) max(rle(x[1:14])$lengths)

data$LS  <- apply(data[mc_items], 1, longstring)   # longstring per memory
data$IRV <- apply(data[mc_items], 1, sd)           # response variability per memory

ls_id <- data %>%
  group_by(ID) %>%
  summarise(long_runs = mean(LS >= 8),   # share of memories with a run of 8+
            LS  = mean(LS),
            IRV = mean(IRV))

describe_var(ls_id$LS)                   # M = 3.07, SD = 0.96

# Flag: Tukey outlier (Q3 + 1.5 * IQR) OR a run of 8+ in at least half of memories
ls_cutoff <- quantile(ls_id$LS, .75) + 1.5 * IQR(ls_id$LS)
ls_cutoff                                # 4.67
ls_id %>% filter(LS > ls_cutoff | long_runs >= .5) %>% arrange(desc(LS))

# Excluded after inspecting the raw responses of the 5 flagged participants
careless_ids <- c("20-HAN08.docx", "112-nil11.docx", "18-sel04.docx")
data <- data %>% filter(!ID %in% careless_ids)

nrow(data)              # 652 memories
n_distinct(data$ID)     # 112 participants


###############################################################################
# 2. VARIABLE PREPARATION
###############################################################################

num_cols <- c("Age", "MRT_TOTAL", "VVIQ_TOTAL", "ZMRT_TOTAL", "ZVVIQ_TOTAL",
              "Object_OSIQ_Factor12item", "Spatial_OSIQ_Factor14item",
              "ZObject_OSIQ", "ZSpatial_OSIQ")
data[num_cols] <- lapply(data[num_cols], to_num)

# Valence: dummy coded, neutral memories as the reference condition
data$Emotion_Condition <- relevel(factor(data$Emotion_Condition), ref = "Neutral")
data$ID <- factor(data$ID)

# Standardized OSIQ scores (12-item object and 14-item spatial PCA solutions).
# MRT and VVIQ are used as the pre-standardized columns in the data file.
data$Object_OSIQ_z  <- as.numeric(scale(data$Object_OSIQ_Factor12item))
data$Spatial_OSIQ_z <- as.numeric(scale(data$Spatial_OSIQ_Factor14item))

# One row per participant for participant-level descriptives
participants <- data[!duplicated(data$ID), ]



##### ZERO-EXCLUDED vs. INCLUDED COMPARISON
# Try things with filtered for excluding people who do not have internal details
# This part is reported to reviewers
data_withoutzero <- data[data$Internal_total != 0, ]


# Check vividness of data with zero internal details
data_zero_only <- data[data$Internal_total == 0, ]

df_summary <- data_zero_only %>%
  summarise(
    n = n(),
    mean_ = mean(Vividness, na.rm = TRUE),
    median = median(Vividness, na.rm = TRUE),
    sd_ = sd(Vividness, na.rm = TRUE),
    min = min(Vividness, na.rm = TRUE),
    max = max(Vividness, na.rm = TRUE)
  )

print(df_summary)

df_summary <- data_withoutzero %>%
  summarise(
    n = n(),
    mean_ = mean(Vividness, na.rm = TRUE),
    median = median(Vividness, na.rm = TRUE),
    sd_ = sd(Vividness, na.rm = TRUE),
    min = min(Vividness, na.rm = TRUE),
    max = max(Vividness, na.rm = TRUE)
  )

print(df_summary)

###############################################################################
# 3. PRELIMINARY ANALYSES
###############################################################################

# --- 3.1 Table 1: descriptives ----------------------------------------------
# Memory level (N = 652)
describe_var(data$Internal_total)   # Episodic details
describe_var(data$Vividness)        # Vividness [1-7]

# Participant level (N = 112)
table(participants$Gender)
describe_var(participants$Age)
describe_var(participants$MRT_TOTAL)
describe_var(participants$VVIQ_TOTAL)
describe_var(participants$Object_OSIQ_Factor12item)
describe_var(participants$Spatial_OSIQ_Factor14item)

# --- 3.2 Appendix Table 1: correlations among imagery measures --------------
imagery_vars <- c("ZObject_OSIQ", "ZSpatial_OSIQ", "ZMRT_TOTAL", "ZVVIQ_TOTAL")
cor_res <- corr.test(participants[imagery_vars], method = "pearson", adjust = "none")
round(cor_res$r, 2)
round(cor_res$p, 3)     # unadjusted p values (as reported in the text)

# --- 3.3 Manipulation checks (Wilcoxon signed-rank, Bonferroni) -------------
# Participant means per condition, then paired comparisons.
manipulation_check <- function(variable) {
  wide <- data %>%
    group_by(ID, Emotion_Condition) %>%
    summarise(value = mean(.data[[variable]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Emotion_Condition, values_from = value)

  # Same pair order as in the reported W statistics (W depends on the order)
  pairs <- list(c("Neutral", "Negative"), c("Neutral", "Positive"), c("Negative", "Positive"))
  res <- lapply(pairs, function(pr) {
    w <- wilcox.test(wide[[pr[1]]], wide[[pr[2]]], paired = TRUE, exact = FALSE)
    z <- statistic(wilcoxsign_test(as.formula(paste(pr[1], "~", pr[2])),
                                   data = wide, distribution = "asymptotic"),
                   type = "standardized")
    data.frame(comparison = paste(pr, collapse = " vs "),
               W = unname(w$statistic), z = round(as.numeric(z), 2),
               p = w$p.value,
               n = sum(complete.cases(wide[[pr[1]]], wide[[pr[2]]])))
  })
  res <- bind_rows(res)
  res$p_bonf <- p.adjust(res$p, method = "bonferroni")
  res
}

# Medians per condition
data %>%
  group_by(Emotion_Condition) %>%
  summarise(valence_mdn = median(Valence_Check), arousal_mdn = median(Arousal_Check))

manipulation_check("Valence_Check")
manipulation_check("Arousal_Check")


###############################################################################
# 4. EPISODIC DETAILS (negative binomial GLMMs)
###############################################################################

# Descriptives by condition (Results text)
data %>%
  group_by(Emotion_Condition) %>%
  summarise(M = mean(Internal_total), SD = sd(Internal_total), n = n())

# Figure 1
plot_by_condition <- function(df, outcome, title, ylab) {
  summ <- df %>%
    group_by(Emotion_Condition) %>%
    summarise(m = mean(.data[[outcome]]), se = sd(.data[[outcome]]) / sqrt(n()))
  ggplot(df, aes(x = Emotion_Condition, y = .data[[outcome]], color = Emotion_Condition)) +
    geom_jitter(width = 0.15, alpha = 0.3, size = 1) +
    geom_point(data = summ, aes(x = Emotion_Condition, y = m),
               inherit.aes = FALSE, color = "black", size = 3) +
    geom_errorbar(data = summ, aes(x = Emotion_Condition, ymin = m - 1.96 * se, ymax = m + 1.96 * se),
                  inherit.aes = FALSE, width = 0.1, color = "black") +
    labs(title = title, x = "Emotion Condition", y = ylab) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none")
}

plot_by_condition(data, "Internal_total",
                  "Mean Number of Episodic Details by Emotion Condition",
                  "Episodic Details (Mean ± 95% CI)")

# --- 4.1 Baseline model (Appendix Table 2) ----------------------------------
ep_base <- glmmTMB(Internal_total ~ Emotion_Condition + (1 | ID),
                   family = nbinom2, data = data)
coef_table(ep_base)

# --- 4.2 Model 1: + spatially oriented measures (Appendix Table 2a) ---------
ep_m1 <- glmmTMB(Internal_total ~ Emotion_Condition + ZMRT_TOTAL + Spatial_OSIQ_z + (1 | ID),
                 family = nbinom2, data = data)
coef_table(ep_m1)
anova(ep_base, ep_m1)          # chi2(2) = 9.70, p = .008

# --- 4.3 Model 2: + object/pictorially oriented measures (Table 2) ----------
ep_m2 <- glmmTMB(Internal_total ~ Emotion_Condition + ZMRT_TOTAL + Spatial_OSIQ_z +
                   ZVVIQ_TOTAL + Object_OSIQ_z + (1 | ID),
                 family = nbinom2, data = data)
coef_table(ep_m2)
anova(ep_m1, ep_m2)            # chi2(2) = 1.77, p = .41
akaike_weights(ep_m1, ep_m2)   # weight for Model 1 = .75
check_collinearity(ep_m2)

# --- 4.4 Hurdle negative binomial model (Appendix Table 3x) -----------------
# Zero part: odds of producing NO episodic details (logit).
# Count part: number of details when at least one was produced (truncated NB).
ep_hurdle <- glmmTMB(Internal_total ~ Emotion_Condition + ZMRT_TOTAL + Spatial_OSIQ_z + (1 | ID),
                     ziformula = ~ Emotion_Condition + ZMRT_TOTAL + Spatial_OSIQ_z + (1 | ID),
                     family = truncated_nbinom2, data = data)
coef_table(ep_hurdle, "cond")                 # count component
coef_table(ep_hurdle, "zi")                   # zero component
round(exp(fixef(ep_hurdle)$zi), 2)            # odds ratios (MRT = 0.69, OSIQ-Spatial = 1.47)
c(AIC = AIC(ep_hurdle), BIC = BIC(ep_hurdle), logLik = as.numeric(logLik(ep_hurdle)))
VarCorr(ep_hurdle)                            # random-intercept SDs (square for variances)
sigma(ep_hurdle)                              # dispersion parameter


###############################################################################
# 5. SUBJECTIVE VIVIDNESS (cumulative link mixed models)
###############################################################################

# Descriptives by condition (Results text)
data %>%
  group_by(Emotion_Condition) %>%
  summarise(M = mean(Vividness), SD = sd(Vividness), n = n())

# Figure 2
plot_by_condition(data, "Vividness",
                  "Mean Vividness Ratings by Emotion Condition",
                  "Vividness Ratings (Mean ± 95% CI)")

# --- 5.1 Baseline model (Appendix Table 3) + pairwise comparisons -----------
viv_base <- clmm(factor(Vividness) ~ Emotion_Condition + (1 | ID), data = data)
coef_table(viv_base)
pairs(emmeans(viv_base, ~ Emotion_Condition), adjust = "bonferroni")

# --- 5.2 Model 3: + object/pictorially oriented measures (Table 3) ----------
viv_m3 <- clmm(factor(Vividness) ~ Emotion_Condition + ZVVIQ_TOTAL + Object_OSIQ_z + (1 | ID),
               data = data)
coef_table(viv_m3)
anova(viv_base, viv_m3)        # chi2(2) = 11.67, p = .003

# --- 5.3 Model 4: + spatially oriented measures (Appendix Table 3a) ---------
viv_m4 <- clmm(factor(Vividness) ~ Emotion_Condition + ZVVIQ_TOTAL + Object_OSIQ_z +
                 ZMRT_TOTAL + Spatial_OSIQ_z + (1 | ID),
               data = data)
coef_table(viv_m4)
anova(viv_m3, viv_m4)          # chi2(2) = 3.92, p = .14
akaike_weights(viv_m3, viv_m4) # weight for Model 3 = .51

# --- 5.4 Model 5: valence x object measures (Appendix Table 4, exploratory) -
viv_m5 <- clmm(factor(Vividness) ~ Emotion_Condition * ZVVIQ_TOTAL +
                 Emotion_Condition * Object_OSIQ_z + (1 | ID),
               data = data)
coef_table(viv_m5)
anova(viv_m3, viv_m5)          # chi2(4) = 5.89, p = .21


###############################################################################
# 6. SENSITIVITY ANALYSIS: EXCLUDING MEMORIES WITHOUT EPISODIC DETAILS
###############################################################################

data_nz <- data %>% filter(Internal_total != 0)
nrow(data) - nrow(data_nz)     # 139 memories removed (21.3%)
nrow(data_nz)                  # 513 memories
n_distinct(data_nz$ID)         # 111 participants

# OSIQ scores re-standardized within this subset
data_nz$Object_OSIQ_z  <- as.numeric(scale(data_nz$Object_OSIQ_Factor12item))
data_nz$Spatial_OSIQ_z <- as.numeric(scale(data_nz$Spatial_OSIQ_Factor14item))
participants_nz <- data_nz[!duplicated(data_nz$ID), ]

# --- 6.1 Appendix Table 1x: descriptives ------------------------------------
describe_var(data_nz$Internal_total)
describe_var(data_nz$Vividness)
table(participants_nz$Gender)
describe_var(participants_nz$Age)
describe_var(participants_nz$MRT_TOTAL)
describe_var(participants_nz$VVIQ_TOTAL)
describe_var(participants_nz$Object_OSIQ_Factor12item)
describe_var(participants_nz$Spatial_OSIQ_Factor14item)

# --- 6.2 Appendix Table 2x: correlations ------------------------------------
cor_nz <- corr.test(participants_nz[imagery_vars], method = "pearson", adjust = "none")
round(cor_nz$r, 2)
round(cor_nz$p, 3)

# --- 6.3 Appendix Figures 1x and 2x -----------------------------------------
plot_by_condition(data_nz, "Internal_total",
                  "Mean Number of Episodic Details by Emotion Condition",
                  "Episodic Details (Mean ± 95% CI)")
plot_by_condition(data_nz, "Vividness",
                  "Mean Vividness Ratings by Emotion Condition",
                  "Vividness Ratings (Mean ± 95% CI)")

# --- 6.4 Episodic details (Appendix Tables 4x-6x) ---------------------------
nz_ep_base <- glmmTMB(Internal_total ~ Emotion_Condition + (1 | ID),
                      family = nbinom2, data = data_nz)
nz_ep_m1 <- glmmTMB(Internal_total ~ Emotion_Condition + ZMRT_TOTAL + Spatial_OSIQ_z + (1 | ID),
                    family = nbinom2, data = data_nz)
nz_ep_m2 <- glmmTMB(Internal_total ~ Emotion_Condition + ZMRT_TOTAL + Spatial_OSIQ_z +
                      ZVVIQ_TOTAL + Object_OSIQ_z + (1 | ID),
                    family = nbinom2, data = data_nz)
coef_table(nz_ep_base)         # Table 4x
coef_table(nz_ep_m1)           # Table 5x
coef_table(nz_ep_m2)           # Table 6x
anova(nz_ep_base, nz_ep_m1)
anova(nz_ep_m1, nz_ep_m2)
akaike_weights(nz_ep_m1, nz_ep_m2)   # AIC 3069.10 vs 3071.93

# --- 6.5 Vividness (Appendix Tables 7x-10x) ---------------------------------
nz_viv_base <- clmm(factor(Vividness) ~ Emotion_Condition + (1 | ID), data = data_nz)
nz_viv_m3 <- clmm(factor(Vividness) ~ Emotion_Condition + ZVVIQ_TOTAL + Object_OSIQ_z + (1 | ID),
                  data = data_nz)
nz_viv_m4 <- clmm(factor(Vividness) ~ Emotion_Condition + ZVVIQ_TOTAL + Object_OSIQ_z +
                    ZMRT_TOTAL + Spatial_OSIQ_z + (1 | ID), data = data_nz)
nz_viv_m5 <- clmm(factor(Vividness) ~ Emotion_Condition * ZVVIQ_TOTAL +
                    Emotion_Condition * Object_OSIQ_z + (1 | ID), data = data_nz)
coef_table(nz_viv_base)        # Table 7x
coef_table(nz_viv_m3)          # Table 8x
coef_table(nz_viv_m4)          # Table 9x
coef_table(nz_viv_m5)          # Table 10x
anova(nz_viv_base, nz_viv_m3)
anova(nz_viv_m3, nz_viv_m4)
anova(nz_viv_m3, nz_viv_m5)
