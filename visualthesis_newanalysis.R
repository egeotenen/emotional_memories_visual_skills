# INSTALL PACKAGES
library(readxl)
data <- read_excel("THESISDATA.xlsx")
install.packages("ordinal")
install.packages("emmeans") # Estimated marginal means
install.packages("multcomp") # Multiple comparisons
install.packages("glmmTMB")
library(emmeans)
library(multcomp)
library(lme4)
library(nlme)
library(MASS)
library(survival)
library(fitdistrplus)
library(glmmTMB)
install.packages("nlme")
install.packages("ordinal")
library(ordinal)
library(performance)
# install.packages(c("fitdistrplus","MASS")) # if needed
library(fitdistrplus)
library(MASS)
library(fitdistrplus)
install.packages("psych")
library(psych)

library(ordinal); library(ggplot2); library(dplyr); library(tidyr)

# 1. FIT THE DISTRIBUTIONS TO FIND BEST GLM MODELS 
install.packages(dplyr)
# Dependent Variables
data$Internal_total
descdist(data$Internal_total)



# Remove outliers in internal 
# Compute IQR boundaries
#Q1 <- quantile(data$Internal_total, 0.25, na.rm = TRUE)
#Q3 <- quantile(data$Internal_total, 0.75, na.rm = TRUE)
#IQR_value <- Q3 - Q1

# Define lower and upper bounds
#lower_bound <- Q1 - 1.5 * IQR_value
#upper_bound <- Q3 + 1.5 * IQR_value

# Filter out the outliers
#data <- data[data$Internal_total >= lower_bound & data$Internal_total <= upper_bound, ]
#data$Internal_total



# Gender Check
data %>%
  group_by(ID, Gender) %>%        # group by both ID and Gender
  summarise(n = n(), .groups="drop") %>%  # count occurrences
  group_by(Gender) %>%            # now group just by Gender
  summarise(total = sum(n))


data %>%
  group_by(ID, Emotion_Condition) %>%        # group by both ID and Gender
  summarise(n = n(), .groups="drop") %>%  # count occurrences
  group_by(Emotion_Condition) %>%            # now group just by Gender
  summarise(total = sum(n))


#Manipulation checks
df_summary <- data %>%
  group_by(Emotion_Condition) %>%
  summarise(
    n = n(),
    mean_valence = mean(Valence_Check, na.rm = TRUE),
    median_valence = median(Valence_Check, na.rm = TRUE),
    sd_valence = sd(Valence_Check, na.rm = TRUE),
    min_valence = min(Valence_Check, na.rm = TRUE),
    max_valence = max(Valence_Check, na.rm = TRUE)
  )

print(df_summary)


#Arousal
df_summary <- data %>%
  group_by(Emotion_Condition) %>%
  summarise(
    n = n(),
    mean_valence = mean(Arousal_Check, na.rm = TRUE),
    median_valence = median(Arousal_Check, na.rm = TRUE),
    sd_valence = sd(Arousal_Check, na.rm = TRUE),
    min_valence = min(Arousal_Check, na.rm = TRUE),
    max_valence = max(Arousal_Check, na.rm = TRUE)
  )

print(df_summary)


# Correlation of Metrics

var_list <- c("ZMRT_TOTAL", "ZVVIQ_TOTAL", "ZSpatial_OSIQ", "ZObject_OSIQ")

df_num <- data[, var_list]  # subset first

df_num[] <- lapply(df_num, function(x) {
  if (is.character(x)) {
    # Replace comma with dot, then convert
    as.numeric(gsub(",", ".", x))
  } else {
    x  # keep as is
  }
})


df_num[] <- lapply(df_num, function(x) {
  if (is.factor(x) || is.character(x)) {
    as.numeric(gsub(",", ".", as.character(x)))
  } else {
    x
  }
})

df_num <- data[!duplicated(data$ID), var_list]


# Correlations
corr <- rcorr.test(df_num, method = "pearson")

# --- 4️⃣ Extract and round results ---
r_values <- round(corr$r, 3)      # Kendall tau correlations
p_values   <- round(corr$p, 4)      # p-values
n_values   <- corr$n                # sample sizes for each pair

# --- 5️⃣ Print neatly ---
cat("\n✅ Pearson's r Correlation Matrix:\n")
print(r_values)

cat("\n📊 Corresponding p-values:\n")
print(p_values)

print(n_values)

# Test vividness difference among emotional categories
# --- 1. Descriptive statistics ---
data %>%
  group_by(Emotion_Condition) %>%
  summarise(
    n = n(),
    mean_vividness = mean(Vividness, na.rm = TRUE),
    sd_vividness = sd(Vividness, na.rm = TRUE)
  )

# --- 2. Check ANOVA assumptions ---
# Normality of residuals
anova_model <- aov(Vividness ~ Emotion_Condition, data = data)
shapiro.test(residuals(anova_model))  # p > .05 means normality is okay

# --- 3. Run one-way ANOVA ---
anova_result <- aov(Vividness ~ Emotion_Condition, data = data)
summary(anova_result)
TukeyHSD(anova_result)




best_fit_distribution <- function(x,
                                  bounds = NULL,
                                  candidates = NULL,
                                  plot = TRUE) {
  x <- x[is.finite(x)]  # remove NA/Inf
  n <- length(x)
  if (n < 20) stop("Not enough observations (need >= 20).")
  
  is_integer <- all(abs(x - round(x)) < .Machine$double.eps^0.5)
  unique_vals <- sort(unique(x))
  is_likert <- is_integer && length(unique_vals) <= 7 &&
    all(unique_vals %in% seq(min(unique_vals), max(unique_vals))) &&
    max(unique_vals) <= 7 && min(unique_vals) >= 1
  
  if (is_likert) {
    message("Likert-style ordinal data detected. Use an ordinal model (e.g., ordinal::clmm).")
    return(invisible(NULL))
  }
  
  # Default candidate sets
  if (is.null(candidates)) {
    if (is_integer) {
      candidates <- c("poisson", "nbinom")
    } else {
      candidates <- c("norm", "lnorm", "gamma")
      if (!is.null(bounds) && bounds[1] >= 0 && bounds[2] > bounds[1]) {
        candidates <- c(candidates, "beta_scaled")
      }
    }
  }
  
  # Safe fitting
  safe_fit <- function(expr) {
    tryCatch(suppressWarnings(expr), error = function(e) NULL)
  }
  
  fits <- list()
  aics <- bic <- numeric(0)
  
  for (fam in candidates) {
    fit <- NULL
    if (fam %in% c("norm","lnorm","gamma","poisson","nbinom")) {
      fit <- safe_fit(fitdist(x, distr = fam))
    } else if (fam == "beta_scaled" && !is.null(bounds)) {
      a <- bounds[1]; b <- bounds[2]
      eps <- 1e-6
      x01 <- pmin(pmax((x - a)/(b - a), eps), 1 - eps)
      fit <- safe_fit(fitdist(x01, distr = "beta"))
      if (!is.null(fit)) attr(fit, "beta_scaled_bounds") <- c(a, b)
    }
    
    if (!is.null(fit)) {
      fits[[fam]] <- fit
      k <- length(fit$estimate)
      ll <- fit$loglik
      aics[fam] <- AIC(fit)
      bic[fam] <- -2*ll + k * log(n)
    }
  }
  
  if (length(fits) == 0) stop("No distributions could be fitted.")
  
  tab <- data.frame(
    distribution = names(fits),
    AIC = as.numeric(aics[names(fits)]),
    BIC = as.numeric(bic[names(fits)]),
    row.names = NULL
  )
  tab <- tab[order(tab$AIC), ]
  
  if (plot) {
    op <- par(mfrow = c(2,2))
    on.exit(par(op))
    denscomp(fits[[tab$distribution[1]]], main = paste("Best fit =", tab$distribution[1]))
    cdfcomp(fits[[tab$distribution[1]]])
    qqcomp(fits[[tab$distribution[1]]])
    ppcomp(fits[[tab$distribution[1]]])
  }
  
  return(list(recommended = tab$distribution[1], table = tab, fits = fits))
}

# Test the distribution fits for memory details--> all fit to nbinom
res_cont <- best_fit_distribution(data$Internal_total, bounds = c(0,23))
res_cont$recommended
res_cont$table

#Test the distribution for memory ratings
descdist(data$Internal_total)

###############################################################################

#2. MEMORY DETAIL MODELS --Rely on Nbinom distributions

# 2.1. Prepare the variables:
num_from_comma <- function(x) as.numeric(gsub(",", ".", x, fixed = TRUE))

data$Object_OSIQ_Factor12item  <- num_from_comma(data$Object_OSIQ_Factor12item)
data$Spatial_OSIQ_Factor14item <- num_from_comma(data$Spatial_OSIQ_Factor14item)
data$ZVVIQ_TOTAL  <- num_from_comma(data$ZVVIQ_TOTAL)
data$ZMRT_TOTAL <- num_from_comma(data$ZMRT_TOTAL)

data$Emotion_Condition <- as.factor(data$Emotion_Condition)
# Reference: Neutral condition
data$Emotion_Condition <- relevel(data$Emotion_Condition, ref='Neutral')
data$ID <- as.factor(data$ID)

# Optional but helps convergence/interpretation
data$Object_OSIQ_z  <- scale(data$Object_OSIQ_Factor12item)
data$Spatial_OSIQ_z <- scale(data$Spatial_OSIQ_Factor14item)


# Descriptives
summary(data$Internal_total)
sd(data$Internal_total)

summary(data$Vividness)
sd(data$Vividness)


#Based on Emotion Condition


data %>%
  group_by(Emotion_Condition) %>%
  summarise(
    mean_x = mean(Internal_total, na.rm = TRUE),
    sd_x   = sd(Internal_total, na.rm = TRUE),
    n      = n()
  )


data %>%
  group_by(Emotion_Condition) %>%
  summarise(
    mean_x = mean(Vividness, na.rm = TRUE),
    sd_x   = sd(Vividness, na.rm = TRUE),
    n      = n()
  )


# Models for Memory Details 

# Visualize Internal Details
library(ggplot2)
library(dplyr)

summary_df <- data %>%
  group_by(Emotion_Condition) %>%
  summarise(
    mean_internal = mean(Internal_total, na.rm = TRUE),
    se = sd(Internal_total, na.rm = TRUE) / sqrt(n()),
    n = n()
  )

ggplot(data, aes(x = Emotion_Condition, y = Internal_total, color = Emotion_Condition)) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1) +
  geom_point(data = summary_df,
             aes(x = Emotion_Condition, y = mean_internal),
             inherit.aes = FALSE, color = "black", size = 3) +
  geom_errorbar(data = summary_df,
                aes(x = Emotion_Condition,
                    ymin = mean_internal - 1.96 * se,
                    ymax = mean_internal + 1.96 * se),
                inherit.aes = FALSE, width = 0.1, color = "black") +
  labs(title = "Mean Number of Internal Details by Emotion Condition",
       x = "Emotion Condition", y = "Internal Details (Mean ± 95% CI)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")


# Visualize Vividness 

summary_df <- data %>%
  group_by(Emotion_Condition) %>%
  summarise(
    mean_vividness = mean(Vividness, na.rm = TRUE),
    se = sd(Vividness, na.rm = TRUE) / sqrt(n()),
    n = n()
  )

ggplot(data, aes(x = Emotion_Condition, y = Vividness, color = Emotion_Condition)) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1) +
  geom_point(data = summary_df,
             aes(x = Emotion_Condition, y = mean_vividness),
             inherit.aes = FALSE, color = "black", size = 3) +
  geom_errorbar(data = summary_df,
                aes(x = Emotion_Condition,
                    ymin = mean_vividness - 1.96 * se,
                    ymax = mean_vividness + 1.96 * se),
                inherit.aes = FALSE, width = 0.1, color = "black") +
  labs(title = "Mean Vividness Ratings by Emotion Condition",
       x = "Emotion Condition", y = "Vividness Ratings (Mean ± 95% CI)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")


# 1. Internal Details

model_null <- glmmTMB(
  Internal_total ~ Emotion_Condition  + (1 | ID),
  family = nbinom2,
  data = data,
  na.action = na.exclude
)
summary(model_null)



model1 <- glmmTMB(
  Internal_total ~ Emotion_Condition + ZMRT_TOTAL + Spatial_OSIQ_z + (1 | ID),
  family = nbinom2,
  data = data,
  na.action = na.exclude
)
summary(model1)

anova(model_null,model1)
compare_performance(model_null,model1)


model2 <- glmmTMB(

  Internal_total ~ Emotion_Condition + Object_OSIQ_z + Spatial_OSIQ_z + ZVVIQ_TOTAL + ZMRT_TOTAL+(1 | ID),
  family = nbinom2,
  data = data,
  na.action = na.exclude
)
summary(model2)


# Compare the models:

# Now compare using anova()
anova(model1, model2)
AIC(model1, model2)

compare_performance(model1, model2)



# 2. Event Details 

model1 <- glmmTMB(
  edint ~ Emotion_Condition + ZMRT_TOTAL + Spatial_OSIQ_z + (1 | ID),
  family = nbinom2,
  data = data,
  na.action = na.exclude
)
summary(model1)


model2 <- glmmTMB(
  edint ~ Emotion_Condition + Object_OSIQ_z + Spatial_OSIQ_z + ZVVIQ_TOTAL + ZMRT_TOTAL+(1 | ID),
  family = nbinom2,
  data = data,
  na.action = na.exclude
)
summary(model2)

# 3. Perceptual Details 

model1 <- glmmTMB(
  percint ~ Emotion_Condition + ZMRT_TOTAL + Spatial_OSIQ_z + (1 | ID),
  family = nbinom2,
  data = data,
  na.action = na.exclude
)
summary(model1)


model2 <- glmmTMB(
  percint ~ Emotion_Condition + Object_OSIQ_z + Spatial_OSIQ_z + ZVVIQ_TOTAL + ZMRT_TOTAL+(1 | ID),
  family = nbinom2,
  data = data,
  na.action = na.exclude
)
summary(model2)


# Now compare using anova()
anova(model1, model2)
AIC(model1, model2)


# 4. Time
model1 <- glmmTMB(
  tint ~ Emotion_Condition + ZMRT_TOTAL + Spatial_OSIQ_z + (1 | ID),
  family = nbinom2,
  data = data,
  na.action = na.exclude
)
summary(model1)


model2 <- glmmTMB(
  tint ~ Emotion_Condition + Object_OSIQ_z + Spatial_OSIQ_z + ZVVIQ_TOTAL + ZMRT_TOTAL+(1 | ID),
  family = nbinom2,
  data = data,
  na.action = na.exclude
)
summary(model2)

# 4.1. Time with interactions
model1 <- glmmTMB(
  tint ~ Emotion_Condition + ZMRT_TOTAL + Spatial_OSIQ_z + (1 | ID),
  family = nbinom2,
  data = data,
  na.action = na.exclude
)
summary(model1)

model2 <- glmmTMB(
  tint ~ Emotion_Condition *ZMRT_TOTAL + Emotion_Condition *Spatial_OSIQ_z + (1 | ID),
  family = nbinom2,
  data = data,
  na.action = na.exclude
)
summary(model2)

# 5. Place
model1 <- glmmTMB(
  plint ~ Emotion_Condition + ZMRT_TOTAL + Spatial_OSIQ_z + (1 | ID),
  family = nbinom2,
  data = data,
  na.action = na.exclude
)
summary(model1)


model2 <- glmmTMB(
  plint ~ Emotion_Condition + Object_OSIQ_z + Spatial_OSIQ_z + ZVVIQ_TOTAL + ZMRT_TOTAL+(1 | ID),
  family = nbinom2,
  data = data,
  na.action = na.exclude
)
summary(model2)

# Now compare using anova()
anova(model1, model2)
AIC(model1, model2)

# 5. Place with Interactions
model1 <- glmmTMB(
  plint ~ Emotion_Condition + ZMRT_TOTAL + Spatial_OSIQ_z + (1 | ID),
  family = nbinom2,
  data = data,
  na.action = na.exclude
)
summary(model1)


model2 <- glmmTMB(
  plint ~ Emotion_Condition*Spatial_OSIQ_z + Emotion_Condition*ZMRT_TOTAL (1 | ID),
  family = nbinom2,
  data = data,
  na.action = na.exclude
)
summary(model2)


# Now compare using anova()
anova(model1, model2)
AIC(model1, model2)

# 6. Thought
model1 <- glmmTMB(
  thoint ~ Emotion_Condition + ZMRT_TOTAL + Spatial_OSIQ_z + (1 | ID),
  family = nbinom2,
  data = data,
  na.action = na.exclude
)
summary(model1)


model2 <- glmmTMB(
  thoint ~ Emotion_Condition + Object_OSIQ_z + Spatial_OSIQ_z + ZVVIQ_TOTAL + ZMRT_TOTAL+(1 | ID),
  family = nbinom2,
  data = data,
  na.action = na.exclude
)
summary(model2)
###############################################################################

# PHENOMENOLOGY MODELS with CLMM:


# Random intercept for ID
# 1. Vividness

model_likert_null <- clmm(
  as.factor(Vividness) ~ Emotion_Condition +   (1|ID),
  data = data
)
summary(model_likert_null)

anova(model_likert, model_likert_null)
compare_performance


model_likert <- clmm(
  as.factor(Vividness) ~ Emotion_Condition +  Object_OSIQ_z + ZVVIQ_TOTAL + (1|ID),
  data = data
)
summary(model_likert)


model_likert2 <- clmm(
  as.factor(Vividness) ~ Emotion_Condition +  Object_OSIQ_z + Spatial_OSIQ_z + ZVVIQ_TOTAL + ZMRT_TOTAL+(1|ID),
  data = data
)
summary(model_likert2)

anova(model_likert, model_likert2)
AIC(model_likert, model_likert2)
compare_performance(model_likert,model_likert2)


#1.1. Vividness with interactions
model_likert2 <- clmm(
  as.factor(Vividness) ~  Emotion_Condition*ZVVIQ_TOTAL +Emotion_Condition*Object_OSIQ_z  +(1|ID),
  data = data
)
summary(model_likert2)

anova(model_likert, model_likert2)
AIC(model_likert, model_likert2)


#Visualization of Interaction Effect


# 2. Reliving
model_likert <- clmm(
  as.factor(Reliving) ~ Emotion_Condition +  Object_OSIQ_z + ZVVIQ_TOTAL + (1|ID),
  data = data
)
summary(model_likert)

model_likert2 <- clmm(
  as.factor(Reliving) ~ Emotion_Condition +  Object_OSIQ_z + Spatial_OSIQ_z + ZVVIQ_TOTAL + ZMRT_TOTAL+(1|ID),
  data = data
)
summary(model_likert2)

anova(model_likert, model_likert2)
AIC(model_likert, model_likert2)

# 2.1. Reliving with Interactions
model_likert <- clmm(
  as.factor(Reliving) ~ Emotion_Condition +  Object_OSIQ_z + ZVVIQ_TOTAL + (1|ID),
  data = data
)
summary(model_likert)

model_likert2 <- clmm(
  as.factor(Reliving) ~ Emotion_Condition*Object_OSIQ_z + ZVVIQ_TOTAL*Emotion_Condition +(1|ID),
  data = data
)
summary(model_likert2)

anova(model_likert, model_likert2)
AIC(model_likert, model_likert2)

# 3. Intensity
model_likert <- clmm(
  as.factor(Intensity) ~ Emotion_Condition +  Object_OSIQ_z + ZVVIQ_TOTAL + (1|ID),
  data = data
)
summary(model_likert)

model_likert2 <- clmm(
  as.factor(Intensity) ~ Emotion_Condition +  Object_OSIQ_z + Spatial_OSIQ_z + ZVVIQ_TOTAL + ZMRT_TOTAL+(1|ID),
  data = data
)
summary(model_likert2)

anova(model_likert, model_likert2)
AIC(model_likert, model_likert2)

# 3.1. Intensity with Interactions
model_likert <- clmm(
  as.factor(Intensity) ~ Emotion_Condition +  Object_OSIQ_z + ZVVIQ_TOTAL + (1|ID),
  data = data
)
summary(model_likert)

model_likert2 <- clmm(
  as.factor(Intensity) ~ Emotion_Condition*Object_OSIQ_z + ZVVIQ_TOTAL*Emotion_Condition +(1|ID),
  data = data
)
summary(model_likert2)

anova(model_likert, model_likert2)
AIC(model_likert, model_likert2)

# 4. Importance 
model_likert <- clmm(
  as.factor(Importance) ~ Emotion_Condition +  Object_OSIQ_z + ZVVIQ_TOTAL + (1|ID),
  data = data
)
summary(model_likert)

model_likert2 <- clmm(
  as.factor(Importance) ~ Emotion_Condition +  Object_OSIQ_z + Spatial_OSIQ_z + ZVVIQ_TOTAL + ZMRT_TOTAL+(1|ID),
  data = data
)
summary(model_likert2)

anova(model_likert, model_likert2)
AIC(model_likert, model_likert2)

# 4.1. Importance with Interactions
model_likert <- clmm(
  as.factor(Importance) ~ Emotion_Condition +  Object_OSIQ_z + ZVVIQ_TOTAL + (1|ID),
  data = data
)
summary(model_likert)

model_likert2 <- clmm(
  as.factor(Importance) ~ Emotion_Condition*Object_OSIQ_z + ZVVIQ_TOTAL*Emotion_Condition +(1|ID),
  data = data
)
summary(model_likert2)

anova(model_likert, model_likert2)
AIC(model_likert, model_likert2)
# maybe i should add spatial imagery as well 


# 5. MTT  
model_likert <- clmm(
  as.factor(MTT) ~ Emotion_Condition +  Object_OSIQ_z + ZVVIQ_TOTAL + (1|ID),
  data = data
)
summary(model_likert)

model_likert2 <- clmm(
  as.factor(MTT) ~ Emotion_Condition +  Object_OSIQ_z + Spatial_OSIQ_z + ZVVIQ_TOTAL + ZMRT_TOTAL+(1|ID),
  data = data
)
summary(model_likert2)

anova(model_likert, model_likert2)
AIC(model_likert, model_likert2)

# 6. Visual  
model_likert <- clmm(
  as.factor(Visual) ~ Emotion_Condition +  Object_OSIQ_z + ZVVIQ_TOTAL + (1|ID),
  data = data
)
summary(model_likert)

model_likert2 <- clmm(
  as.factor(Visual) ~ Emotion_Condition +  Object_OSIQ_z + Spatial_OSIQ_z + ZVVIQ_TOTAL + ZMRT_TOTAL+(1|ID),
  data = data
)
summary(model_likert2)

anova(model_likert, model_likert2)
AIC(model_likert, model_likert2)


# 7. Auditory  
model_likert <- clmm(
  as.factor(Auditory) ~ Emotion_Condition +  Object_OSIQ_z + ZVVIQ_TOTAL + (1|ID),
  data = data
)
summary(model_likert)

model_likert2 <- clmm(
  as.factor(Auditory) ~ Emotion_Condition +  Object_OSIQ_z + Spatial_OSIQ_z + ZVVIQ_TOTAL + ZMRT_TOTAL+(1|ID),
  data = data
)
summary(model_likert2)

anova(model_likert, model_likert2)
AIC(model_likert, model_likert2)

# 7.1. Auditory with Interactions
model_likert <- clmm(
  as.factor(Auditory) ~ Emotion_Condition +  Object_OSIQ_z + ZVVIQ_TOTAL + (1|ID),
  data = data
)
summary(model_likert)

model_likert2 <- clmm(
  as.factor(Auditory) ~ Emotion_Condition*Object_OSIQ_z + ZVVIQ_TOTAL*Emotion_Condition +(1|ID),
  data = data
)
summary(model_likert2)

anova(model_likert, model_likert2)
AIC(model_likert, model_likert2)

# 8. Odor-taste  
model_likert <- clmm(
  as.factor(Odortaste) ~ Emotion_Condition +  Object_OSIQ_z + ZVVIQ_TOTAL + (1|ID),
  data = data
)
summary(model_likert)

model_likert2 <- clmm(
  as.factor(Odortaste) ~ Emotion_Condition +  Object_OSIQ_z + Spatial_OSIQ_z + ZVVIQ_TOTAL + ZMRT_TOTAL+(1|ID),
  data = data
)
summary(model_likert2)

anova(model_likert, model_likert2)
AIC(model_likert, model_likert2)

# 8.1. Odortaste with Interactions
model_likert <- clmm(
  as.factor(Odortaste) ~ Emotion_Condition +  Object_OSIQ_z + ZVVIQ_TOTAL + (1|ID),
  data = data
)
summary(model_likert)

model_likert2 <- clmm(
  as.factor(Odortaste) ~  ZVVIQ_TOTAL*Emotion_Condition +(1|ID),
  data = data
)
summary(model_likert2)

anova(model_likert, model_likert2)
AIC(model_likert, model_likert2)

# 9. Tactile
model_likert <- clmm(
  as.factor(Tactile) ~ Emotion_Condition +  Object_OSIQ_z + ZVVIQ_TOTAL + (1|ID),
  data = data
)
summary(model_likert)

model_likert2 <- clmm(
  as.factor(Tactile) ~ Emotion_Condition +  Object_OSIQ_z + Spatial_OSIQ_z + ZVVIQ_TOTAL + ZMRT_TOTAL+(1|ID),
  data = data
)
summary(model_likert2)

anova(model_likert, model_likert2)
AIC(model_likert, model_likert2)

# 9. Verbal
model_likert <- clmm(
  as.factor(Verbal) ~ Emotion_Condition +  Object_OSIQ_z + ZVVIQ_TOTAL + (1|ID),
  data = data
)
summary(model_likert)

model_likert2 <- clmm(
  as.factor(Verbal) ~ Emotion_Condition +  Object_OSIQ_z + Spatial_OSIQ_z + ZVVIQ_TOTAL + ZMRT_TOTAL+(1|ID),
  data = data
)
summary(model_likert2)

anova(model_likert, model_likert2)
AIC(model_likert, model_likert2)
#############################################################################

# POST HOC ANALYSIS
# Compute estimated marginal means
em <- emmeans(model, ~  part)
pairwise_part <- contrast(em, method = "pairwise")
print(pairwise_part)

# Pairwise comparisons for a specific factor within another factor's levels
emmeans(model, pairwise ~ part | Group |round)

# Pairwise comparisons
pairs(em)

# Bonferroni-adjusted comparisons
pairs(em, adjust = "bonferroni")


#PLOT THEM
plot(em, comparisons = TRUE)
interaction.plot(dataforSentimentParts_allrounds$part, dataforSentimentParts_allrounds$Group, dataforSentimentParts_allrounds$value)

