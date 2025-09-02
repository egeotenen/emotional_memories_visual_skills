# INSTALL PACKAGES
library(readxl)
data <- read_excel("THESISDATA.xlsx")

install.packages("emmeans") # Estimated marginal means
install.packages("multcomp") # Multiple comparisons
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

# 1. FIT THE DISTRIBUTIONS TO FIND BEST GLM MODELS 

# Dependent Variables
data$Internal_total
descdist(data$Internal_total)

# install.packages(c("fitdistrplus","MASS")) # if needed
library(fitdistrplus)
library(MASS)
library(fitdistrplus)

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
res_cont <- best_fit_distribution(data$External, bounds = c(0,40))
res_cont$recommended
res_cont$table

#Test the distribution for memory ratings
descdist(data$edint)

###############################################################################

#2. MEMORY DETAIL MODELS --Rely on Nbinom distributions

# 2.1. Prepare the variables:
num_from_comma <- function(x) as.numeric(gsub(",", ".", x, fixed = TRUE))

data$Object_OSIQ_Factor12item  <- num_from_comma(data$Object_OSIQ_Factor12item)
data$Spatial_OSIQ_Factor14item <- num_from_comma(data$Spatial_OSIQ_Factor14item)
data$ZVVIQ_TOTAL  <- num_from_comma(data$ZVVIQ_TOTAL)
data$ZMRT_TOTAL <- num_from_comma(data$ZMRT_TOTAL)

data$Emotion_Condition <- as.factor(data$Emotion_Condition)
data$ID <- as.factor(data$ID)

# Optional but helps convergence/interpretation
data$Object_OSIQ_z  <- scale(data$Object_OSIQ_Factor12item)
data$Spatial_OSIQ_z <- scale(data$Spatial_OSIQ_Factor14item)





# Models for Memory Details 

# 1. Internal Details

model1 <- glmmTMB(
  Internal_total ~ Emotion_Condition + ZMRT_TOTAL + Spatial_OSIQ_z + (1 | ID),
  family = nbinom2,
  data = data,
  na.action = na.exclude
)
summary(model1)


model2 <- glmmTMB(
  Internal_total ~ Emotion_Condition + Object_OSIQ_z + Spatial_OSIQ_z + ZVVIQ_TOTAL + ZMRT_TOTAL+(1 | ID),
  family = nbinom2,
  data = data,
  na.action = na.exclude
)
summary(model2)

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

#1.1. Vividness with interactions

model_likert <- clmm(
  as.factor(Vividness) ~ Emotion_Condition +  Object_OSIQ_z + ZVVIQ_TOTAL + (1|ID),
  data = data
)
summary(model_likert)

model_likert2 <- clmm(
  as.factor(Vividness) ~  Emotion_Condition*ZVVIQ_TOTAL +(1|ID),
  data = data
)
summary(model_likert2)

anova(model_likert, model_likert2)
AIC(model_likert, model_likert2)


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

