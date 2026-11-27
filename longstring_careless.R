# LONGSTRING ANALYSIS WITH THE careless PACKAGE (Yentes & Wilhelm, 2021)
# install.packages("careless")   # run once
library(careless)
library(dplyr)
library(readxl)

data <- read_excel("data_osf.xlsx")
names(data) <- trimws(names(data))

# Memory characteristics items (1-7), in survey order (Month_ago is not a 1-7 rating and is left out)
items <- c("Vividness", "Reliving", "Visual", "Auditory", "Odortaste", "Tactile",
           "Verbal", "Access", "MTT", "Pos_Emo", "Neg_Emo", "Intensity",
           "Importance", "Talk", "Observer")

# Per memory: longest run of identical answers and response variability
data$LS  <- careless::longstring(as.data.frame(data[mc_items]))
data$IRV <- careless::irv(as.data.frame(data[mc_items]))

# Per participant
ls_id <- data %>%
  group_by(ID) %>%
  summarise(long_runs = mean(LS >= 8),   # share of memories with a run of 8+
            LS  = mean(LS),
            IRV = mean(IRV))

mean(ls_id$LS); sd(ls_id$LS)                         # M = 3.12, SD = 1.02

# Flag: Tukey outlier (Q3 + 1.5 * IQR)
cutoff <- quantile(ls_id$LS, .75) + 1.5 * IQR(ls_id$LS)
cutoff                                               # 4.50
ls_id %>% filter(LS > cutoff) %>% arrange(desc(LS))  # 7 participants


# EXCLUDE (chosen after inspecting raw responses of flagged participants)
careless_ids <- c("20-HAN08.docx", "112-nil11.docx", "18-sel04.docx")
data <- data %>% filter(!ID %in% careless_ids)
nrow(data)              # 652
length(unique(data$ID)) # 112
