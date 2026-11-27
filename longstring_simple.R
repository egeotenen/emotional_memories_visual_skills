# LONGSTRING ANALYSIS (careless responding)
library(dplyr)
library(readxl)

data <- read_excel("data_osf.xlsx")
names(data) <- trimws(names(data))

# Memory characteristics items (1-7), in survey order
items <- c("Vividness", "Reliving", "Visual", "Auditory", "Odortaste", "Tactile",
           "Verbal", "Access", "MTT", "Pos_Emo", "Neg_Emo", "Intensity",
           "Importance", "Talk", "Observer")

# Longest run of identical answers (Month_ago, between Talk and Observer, breaks the run)
longstring <- function(x) max(rle(x[1:14])$lengths)

data$LS  <- apply(data[items], 1, longstring)   # per memory
data$IRV <- apply(data[items], 1, sd)           # response variability per memory

# Per participant
ls_id <- data %>%
  group_by(ID) %>%
  summarise(long_runs = mean(LS >= 8), LS = mean(LS), IRV = mean(IRV))

summary(ls_id$LS)

# Flag: Tukey outlier OR run of 8+ in at least half of memories
cutoff <- quantile(ls_id$LS, .75) + 1.5 * IQR(ls_id$LS)
ls_id %>% filter(LS > cutoff | long_runs >= .5) %>% arrange(desc(LS))

# EXCLUDE (chosen after inspecting raw responses of flagged participants)
careless_ids <- c("20-HAN08.docx", "112-nil11.docx", "18-sel04.docx")
data <- data %>% filter(!ID %in% careless_ids)
