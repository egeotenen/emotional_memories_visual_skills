###############################################################
#### LONGSTRING (STRAIGHTLINING) ANALYSIS ######################
###############################################################
# Careless-responding check for the memory characteristic ratings
# (memory characteristics questionnaire items only).
#
# Data: data_osf.xlsx, Sheet1 (long format: one row = one memory,
#       up to 6 memories per participant).
#
# Longstring = length of the longest run of identical consecutive
# answers (Johnson, 2005; Meade & Craig, 2012). It is computed
# within each memory, over the 1-7 rating items in the order they
# were shown in the survey (same order as the Qualtrics export):
#
#   Vividness, Reliving, Visual, Auditory, Odortaste, Tactile,
#   Verbal, Access, MTT, Pos_Emo, Neg_Emo, Intensity, Importance,
#   Talk, Observer                                     (15 items)
#
# The manipulation checks (Valence_Check, Arousal_Check) and the
# scene items (S_*) are NOT included.
#
# Month_ago is left out because it is a free number of months,
# not a 1-7 rating. Because it sat between Talk and Observer,
# it is treated as a break in the string (a run cannot continue
# across it).
#
# The OSIQ / VVIQ / MRT / SRIT are only available as totals,
# so longstring cannot be computed for them here.
#
# Output (written to the working directory):
#   longstring_by_memory.csv       one row per memory
#   longstring_by_participant.csv  one row per participant + flags
###############################################################

library(dplyr)
library(tidyr)
library(readxl)


# IMPORT DATA
data <- read_excel("data_osf.xlsx", sheet = "Sheet1")
names(data) <- trimws(names(data))   # some column names have spaces


# ITEMS (in presentation order)
mc_part1   <- c("Vividness", "Reliving", "Visual", "Auditory", "Odortaste",
                "Tactile", "Verbal", "Access", "MTT", "Pos_Emo", "Neg_Emo",
                "Intensity", "Importance", "Talk")
mc_part2   <- c("Observer")            # comes after Month_ago

mc_items <- c(mc_part1, mc_part2)


# HELPER FUNCTIONS
# Longest run of identical consecutive values (NAs break the run)
longstring_max <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(NA_real_)
  r <- rle(x)
  max(r$lengths[!is.na(r$values)])
}

# Average run length (Meade & Craig's "average longstring")
longstring_avg <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(NA_real_)
  r <- rle(x)
  mean(r$lengths[!is.na(r$values)])
}

# Longstring over several segments: runs cannot cross segment borders
longstring_max_segments <- function(row, segments) {
  max(sapply(segments, function(s) longstring_max(unlist(row[s]))))
}
longstring_avg_segments <- function(row, segments) {
  runs <- unlist(lapply(segments, function(s) {
    r <- rle(as.numeric(unlist(row[s])))
    r$lengths[!is.na(r$values)]
  }))
  mean(runs)
}


###############################################################
##### 1. LONGSTRING PER MEMORY ################################
###############################################################

# The sequence is split at Month_ago (see header)
mc_segments <- list(mc_part1, mc_part2)

ls_memory <- data %>%
  dplyr::select(ID, Cue, Emotion_Condition, all_of(mc_items)) %>%
  mutate(row_id = row_number())

ls_memory$LS_max <- apply(ls_memory[mc_items], 1, function(r)
  longstring_max_segments(as.list(r), mc_segments))
ls_memory$LS_avg <- apply(ls_memory[mc_items], 1, function(r)
  longstring_avg_segments(as.list(r), mc_segments))

# Within-memory SD across the 15 items (intra-individual response
# variability, IRV; Dunn et al., 2018). Used as a companion index:
# straightliners have IRV close to 0.
ls_memory$IRV <- apply(ls_memory[mc_items], 1, function(r)
  sd(as.numeric(r), na.rm = TRUE))

# Memories where all 14 items before Month_ago got one single answer
ls_memory <- ls_memory %>%
  mutate(all_same = LS_max >= length(mc_part1))

# Descriptives at memory level
ls_memory %>%
  summarise(
    n_memories    = n(),
    mean_LS   = mean(LS_max),
    sd_LS     = sd(LS_max),
    median_LS = median(LS_max),
    max_LS    = max(LS_max),
    mean_LS_avg = mean(LS_avg),
    n_all_same  = sum(all_same)
  ) %>%
  print(width = Inf)

table(ls_memory$LS_max)


###############################################################
##### 2. LONGSTRING PER PARTICIPANT ###########################
###############################################################

ls_participant <- ls_memory %>%
  group_by(ID) %>%
  summarise(
    n_memories        = n(),
    LS_max_mean = mean(LS_max),   # typical longest run per memory
    LS_max_max  = max(LS_max),    # single worst memory
    LS_avg_mean = mean(LS_avg),   # average run length
    n_all_same  = sum(all_same),
    IRV_mean    = mean(IRV),
    .groups = "drop"
  )

describe_cols <- c("LS_max_mean", "LS_max_max", "LS_avg_mean", "IRV_mean")
ls_participant %>%
  summarise(across(all_of(describe_cols),
                   list(mean = mean, sd = sd, median = median,
                        min = min, max = max))) %>%
  pivot_longer(everything(),
               names_to = c("index", ".value"),
               names_pattern = "(.*)_(mean|sd|median|min|max)$") %>%
  print(width = Inf)


###############################################################
##### 3. FLAGGING LIKELY CARELESS RESPONDERS ##################
###############################################################
# There is no universal longstring cutoff, so two rules are used
# and reported together:
#
# (a) Data-driven (Curran, 2016; Meade & Craig, 2012): a participant
#     is an outlier if their mean per-memory longstring is above
#     Q3 + 1.5*IQR of the sample (Tukey's rule).
#
# (b) Rational: a single answer for at least half of the 15 items
#     in a row (longstring >= 8) is very unlikely for attentive
#     responding, because the items cover different senses and
#     both positive and negative emotion. A participant is flagged
#     if this happens in at least half of their memories.
#
# A participant flagged by either rule is "flagged"; flagged by
# both is a "strong" flag. Flags are for inspection and sensitivity
# analyses, not for automatic exclusion.

half_items <- ceiling(length(mc_items) / 2)   # 8 of 15

q  <- quantile(ls_participant$LS_max_mean, c(.25, .75))
tukey_cutoff <- q[2] + 1.5 * (q[2] - q[1])
cat("\nTukey cutoff for mean per-memory longstring:", round(tukey_cutoff, 2), "\n")
cat("Rational cutoff: longstring >=", half_items, "in at least half of memories\n")

prop_long <- ls_memory %>%
  group_by(ID) %>%
  summarise(prop_long_memories = mean(LS_max >= half_items), .groups = "drop")

ls_participant <- ls_participant %>%
  left_join(prop_long, by = "ID") %>%
  mutate(
    flag_tukey    = LS_max_mean > tukey_cutoff,
    flag_rational = prop_long_memories >= 0.5,
    flag_any      = flag_tukey | flag_rational,
    flag_strong   = flag_tukey & flag_rational
  ) %>%
  arrange(desc(LS_max_mean))

cat("\nParticipants flagged (Tukey rule):  ", sum(ls_participant$flag_tukey), "\n")
cat("Participants flagged (rational rule):", sum(ls_participant$flag_rational), "\n")
cat("Participants flagged (either):       ", sum(ls_participant$flag_any), "\n")
cat("Participants flagged (both):         ", sum(ls_participant$flag_strong), "\n\n")

ls_participant %>%
  filter(flag_any) %>%
  dplyr::select(ID, n_memories, LS_max_mean, LS_max_max,
                prop_long_memories, n_all_same,
                IRV_mean, flag_tukey, flag_rational) %>%
  print(n = Inf, width = Inf)

# How well does longstring agree with low variability?
cor.test(ls_participant$LS_max_mean, ls_participant$IRV_mean,
         method = "spearman", exact = FALSE)


###############################################################
##### 4. SAVE RESULTS #########################################
###############################################################

write.csv(ls_memory, "longstring_by_memory.csv", row.names = FALSE)
write.csv(ls_participant, "longstring_by_participant.csv", row.names = FALSE)

# Vector of flagged IDs for sensitivity analyses in other scripts, e.g.
#   data_clean <- data %>% filter(!ID %in% flagged_ids)
flagged_ids <- ls_participant$ID[ls_participant$flag_any]


###############################################################
##### 5. PLOT ################################################
###############################################################

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  p <- ggplot(ls_participant, aes(x = LS_max_mean, fill = flag_any)) +
    geom_histogram(binwidth = 0.5, color = "white") +
    geom_vline(xintercept = tukey_cutoff, linetype = "dashed") +
    scale_fill_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick"),
                      labels = c("Not flagged", "Flagged")) +
    labs(title = "Mean per-memory longstring by participant",
         x = "Mean longest run of identical answers (15 memory items)",
         y = "Participants", fill = NULL) +
    theme_minimal()
  print(p)
  ggsave("longstring_histogram.png", p, width = 7, height = 4.5, dpi = 150)
}
