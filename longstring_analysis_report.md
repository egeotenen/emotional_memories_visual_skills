# Longstring Analysis: Careless Responding Check (Memory Characteristics Items)

**Script:** `longstring_analysis.R`
**Data:** `data_osf.xlsx`, Sheet1 (670 memories from 115 participants)
**Outputs:** `longstring_by_memory.csv`, `longstring_by_participant.csv`, `longstring_histogram.png`

---

## 1. What longstring analysis is

Longstring analysis is a standard way to screen survey data for careless responding (Johnson, 2005; Meade & Craig, 2012; Curran, 2016). It looks for **straightlining**, where a participant gives the same answer to many items in a row without reading them, for example answering "4" to everything.

The **longstring index** is the length of the longest unbroken run of identical answers. If a person answers `7 7 6 7 7 7 7 3`, their longstring is 4.

Attentive participants also give the same answer several times in a row now and then. That is normal, especially when nearby items measure similar things (here, for example, Vividness, Reliving, Visual and Auditory are often all rated 7 for a vivid memory). Longstring therefore isn't a yes/no test. It tells you which participants are unusual compared with the rest of the sample, and those participants are worth a closer look.

---

## 2. How the analysis was done

### 2.1 Unit of analysis

The data are in long format: one row is one recalled memory, and each participant recalled up to 6 memories (101 participants have 6, and 14 have 3 to 5). Each memory was rated on the same items, so longstring is computed **within each memory** and then summarised **per participant**.

Runs across memories are not counted, because the order in which memories were presented is not in the file.

### 2.2 Items included

Only the **memory characteristics questions** are used: 15 items, all on a 1–7 scale, in survey order.

> Vividness, Reliving, Visual, Auditory, Odortaste, Tactile, Verbal, Access, MTT, Pos_Emo, Neg_Emo, Intensity, Importance, Talk, Observer

The item order follows the Qualtrics export columns (see the "Deniz anıları" and "Lunapark" sheets).

**Not included:**
- the manipulation checks (Valence_Check, Arousal_Check)
- the scene items (S_Perspective, S_Content, S_Layout and their repeats)
- the imagery questionnaires (OSIQ, VVIQ, MRT, SRIT), which are only available as totals anyway

**Month_ago.** This item is a number of months, not a 1–7 rating, so it is left out. It appeared between *Talk* and *Observer*, so it is treated as a **break**: a run cannot continue across it. As a result, the longest possible run is 14 (Vividness through Talk).

### 2.3 Indices computed

**Per memory** (`longstring_by_memory.csv`)

| Column | Meaning |
|---|---|
| `LS_max` | Longest run of identical answers across the 15 items |
| `LS_avg` | Average length of all runs (higher = fewer answer changes) |
| `IRV` | Intra-individual response variability: the SD of the 15 answers (Dunn et al., 2018). Close to 0 means nearly the same answer everywhere. |
| `all_same` | TRUE if all 14 items before Month_ago got the same answer |

**Per participant** (`longstring_by_participant.csv`)

| Column | Meaning |
|---|---|
| `LS_max_mean` | Mean of `LS_max` over the participant's memories. **This is the main index.** |
| `LS_max_max` | The participant's single worst memory |
| `LS_avg_mean`, `IRV_mean` | Participant means of the per-memory indices |
| `n_all_same` | Number of memories with all 14 items answered the same |
| `prop_long_memories` | Share of memories with a run of 8 or more |
| `flag_tukey`, `flag_rational`, `flag_any`, `flag_strong` | Flags (see below) |

### 2.4 Flagging rules

The literature has no single agreed cutoff for longstring, so two complementary rules are used.

1. **Outlier rule (data-driven).** A participant is flagged if their `LS_max_mean` is above **Q3 + 1.5 × IQR** (Tukey's rule). In this sample Q1 = 2.50, Q3 = 3.37, so the **cutoff is 4.67**.
2. **Half-the-items rule (rational).** A participant is flagged if, in **at least half of their memories**, they gave the same answer to **8 or more of the 15 items in a row**. The items cover different senses and both positive and negative emotion, so such long runs are very unlikely from attentive responding.

- `flag_any`: flagged by either rule
- `flag_strong`: flagged by both rules

---

## 3. Findings

### 3.1 Sample overview

| | Mean | SD | Median | Min | Max |
|---|---|---|---|---|---|
| Longest run per memory, participant mean (`LS_max_mean`) | 3.07 | 0.96 | 2.83 | 2.00 | 9.83 |
| Worst memory per participant (`LS_max_max`) | 4.58 | 1.84 | 4 | 2 | 14 |
| Average run length (`LS_avg_mean`) | 1.40 | 0.36 | 1.31 | 1.12 | 4.36 |
| Response SD (`IRV_mean`) | 2.04 | 0.36 | 2.10 | 0.60 | 2.84 |

- A typical participant's longest run is about **3 identical answers** per memory, which is normal.
- Of 670 memories, only **9** have a run of 8 or more, and only **3** have the same answer for all 14 items before Month_ago.
- The distribution is narrow (most participants are between 2.5 and 3.4), so the Tukey cutoff (4.67) is fairly strict and also catches some participants whose runs have an innocent explanation (see 3.2).

### 3.2 Flagged participants

| ID | Memories | Mean longstring | Worst memory | Memories with run ≥ 8 | Response SD | Rules met | Assessment |
|---|---|---|---|---|---|---|---|
| **20-HAN08** | 6 | 9.83 | 14 | 4 of 6 | 0.72 | Both | **Very likely careless** |
| **112-nil11** | 6 | 5.50 | 14 | 1 of 6 | 0.60 | Outlier | **Suspicious, inspect** |
| **18-sel04** | 6 | 5.33 | 9 | 2 of 6 | 2.19 | Outlier | **Partly suspicious, inspect** |
| **51-GAM09** | 6 | 5.17 | 6 | 0 of 6 | 2.29 | Outlier | **Probably a false alarm** |
| **90-ZEY11** | 6 | 5.17 | 7 | 0 of 6 | 2.21 | Outlier | **Probably a false alarm** |

**20-HAN08.** For two memories (Mezarlik, Hastane) the participant answered **4 to every item**, and Lunapark is 4 throughout apart from Vividness. Deniz and Koridor also have long runs of 4s. Only Dag shows clearly varied answers. Answering 4 (the midpoint) for both positive and negative emotion, and for every sense, is the classic careless pattern. This is the clearest case in the data.

**112-nil11.** Almost every answer is 6 or 7, and this participant has the lowest response SD in the sample (0.60). The cemetery memory (Mezarlik) got **7 on every item, including both Pos_Emo and Neg_Emo**, which is contradictory. This could be careless responding or an extreme "agree" response style. Reading the memory narratives would help decide.

**18-sel04.** Most memories look attentive, with varied answers and sensible emotion ratings. However, two memories look careless: Lunapark is almost all 4s (a run of 8), and Hastane is almost all 7s, including **7 for both Pos_Emo and Neg_Emo** (a run of 9). Consider dropping just those two memories.

**51-GAM09 and 90-ZEY11.** Their runs come mainly from rating the first items (Vividness, Reliving, Visual, Auditory, and sometimes Odortaste and Tactile) all 7, which is plausible for vivid memories. Their other answers vary, and the emotion items match the memory type (e.g. low Neg_Emo for positive memories). Their response SD is above the sample average. **These do not look like careless responding.**

### 3.3 Longstring and response variability

Across participants, mean longstring and response SD were weakly related (Spearman ρ = 0.19, p = .04). The two indices pick up different things: someone can use the whole 1–7 range (high SD) and still have long runs of 7s on related items. The clearest careless cases are those that are **extreme on both** (long runs **and** low SD), which describes 20-HAN08 and 112-nil11.

---

## 4. How to interpret and use the results

1. **Flags call for a closer look, not automatic exclusion.** For each flagged participant, check the raw answers (in `longstring_by_memory.csv`) and, ideally, the memory narratives. Do they describe a real memory with detail? Do the ratings contradict the content?
2. **Runs on related items are expected.** High ratings on Vividness, Reliving, Visual and Auditory often go together. Runs that span unrelated items (e.g. Odortaste and Tactile at the same level as Verbal, Access and Importance) are more telling.
3. **Contradictions are strong evidence.** The same high or middle score on both Pos_Emo and Neg_Emo, or on every sense, is hard to explain by real memory characteristics.
4. **Report it as a sensitivity analysis.** A common approach is to run the main models (a) on the full sample and (b) without flagged participants, and report whether conclusions change. The script creates `flagged_ids` (all 5 flags) for this:

   ```r
   source("longstring_analysis.R")
   data_clean <- data %>% filter(!ID %in% flagged_ids)
   # only the strongest case (20-HAN08):
   strong_ids <- ls_participant$ID[ls_participant$flag_strong]
   ```

   Based on the inspection above, a reasonable middle option is to exclude **20-HAN08** and **112-nil11**, and the two suspicious memories of **18-sel04** (Lunapark, Hastane). Dropping single memories works because the models already use memories as rows with `(1|ID)`.
5. **Suggested write-up** (adapt as needed):

   > *Careless responding was assessed with the longstring index (Johnson, 2005), computed within each memory as the longest run of identical consecutive responses across the 15 memory characteristics items (1–7 scale). Participants were flagged if their mean longstring exceeded Q3 + 1.5 × IQR of the sample (4.67) or if at least half of their memories contained a run of 8 or more identical responses. Five of 115 participants were flagged. Inspection of the raw responses showed one clear case of straightlining (identical midpoint responses on all items for several memories), one participant with uniformly high ratings including contradictory emotion ratings, one participant with two straightlined memories, and two flags driven by consistently high ratings on the related vividness and sensory items. Main analyses were repeated excluding [participants/memories]; [results were unchanged / describe changes].*

---

## 5. Limitations

- Longstring only detects **straightlining**. It does not catch random or patterned careless answers (e.g. 1-7-1-7), which would need other indices such as Mahalanobis distance or psychometric synonyms and antonyms.
- Runs are counted within each memory using the export column order, which is assumed to match the on-screen order.
- With 15 items per memory, runs are short, so the Tukey cutoff is sensitive to small differences. It should be recomputed if the sample changes (the script does this automatically).

---

## References

- Curran, P. G. (2016). Methods for the detection of carelessly invalid responses in survey data. *Journal of Experimental Social Psychology, 66*, 4–19.
- Dunn, A. M., Heggestad, E. D., Shanock, L. R., & Theilgard, N. (2018). Intra-individual response variability as an indicator of insufficient effort responding. *Journal of Business and Psychology, 33*, 105–121.
- Johnson, J. A. (2005). Ascertaining the validity of individual protocols from Web-based personality inventories. *Journal of Research in Personality, 39*, 103–129.
- Meade, A. W., & Craig, S. B. (2012). Identifying careless responses in survey data. *Psychological Methods, 17*, 437–455.
