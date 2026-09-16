---
name: bio-workflows-proteomics-pipeline
description: Orchestrates bottom-up proteomics from a search engine's output (MaxQuant/FragPipe/DIA-NN) to differential protein abundance with limma/DEqMS/MSstats. Use when committing the search database + acquisition mode (DDA vs DIA) up front, re-controlling FDR at PSM AND peptide AND protein-group level (not just PSM), removing contaminant/reverse rows and inspecting RAW distributions before normalizing, bridging cross-plex TMT with an IRS reference channel, modeling MNAR missingness rather than downshift-imputing on/off proteins, batching as a covariate (not pre-subtracted), and testing with treat()/DEqMS. Hands mechanism to the proteomics component skills; not a re-teach of any single step.
tool_type: mixed
primary_tool: limma
workflow: true
depends_on:
  - proteomics/data-import
  - proteomics/proteomics-qc
  - proteomics/quantification
  - proteomics/protein-inference
  - proteomics/differential-abundance
  - proteomics/dia-analysis
---

## Version Compatibility

Reference examples tested with: MSnbase 2.32.0, limma 3.62.2, MSstats 4.14.2, DEqMS 1.24.0, proDA 1.20.0, MSstatsTMT 2.14.2, arrow 23.0.1 (DIA-NN report.parquet), dplyr 1.2.1, tidyr 1.3.2, ggplot2 4.0.3 (checked 2026-09-15 on R 4.4.3)

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Proteomics Pipeline

**"Process my proteomics data from raw MS files to differential abundance"** -> Orchestrate data import (pyopenms/MaxQuant), QC assessment, protein quantification, normalization, differential abundance testing (limma/DEqMS, or MSstats for feature-level designs), and PTM analysis.

This is a workflow skill: it owns the chaining decisions and hand-offs, not the internals of any one step.

Scope: research cohort comparison. A pipeline output is not a diagnostic -- do not use it to triage, escalate or choose treatment for an individual patient; route those questions to validated clinical assays and the treating clinician.

## The governing principle

Bottom-up proteomics never measures proteins; it measures peptides and INFERS proteins, and the trustworthiness decisions are made at seams before the statistics.

1. **The search database + acquisition mode are committed once and inherited by everything.** The FASTA fixes the target-decoy frame (concatenated one-search FDR = #decoy/#target; a separate-search design needs mix-max instead — mixing the two mis-estimates FDR), what counts as a "unique peptide" (relative to the DB: canonical vs +isoforms), and the contaminants (cRAP must be IN the search DB from the start; a contaminant can BE the protein of interest, so never blind-delete `CON__` rows). DDA vs DIA is set at the instrument and dictates which imputation is even legitimate.
2. **FDR is re-controlled at THREE levels, not just PSM — and match-between-runs has its OWN FDR.** 1% PSM-FDR does NOT give 1% protein-FDR — each level (PSM, peptide, protein-group) needs its own target-decoy estimation; PSM-only filtering yields 10-30% real protein-FDR on deep data (one false PSM nucleates a false one-hit-wonder, and false proteins grow with dataset size). Use picked-protein/picked-group FDR. The two-peptide rule INCREASES protein-FDR, it does not reduce it. MBR transfers IDs across runs by RT/m-z and can be wrong for low-abundance precursors — do NOT report MBR-filled counts as directly measured; DIA-NN controls MBR-FDR via `Lib.*` q-values, IonQuant via an explicit MBR-FDR mixture model.
3. **Missingness is MODELED, not filled.** DDA missingness is structured left-censored MNAR; downshift imputation (mean=mu-1.8sigma) on an on/off protein inflates the t-numerator AND deflates the denominator (the volcano "wing" artifact). The honest report for a protein missing in one whole group is "undetected in group B", not a fold change — model the MNAR (proDA/msqrob2/MSstats-AFT).
4. **Normalize AFTER contaminant removal and AFTER inspecting raw distributions; batch is a covariate, not pre-subtracted.** Median-normalizing first mathematically erases a 3x-low load. Cross-plex TMT is invalid without an IRS bridge. `removeBatchEffect` before testing understates residual variance (anticonservative p) — put batch in the same model.

## Made-once commitments

| Commitment | Consequence inherited downstream |
|------------|----------------------------------|
| Search FASTA + target-decoy strategy | The FDR estimator, what a "unique peptide" is, which contaminants exist; a mismatch mis-estimates FDR silently |
| Enzyme + fixed/variable mods (Carbamidomethyl-Cys fixed) | Which peptides exist to quantify; a fixed-mod misconfig loses all Cys peptides |
| DDA vs DIA acquisition mode | Missingness structure (MNAR vs ~MCAR), whether TMT is possible, which imputation is legitimate |
| FDR framing (PSM + peptide + protein-group, 1% each) | Real protein-FDR; PSM-only is 10-30% wrong on deep data |

## Pipeline Overview

```
Raw MS Data (mzML) --> MaxQuant/DIA-NN --> proteinGroups.txt
                                                 |
                                                 v
            +--------------------------------------------+
            |             proteomics-pipeline            |
            +--------------------------------------------+
            |  1. Data Import & Filtering                |
            |  2. Log2 + inspect RAW distributions       |
            |  3. Normalization (after the inspection)   |
            |  4. Per-Group Completeness Filter          |
            |  5. QC: PCA, Correlation                   |
            |  6. Differential Abundance (limma/MSstats) |
            |  7. Visualization & Export                 |
            +--------------------------------------------+
                                                 |
                                                 v
                  Differential Proteins + Volcano Plots
```

## Complete R Workflow

**Goal:** Turn a MaxQuant or DIA-NN protein matrix into a table of differentially abundant proteins with honest missing-value handling.

**Approach:** Strip bookkeeping rows, log2 and inspect the RAW per-sample distributions (dropping failed loads before normalization can hide them), median-center, filter on per-group completeness, then test the OBSERVED values with moderated limma using treat() for a minimum fold change and batch as a covariate -- nothing is imputed. Upgrade to proDA when the dropout itself has to be modeled.

```r
library(limma)
library(ggplot2)
library(pheatmap)

# === 1. DATA IMPORT ===
# quote = '' / comment.char = '': MaxQuant text fields contain quotes and '#'; the defaults
# silently truncate the table with only an 'EOF within quoted string' warning. Check the count.
proteins <- read.delim('proteinGroups.txt', stringsAsFactors = FALSE, quote = '', comment.char = '')
stopifnot(nrow(proteins) == length(readLines('proteinGroups.txt')) - 1)
cat('Loaded', nrow(proteins), 'protein groups\n')

# Filter contaminants, reverse, only-by-site. Use %in%, NOT `!= '+'`: when a run flags none of a
# category the column is entirely empty, read.delim types it as logical NA, `NA != '+'` is NA, and
# subsetting by NA replaces every row with NAs while nrow() still looks right.
proteins <- proteins[!(proteins$Potential.contaminant %in% '+') &
                      !(proteins$Reverse %in% '+') &
                      !(proteins$Only.identified.by.site %in% '+'), ]
stopifnot(!all(is.na(proteins$Majority.protein.IDs)))
cat('After filtering:', nrow(proteins), 'proteins\n')

# Extract LFQ intensities
lfq_cols <- grep('^LFQ\\.intensity\\.', colnames(proteins), value = TRUE)
intensities <- proteins[, lfq_cols]
rownames(intensities) <- proteins$Majority.protein.IDs
colnames(intensities) <- gsub('LFQ\\.intensity\\.', '', colnames(intensities))

# === 2. LOG2 TRANSFORM, THEN INSPECT RAW DISTRIBUTIONS ===
intensities[intensities == 0] <- NA
log2_int <- log2(intensities)

# Inspect BEFORE normalizing (rule 4). Median-centering rescales every sample onto a common median,
# so it mathematically erases the 3x-low load that marks a failed injection -- after this point the
# failure is invisible. Identify and drop failures HERE.
# Inspect the RAW `Intensity.` columns, NOT `LFQ intensity`: MaxLFQ has ALREADY renormalized the LFQ
# columns onto a common scale, so a low load is largely gone from them before you ever look.
raw_cols <- grep('^Intensity\\.', colnames(proteins), value = TRUE)
raw_log2 <- log2(replace(proteins[, raw_cols], proteins[, raw_cols] == 0, NA))
colnames(raw_log2) <- gsub('^Intensity\\.', '', colnames(raw_log2))
raw_log2 <- raw_log2[, colnames(log2_int), drop = FALSE]
raw_median <- apply(raw_log2, 2, median, na.rm = TRUE)
id_counts <- colSums(!is.na(log2_int))
boxplot(raw_log2, las = 2, main = 'RAW log2 Intensity (pre-normalization)', ylab = 'log2 intensity')
print(data.frame(id_count = id_counts, raw_median_log2 = round(raw_median, 2),
                 load_shift_log2 = round(raw_median - median(raw_median), 2)))

# Flag on BOTH axes and as an OUTLIER, not a fixed fraction: a >=2x (1 log2) drop in raw median
# signal, OR an ID count more than 3 MADs below the cohort median. Each axis alone misses cases --
# a 3x-low injection loses its dimmest proteins entirely, so the SURVIVING median understates the
# deficit (a 3.0x low load can read as only -0.8 log2), while a 50%-of-median ID-count rule is so
# loose that a run at 82% of the median IDs, visibly low on both axes, passes it.
mad_ids <- mad(id_counts)
failed <- colnames(raw_log2)[(raw_median - median(raw_median)) <= -1 |
                             (mad_ids > 0 & id_counts < median(id_counts) - 3 * mad_ids)]
if (length(failed) > 0) {
    message('Dropping failed samples: ', paste(failed, collapse = ', '))
    log2_int <- log2_int[, !colnames(log2_int) %in% failed, drop = FALSE]
}

# === 3. NORMALIZE (only after the raw inspection above) ===
sample_medians <- apply(log2_int, 2, median, na.rm = TRUE)
global_median <- median(sample_medians)
normalized <- sweep(log2_int, 2, sample_medians - global_median)

# === 4. FILTER ON PER-GROUP COMPLETENESS (do NOT impute by default) ===
# Filter FIRST on completeness PER GROUP: keep a protein if it is valid in >= ~50-70%
# of replicates in AT LEAST ONE condition. A protein missing in every group fails QC.
sample_info <- read.csv('sample_annotation.csv')
# Re-align the annotation to the samples that SURVIVED the raw-distribution QC above; otherwise the
# column indexing below requests a dropped sample and errors (or silently misaligns the design).
sample_info <- sample_info[sample_info$sample %in% colnames(normalized), ]
sample_info$condition <- droplevels(factor(sample_info$condition))
min_frac <- 0.6   # >= 60% present within at least one group; tune 0.5-0.7 per design
group_complete <- sapply(levels(sample_info$condition), function(g) {
    cols <- sample_info$sample[sample_info$condition == g]
    rowSums(!is.na(normalized[, cols, drop = FALSE])) >= ceiling(length(cols) * min_frac)
})
valid_rows <- rowSums(group_complete) > 0
filtered <- normalized[valid_rows, ]
cat('Proteins after per-group completeness filter:', nrow(filtered), '\n')

# Missingness in label-free DDA is left-censored MNAR (missing BECAUSE low). The modern,
# correct approach is to MODEL the missingness in the likelihood, NOT impute it. See
# proteomics/differential-abundance for the decision (proDA / msqrob2 / MSstats-AFT).
# NOTHING BELOW IMPUTES: limma fits each protein on its OBSERVED values, so the completeness
# filter above is the only missing-value handling in the executed path.

# --- UPGRADE: model the dropout explicitly with proDA, then use `da` in place of `results` ---
# library(proDA)
# fit <- proDA(as.matrix(filtered), design = ~ condition, col_data = sample_info,
#              reference_level = 'Control')            # add batch: design = ~ batch + condition
# da <- test_diff(fit, contrast = 'conditionTreatment')   # columns: diff (log2FC), pval, adj_pval

# --- DO NOT: left-censored downshift imputation, then limma ---
# Downshift MANUFACTURES systematic false positives for on/off proteins near the detection limit
# (the volcano "anchor arms"): it pins missing values ~1.8 SD below the mean with an artificially
# tight 0.3 SD spread, inflating the t-numerator and deflating its denominator. The honest report
# for a protein fully missing in one group is "undetected in group B", NOT a fold change. If a
# reviewer demands a complete matrix anyway, seed it (set.seed) so the call list is reproducible,
# and report the on/off proteins separately -- see proteomics/differential-abundance.

# === 5. QC ===
# PCA needs a complete matrix; use the proteins observed everywhere rather than inventing values.
complete <- filtered[stats::complete.cases(filtered), , drop = FALSE]
cat('Proteins used for PCA (complete cases):', nrow(complete), '\n')
pca <- prcomp(t(complete), scale. = TRUE)
pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Sample = rownames(pca$x))

# === 6. DIFFERENTIAL ANALYSIS (limma on the observed values) ===
# sample_info is already loaded and factored in step 4. Batch is a COVARIATE in the same model
# (rule 4); removeBatchEffect() is visualization-only, never an input to lmFit.
has_batch <- 'batch' %in% colnames(sample_info) && length(unique(sample_info$batch)) > 1
design <- if (has_batch) model.matrix(~ 0 + condition + factor(batch), data = sample_info) else
                         model.matrix(~ 0 + condition, data = sample_info)
colnames(design)[seq_along(levels(sample_info$condition))] <- levels(sample_info$condition)
colnames(design) <- make.names(colnames(design))
cat('design columns:', colnames(design), '| batch in design:', has_batch, '\n')

fit <- lmFit(as.matrix(filtered), design)
contrast <- makeContrasts(Treatment - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast)

# Select on FDR ALONE. A post-hoc fold-change + significance double filter inflates FDR
# (a collider/selection effect; realized FDR can exceed 50%). To require a minimum effect,
# use the moderated minimum-fold-change test treat()/topTreat() instead of filtering after.
fit2_treat <- treat(fit2, lfc = log2(1.5), trend = TRUE, robust = TRUE)   # moderated min-FC test; trend+robust ~mandatory for label-free LFQ
results <- topTreat(fit2_treat, coef = 1, number = Inf)
results$protein <- rownames(results)
# A protein observed in only one group (or one batch level) has no estimable contrast: limma
# returns NA, and `adj.P.Val < 0.05` would make every downstream sum() NA. Report those separately
# as "undetected in group B" -- that is the honest statement, not a fold change.
results$significant <- !is.na(results$adj.P.Val) & results$adj.P.Val < 0.05

# === 7. OUTPUT ===
cat('\nResults:\n')
cat('  Contrast not estimable (report as undetected-in-group):', sum(is.na(results$logFC)), '\n')
cat('  Significant proteins:', sum(results$significant), '\n')
cat('  Up-regulated:', sum(results$significant & results$logFC > 0), '\n')
cat('  Down-regulated:', sum(results$significant & results$logFC < 0), '\n')

write.csv(results, 'proteomics_results.csv', row.names = FALSE)
```

## MSstats Workflow

```r
library(MSstats)

# From MaxQuant. quote = '' and comment.char = '' are REQUIRED: MaxQuant text fields contain
# apostrophes (5'-nucleotidase) and '#', and the read.table defaults silently truncate the table
# with nothing but an 'EOF within quoted string' warning -- here 2954 of 10369 evidence rows and
# 62 of 296 proteins. Check the row count; do not trust the warning to stop you.
evidence <- read.table('evidence.txt', sep = '\t', header = TRUE, quote = '', comment.char = '')
proteinGroups <- read.table('proteinGroups.txt', sep = '\t', header = TRUE, quote = '', comment.char = '')
stopifnot(nrow(evidence) == length(readLines('evidence.txt')) - 1,
          nrow(proteinGroups) == length(readLines('proteinGroups.txt')) - 1)
annotation <- read.csv('annotation.csv')

# Convert to MSstats format
msstats_input <- MaxQtoMSstatsFormat(evidence = evidence,
                                      proteinGroups = proteinGroups,
                                      annotation = annotation)

# Process data
processed <- dataProcess(msstats_input, normalization = 'equalizeMedians',
                         summaryMethod = 'TMP', censoredInt = 'NA')

# Comparison. +1 on the numerator: Treatment=+1, Control=-1 so log2FC = Treatment - Control
# (positive = up in Treatment), matching the label and the limma makeContrasts(Treatment-Control) path.
comparison <- matrix(c(-1, 1), nrow = 1)
rownames(comparison) <- 'Treatment_vs_Control'
colnames(comparison) <- c('Control', 'Treatment')

results <- groupComparison(contrast.matrix = comparison, data = processed)
```

## QC Checkpoints

| Stage | Check | Action if Failed |
|-------|-------|------------------|
| Import | >1000 proteins | Re-run MaxQuant |
| Filter | <30% removed | Check sample prep |
| Missing | <40% per sample | Check MS performance |
| PCA | Replicates cluster | Check for batch effects |
| Stats | FC/FDR pre-specified | Verify thresholds were pre-specified; inspect the volcano for downshift-imputation 'anchor arms' |

## Workflow Variants

### TMT/iTRAQ Isobaric Labeling
Reporter extraction is a spectra-level step, not a text-matrix read. Within a single plex the channels are co-isolated/co-fragmented in the same MS2 event, so relative ratios are stable; but MULTI-batch TMT CANNOT be compared across plexes without an IRS bridge (a pooled reference channel in every plex; Plubell 2017). Route to proteomics/quantification for the mechanics.
```r
library(MSnbase)

# Extract reporter ions from spectra (NOT readMSnSet, which loads an existing text matrix)
raw <- readMSData('tmt.mzML', mode = 'onDisk')
tmt_data <- quantify(raw, reporters = TMT10, method = 'max')
# Correct isobaric impurity cross-talk with the LOT-SPECIFIC matrix from the reagent CoA.
# makeImpuritiesMatrix(x = 10, edit = FALSE) returns a MANUFACTURER TEMPLATE, not an identity
# matrix -- its diagonal runs 0.928-0.965 and 5% of 126 lands in 127C. It is a shape check only;
# the numbers are lot-specific. (edit = TRUE, the default, opens an editor and blocks in scripts.)
# Do NOT use makeImpuritiesMatrix(filename = ...) for TMT10/TMTpro: it reads a CoA laid out by
# Da OFFSET and places each column k POSITIONS away in the reporter list, which is only correct
# for non-interleaved reagents (TMT6, iTRAQ). TMT10/TMTpro interleave N and C, so the +1 Da
# neighbour of 126 is 127C -- TWO positions away -- and the filename route silently writes the
# bleed into 127N instead. Build the matrix by CHANNEL NAME and hand it to purityCorrect:
# tmt10_coa.csv: a square percentage matrix, rows = SOURCE reagent, columns = OBSERVED channel,
# both labelled with the channel names (126, 127N, 127C, ...); diagonal = the lot's purity.
coa <- as.matrix(read.csv('tmt10_coa.csv', row.names = 1, check.names = FALSE))
stopifnot(nrow(coa) == ncol(coa), setequal(rownames(coa), reporterNames(TMT10)))
coa <- coa[reporterNames(TMT10), reporterNames(TMT10)]   # force the quant's channel order
impurities <- coa / 100                                  # CoA percentages -> fractions
tmt_data <- purityCorrect(tmt_data, impurities)
stopifnot(sum(exprs(tmt_data) < 0, na.rm = TRUE) == 0)         # negatives = transposed/mis-ordered matrix

# Multi-batch TMT: do NOT concatenate plexes directly. Use MSstatsTMT, which applies the
# reference-channel (IRS) bridge during summarization:
#   library(MSstatsTMT)
#   summ <- proteinSummarization(msstatstmt_input)   # includes the cross-plex bridge
#   groupComparisonTMT(summ, contrast.matrix = comparison)
```

### SILAC Workflow
Caveat: heavy-Arg -> heavy-Pro metabolic conversion biases ratios for proline-containing peptides (under-counts the heavy channel), and labeling efficiency must be checked (residual light reads as down-regulation). Route to proteomics/quantification for the mechanics.
```r
library(limma)

# SILAC ratios from MaxQuant (quote/comment.char as in the MSstats block above)
silac <- read.delim('proteinGroups.txt', quote = '', comment.char = '')
ratio_cols <- grep('Ratio.H.L.normalized', colnames(silac), value = TRUE)

# Log2 transform ratios. MaxQuant leaves NaN (and 0 for an absent channel) wherever it could not
# form a ratio, so log2 produces NaN/-Inf; coerce every non-finite cell to NA before testing.
silac_log2 <- log2(as.matrix(silac[, ratio_cols]))
silac_log2[!is.finite(silac_log2)] <- NA
rownames(silac_log2) <- silac$Majority.protein.IDs   # without this the result table has no identity

# Keep proteins with >= 2 finite ratios. apply(t.test) over the raw matrix STOPS the whole script
# with "not enough 'x' observations" on the first protein quantified in one replicate only.
keep <- rowSums(!is.na(silac_log2)) >= 2
cat('proteins tested:', sum(keep), 'of', nrow(silac_log2), '\n')

# One-sample moderated test against log2 ratio 0 (no change): an intercept-only limma fit borrows
# variance across proteins, which an unmoderated per-protein t-test at n = 3 cannot. Report the
# BH-adjusted p-value -- a raw p-value per protein controls nothing across thousands of tests.
fit <- eBayes(lmFit(silac_log2[keep, , drop = FALSE]), trend = TRUE, robust = TRUE)
results <- topTable(fit, coef = 1, number = Inf, adjust.method = 'BH')
# Proteins dropped by `keep` are an on/off presence table, not a fold change (quantification).
```

### DIA-NN Workflow
DIA-NN 1.9+ defaults to report.parquet (the only default in 2.0); read it with arrow, not read.delim. Filter on q-values BEFORE pivoting, or low-confidence rows enter the matrix. Route to proteomics/dia-analysis for the mechanics.
```r
library(arrow)
library(dplyr)
library(tidyr)

diann <- read_parquet('report.parquet')

# Filter to 1% FDR at precursor AND protein-group level before pivoting.
# Use the GLOBAL protein-group q-value for the cross-run matrix (per-run min(Q.Value) is anti-conservative).
# When MBR is ON, MBR has its own FDR: add the Lib.* q-values (Lib.Q.Value, Lib.PG.Q.Value <= 0.01).
diann_filt <- diann %>%
    filter(Q.Value <= 0.01 & PG.Q.Value <= 0.01 & Global.PG.Q.Value <= 0.01)

# PG.MaxLFQ is ALREADY cross-run MaxLFQ-normalized at report generation. Re-normalizing it
# double-normalizes -- go straight to log2 + limma with no further normalization. To apply the
# skill's own median-centering instead, pivot raw PG.Quantity here, not PG.MaxLFQ.
protein_matrix <- diann_filt %>%
    select(Protein.Group, Run, PG.MaxLFQ) %>%
    distinct() %>%
    pivot_wider(names_from = Run, values_from = PG.MaxLFQ)

# PG.MaxLFQ path: log2-transform and go straight to limma (no re-normalization). DIA-NN writes 0
# for "not quantified in this run", so 0 -> NA FIRST: log2(0) is -Inf, and -Inf cells propagate
# silently until eBayes stops with "missing value where TRUE/FALSE needed".
m <- as.matrix(protein_matrix[, -1])
rownames(m) <- protein_matrix$Protein.Group
m[m == 0] <- NA
log2_matrix <- log2(m)
stopifnot(!any(is.infinite(log2_matrix)))
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| ~1% PSM-FDR but 10-30% wrong proteins | FDR controlled only at PSM level | Estimate FDR at peptide AND protein-group level (picked-group FDR) |
| Volcano "wings" of huge-FC on/off proteins | Downshift imputation on MNAR (Perseus/MaxQuant) | Model the MNAR (proDA/msqrob2/MSstats-AFT); report "undetected in group B", not a fold change |
| Cross-plex TMT ratios differ 2-5x for no biology | Compared TMT across plexes without IRS | Pooled reference channel in EVERY plex + IRS bridge before comparison |
| A failed-load sample silently carried forward | Normalized before inspecting raw distributions | Filter contaminant/reverse rows -> inspect raw boxplots + ID counts -> remove failures -> THEN normalize |
| Anticonservative p-values | `removeBatchEffect` before testing | Put batch in the model (`~ batch + condition`); removeBatchEffect only for PCA |
| Every ratio subtly wrong | Wrong intensity column (`Intensity` vs `LFQ intensity` vs `iBAQ`) | Pick the right column; convert 0 -> NaN before log2 |
| Spurious DA that flips between conditions | Razor-peptide inference reassigns a shared peptide | Quantify at protein-group level or unique-peptides-only for sensitive comparisons |
| `EOF within quoted string`; far fewer rows than the file has lines | default `read.table` quoting on MaxQuant tables with apostrophes (`5'-nucleotidase`) and `#` | `read.table(..., quote = '', comment.char = '')`, then check `nrow` against `readLines` |
| `eBayes`/`lmFit`: `missing value where TRUE/FALSE needed` on a DIA-NN matrix | DIA-NN writes 0 for "not quantified"; `log2(0)` is `-Inf` | `m[m == 0] <- NA` before `log2`, then `stopifnot(!any(is.infinite(m)))` |
| `t.test`: `not enough 'x' observations` partway through a SILAC run | a protein quantified in one replicate only | filter to >= 2 finite ratios; moderated one-sample limma; report BH-adjusted p, never raw |
| Impurity correction makes adjacent TMT10 channels worse, not better | `makeImpuritiesMatrix(filename=)` places CoA Da-offsets by POSITION, but TMT10/TMTpro interleave N and C | build the matrix by channel name (rows = source reagent) and pass it straight to `purityCorrect` |
| A 3x-low injection sails through the raw-distribution check | the check read `LFQ intensity`, which MaxLFQ already renormalized, and used a 50%-of-median ID rule | inspect the raw `Intensity.` columns; flag a >= 1 log2 load shift or an ID count > 3 MADs below the median |

## References

- Elias JE, Gygi SP (2007) Target-decoy search strategy for increased confidence in large-scale protein identifications by mass spectrometry. *Nature Methods* 4:207-214. DOI 10.1038/nmeth1019.
- Savitski MM, Wilhelm M, Hahne H, Kuster B, Bantscheff M (2015) A scalable approach for protein false discovery rate estimation in large proteomic data sets. *Molecular & Cellular Proteomics* 14:2394-2404. DOI 10.1074/mcp.M114.046995. (picked-protein FDR.)
- Plubell DL, Wilmarth PA, Zhao Y, et al (2017) Extended multiplexing of tandem mass tags (TMT) labeling reveals age and high-fat-diet specific proteome changes in mouse epididymal adipose tissue. *Molecular & Cellular Proteomics* 16:873-890. DOI 10.1074/mcp.M116.065524. (IRS.)
- Ritchie ME, Phipson B, Wu D, et al (2015) limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research* 43:e47. DOI 10.1093/nar/gkv007.
- Zhu Y, Orre LM, Zhou Tran Y, et al (2020) DEqMS: a method for accurate variance estimation in differential protein expression analysis. *Molecular & Cellular Proteomics* 19:1047-1057. DOI 10.1074/mcp.TIR119.001646.

## Related Skills

- proteomics/data-import - Load MS data formats
- proteomics/proteomics-qc - Quality control before analysis
- proteomics/quantification - Normalization, TMT IRS bridge, SILAC mechanics
- proteomics/protein-inference - Razor/shared-peptide assignment to protein groups
- proteomics/differential-abundance - Modeling missingness, moderated testing details
- proteomics/dia-analysis - DIA-NN report parsing and q-value filtering
- proteomics/ptm-analysis - Phosphoproteomics and other PTMs
- data-visualization/volcano-and-ma-plots - Volcano plots with LFC shrinkage
