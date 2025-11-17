# CCGA Cancer Detection Study

## Overview

This repository contains analysis and replication materials for the first CCGA (Circulating Cell-free Genome Atlas) substudy, which compares different genomic approaches for detecting cancer from blood samples.

**Paper:** [Jamshidi et al., Cancer Cell (2022)](https://www.cell.com/cancer-cell/fulltext/S1535-6108(22)00513-X)

**Study Design:** Multi-assay comparison evaluating which cell-free DNA (cfDNA) sequencing methods are most effective for early cancer detection across diverse cancer types and stages.

**Cohort Size:** 2,261 participants

- Training set: 1,414 participants
- Validation set: 847 participants

---

## Table of Contents

1. [Data Overview](#data-overview)
2. [Main Analyses](#main-analyses)
3. [Statistical Methods](#statistical-methods)
4. [Replication Guide](#replication-guide)
5. [Key Findings](#key-findings)

---

## Data Overview

### 1. Clinical Data

**Description:** Participant demographics and cancer characteristics

**Contents:**

- Age, sex, smoking status
- Cancer type and clinical stage (I-IV)
- Training/validation assignment

**Usage:** Stratification, classifier training, outcome analysis

### 2. Cell-Free DNA (cfDNA) Sequencing Data

All three assays below were performed on **the same blood samples** from each participant:

#### a) Whole-Genome Bisulfite Sequencing (WGBS)

**Measurement:** DNA methylation patterns across ~30 million CpG sites
**Coverage:** 30X
**File Format:** Methylation beta values or fragment-level methylation data
**Classifiers Used:** WG methylation (best-performing single assay!)
**Key Feature:** Methylation patterns are tissue-specific, enabling cancer type prediction

#### b) Targeted Sequencing (TS)

**Measurement:** Mutations in 507 cancer-related genes
**Coverage:** 60,000X raw depth (3,000X unique molecules with UMI correction)
**File Format:** Variant calls (SNVs, indels), gene-level features
**Classifiers Used:** SNV, SNV-WBC
**Key Feature:** Detects rare cancer-specific somatic mutations

#### c) Whole-Genome Sequencing (WGS)

**Measurement:** Copy number alterations and fragment characteristics
**Coverage:** 30X
**File Format:** Read counts in genomic bins, fragment length distributions
**Classifiers Used:**

- SCNA (somatic copy number alterations)
- SCNA-WBC
- Fragment endpoints
- Fragment lengths
- Allelic imbalance

### 3. Matched White Blood Cell (WBC) Data

**Source:** Buffy coat fraction from the same blood draw
**Sequencing:** TS and WGS protocols (matching cfDNA assays)
**Purpose:** Filter out clonal hematopoiesis (CH) variants
**Why Critical:** Approximately 70% of variants detected in cancer patient cfDNA originate from blood cells (CH), not tumors

### 4. Tumor Tissue Data

**Source:** FFPE tumor biopsies
**Sequencing:** WGS at 60X coverage
**Purpose:** Calculate circulating tumor allele fraction (cTAF)
**Availability:** 409/1,297 cancer patients (31.5%)
**Used For:** Clinical limit of detection (LOD) analyses (Figures 3, 5, S2, S6, S7)

### Data Processing Flow

```
Blood Sample
    ├─> Plasma → cfDNA extraction → WGBS + TS + WGS (3 assays)
    └─> Buffy coat → WBC DNA → TS + WGS (noise filtering)

Tumor Tissue (subset only)
    └─> FFPE biopsy → WGS → Somatic variant calling → cTAF estimation
```

---

## Main Analyses

### Analysis 1: Cancer Signal Detection (Table 3, Figure 2)

**Research Question:** Can we distinguish cancer from non-cancer samples using cfDNA?

**Data Required:**

- cfDNA sequencing from all three assays (WGBS, TS, WGS)
- Matched WBC sequencing (for WBC-filtered classifiers)
- Clinical labels (cancer vs. non-cancer)

**The 10 Classifiers:**

#### 1. Allelic Imbalance

- **Data Source:** WGS
- **What it detects:** Imbalances in allele frequencies at heterozygous SNP sites
- **How it works:** In normal DNA, heterozygous SNPs show ~50/50 ratio between alleles. Cancer DNA can show skewed ratios due to copy number changes or loss of heterozygosity (LOH)
- **Features:** Deviation from expected 0.5 allele frequency across genomic regions
- **Key insight:** Detects chromosomal-level abnormalities without needing matched normal tissue

#### 2. SNV (Single Nucleotide Variants)

- **Data Source:** Targeted sequencing (507 cancer genes)
- **What it detects:** Somatic mutations in cancer-related genes
- **How it works:** Identifies low-frequency variants in cfDNA that differ from germline
- **Features:** Gene-level maximum variant allele frequency (VAF)
- **Limitation:** Cannot distinguish tumor mutations from clonal hematopoiesis (CH)

#### 3. SNV-WBC (SNV with White Blood Cell filtering)

- **Data Source:** Targeted sequencing + matched WBC sequencing
- **What it detects:** Tumor-specific mutations after removing CH variants
- **How it works:** Same as SNV, but subtracts variants also found in WBC DNA
- **Features:** Gene-level maximum VAF (CH-filtered)
- **Key advantage:** Removes the ~70% of variants that come from blood cells, not tumors

#### 4. SCNA (Somatic Copy Number Alterations)

- **Data Source:** WGS
- **What it detects:** Genome-wide copy number changes (gains/losses of chromosomal regions)
- **How it works:** Counts reads in 100kb genomic bins and detects deviations from diploid coverage
- **Features:** Read depth across genomic bins (GC-corrected, PCA-normalized)
- **Model:** Convolutional neural network to learn spatial patterns
- **Limitation:** CH can also cause copy number changes

#### 5. SCNA-WBC (SCNA with WBC filtering)

- **Data Source:** WGS + matched WBC WGS
- **What it detects:** Tumor-specific copy number alterations
- **How it works:** Uses WBC data to mask clonal hematopoiesis-derived CNAs
- **Features:** Same as SCNA but with WBC-masking applied
- **Key advantage:** Distinguishes tumor CNAs from blood cell CNAs

#### 6. Fragment Endpoints

- **Data Source:** WGS
- **What it detects:** Preferential cfDNA fragmentation patterns at specific genomic locations
- **How it works:** Cancer and normal cells fragment DNA differently due to nucleosome positioning and chromatin accessibility
- **Features:** Frequency of fragment start/end positions across the genome
- **Biological basis:** Tumor-derived cfDNA has different fragmentation profiles than normal cfDNA

#### 7. Fragment Lengths

- **Data Source:** WGS
- **What it detects:** Shifts in the size distribution of cfDNA fragments
- **How it works:** Measures distribution of DNA fragment sizes (typically 150-200bp for normal cfDNA)
- **Features:** Fragment length distributions and summary statistics
- **Biological basis:** Tumor cfDNA often shows shorter fragments and different periodicity patterns

#### 8. Clinical Data

- **Data Source:** Patient demographics and clinical history
- **What it detects:** Cancer risk based on non-genomic factors
- **Features:** Age, sex, smoking status, family history
- **Purpose:** Baseline comparison to show added value of genomic classifiers
- **Note:** This is the only non-molecular classifier

#### 9. WG Methylation (Whole-Genome Methylation)

- **Data Source:** WGBS (Whole-Genome Bisulfite Sequencing)
- **What it detects:** Cancer-specific DNA methylation patterns
- **How it works:**
  - Identifies cfDNA fragments with extreme methylation (>90% or <10%)
  - Scores fragments by cancer-specificity using training data
  - Selects top-ranked cancer-associated methylation patterns
- **Features:** Extreme methylation fragments at cancer-specific genomic locations
- **Model:** Kernel logistic regression (RBF kernel)
- **Key advantages:**
  - ~30 million CpG sites provide massive feature space
  - Tissue-specific patterns enable cancer type prediction
  - Not affected by clonal hematopoiesis
- **Best performer:** 34% sensitivity at 98% specificity

#### 10. Pan-Feature (Multi-Modal Ensemble)

- **Data Source:** Combines all genomic assays
- **What it detects:** Cancer signal by integrating multiple data types
- **How it works:** Meta-classifier that combines predictions from multiple individual classifiers
- **Features:** Scores or predictions from other genomic classifiers (1-7, 9)
- **Purpose:** Test whether multi-assay integration outperforms single assays
- **Note:** Excludes clinical-only classifier (#8)

<p align="center">
  <img src="figure2_roc_curves.png" alt="Performance Graph " width="700" />
</p>

---

**Methodology:**

1. Extract genomic features from each sequencing assay
2. Train 10 classifiers (8 single-feature genomic + 1 multi-feature genomic + 1 clinical-only)
3. Set decision threshold at 98% specificity (on training set)
4. Evaluate sensitivity on validation set

**Key Result:** WG methylation achieved 34% sensitivity at 98% specificity (best single-assay performance)

---

### Analysis 2: Clinical Limit of Detection (Figures 3, S2)

**Research Question:** What is the minimum tumor burden each classifier can detect?

**Data Required:**

- Subset with available tumor tissue (n=409 training, n=113 validation)
- Tumor WGS data for cTAF calculation
- cfDNA classifier scores

**Methodology:**

1. Call somatic variants in tumor tissue
2. Quantify those variants in matched cfDNA to estimate cTAF
3. Fit logistic regression: P(detection) ~ log₁₀(cTAF)
4. Define clinical LOD as cTAF where P(detection) = 50%

**Key Result:** WG methylation LOD ≈ 0.001 (detectable at 0.1% tumor fraction)

---

### Analysis 3: Tumor Shedding Dynamics (Figures 5, S6, S7)

**Research Question:** How does circulating tumor DNA vary by cancer type and stage?

**Data Required:**

- Tumor tissue WGS
- cfDNA targeted sequencing data
- Clinical stage annotations

**Methodology:**

1. Identify somatic variants in tumor (tumor WGS minus WBC WGS)
2. Quantify variant allele frequencies in cfDNA TS data
3. Calculate cTAF = (tumor fraction) × (median tumor VAF)
4. Stratify by cancer type and stage

**Key Result:** cTAF explains 72% of variance in detectability (stage alone is a poor predictor)

---

### Analysis 4: Cancer Type Prediction (Figures 4, S5)

**Research Question:** Can we predict the tissue of origin for detected cancers?

**Data Required:**

- Samples with positive cancer signal from all three assay types (n=127)
- Multi-modal classifier scores (WGBS, TS, WGS)
- True cancer type labels

**Methodology:**

1. Train multinomial classifiers on genomic features
2. Predict cancer type for each detected sample
3. Calculate top-1 and top-3 accuracy

**Key Result:** WG methylation achieved 75% accuracy (substantially better than mutation-based methods)

---

## Statistical Methods

### McNemar's Test

**Purpose:** Compare paired classifier performance on the same samples

**Application:** Test whether WG methylation significantly outperforms SNV classifiers

**Why Paired:** Each sample receives scores from multiple classifiers, creating dependent observations

---
