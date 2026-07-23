# In Silico Crop Breeding Optimization using AlphaSimR & Deep Learning (AI) assisted crossing concept

An advanced stochastic simulation of a recurrent Doubled Haploid (DH) crop breeding program coupled with a deep learning predictive mating engine. This project utilizes the `AlphaSimR` package in R to model and evaluate mating strategies, balancing short-term genetic yield with the long-term preservation of genetic diversity. Furthermore, it integrates a neural network (via `keras3` and `reticulate`) to intelligently predict and optimize future crossing blocks.

## Overview

Modern commercial breeding programs face a critical trade-off: aggressive truncation selection maximizes short-term yield but rapidly depletes additive genetic variance, leading to genetic plateaus and the extinction of rare alleles.

This project provides a fully modular R script to test specialized crossing block interventions designed to rescue these rare alleles without crashing the competitive yield of the program. Building on this stochastic foundation, the framework now includes a **machine learning module** that harvests 60 cycles of historical genomic and phenotypic data to train a neural network. This AI engine evaluates thousands of synthetic F1 combinations to recommend elite crosses that maximize genetic gain.

## Prerequisites

To run this hybrid simulation, you will need R installed on your machine, alongside a Python environment configured for deep learning.

**R Packages Required:**

* `AlphaSimR`: For stochastic breeding program simulation.
* `ggplot2`: For visualizing genetic gain and diversity metrics.
* `reticulate` & `keras3`: To bridge R to the Python TensorFlow backend.

You can install these dependencies in R using:

```R
install.packages(c("AlphaSimR", "ggplot2", "reticulate", "keras3"))

# Note: On first run, you may need to initialize the Keras backend:
# keras3::install_keras()

```

## Usage

The entire project is consolidated into a master script, designed to be highly modular.

1. Clone this repository or download the master `.R` script.
2. Open the script in RStudio or your preferred R environment.
3. Run the script. The script is structured into modular phases, allowing you to isolate specific simulations (e.g., Structural Schemes, Rare Allele Schemes, or AI Predictive Mating).
4. Visualizations will automatically generate in your plot viewer, and AI mating recommendations will print to the console at the end of the execution block.

## Project Methodology

### 1. Baseline Pipeline Architecture

The simulation models a complex quantitative trait controlled by 20 QTLs with a baseline heritability of 0.35. The annual pipeline mirrors a real-world multi-stage DH program:

* **Year 1:** Parents are crossed (P1 x P2) to create 200 F1 families.
* **Year 1-2:** 100 Doubled Haploid (DH) lines are derived per cross (20,000 lines total).
* **Year 3 (HDRW):** Visual selection in Head Rows; the top 5 lines per family are advanced.
* **Year 4 (PYT):** Preliminary Yield Trial (1 rep/location). Top lines are recycled as parents for the next recurrent cycle.
* **Year 5-6 (AYT & EYT):** Advanced and Elite Yield Trials with increasing spatial replication (4 and 16 locations) to simulate standard commercial advancement.

### 2. Diversity Benchmarking

The script evaluates how foundational genetic diversity impacts the variance depletion rate. It benchmarks three starting populations:

* **High Diversity:** 50 randomly selected, unrelated founders.
* **Medium Diversity:** 50 founders derived from 5 full-sib families.
* **Low Diversity:** 50 founders derived from a single full-sib cross.

### 3. Structural Mating Designs

To delay the genetic plateau, the crossing block structure is dynamically modified across 60 cycles:

* **Base Scheme:** Random mating among the top 50 strictly selected yielding lines.
* **Scheme A (High x Medium):** Forced crosses between absolute elite lines (Top 10) and moderately yielding lines (Ranks 11-50).
* **Scheme B (High x Diverse Medium):** Genomic distance integration. The top 10 elites are crossed with 40 lines selected from a wider pool (Ranks 11-150) based strictly on their Euclidean genetic distance from the elite group.

### 4. Rare Allele Preservation (RAS)

The core feature of this project is the Rare Allele Score (RAS) algorithm. It tracks population-wide minor allele frequencies (MAF) in real-time, isolating strictly rare polymorphic loci (MAF < 0.05).

To prevent severe yield penalties ("ambulance chasing"), a **30% Flexibility Buffer** is enforced:

* **30% of Parent 2** (12 lines) are selected strictly based on the RAS preservation scheme (e.g., Unrestricted RAS, Yield-Gated RAS, Minimum Thresholds, or Standardized Indices).
* **70% of Parent 2** (28 lines) are filled by the absolute highest-yielding candidates available.
* **Parent 1** remains the top 10 elite yielders.

### 5. AI-Driven Predictive Mating

To optimize future generations, the pipeline embeds a deep neural network to predict the yielding potential of unseen crosses.

* **Data Generation:** Throughout the 60-cycle simulation, segregating site genotypes and phenotypic values are aggregated, building a massive training matrix.
* **Deep Learning Architecture:** A feed-forward Keras neural network (256-128-64 dense layers) uses ReLU activations and dropout regularization to learn non-linear genetic interactions mapping to yield.
* **Combinatoric Simulation:** For Cycle 61, the framework generates all unique possible F1 crosses (150 choose 2 = 11,175 combinations) by calculating mid-parent genomic profiles.
* **Elite Recommendation:** The model scores all 11,175 synthetic crosses in a fraction of a second, ranking them to bypass standard random mating and maximize targeted genetic gain.

## Outputs and Visualization

Upon completion, the script generates two distinct sets of outputs:

**1. ggplot2 Visualizations (Historical Metrics):**

* **Structural Schemes - Genetic Gain:** Maps the mean genetic value over time to show how quickly truncation selection plateaus compared to distance-based mating.
* **Rare Allele Schemes - Genetic Gain:** Evaluates the yield drag associated with different rare-allele rescue interventions.
* **Rare Allele Schemes - RAS Preservation:** Tracks the average number of rare alleles successfully maintained per individual across the breeding cycles.

**2. AI Mating Console Reports (Predictive Metrics):**

* **Cycle 61 Projections:** Calculates the baseline Cycle 60 Mean Genetic Value versus the AI's predicted Mean Genetic Value, computing the potential percentage of genetic gain.
* **Top 10 Crosses:** Outputs a ranked list of the specific Parent 1 x Parent 2 combinations expected to yield the highest-performing offspring.
