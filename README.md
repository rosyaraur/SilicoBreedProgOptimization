# In Silico Crop Breeding Optimization using AlphaSimR

An advanced stochastic simulation of a recurrent Doubled Haploid (DH) crop breeding program. This project utilizes the `AlphaSimR` package in R to model, evaluate, and optimize mating and selection strategies. The primary goal is to balance short-term genetic gain (yield) with the long-term preservation of genetic diversity and rare beneficial alleles.

## Table of Contents

* [Overview](https://www.google.com/search?q=%23overview)
* [Prerequisites](https://www.google.com/search?q=%23prerequisites)
* [Usage](https://www.google.com/search?q=%23usage)
* [Project Methodology](https://www.google.com/search?q=%23project-methodology)
* [1. Baseline Pipeline Architecture](https://www.google.com/search?q=%231-baseline-pipeline-architecture)
* [2. Diversity Benchmarking](https://www.google.com/search?q=%232-diversity-benchmarking)
* [3. Structural Mating Designs](https://www.google.com/search?q=%233-structural-mating-designs)
* [4. Rare Allele Preservation (RAS)](https://www.google.com/search?q=%234-rare-allele-preservation-ras)


* [Outputs and Visualization](https://www.google.com/search?q=%23outputs-and-visualization)

---

## Overview

Modern commercial breeding programs face a critical trade-off: aggressive truncation selection maximizes short-term yield but rapidly depletes additive genetic variance, leading to genetic plateaus and the extinction of rare alleles. This project provides a fully modular `AlphaSimR` script to test specialized crossing block interventions designed to rescue these rare alleles without crashing the competitive yield of the program.

## Prerequisites

To run this simulation, you will need **R** installed on your machine along with the following packages:

* `AlphaSimR`: For stochastic breeding program simulation.
* `ggplot2`: For visualizing genetic gain and diversity metrics.

You can install these dependencies in R using:

```R
install.packages(c("AlphaSimR", "ggplot2"))

```

## Usage

The entire project is consolidated into a single master script.

1. Clone this repository or download the master `.R` script.
2. Open the script in RStudio or your preferred R environment.
3. Run the script from top to bottom. The script is structured into modular phases, allowing you to run specific simulations (e.g., just the Structural Schemes or just the Rare Allele Schemes) as needed.
4. Visualizations will automatically generate in your plot viewer at the end of the execution block.

---

## Project Methodology

### 1. Baseline Pipeline Architecture

The simulation models a complex quantitative trait controlled by 20 QTLs with a baseline heritability of 0.35. The annual pipeline mirrors a real-world multi-stage DH program:

* **Year 1:** Parents are crossed (P1 x P2) to create 200 F1 families.
* **Year 1-2:** 100 Doubled Haploid (DH) lines are derived per cross (20,000 lines total).
* **Year 3 (HDRW):** Visual selection in Head Rows; the top 5 lines per family are advanced.
* **Year 4 (PYT):** Preliminary Yield Trial (1 rep/location). Top lines are recycled as parents for the next recurrent cycle.
* **Year 5-6 (AYT & EYT):** Advanced and Elite Yield Trials with increasing spatial replication (4 and 16 locations, respectively) to simulate standard commercial advancement.

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

The core feature of this project is the **Rare Allele Score (RAS)** algorithm. It tracks population-wide minor allele frequencies (MAF) in real-time, isolating strictly rare polymorphic loci (MAF < 0.05).

Six targeted selection schemes are evaluated to rescue these alleles. To prevent severe yield penalties ("ambulance chasing"), a **30% Flexibility Buffer** is enforced:

* **30% of Parent 2 (12 lines)** are selected strictly based on the RAS preservation scheme (e.g., Unrestricted RAS, Yield-Gated RAS, Minimum Thresholds, or Standardized Indices).
* **70% of Parent 2 (28 lines)** are filled by the absolute highest-yielding candidates available.
* Parent 1 remains the top 10 elite yielders.

---

## Outputs and Visualization

Upon completion, the script generates three primary `ggplot2` visualizations comparing the 60-cycle trajectories:

1. **Structural Schemes - Genetic Gain:** Maps the mean genetic value over time to show how quickly truncation selection plateaus compared to distance-based mating.
2. **Rare Allele Schemes - Genetic Gain:** Evaluates the yield drag associated with the different rare-allele rescue interventions.
3. **Rare Allele Schemes - RAS Preservation:** Tracks the average number of rare alleles successfully maintained per individual across the breeding cycles.
