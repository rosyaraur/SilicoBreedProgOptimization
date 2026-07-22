## In Silico Breeding Program Optimization

The primary objective of this project is to use the stochastic simulation software `AlphaSimR` to design, model, and evaluate a doubled haploid (DH) crop breeding program. The project aims to identify mating and selection strategies that maximize long-term genetic gain while mitigating the rapid loss of genetic diversity and the extinction of rare beneficial alleles.

---

### Phase 1: Baseline Pipeline Architecture

The project began by establishing a computational baseline that mirrors a standard, real-world DH breeding pipeline.

* **Founder Population Simulation:** Generated a base population using the coalescent simulator (Macs) to mimic historical linkage disequilibrium.
* **Trait Architecture:** Modeled a complex quantitative trait controlled by 20 QTLs (2 per chromosome) with a moderate baseline heritability ($h^2 = 0.35$).
* **Pipeline Execution:** Replicated a multi-year cyclic pipeline: Crossing (P1 $\times$ P2) $\rightarrow$ F1 generation $\rightarrow$ DH line derivation $\rightarrow$ Head Row (HDRW) visual selection $\rightarrow$ Multi-environment yield trials (PYT, AYT, EYT).
* **Recurrent Selection Loop:** Implemented a continuous looping structure to simulate 60 to 100 cycles of recurrent truncation selection, feeding the top Preliminary Yield Trial (PYT) performers back into the crossing block.

---

### Phase 2: Diversity Benchmarking

To understand the initial genetic boundaries of the population, we established baseline benchmarks for genetic diversity.

* **Diversity Stratification:** Constructed three distinct founder populations: High (unrelated individuals), Medium (5 full-sib families), and Low (1 full-sib family).
* **Variance Tracking:** Measured the additive genetic variance ($\sigma^2_A$) at each specific stage of the pipeline to quantify the bottleneck effect of standard truncation selection.
* **Distance Metrics:** Utilized Euclidean distance on marker genotypes to establish quantitative baselines for founder relatedness.

---

### Phase 3: Structural Mating Design Evaluation

Recognizing the rapid allele fixation and genetic plateau inherent in the baseline approach, we modified the crossing block structure to evaluate alternative parental pairings.

* **Base Scheme:** Random mating among the top 50 strictly selected yielding lines.
* **Scheme A (High $\times$ Medium):** Stratified the crossing block to force matings between the absolute elite lines (Top 10) and moderately yielding lines (Ranks 11-50) to introduce slight phenotypic variance.
* **Scheme B (High $\times$ Diverse Medium):** Integrated genomic distance into the crossing block. Evaluated a wide pool of medium-yielding lines (Ranks 11-150) and selected the 40 most genetically distant individuals to cross with the Top 10 elite lines.

---

### Phase 4: Rare Allele Preservation Strategy

The final and most complex phase involved designing active interventions to rescue rare minor alleles from extinction without crashing the short-term yield potential of the program.

* **Dynamic Rare Allele Scoring (RAS):** Developed a custom algorithmic function to track population-wide allele frequencies in real-time, isolating strictly rare polymorphic loci ($MAF < 0.05$) and scoring individuals based on their rare allele load.
* **Hybrid Selection Schemes:** Designed six targeted selection schemes for Parent 2 (P2) ranging from unrestricted RAS maximization to sophisticated, standardized Yield + RAS indices.
* **The 30% Flexibility Buffer:** To prevent "ambulance chasing" (where the program severely penalizes yield to rescue a single allele), a split-selection protocol was implemented. Exactly 30% of the P2 crossing block is selected to satisfy the strict rare allele criteria, while the remaining 70% is filled by the highest-yielding candidates, ensuring a pragmatic balance between introgression and genetic gain.

---

### Phase 5: Software Consolidation

The conceptual methodologies were unified into a single, modular R script relying on `AlphaSimR` and `ggplot2` to allow for rapid, reproducible execution and visual comparison of genetic gain and allele preservation trajectories across all schemes.

Here is the comprehensive function-wise documentation for the AlphaSimR simulation script, detailing the objectives and underlying methodologies for each custom function.


# R Functions 

### 1. `calc_ras(pop)`

**Objective**
To quantify the number of strictly rare alleles carried by each individual within a given population, providing a quantitative metric (Rare Allele Score, or RAS) to guide diversity-preserving selection.

**Methodology**

* **Genotype Extraction:** Pulls the segregating site genotype matrix (values of 0, 1, or 2) for the entire provided population.
* **MAF Calculation:** Calculates the allele frequency of the '1' allele at each locus, then converts this to the Minor Allele Frequency (MAF).
* **Threshold Filtering:** Identifies loci where the allele is both segregating and strictly rare ($0 < \text{MAF} < 0.05$). If no rare alleles are present (due to fixation or loss), the function immediately returns a score of 0 for all individuals.
* **Scoring:** For the identified rare loci, the function determines the dominant (major) allele dosage. It then calculates the absolute deviation from this major dosage for each individual, effectively counting how many rare (minor) alleles they possess.

---

### 2. `simulate_single_cycle(parents, scenario_name)`

**Objective**
To benchmark initial genetic diversity and track the phase-by-phase depletion of additive genetic variance during a single, standard execution of a doubled haploid (DH) breeding pipeline.

**Methodology**

* **Initial Diversity Measurement:** Computes the mean Euclidean genetic distance among the founder `parents` based on their marker genotypes.
* **Pipeline Execution:** Runs the classic DH pipeline simulating Years 1 through 6:
* Creates 200 random F1 crosses.
* Generates 100 DH lines per cross (20,000 total).
* Applies single-rep phenotyping (HDRW) and strictly selects the top 5 lines within each family.
* Simulates multi-environment yield trials (PYT, AYT, EYT) with increasing replication accuracy, applying strict truncation selection at each stage (100 $\rightarrow$ 10 $\rightarrow$ 1).


* **Variance Tracking:** Calculates and logs the additive genetic variance (`VarG`) of the population at the parental, DH, PYT, AYT, and EYT stages to map the bottleneck effect of selection.

---

### 3. `simulate_structural_schemes(scheme_name, start_parents, n_cycles)`

**Objective**
To compare the long-term genetic gain (yield) and variance preservation of standard truncation selection against structurally modified mating designs across multiple cycles of recurrent selection.

**Methodology**

* **Mating Design Logic:** Implements three distinct parental crossing strategies:
* **Base:** Random mating among the top 50 yielding lines.
* **Scheme A (High x Med):** Splits the top 50 yielders into a "High" tier (1-10) and a "Medium" tier (11-50), forcing crosses strictly between the two tiers.
* **Scheme B (High x Diverse Med):** Selects the top 10 as "High". For the second parent pool, it evaluates a wider pool of medium yielders (ranks 11-150), calculates their Euclidean genetic distance to the High tier, and selects the 40 most genetically distinct lines to use as crossing partners.


* **Cyclic Execution:** Executes the DH pipeline (F1 $\rightarrow$ DH $\rightarrow$ HDRW $\rightarrow$ PYT) over the specified `n_cycles`.
* **Population Update:** Phenotypically selects the parent pool from the PYT stage to drive the next cycle, logging the mean genetic value (`MeanG`) and genetic variance (`VarG`) at the start of each cycle.

---

### 4. `simulate_rare_allele_schemes(scheme_num, start_parents, n_cycles)`

**Objective**
To evaluate specialized selection indices aimed at rescuing rare alleles from extinction without severely compromising the program's overall genetic gain, utilizing a flexible 30% inclusion buffer.

**Methodology**

* **Fixed Elite Parent:** Parent 1 (P1) is rigidly set as the top 10 highest-yielding lines from the available candidate pool.
* **Targeted Parent 2 (P2) Selection (30% Flexibility):** P2 consists of 40 lines, selected in two distinct steps to balance yield and diversity:
* **Step A (Rare Allele Target):** Exactly 12 lines (30%) are selected based on the specific rules of the chosen scheme (Schemes 1-6). These schemes variously filter the population using raw RAS maximization, yield-gated RAS (e.g., must be in the top 50 or 150 for yield), minimum RAS thresholds, distance hybrids, or a 50/50 standardized Yield + RAS index.
* **Step B (Yield Backfill):** The remaining 28 lines (70%) are filled strictly by the highest-yielding candidates left in the pool, ensuring strong forward momentum for the trait.


* **Cyclic Execution & Tracking:** Crosses P1 with the dynamically constructed P2, executes the pipeline to the PYT stage, and selects the top 150 lines as the candidate pool for the next cycle. Logs `MeanG`, `VarG`, and the `MeanRAS` to track the successful preservation of the targeted minor alleles.