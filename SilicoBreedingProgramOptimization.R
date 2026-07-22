# =========================================================================
# COMPREHENSIVE ALPHASIMR MASTER SCRIPT
# Includes: Diversity tracking, Structural Schemes, and Rare Allele Schemes
# =========================================================================

# Install packages if necessary: install.packages(c("AlphaSimR", "ggplot2"))
library(AlphaSimR)
library(ggplot2)

# =========================================================================
# PART 1: GLOBAL SETUP & POPULATION INITIALIZATION
# =========================================================================
set.seed(42) # For reproducible results

# Generate a global base population (200 individuals)
cat("Generating founder population...\n")
founderPop = runMacs(nInd = 200, nChr = 10, segSites = 100)

# Define simulation parameters: 20 QTLs total, Heritability = 0.35
SP = SimParam$new(founderPop)
SP$addTraitA(nQtlPerChr = 2, mean = 0, var = 1)
SP$setVarE(h2 = 0.35) 

pop_base = newPop(founderPop)

# Create the three starting diversity scenarios
# 1. HIGH: 50 completely random unrelated individuals
parents_high = pop_base[1:50]

# 2. MEDIUM: 5 crosses, 10 DH lines each (50 total)
cross_med = randCross(pop_base[51:60], nCrosses = 5)
parents_med = makeDH(cross_med, nDH = 10)

# 3. LOW: 1 cross, 50 DH siblings (50 total)
cross_low = randCross(pop_base[61:62], nCrosses = 1)
parents_low = makeDH(cross_low, nDH = 50)


# =========================================================================
# PART 2: HELPER FUNCTIONS
# =========================================================================

# Calculate Rare Allele Score (RAS) with a strict < 0.05 MAF threshold
calc_ras <- function(pop) {
  geno = pullSegSiteGeno(pop)
  af = colMeans(geno) / 2
  maf = ifelse(af > 0.5, 1 - af, af)
  
  rare_loci = which(maf > 0 & maf < 0.05)
  
  if(length(rare_loci) == 0) return(rep(0, nInd(pop)))
  
  major_allele_dosage = ifelse(af[rare_loci] > 0.5, 2, 0)
  
  ras = apply(geno[, rare_loci, drop=FALSE], 1, function(ind) {
    sum(abs(ind - major_allele_dosage))
  })
  return(ras)
}


# =========================================================================
# PART 3: SIMULATION WRAPPER FUNCTIONS
# =========================================================================

# -------------------------------------------------------------------------
# WRAPPER 1: Single-Cycle Pipeline (Tracks Genetic Variance Loss)
# -------------------------------------------------------------------------
simulate_single_cycle <- function(parents, scenario_name) {
  cat("Running 1-Cycle Pipeline for:", scenario_name, "\n")
  geno_matrix = pullSegSiteGeno(parents)
  mean_dist = mean(dist(geno_matrix))
  
  f1 = randCross(parents, nCrosses = 200)
  dhLines = makeDH(f1, nDH = 100)
  
  dhLines = setPheno(dhLines, reps = 1)
  pytLines = selectWithinFam(dhLines, nInd = 5, use = "pheno")
  
  pytLines = setPheno(pytLines, reps = 1)
  aytLines = selectInd(pytLines, nInd = 100, use = "pheno")
  
  aytLines = setPheno(aytLines, reps = 4)
  eytLines = selectInd(aytLines, nInd = 10, use = "pheno")
  
  results = data.frame(
    Scenario = scenario_name,
    Founder_Dist = round(mean_dist, 2),
    VarG_Parents = round(varG(parents), 4),
    VarG_DH = round(varG(dhLines), 4),
    VarG_PYT = round(varG(pytLines), 4),
    VarG_AYT = round(varG(aytLines), 4),
    VarG_EYT = round(varG(eytLines), 4)
  )
  return(results)
}

# -------------------------------------------------------------------------
# WRAPPER 2: Structural Mating Schemes (Base, A, B)
# -------------------------------------------------------------------------
simulate_structural_schemes <- function(scheme_name, start_parents, n_cycles = 60) {
  parents = start_parents
  results = data.frame(Scheme = character(), Cycle = integer(), MeanG = numeric(), VarG = numeric())
  
  cat("Running Structural:", scheme_name, "...\n")
  
  suppressWarnings({
    for (i in 1:n_cycles) {
      if (scheme_name == "Base") {
        f1 = randCross(parents, nCrosses = 200)
        
      } else if (scheme_name == "Scheme A (High x Med)") {
        high_parents = parents[1:10]
        med_parents = parents[11:50]
        f1 = randCross2(high_parents, med_parents, nCrosses = 200)
        
      } else if (scheme_name == "Scheme B (High x Diverse Med)") {
        high_parents = parents[1:10]
        if (nInd(parents) > 50) {
          med_pool = parents[11:nInd(parents)]
          geno_high = pullSegSiteGeno(high_parents)
          geno_med = pullSegSiteGeno(med_pool)
          dists = apply(geno_med, 1, function(x) mean(sqrt(rowSums(sweep(geno_high, 2, x)^2))))
          diverse_idx = order(dists, decreasing = TRUE)[1:40]
          med_parents = med_pool[diverse_idx]
        } else {
          med_parents = parents[11:50]
        }
        f1 = randCross2(high_parents, med_parents, nCrosses = 200)
      }
      
      dhLines = makeDH(f1, nDH = 100)
      dhLines = setPheno(dhLines, reps = 1)
      pytLines = selectWithinFam(dhLines, nInd = 5, use = "pheno")
      pytLines = setPheno(pytLines, reps = 1)
      
      if (scheme_name == "Scheme B (High x Diverse Med)") {
        next_parents = selectInd(pytLines, nInd = 150, use = "pheno")
      } else {
        next_parents = selectInd(pytLines, nInd = 50, use = "pheno")
      }
      
      results = rbind(results, data.frame(Scheme = scheme_name, Cycle = i, MeanG = meanG(parents), VarG = varG(parents)))
      parents = next_parents
      rm(f1, dhLines, pytLines); gc(verbose = FALSE)
    }
  })
  return(results)
}

# -------------------------------------------------------------------------
# WRAPPER 3: Rare Allele Targeting Schemes (30% Flexibility)
# -------------------------------------------------------------------------
simulate_rare_allele_schemes <- function(scheme_num, start_parents, n_cycles = 60) {
  parents = start_parents
  scheme_label = paste("Rare Scheme", scheme_num)
  results = data.frame(Scheme = character(), Cycle = integer(), MeanG = numeric(), VarG = numeric(), MeanRAS = numeric())
  cat("Running", scheme_label, "(30% Flexibility)...\n")
  
  n_p2_total = 40
  n_ras_target = round(n_p2_total * 0.30) # 12 slots for RAS
  
  suppressWarnings({
    for (i in 1:n_cycles) {
      p1 = parents[1:10] # Top 10 yielders
      
      ras_scores = calc_ras(parents)
      yield_pheno = pheno(parents)
      
      # Step A: Select 12 for RAS
      if (scheme_num == 1) {
        ras_picks = order(ras_scores, decreasing = TRUE)[1:n_ras_target]
      } else if (scheme_num == 2) {
        pool_idx = 1:min(150, nInd(parents))
        best_ras = order(ras_scores[pool_idx], decreasing = TRUE)[1:n_ras_target]
        ras_picks = pool_idx[best_ras]
      } else if (scheme_num == 3) {
        pool_idx = 1:min(50, nInd(parents))
        best_ras = order(ras_scores[pool_idx], decreasing = TRUE)[1:n_ras_target]
        ras_picks = pool_idx[best_ras]
      } else if (scheme_num == 4) {
        mean_ras = mean(ras_scores)
        meets_target = which(ras_scores > mean_ras)
        if (length(meets_target) >= n_ras_target) {
          target_yields = yield_pheno[meets_target]
          best_yields = order(target_yields, decreasing = TRUE)[1:n_ras_target]
          ras_picks = meets_target[best_yields]
        } else {
          ras_picks = meets_target 
        }
      } else if (scheme_num == 5) {
        geno_p1 = pullSegSiteGeno(p1)
        geno_all = pullSegSiteGeno(parents)
        dists = apply(geno_all, 1, function(x) mean(sqrt(rowSums(sweep(geno_p1, 2, x)^2))))
        hybrid_score = rank(ras_scores) + rank(dists)
        ras_picks = order(hybrid_score, decreasing = TRUE)[1:n_ras_target]
      } else if (scheme_num == 6) {
        std_yield = scale(yield_pheno)
        std_ras = scale(ras_scores)
        index_score = (0.5 * std_yield) + (0.5 * std_ras)
        ras_picks = order(index_score, decreasing = TRUE)[1:n_ras_target]
      }
      
      # Step B: Fill the rest (28 slots) with highest yielders
      remaining_idx = setdiff(1:nInd(parents), ras_picks)
      n_to_fill = n_p2_total - length(ras_picks)
      yield_picks = remaining_idx[1:n_to_fill]
      
      p2_idx = c(ras_picks, yield_picks)
      p2 = parents[p2_idx]
      
      # Pipeline
      f1 = randCross2(p1, p2, nCrosses = 200)
      dhLines = makeDH(f1, nDH = 100)
      dhLines = setPheno(dhLines, reps = 1)
      pytLines = selectWithinFam(dhLines, nInd = 5, use = "pheno")
      pytLines = setPheno(pytLines, reps = 1)
      
      next_parents = selectInd(pytLines, nInd = 150, use = "pheno")
      
      results = rbind(results, data.frame(Scheme = scheme_label, Cycle = i, MeanG = meanG(parents), VarG = varG(parents), MeanRAS = mean(ras_scores)))
      parents = next_parents
      rm(f1, dhLines, pytLines); gc(verbose = FALSE)
    }
  })
  return(results)
}


# =========================================================================
# PART 4: EXECUTION BLOCK
# =========================================================================

# 1. Run Initial Diversity Benchmarks (1 Cycle)
cat("\n--- Phase 1: Diversity Base Benchmarks ---\n")
div_results = rbind(
  simulate_single_cycle(parents_high, "High Diversity"),
  simulate_single_cycle(parents_med,  "Medium Diversity"),
  simulate_single_cycle(parents_low,  "Low Diversity")
)
print(div_results)

# 2. Run Structural Schemes (60 Cycles)
cat("\n--- Phase 2: Structural Mating Schemes ---\n")
struct_results = data.frame()
for(scheme in c("Base", "Scheme A (High x Med)", "Scheme B (High x Diverse Med)")) {
  struct_results = rbind(struct_results, simulate_structural_schemes(scheme, parents_med, n_cycles = 60))
}

# 3. Run Rare Allele Targeting Schemes (60 Cycles)
cat("\n--- Phase 3: Rare Allele Schemes (1-6) ---\n")
rare_results = data.frame()
for(s in 1:6) {
  rare_results = rbind(rare_results, simulate_rare_allele_schemes(s, parents_med, n_cycles = 60))
}

# =========================================================================
# PART 5: PLOTTING
# =========================================================================

# Plot Structural Schemes (Genetic Gain)
p_struct = ggplot(struct_results, aes(x = Cycle, y = MeanG, color = Scheme)) +
  geom_line(linewidth = 1) +
  theme_minimal() +
  labs(title = "Structural Schemes: Genetic Gain", y = "Mean Genetic Value", x = "Cycle") +
  theme(legend.position = "bottom")

# Plot Rare Allele Schemes (Genetic Gain)
p_rare_gain = ggplot(rare_results, aes(x = Cycle, y = MeanG, color = Scheme)) +
  geom_line(linewidth = 1) +
  theme_minimal() +
  labs(title = "Rare Allele Schemes: Genetic Gain", y = "Mean Genetic Value", x = "Cycle") +
  theme(legend.position = "bottom")

# Plot Rare Allele Schemes (RAS Preservation)
p_rare_ras = ggplot(rare_results, aes(x = Cycle, y = MeanRAS, color = Scheme)) +
  geom_line(linewidth = 1) +
  theme_minimal() +
  labs(title = "Rare Allele Schemes: RAS Preservation", y = "Mean RAS", x = "Cycle") +
  theme(legend.position = "bottom")

# Display plots
print(p_struct)
print(p_rare_gain)
print(p_rare_ras)