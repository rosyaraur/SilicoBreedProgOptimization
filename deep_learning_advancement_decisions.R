install.packages("reticulate")
install.packages("keras3")
# Run this once to setup the Python backend:
# install_keras()

library(AlphaSimR)
library(keras3) # Updated library
library(reticulate)

set.seed(42)

# =========================================================================
# PHASE 1: DATA SIMULATION (Targeted Selection Training Data)
# =========================================================================
cat("Initializing AlphaSimR and running 60 cycles with targeted data logging...\n")

founderPop = runMacs(nInd = 200, nChr = 10, segSites = 100)
SP = SimParam$new(founderPop)
SP$addTraitA(nQtlPerChr = 2, mean = 0, var = 1)
SP$setVarE(h2 = 0.35) 

pop_base = newPop(founderPop)
parents = pop_base[1:50] 

X_list = list()
Y_list = list()

for (i in 1:60) {
  f1 = randCross(parents, nCrosses = 200)
  dhLines = makeDH(f1, nDH = 100)
  dhLines = setPheno(dhLines, reps = 1)
  
  # Stage 1 (PYT)
  pytLines = selectWithinFam(dhLines, nInd = 5, use = "pheno")
  pytLines = setPheno(pytLines, reps = 1)
  
  # Stage 2 (AYT)
  aytLines = selectInd(pytLines, nInd = 150, use = "pheno")
  
  # Data Logging Logic
  if (i == 1) {
    X_list[[i]] = pullSegSiteGeno(dhLines)
    Y_list[[i]] = pheno(dhLines)
  } else {
    X_list[[i]] = pullSegSiteGeno(pytLines)
    Y_list[[i]] = pheno(pytLines)
  }
  
  parents = aytLines
}

# Aggregate full training dataset natively in R
X_train = do.call(rbind, X_list)
Y_train = do.call(rbind, Y_list)

X_candidates = pullSegSiteGeno(parents)
mean_G_C60 = mean(gv(parents))

cat(sprintf("Simulation Complete. Extracted %d targeted records.\n", nrow(X_train)))
cat(sprintf("Cycle 60 Mean Genetic Value: %.4f\n\n", mean_G_C60))

# =========================================================================
# PHASE 2: DEEP LEARNING MODEL (Via Keras 3)
# =========================================================================
cat("Training Deep Learning Model via embedded Python (keras3)...\n")

# Define the sequential model 
model <- keras_model_sequential(input_shape = c(ncol(X_train))) %>%
  layer_dense(units = 256, activation = 'relu') %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 128, activation = 'relu') %>%
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 64, activation = 'relu') %>%
  layer_dense(units = 1, activation = 'linear')

# Compile the model
model %>% compile(
  optimizer = optimizer_adam(learning_rate = 0.001),
  loss = 'mse',
  metrics = c('mae')
)

# Train the model 
history <- model %>% fit(
  X_train, Y_train,
  epochs = 15,
  batch_size = 128,
  validation_split = 0.1,
  verbose = 1
)

cat("Model Training Complete.\n\n")

# =========================================================================
# PHASE 3: CYCLE 61 PREDICTIVE MATING
# =========================================================================
cat("Evaluating all possible Cycle 61 combinations...\n")

num_candidates = nrow(X_candidates)
num_markers = ncol(X_candidates)

# Generate all unique pair combinations (150 choose 2)
# combn returns a 2 x 11175 matrix of indices
combos = combn(1:num_candidates, 2)
num_combos = ncol(combos)

cat(sprintf("Simulating %d potential synthetic F1 genotypes...\n", num_combos))

# Pre-allocate matrix for expected genotypes
expected_genos = matrix(0, nrow = num_combos, ncol = num_markers)

for(i in 1:num_combos) {
  p1 = combos[1, i]
  p2 = combos[2, i]
  expected_genos[i, ] = (X_candidates[p1, ] + X_candidates[p2, ]) / 2.0
}

cat("Predicting offspring yield for all combinations...\n")
predictions = model %>% predict(expected_genos, batch_size = 512)

# Rank the crosses
top_indices = order(predictions, decreasing = TRUE)
top_200_idx = top_indices[1:200]

predicted_yields = predictions[top_200_idx]
mean_predicted_yield = mean(predicted_yields)

potential_gain = mean_predicted_yield - mean_G_C60
percent_gain = (potential_gain / abs(mean_G_C60)) * 100

cat("\n==================================================\n")
cat("CYCLE 61 AI MATING RECOMMENDATIONS\n")
cat("==================================================\n")
cat(sprintf("Cycle 60 Baseline Mean G:  %.4f\n", mean_G_C60))
cat(sprintf("Cycle 61 Predicted Mean G: %.4f\n", mean_predicted_yield))
cat(sprintf("Potential Genetic Gain:    +%.4f (+%.2f%%)\n\n", potential_gain, percent_gain))

cat("Top 10 Recommended Crosses:\n")
for (i in 1:10) {
  idx = top_200_idx[i]
  p1 = combos[1, idx]
  p2 = combos[2, idx]
  cat(sprintf("Rank %d: Parent %d x Parent %d | Predicted Yield: %.4f\n", 
              i, p1, p2, predictions[idx]))
}

