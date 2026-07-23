# title: "AlphaSimR & Keras Deep Learning Breeding Pipeline"

### **1. Environment Setup & Initialization**
Before beginning the simulation, we must load the necessary libraries and set a random seed to ensure our results are reproducible.

*   **`AlphaSimR`**: The core engine driving the quantitative genetics and breeding program simulation.
*   **`keras3` & `reticulate`**: The bridge that allows R to seamlessly communicate with the Python-based TensorFlow backend to build and train our deep learning model. 

```{r setup, message=FALSE, warning=FALSE}
library(AlphaSimR)
library(keras3)
library(reticulate)

# Set seed for reproducibility
set.seed(42)

```

### **2. Phase 1: Simulation Setup & Trait Architecture**

In this step, we initialize the biological foundation of our virtual breeding program.

* **Founder Population**: We use the `runMacs` function to generate a historical base population of 200 individuals across 10 chromosomes, with 100 segregating sites per chromosome.
* **Trait Architecture**: We define a single quantitative trait (e.g., yield) controlled by 2 Quantitative Trait Loci (QTLs) per chromosome. We apply an environmental variance to achieve a narrow-sense heritability of 0.35, simulating a moderately complex trait.
* **Initial Parents**: We establish the `pop_base` and select the first 50 individuals to serve as the initial crossing block.

```{r simulation_setup}
cat("Initializing AlphaSimR...\n")

founderPop = runMacs(nInd = 200, nChr = 10, segSites = 100)
SP = SimParam$new(founderPop)
SP$addTraitA(nQtlPerChr = 2, mean = 0, var = 1)
SP$setVarE(h2 = 0.35) 

pop_base = newPop(founderPop)
parents = pop_base 

# Initialize lists to store historical data
X_list = list()
Y_list = list()

```

### **3. Phase 1: The 60-Cycle Breeding Loop & Data Logging**

This is the workhorse of the simulation. We loop through 60 cycles of a Doubled Haploid (DH) breeding pipeline, generating historical training data for our AI model.

* **Crossing & DH Creation**: In each cycle, we make 200 random crosses from the parent pool and generate 100 DH lines from those F1s.
* **Multi-Stage Selection**: The DH lines undergo a Preliminary Yield Trial (PYT), where the top 5 individuals per family are selected based on phenotypic performance. These PYT lines are then advanced to an Advanced Yield Trial (AYT), where the top 150 individuals overall are selected to become the parents for the next cycle.
* **Genomic & Phenotypic Extraction**: At each cycle, we extract the segregating site genotypes (`X`) and their corresponding phenotypic values (`Y`), appending them to a growing list. By the end of 60 cycles, this forms a massive, multi-generational dataset.
* **Aggregation**: We use `do.call(rbind, ...)` to collapse these lists into large, unified matrices (`X_train` and `Y_train`) ready for machine learning.

```{r breeding_loop}
cat("Running 60 cycles with targeted data logging...\n")

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

```

### **4. Phase 2: Deep Learning Model Architecture**

With our 60-cycle dataset assembled, we construct a deep neural network to learn the complex, non-linear mapping between a plant's genomic profile and its final yield.

* **Sequential API**: We build a feed-forward neural network using `keras_model_sequential()`.
* **Hidden Layers**: The network consists of three dense (fully connected) layers with 256, 128, and 64 neurons, respectively. They use the `relu` (Rectified Linear Unit) activation function to capture non-linear genetic interactions (like epistasis).
* **Regularization**: We inject `layer_dropout` steps (dropping 30% and 20% of connections) to prevent the model from overfitting to the historical training data.
* **Output Layer**: A final dense layer with a single neuron and a `linear` activation function outputs our continuous yield prediction.

```{r build_model}
cat("Constructing Deep Learning Model via embedded Python (keras3)...\n")

# Define the sequential model 
model <- keras_model_sequential(input_shape = c(ncol(X_train))) %>%
  layer_dense(units = 256, activation = 'relu') %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 128, activation = 'relu') %>%
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 64, activation = 'relu') %>%
  layer_dense(units = 1, activation = 'linear')

```

### **5. Phase 2: Model Compilation & Training**

Before the model can learn, it must be compiled with a strategy for measuring and minimizing its errors.

* **Optimizer & Loss**: We use the `adam` optimizer (with a learning rate of 0.001) to update the network weights, and Mean Squared Error (`mse`) as our loss function to heavily penalize large prediction errors.
* **Model Fitting**: The `fit` function trains the model for 15 epochs. We pass it our historical genotypes (`X_train`) and phenotypes (`Y_train`), processing them in batches of 128. We withhold 10% of the data as a validation set to monitor performance during training.

```{r train_model}
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

```

### **6. Phase 3: Generating Cycle 61 Synthetic Genotypes**

Instead of selecting parents randomly for Cycle 61, we use our trained model to perform predictive mating.

* **Combinatorics**: We use the `combn` function to generate every possible unique pairing of our 150 Cycle 60 candidate parents (resulting in 11,175 potential crosses).
* **Mid-Parent Genotypes**: We pre-allocate an empty matrix and calculate the expected genotype of every possible F1 offspring. Because these are inbred DH lines, the expected F1 genotype is simply the mathematical average of the two parents' genomic matrices.

```{r generate_crosses}
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

```

### **7. Phase 3: AI Yield Prediction & Elite Selection**

In the final step, we unleash the trained Keras model on our synthetic F1 combinations.

* **Scoring**: The `predict` function feeds all 11,175 expected genotypes into the neural network, generating a predicted yield for every possible cross in a fraction of a second.
* **Ranking & Selection**: We use the `order` function to rank these predictions from highest to lowest. We extract the indices of the top 200 crosses, representing the AI's mathematically optimized crossing recommendations.
* **Output Metrics**: Finally, the script calculates the projected genetic gain over the Cycle 60 baseline and prints the top 10 recommended parent pairings to the console.

```{r predict_and_select}
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

```

```

```