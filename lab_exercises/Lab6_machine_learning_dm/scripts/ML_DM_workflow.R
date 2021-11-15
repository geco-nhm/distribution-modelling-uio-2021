###################################################################
# BIOS 5211/9211 - Autumn, 2021
# Lasse T. Keetz, NHM, University of Oslo, Nov 2021
#
# Lab 6: Building a Machine Learning workflow
###################################################################

##########################
### 1. Library imports ###
##########################

package_names = c("caret",        # ML tools
                  "ranger",       # Random Forest classifier
                  "randomForest", # Diff. RF classifier (needed for caret)
                  "raster",       # Spatial data, raster operations
                  "rstudioapi",   # Automatically set work dir.
                  "rgbif",        # GBIF R interface, species data
                  "MLmetrics",    # Machine learning metrics (AUC, etc.)
                  "RColorBrewer", # Color palettes
                  "ggplot2"       # Nicer plots
                  )

### Install and load packages
for(p_name in package_names) {
  if(p_name %in% rownames(installed.packages())){
    require(p_name, character.only = TRUE)
  }
  else {
    install.packages(p_name)
    require(p_name, character.only = TRUE)
  }
}

### Set working directory to current document's location
current_path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_path)

### Set random seed, ensures that steps involving randomness
# will always give exactly the same outcome (--> reproducibility)
set.seed(112)

###################################################################
###################################################################

################################
####### 2. Prepare data ########
################################

############### Create stack of predictor rasters #################

# Path to raster (.tif) files
pred_path <- paste0(dirname(dirname(dirname(current_path))),
                    "/lab_data/nor_predictor_rasters/2km")

# Store all paths to files ending in ".tif" in a vector
raster_fnames = list.files(path = pred_path, pattern = "\\.tif$",
                           ignore.case = TRUE, full.names = TRUE)
# Create a stack of all individual rasters
pred_stack = raster::stack(raster_fnames)
pred_stack

# Plot DEM to test if everything works as expected
raster::plot(pred_stack[["dem100"]])

### Factorize categorical predictors
cat_var_names = c("ar50_artype",  # Land cover type
                  "geo_norge123"  # Nutrient content in bedrock
                  )

# Convert the corresponding raster values, see documentation
for(cat_var in cat_var_names) {
  pred_stack[[cat_var]] = as.factor(pred_stack[[cat_var]])
}

# Quick test if it worked
# Class of pred. column 'factor'? True/False
for(ras_idx in 1:nlayers(pred_stack)) {
  print(paste0(names(pred_stack)[ras_idx],": ",
               is.factor(pred_stack[[ras_idx]])))
}


################## Load target (species) data #####################

### Pick your species here!
species_name = "Picea sitchensis"  # Group 1
species_name = "Silene acaulis" # Group 2

key <- rgbif::name_backbone(name = species_name,
                            kingdom = "Plantae")$speciesKey

sp <- rgbif::occ_search(taxonKey = key,
                        hasCoordinate = TRUE,
                        country = "NO",
                        limit = 6000)
### Print info
print(sp)

# Convert lat-long coordinates to coordinate system of environmental raster data
occ_points <- data.frame(x = sp$data$decimalLongitude,
                         y = sp$data$decimalLatitude)

# Convert data.frame to SpatialPoints object
occ_points <- sp::SpatialPoints(occ_points,
                                proj4string = CRS("+proj=longlat +datum=WGS84"))
occ_UTM33 <- sp::spTransform(occ_points,
                             CRS("+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs"))
sp$data$x <- occ_UTM33$x
sp$data$y <- occ_UTM33$y

# Make background raster for Norway 
# Value=1 where we have DEM observations, NA outside
norway_mask <- pred_stack[["dem100"]]/pred_stack[["dem100"]]
plot(norway_mask)

# raster with counts of occurrences in each cell of Norway
occ_ras <- raster::rasterize(x = occ_UTM33,
                             y = norway_mask,
                             fun = 'count',
                             background = 0)
# Mask raster by setting values outside of Norway bounds to NA
occ_ras <- occ_ras*norway_mask

# Zoom in to visualize occurrences, random location
plot(occ_ras,
     xlim = c(-40000,-10000),
     ylim = c(6630000, 6650000))
summary(occ_ras)


################## Create feature matrix (df) #####################

### Presences ###
# Extract indices of raster matrix where species occurs
presence_indices <- which(values(occ_ras)==1)

### Absences ###
# Extract indices of raster matrix where species does not occur
absence_indices <- which(values(occ_ras)==0)
# Create random sub-sample with the same amount of absence points as there are
# presence points (OBS! Think about potential problems with this approach!)
absences_sample_indices <- sample(absence_indices,
                                  size = length(presence_indices))

# Combine presences, absences and environmental data in a data.frame
training_indices <- c(presence_indices, absences_sample_indices)
training_xy <- coordinates(occ_ras)[training_indices, ]
training_data <- data.frame(training_xy, 
                            presence = values(occ_ras)[training_indices],
                            extract(pred_stack, training_xy))

# Remove rows that contain NA values in any of the columns
# OBS! Also debatable! Other techniques (e.g., imputation) exist
training_data <- training_data[complete.cases(training_data),]

# Convert factor variables
for(cat_var in cat_var_names) {
  training_data[[cat_var]] <- factor(training_data[[cat_var]])
}
# Print data types
sapply(training_data, class)

# Print first rows to inspect
head(training_data)
# How many samples remained? Balance between pres/abs okay?
print(paste0("Number of presences: ",
             length(which(training_data$presence==1))))
print(paste0("Number of absences: ",
             length(which(training_data$presence==0))))

# Plot norway map
raster::plot(norway_mask, 
             col = rgb(red = 0.9, green = 0.9, blue = 0.9, alpha=1))
# Create a vector with colors based on pres/abs
# --> presences in red, absences in blue
colors <- c(rgb(red = 1, green = 0, blue = 0, alpha = 0.2),
            rgb(red = 0, green = 0, blue = 1, alpha = 0.2))[training_data$presence+1]
# Add points
points(training_data[,c("x","y")],
       col = colors, pch = 16, cex=0.3)



#################################################################
### 3. Splitting data
#################################################################


### 10-fold cross-validation
cv_indices <- caret::createFolds(
  y = training_data$presence,
  k = 10,
  returnTrain = TRUE
)

### Bootstrapping
bs_indices <- caret::createResample(
  y = training_data$presence,
  times = 10
)
# Print table, note that some indices are contained multiple times
table(bs_indices$Resample01)

### State-of-the-art for DMs: spatial cross-validation!
# (Out of the scope for this lab)

#################################################################
#################################################################

### 4. Model training

# We will use customized functions for model training. Note that
# the following steps can be performed using the caret package
# directly (-> see documentation for trainControl() and train()).
# However, this comes at the cost of flexibility as you can only
# tune a subset of the hyperparameters for each included package
# / model. For a list of all included models (packages), see:
# https://topepo.github.io/caret/available-models.html

#################################################################
#################################################################

################################
######### Random Forest ########
################################

### 3.1 Hyperparameter tuning ###
n_features <- nlayers(pred_stack)
sqrt_n_featutres <- as.integer(sqrt(n_features))

# Define Tuning grid
tune_param_grid <- base::expand.grid(
  # Number of randomly drawn vars used for splitting in each node
  # Default: square root of number of features
  mtry = sqrt_n_featutres-1:sqrt_n_featutres+1,
  # Maximum depth of trees, default 0 -> no restriction
  max.depth = seq(0, 10, by = 2)
  ### See "?ranger" for a full list of tunable parameters. E.g.:
  # Number of trees to grow, default: 500
  # num.trees = c(200, 500)
)
# Print result
tune_param_grid

### Instantiate tuning result data frame
# Create vector with names for the folds
fold_names = rep("F", 10)
for(idx in 1:length(fold_names)) {
  fold_names[idx] = paste0(fold_names[idx], as.character(idx))
}

# Retrieve tuneable hyperparameter names
param_names <- colnames(tune_param_grid)

# Define metric names
metric_names <- c("AUC", "Accuracy", "Kappa",
                  "Precision", "Recall", "F1_score")

# Determine number of total columns of result data frame
length_inner_df = length(param_names) + length(metric_names)

# Instantiate the result data frame
hyperparam_df <- data.frame(
  matrix(ncol = length_inner_df, nrow = 0)
  )
# Set the column names to the hyperparameter names and the metrics
colnames(hyperparam_df) <- c(param_names, metric_names)

# Instantiate a list to store the results of each training iteration
rf_tuning_list = sapply(fold_names, function(x) NULL)
# Each list entry will contain one result data frame
for(idx in 1:length(rf_tuning_list)) {
  rf_tuning_list[[idx]] = hyperparam_df
}
rf_tuning_list

##### CV tuning #####

### SET SPLITTING METHOD HERE! One of: cv_indices, bs_indices
split_indices = cv_indices

### Start training
# Index keeping track of the current fold
fold_idx = 1
# Index keeping track of the number of iterations
current_iter = 1
# Determine total number of iterations 
# (-> how often the foor loop "starts over")
total_iterations = length(split_indices) * nrow(tune_param_grid)
# Store predictor variable names in a vector
predictor_names = names(pred_stack)

### Outer for-loop: run through each training sub-set (e.g., fold in cv)
for(fold_indices in split_indices) {

  ### Inner for-loop: run through list of hyperparam. combinations
  for(param_grid_idx in 1:nrow(tune_param_grid)) {
    
    # X -> feature matrix (predictors), y -> target variable (presence/absence)
    X_train = training_data[fold_indices, predictor_names]
    y_train = training_data[fold_indices, "presence"]
    
    X_test = training_data[-fold_indices, predictor_names]
    y_test = training_data[-fold_indices, "presence"]
    
    ### Train random forest on current training set
    current_rf <- ranger::ranger(
      x = X_train,
      y = y_train,
      classification = TRUE,
      probability = TRUE, # Return class probabilities (0-1)
      num.trees = 200, # tune_param_grid[param_grid_idx, "num.trees"],
      mtry = tune_param_grid[param_grid_idx, "mtry"],
      max.depth = tune_param_grid[param_grid_idx, "max.depth"],
      importance = "permutation"
    )
    
    ### Make predictions
    current_predictions = predict(
      current_rf, 
      X_test
    )
    
    ### Transform the probabilities to actual class predictions
    # Instantiate empty vector to store the class predictions
    class_predictions <- vector(
      mode="numeric", 
      length=nrow(current_predictions$predictions)
      )
    
    # Assign class that is more likely (prob. > 0.5)
    cut_off_val = 0.5
    for(idx in 1:nrow(current_predictions$predictions)) {
      if(current_predictions$predictions[idx, 1] < cut_off_val){
        class_predictions[idx] = 0
      } else {
        class_predictions[idx] = 1
      }
    }
    
    ### Generate confusion matrix, metrics like Accuracy, F1-score, ...
    cur_conf_mat <- caret::confusionMatrix(
      as.factor(class_predictions),
      as.factor(y_test)
      )
    
    ### Also calculate the AUC value
    cur_auc = MLmetrics::AUC(
      y_true = as.factor(y_test),
      y_pred = current_predictions$predictions[, 1]
    )
    
    ### Write current hyperparams to data frame
    rf_tuning_list[[fold_idx]][param_grid_idx, "mtry"] <-
      tune_param_grid[param_grid_idx, "mtry"]
    
    #rf_tuning_list[[fold_idx]][param_grid_idx, "num.trees"] <-
     # tune_param_grid[param_grid_idx, "num.trees"]
    
    rf_tuning_list[[fold_idx]][param_grid_idx, "max.depth"] <-
      tune_param_grid[param_grid_idx, "max.depth"]
    
    ### Add classification results to the result data frame
    rf_tuning_list[[fold_idx]][param_grid_idx, "AUC"] <-
      cur_auc
    rf_tuning_list[[fold_idx]][param_grid_idx, "Accuracy"] <-
      cur_conf_mat$overall["Accuracy"]
    rf_tuning_list[[fold_idx]][param_grid_idx, "Kappa"] <-
      cur_conf_mat$overall["Kappa"]
    rf_tuning_list[[fold_idx]][param_grid_idx, "Precision"] <-
      cur_conf_mat$byClass["Precision"]
    rf_tuning_list[[fold_idx]][param_grid_idx, "Recall"] <-
      cur_conf_mat$byClass["Recall"]
    rf_tuning_list[[fold_idx]][param_grid_idx, "F1_score"] <-
      cur_conf_mat$byClass["F1"]
    
    print(paste0("Current iter: ",current_iter," (out of ", total_iterations,")"))
    current_iter = current_iter + 1
    
  }
  
  print(paste0("Done with fold: ", fold_idx))
  fold_idx = fold_idx + 1
  
}


### Calculate mean performance ###
mean_results_df <- rf_tuning_list[[1]]
for(idx in 2:length(rf_tuning_list)) {
  mean_results_df <- mean_results_df + rf_tuning_list[[idx]]
}
mean_results_df <- mean_results_df / length(rf_tuning_list)

### Transform to matrix ###
mean_results_mat <- matrix(
  nrow = nrow(unique(tune_param_grid["mtry"])),
  ncol = nrow(unique(tune_param_grid["max.depth"]))
  )

# Unique hyperparam values from tuning grid
unique_mtry <- unique(tune_param_grid["mtry"])
unique_maxd <- unique(tune_param_grid["max.depth"])

# Name columns and rows of the matrix (with hyperparam values)
rownames(mean_results_mat) <- sapply(unique_mtry, as.character)
colnames(mean_results_mat) <- sapply(unique_maxd, as.character)

### Retrieve AUC values for each hyperparameter
### combination in averaged results data frame
for(i in 1:nrow(mean_results_mat)) {
  for(j in 1:ncol(mean_results_mat)) {
    df_row_idx = intersect(
      which(mean_results_df["mtry"] == rownames(mean_results_mat)[i]),
      which(mean_results_df["max.depth"] == colnames(mean_results_mat)[j])
      )
    mean_results_mat[i, j] <- mean_results_df[df_row_idx, "AUC"]
  }
}

### Plot results ###
plot_heatmap <- function() {
  dev.new()
  # Create heatmap, darker red colors -> higher AUC values
  heatmap(mean_results_mat,
          main = "Random forest - hyperparameter tuning", 
          xlab = "max.depth",
          ylab = "mtry",
          cexRow = 1.5, cexCol = 1.5,
          col = colorRampPalette(brewer.pal(8, "Oranges"))(10),
          Rowv = NA, Colv=NA)
  legend(x="topright", 
         title = "AUC",
         legend=c(
           paste0("min - ", round(min(mean_results_mat),3)),
           paste0("mean - ", round(mean(mean_results_mat),3)),
           paste0("max - ", round(max(mean_results_mat),3))),
         fill = colorRampPalette(brewer.pal(8, "Oranges"))(3)
  )
}

# Call function
plot_heatmap()

#####
# Pick model with best hyperparameters, plot performance across folds ###
#####

# Array indices pointing to maximum AUC in mean results mat
best_idx = which(mean_results_mat == max(mean_results_mat),
                 arr.ind = T)
# Retrieve correspoding mtry and max.depth values
best_mtry = as.integer(rownames(mean_results_mat)[best_idx[1]])
best_maxd = as.integer(colnames(mean_results_mat)[best_idx[2]])

### Instantiate matrix to store metrics for each training set
fold_res_mat = matrix(nrow = length(rf_tuning_list),
                      ncol = length(metric_names))
colnames(fold_res_mat) <- metric_names

### Fill matrix with metric values from best hyperparams
for(idx in 1:length(rf_tuning_list)){
  cur_row <- intersect(
    which(rf_tuning_list[[idx]]$mtry == best_mtry),
    which(rf_tuning_list[[idx]]$max.depth == best_maxd)
    )
  for(metric in metric_names){
    fold_res_mat[idx, metric] <- 
      rf_tuning_list[[idx]][cur_row, metric]
  }
}

### Plot results as a boxplot ###
plot_boxplot <- function() {
  dev.new()
  boxplot(fold_res_mat, 
          main = "Best RF model - 10-fold-cv performance",
          col = "#E6E6FA")
}

plot_boxplot()

#################################################################
### Variable importance ###
#################################################################

### Retrain RF models with best hyperparameters, store objects

rf_list = list()
idx = 1
for(fold_indices in split_indices) {
    # X -> feature matrix (predictors), y -> target variable (presence/absence)
    X_train = training_data[fold_indices, predictor_names]
    y_train = training_data[fold_indices, "presence"]
    
    X_test = training_data[-fold_indices, predictor_names]
    y_test = training_data[-fold_indices, "presence"]
    
    ### Train random forest on current training set
    current_rf <- ranger::ranger(
      x = X_train,
      y = y_train,
      classification = TRUE,
      probability = TRUE,
      num.trees = 500,
      mtry = best_mtry,
      max.depth = best_maxd,
      importance = "permutation"
    )
    
    ### Make predictions
    current_predictions = predict(
      current_rf, 
      X_test
    )
    
    rf_list[[idx]] = current_rf
    idx = idx + 1
}

### Plot variable importance ###
# Note that the permutation importance values contained in the
# ranger RF objects are calculated based on out-of-bag (OOB)
# error.

importance_mat <- matrix(nrow = length(predictor_names),
                         ncol = length(rf_list) + 2)
colnames(importance_mat) <- c(fold_names, "mean_importance", "sd")
rownames(importance_mat) <- predictor_names

for(idx in 1:length(rf_list)) {
  importance_mat[,idx] <- rf_list[[idx]]$variable.importance
}
### Determine mean and st. dev.
for(idx in 1:nrow(importance_mat)){
  importance_mat[idx, "mean_importance"] <-
    mean(importance_mat[idx, 1:length(rf_list)])
  importance_mat[idx, "sd"] <-
    sd(importance_mat[idx, 1:length(rf_list)])
}
print(importance_mat)

# Convert to data frame for ggplot
importance_df <- as.data.frame(importance_mat)

# Plot ranked mean importances
ggplot(data = importance_df, 
       aes(x=reorder(predictor_names, mean_importance), 
           y=mean_importance,
           fill = mean_importance)) + 
  geom_bar(stat="identity", position="dodge") + 
  geom_errorbar(data = importance_df,
                aes(ymin=mean_importance-sd,
                    ymax=mean_importance+sd),
                width=.2,
                position=position_dodge(.9)) +
  coord_flip() +
  ylab("Variable Importance (Permutation)") +
  xlab("") +
  guides(fill="none") +
  scale_fill_gradient(low="#d0efff", high="#1167b1")


#################################################################
### Feature selection ###
#################################################################

# You can ignore this part, but you need to run it :)
rfFuncsEdit <- rfFuncs
rfFuncsEdit$fit <- function(x, y, first, last, ...) {
  loadNamespace("randomForest")
  if(!is.factor(y)) y = as.factor(y)
  randomForest::randomForest(
    x, y,
    ntree = 200,
    importance = TRUE,
    type = "classification",
    ...)
}

### Recursive feature elimination (RFE), backwards selection ###

# RFE control options
rfe_ctrl <- caret::rfeControl(
  functions = rfFuncsEdit, # Use random forest classifier
  method = "cv", # Cross-validation
  number = 5, # 5 folds
  saveDetails = TRUE,
  verbose = TRUE,
  returnResamp = "all"
)

# Perform recursive feature eliminations (RFE)
rfe <- caret::rfe(
  x = as.matrix(training_data[predictor_names]),
  y = as.factor(as.matrix(training_data["presence"])),
  sizes = 1:length(predictor_names),
  metric = "Kappa",
  rfeControl = rfe_ctrl
)

### Plot results
plot(rfe, type=c("g", "o"))
print(rfe$results)
print(rfe$optVariables)

#################################################################
### Make predictions ###
#################################################################

### Train model on full data set with tuned hyperparameters ###
rf <- ranger::ranger(
  x = training_data[predictor_names],
  y = as.matrix(training_data["presence"]),
  classification = TRUE,
  probability = TRUE,
  num.trees = 500,
  mtry = best_mtry,
  max.depth = best_maxd,
  importance = "permutation"
)

### Spatial prediction
print(names(pred_stack))
print(rf$forest$independent.variable.names)
prediction_ras <- raster::predict(
  pred_stack,
  model = rf,
  type = "response",
  progress = "text",
  fun = function(model, ...) predict(model, ...)$predictions
  )

### Plot results, compare to true occurrences
plot(prediction_ras)
points(training_data$x[training_data["presence"]==1],
       training_data$y[training_data["presence"]==1],
       col = rgb(red = 0, green = 0, blue = 0, alpha=0.1),
       cex = 0.3, pch = 2)

# Zoom in to visualize occurrences, use Oslo extent
oslo_extent <- extent(220000, 280000, 6620000, 6680000)

plot(occ_ras, ext = oslo_extent)
plot(prediction_ras, ext = oslo_extent)

### How do you evaluate our model? Do the predictions make sense?
### Which variables are important predictors and why? Do you have 
### ideas for improving the workflow?


### Final exercise (if time allows):
#
# Option 1: Machine Learning and Distribution Modelling
#
# Write your own code (with inspiration from the code above) to test
# a different machine learning classifier on the same data. Read the
# documentation about the new classifier's hyper-parameters and use
# grid search tuning to optimize the model. You can basically copy-
# paste the entire script and adapt the parts where the ranger RF
# objects are involved. You do not need to re-run the data extraction,
# these objects will stay loaded in your R environment.
#
# Suggestions for classifiers (package):
#
# Support Vector Machine - SVM (e1071)
# Artificial Neural Networks (neuralnet)
# Adaptive Boosting / eXtreme gradient boosting (fastAdaboost, xgboost)
#
# Use one of the models in this list: 
# https://topepo.github.io/caret/available-models.html
#
# ATTENTION! Unlike Random Forest, most ML classifiers cannot handle 
# categorical features by themselves without prior transformations
# (e.g., one-hot encoding). Easiest fix: remove the categorical predictors
# from the feature matrix before training your model, and from the raster
# stack before making spatial predictions.



### Option 2: Practice R coding
# Write a function in R that takes an integer number (n) as input and 
# then prints the first n numbers of the Fibonacci series (0, 1, 1, 2, 
# 3, 5, 8, 13, ...).
# See: https://www.mathsisfun.com/numbers/fibonacci-sequence.html)
#
# Expected output:
# print_fibo(n = 5)
# > The first 5 elements of the Fibonacci series are: 0, 1, 1, 2, 3


### The End. ###
