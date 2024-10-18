# AppliedStast
#Read the Initial Data set and displaying the structure and set the dimension
InitialData <- read.csv(file = "gene-expression-invasive-vs-noninvasive-cancer (1).csv")
str(InitialData)
dim(InitialData)


#Creating the random set from largest register number from the group
set.seed(2316529)

#Creating the subset of 2000 columns in data set
indices <- rank(runif(1:4948))[1:2000]
New_Dataset <- InitialData[indices]
New_Dataset$Class <- InitialData$Class
dim(New_Dataset)

#Check null values and fill the null values with median
sum(is.na(New_Dataset))

#Calculate the sum of rows in this dataset and see how many na value or row or column is present
New_Dataset<-New_Dataset[-c(54),]
sum(is.na(New_Dataset))
which (is.na(New_Dataset),arr.ind = T)
New_Dataset[23,983]=median(New_Dataset[,983],na.rm = TRUE)
sum(is.na(New_Dataset))

#Assigning the New Dataset to Remove Null
Remove_Null_values <- New_Dataset
sum(is.na(Remove_Null_values))

#Check null values and fill the null values with knn imputation
# library(DMwR2)
# Remove_Null_values<- knnImputation(New_Dataset)
# sum(is.na(Remove_Null_values))


#Unsupervised using variance method- Dimensional Reduction
VARIANCE_DATASET <- apply(Remove_Null_values, 2, var)
FILTER_VARIANCE_DATASET <- names(VARIANCE_DATASET[VARIANCE_DATASET > 0.24])
FILTER_GENE_SUBSET <- Remove_Null_values[, FILTER_VARIANCE_DATASET]
dim(FILTER_GENE_SUBSET)


# Function for Two sample t-test SUPERVISED
gene_func <- function(gene) {
  results <- t.test(Remove_Null_values[[gene]] ~ Remove_Null_values$Class)
  return(results$p.value)
}

# Applying t-test function for the data set
p_values <- sapply(names(Remove_Null_values)[-ncol(Remove_Null_values)], gene_func)

# Variable to store genes names
col_names = names(p_values)

# Storing t-test results in a data frame
t_test_results_df <- data.frame(Genes = col_names, P_Value = p_values)

# Filtering genes for t-test
filtered_genes <- t_test_results_df[t_test_results_df$P_Value < 0.05, ]

# Transposing the filtered genes
transp_filtered_genes <- t(filtered_genes)

# Filtered column names
filtered_col_names <- colnames(transp_filtered_genes)

reduced_data_T_Test <- Remove_Null_values[filtered_col_names]

reduced_data_T_Test$Class <- Remove_Null_values$Class

dim(reduced_data_T_Test)

# Split the data into training and testing sets(Model Building)
library(caTools)
split_TTest_REDUCED_DATASET <- sample.split(reduced_data_T_Test$Class, SplitRatio = 0.7, set.seed (2316529))
TRAIN_TTest_REDUCED_DATASET  <- subset(reduced_data_T_Test, split_TTest_REDUCED_DATASET  == TRUE)
TEST_TTest_REDUCED_DATASET  <- subset(reduced_data_T_Test, split_TTest_REDUCED_DATASET  == FALSE)

#Spliting the train data in 0 and 1
table(TRAIN_TTest_REDUCED_DATASET$Class)
TRAIN_TTest_REDUCED_DATASET$Class[TRAIN_TTest_REDUCED_DATASET$Class==1]<- 0
TRAIN_TTest_REDUCED_DATASET$Class[TRAIN_TTest_REDUCED_DATASET$Class==2]<- 1
TEST_TTest_REDUCED_DATASET$Class[TEST_TTest_REDUCED_DATASET$Class==1]<- 0
TEST_TTest_REDUCED_DATASET$Class[TEST_TTest_REDUCED_DATASET$Class==2]<-1

# importing Libraries#
library(caret)
library(MASS)
library(randomForest)
library(class)
library(e1071)
# LOGISTICS REGRISSION MODEL T-test ##

# Fit a logistic regression model using training data
LOGISTICS_MODEL <- glm(Class ~ ., data = TRAIN_TTest_REDUCED_DATASET, family = 'binomial')

# Predict probabilities on test data
LOGISTICS_MODEL_PREDICT <- predict(LOGISTICS_MODEL, newdata = TEST_TTest_REDUCED_DATASET, type = "response")

# Convert predicted probabilities to binary predictions (0 or 1)
LOGISTICS_PREDICTED_CLASS <- ifelse(LOGISTICS_MODEL_PREDICT > 0.5, 1, 0)

# Create confusion matrix
LOGISTICS_confusion_Matrix_summary <- confusionMatrix(as.factor(LOGISTICS_PREDICTED_CLASS), as.factor(TEST_TTest_REDUCED_DATASET$Class))
LOGISTICS_confusion_Matrix_summary

# Calculate misclassification error
LOGISTICS_misclassification_error <- 1 - LOGISTICS_confusion_Matrix_summary$overall["Accuracy"]
LOGISTICS_misclassification_error

# LDA MODEL #

# Fit LDA model
LDA_MODEL <- lda(Class ~ ., data = TRAIN_TTest_REDUCED_DATASET)

# Predict on test data
LDA_MODEL_PREDICT <- predict(LDA_MODEL, newdata = TEST_TTest_REDUCED_DATASET)$class

# Create confusion matrix
LDA_confusion_Matrix_summary <- confusionMatrix(as.factor(LDA_MODEL_PREDICT), as.factor(TEST_TTest_REDUCED_DATASET$Class))
LDA_confusion_Matrix_summary

# Calculate misclassification error
misclassification_error_LDA <- 1 - LDA_confusion_Matrix_summary$overall["Accuracy"]
misclassification_error_LDA

# QDA MODEL NOT WORKING #

# Fit a QDA model using training data
QDA_MODEL <- qda(Class ~ ., data = TRAIN_TTest_REDUCED_DATASET)

# Predict using the QDA model
QDA_PREDICTION <- predict(QDA_MODEL, newdata = TEST_TTest_REDUCED_DATASET)$class
QDA_PREDICTION

# Create confusion matrix for table
confusion_Matrix_QDA <- table(QDA_PREDICTION, TEST_TTest_REDUCED_DATASET$Class)
confusion_Matrix_QDA

# Calculate misclassification error for QDA
misclassification_error_QDA <- mean(QDA_PREDICTION != TEST_TTest_REDUCED_DATASET$Class)
misclassification_error_QDA

# RANDOM FOREST MODEL #

# Train the Random Forest model
RF_MODEL <- randomForest(Class ~ ., data = TRAIN_TTest_REDUCED_DATASET,ntree  = 3,set.seed(2316529))

# Predict on test data
RF_PREDICTIONS <- predict(RF_MODEL, newdata = TEST_TTest_REDUCED_DATASET,type= "response")
RF_PREDICTIONS
RF_PREDICTED_CLASS<-ifelse(RF_PREDICTIONS > 0.5, 1, 0)

# Create confusion matrix
RF_confusion_Matrix_summary <- confusionMatrix(as.factor(RF_PREDICTED_CLASS), as.factor(TEST_TTest_REDUCED_DATASET$Class))
RF_confusion_Matrix_summary

# Or if you want to use table function
RF_Conf_matrix <- table(RF_PREDICTED_CLASS, TEST_TTest_REDUCED_DATASET$Class)
RF_Conf_matrix

# Calculate misclassification error
RF_misclassification_error <- 1 - RF_confusion_Matrix_summary$overall["Accuracy"]
RF_misclassification_error

## KNN MODEL ##

# Train the k-NN model
KNN_MODEL <- knn(train = TRAIN_TTest_REDUCED_DATASET[, -which(names(TRAIN_TTest_REDUCED_DATASET) == "Class")], 
                 test = TEST_TTest_REDUCED_DATASET[, -which(names(TEST_TTest_REDUCED_DATASET) == "Class")], 
                 cl = TRAIN_TTest_REDUCED_DATASET$Class, 
                 k = 3)  

# Predict on test data
KNN_PREDICTIONS <- KNN_MODEL

# Create confusion matrix
KNN_confusion_Matrix_summary <- confusionMatrix(as.factor(KNN_PREDICTIONS), as.factor(TEST_TTest_REDUCED_DATASET$Class))
KNN_confusion_Matrix_summary

# Create confusion matrix table
KNN_confusion_Matrix <- table(KNN_PREDICTIONS, TEST_TTest_REDUCED_DATASET$Class)
KNN_confusion_Matrix

# Calculate misclassification error using KNN means.
KNN_misclassification_error <- 1 - KNN_confusion_Matrix_summary$overall["Accuracy"]
KNN_misclassification_error


### SVM MODEL ###
# Train the SVM model
SVM_MODEL <- svm(Class ~ ., data = TRAIN_TTest_REDUCED_DATASET, kernel = "radial")

# Predict on test data
SVM_PREDICTIONS <- predict(SVM_MODEL, newdata = TEST_TTest_REDUCED_DATASET)
SVM_PREDICTED_CLASS<-ifelse(SVM_PREDICTIONS > 0.5, 1, 0)
# Create confusion matrix
SVM_confusion_Matrix_summary <- confusionMatrix(as.factor(SVM_PREDICTED_CLASS), as.factor(TEST_TTest_REDUCED_DATASET$Class))
SVM_confusion_Matrix_summary

# Create confusion matrix
SVM_confusion_Matrix <- table(SVM_PREDICTED_CLASS, TEST_TTest_REDUCED_DATASET$Class)
SVM_confusion_Matrix

# Calculate misclassification error
SVM_misclassification_error <- 1 - SVM_confusion_Matrix_summary$overall["Accuracy"]
SVM_misclassification_error

table(New_Dataset$Class)
New_Dataset$Class[New_Dataset$Class==1]<- 0
New_Dataset$Class[New_Dataset$Class==2]<- 1

# Re sampling Technique of k-fold cross-validation LOGISTICS REGRISSION MODEL for 2000 DATASET # 
# Set the number of folds for cross-validation #
k <- 5

# Create an empty vector to store classification errors
LOGISTICS_misclassification_error_KFOLD <- numeric(k)

# Perform k-fold cross-validation
for (i in 1:(k-1)) {  
  # Split the data into training and testing sets for the current fold
  set.seed (2316529)
  indices <- sample(1:nrow(New_Dataset), replace = FALSE)
  split <- cut(indices, breaks = k, labels = FALSE)
  
  # Fit a logistic regression model using training data
  LOGISTICS_MODEL_KFOLD <- glm(Class ~ ., data = New_Dataset[split != i, ], family = 'binomial')
  
  # Predict probabilities on the testing set
  LOGISTICS_MODEL_PREDICT_KFOLD <- predict(LOGISTICS_MODEL_KFOLD, newdata = New_Dataset[split == i, ], type = "response")
  
  # Convert predicted probabilities to binary predictions (0 or 1)
  LOGISTICS_PREDICTED_CLASS_KFOLD <- ifelse(LOGISTICS_MODEL_PREDICT_KFOLD > 0.5, 1, 0)
  
  # Create confusion matrix
  LOGISTICS_confusion_Matrix_summary_KFOLD <- confusionMatrix(as.factor(LOGISTICS_PREDICTED_CLASS_KFOLD), as.factor(New_Dataset$Class[split == i]))
 
  # Calculate misclassification error for the current fold
  misclassification_error_KFOLD <- 1 - LOGISTICS_confusion_Matrix_summary_KFOLD$overall["Accuracy"]
  LOGISTICS_misclassification_error_KFOLD[i] <- misclassification_error_KFOLD
}

# Calculate the average misclassification error across all folds #
AVERAGE_misclassification_error_logit <- mean(LOGISTICS_misclassification_error_KFOLD)
AVERAGE_misclassification_error_logit

# Re sampling Technique of k-fold cross-validation LDA MODEL for 2000 dataset #
library(MASS)
# Set the number of folds for cross-validation
k <- 5

# Create an empty vector to store classification errors
LDA_misclassification_errors <- numeric(k)

# Perform k-fold cross-validation
for (i in 1:(k-1)) {  
  # Split the data into training and testing sets for the current fold
  set.seed (2316529)
  indices <- sample(1:nrow(New_Dataset), replace = FALSE)
  split <- cut(indices, breaks = k, labels = FALSE)
  
  # Train the LDA model on the training set
  LDA_MODEL_KFOLD <- lda(Class ~ ., data = New_Dataset[split != i, ])
  
  # Predict on the testing set
  LDA_MODEL_PREDICT <- predict(LDA_MODEL_KFOLD, newdata = New_Dataset[split == i, ])$class
  
  # Create confusion matrix
  LDA_confusion_Matrix_summary <- confusionMatrix(as.factor(LDA_MODEL_PREDICT), as.factor(New_Dataset$Class[split == i]))
  
  # Calculate misclassification error for the current fold
  LDA_misclassification_error <- 1 - LDA_confusion_Matrix_summary$overall["Accuracy"]
  LDA_misclassification_errors[i] <- LDA_misclassification_error
}

# Calculate the average misclassification error across all folds
average_misclassification_error_LDA <- mean(LDA_misclassification_errors)
average_misclassification_error_LDA


# Re sampling Technique of k-fold cross-validation Random Forest MODEL for 2000 dataset # 

library(randomForest)

# Set the number of folds for cross-validation
k <- 5
min_errors_rf <- 1
error_rf=list()

# Create an empty vector to store misclassification errors for each fold
RF_misclassification_errors <- numeric(k)

# Perform k-fold cross-validation
for (j in 1:50){
  for (i in 1 :(k-1)){ 
    # Split the data into training and testing sets for the current fold
    set.seed (2316529)
    indices <- sample(1:nrow(New_Dataset), replace = FALSE)
    
    split <- cut(indices, breaks = k, labels = FALSE)
    
    # Train the Random Forest model on the training set
    RF_MODEL_KFOLD <- randomForest(Class ~ ., data = New_Dataset[split != i, ], ntree = j, set.seed (2316529))
    
    # Predict on the testing set
    RF_PREDICTIONS_KFOLD <- predict(RF_MODEL_KFOLD, newdata = New_Dataset[split == i, ], type = "response")
    RF_PREDICTED_CLASS_KFOLD <- ifelse(RF_PREDICTIONS_KFOLD > 0.5, 1, 0)
    
    # Create confusion matrix
    RF_confusion_Matrix_summary <- confusionMatrix(as.factor(RF_PREDICTED_CLASS_KFOLD), as.factor(New_Dataset$Class[split == i]))
    
    # Calculate misclassification error for the current fold
    RF_misclassification_error <- 1 - RF_confusion_Matrix_summary$overall["Accuracy"]
    
    #Append error to list
    error_rf<- append(error_rf,RF_misclassification_error)
  }
  
  error_rf
  aplha1<-unlist(error_rf)
  mean_aplha1<-mean(aplha1)
  mean_aplha1
  if (mean_aplha1<min_errors_rf){
    min_errors_rf<-mean_aplha1
    V<-j
  }
}
min_errors_rf
V

# Re sampling Technique of k-fold cross-validation for KNN MODEL for 2000 dataset #

# Set the number of folds for cross-validation
k <- 5
min_errorsknn <- 1
error_knn=list()

# Create an empty vector to store classification errors
KNN_misclassification_errors <- numeric(k)
for (j in 1:50){
  for (i in 1 :(k-1)){
    
    # Split the data into training and testing sets for the current fold
    set.seed(2316529)
    indices <- sample(1:nrow(New_Dataset), replace = FALSE)
    
    split <- cut(indices, breaks = k, labels = FALSE)
    
    # Train the k-NN model on the training set
    
    KNN_MODEL_KFOLD <- knn(train = New_Dataset[split != i, -which(names(New_Dataset) == "Class")], 
                           test = New_Dataset[split == i, -which(names(New_Dataset) == "Class")], 
                           cl = New_Dataset$Class[split != i],k=j)  
    
    # Predict on the testing set
    KNN_PREDICTIONS_KFOLD <- KNN_MODEL_KFOLD
    
    # Create confusion matrix
    KNN_confusion_Matrix_summary <- confusionMatrix(as.factor(KNN_PREDICTIONS_KFOLD), as.factor(New_Dataset$Class[split == i]))
    KNN_confusion_Matrix_summary
    
    # Calculate misclassification error for the current fold
    KNN_misclassification_error <- 1 - KNN_confusion_Matrix_summary$overall["Accuracy"]
    
    #Append error to list
    error_knn<- append(error_knn,KNN_misclassification_error)
  }
  
  error_knn
  aplha<-unlist(error_knn)
  mean_aplha<-mean(aplha)
  mean_aplha
  if (mean_aplha<min_errorsknn){
    min_errorsknn<-mean_aplha
    K<-j
  }
}
min_errorsknn
K

# Re sampling Technique of k-fold cross-validation for SVM MODEL for 2000 dataset #

# Set the number of folds for cross-validation
k <- 5

# Create an empty vector to store classifications errors
SVM_misclassification_errors <- numeric(k)

# Perform k-fold cross-validation
for (i in 1:(k-1)) {  # Iterate from 1 to 8 instead of 1 to 9
  # Split the data into training and testing sets for the current fold
  set.seed(2316529)
  indices <- sample(1:nrow(New_Dataset), replace = FALSE)
  split <- cut(indices, breaks = k, labels = FALSE)
  
  # Train the SVM model on the training set
  SVM_MODEL_KFOLD <- svm(Class ~ ., data = New_Dataset[split != i, ], kernel = "radial")
  
  # Predict on the testing set
  SVM_PREDICTIONS_KFOLD <- predict(SVM_MODEL_KFOLD, newdata = New_Dataset[split == i, ])
  SVM_PREDICTED_CLASS_KFOLD<-ifelse(SVM_PREDICTIONS_KFOLD > 0.5, 1, 0)
  SVM_confusion_Matrix_summary <- confusionMatrix(as.factor(SVM_PREDICTED_CLASS_KFOLD), as.factor(New_Dataset$Class[split == i]))
  
  # Calculate misclassification error for the current fold
  SVM_misclassification_error <- 1 - SVM_confusion_Matrix_summary$overall["Accuracy"]
  SVM_misclassification_errors[i]<-SVM_misclassification_error
}

# Calculate the average misclassification error across all folds
SVM_average_misclassification_error <- unlist(SVM_misclassification_errors)
SVM_average_misclassification_error
mean(SVM_average_misclassification_error)

# Split the dataset into 0 and 1
table(reduced_data_T_Test$Class)
reduced_data_T_Test$Class[reduced_data_T_Test$Class==1]<- 0
reduced_data_T_Test$Class[reduced_data_T_Test$Class==2]<- 1

# Re sampling Technique of k-fold cross-validation LOGISTICS REGRISSION MODEL # 
#reduced_data_T_Test
# Set the number of folds for cross-validation
k <- 5

# Create an empty vector to store classification errors
LOGISTICS_misclassification_error_KFOLD <- numeric(k)

# Perform k-fold cross-validation
for (i in 1:(k-1)) {  
  # Split the data into training and testing sets for the current fold
  set.seed (2316529)
  indices <- sample(1:nrow(reduced_data_T_Test), replace = FALSE)
  split <- cut(indices, breaks = k, labels = FALSE)
  
  # Fit a logistic regression model using training data
  LOGISTICS_MODEL_KFOLD <- glm(Class ~ ., data = reduced_data_T_Test[split != i, ], family = 'binomial')
  
  # Predict probabilities on the testing set
  LOGISTICS_MODEL_PREDICT_KFOLD <- predict(LOGISTICS_MODEL_KFOLD, newdata = reduced_data_T_Test[split == i, ], type = "response")
  
  # Convert predicted probabilities to binary predictions (0 or 1)
  LOGISTICS_PREDICTED_CLASS_KFOLD <- ifelse(LOGISTICS_MODEL_PREDICT_KFOLD > 0.5, 1, 0)
  
  # Create confusion matrix
  LOGISTICS_confusion_Matrix_summary_KFOLD <- confusionMatrix(as.factor(LOGISTICS_PREDICTED_CLASS_KFOLD), as.factor(reduced_data_T_Test$Class[split == i]))
  
  # Calculate misclassification error for the current fold
  misclassification_error_KFOLD <- 1 - LOGISTICS_confusion_Matrix_summary_KFOLD$overall["Accuracy"]
  LOGISTICS_misclassification_error_KFOLD[i] <- misclassification_error_KFOLD
}

# Calculate the average misclassification error across all folds
AVERAGE_misclassification_error_logit <- mean(LOGISTICS_misclassification_error_KFOLD)
AVERAGE_misclassification_error_logit

# Re sampling Technique of k-fold cross-validation LDA MODEL #
# Set the number of folds for cross-validation

k <- 5

# Create an empty vector to store classification errors
LDA_misclassification_errors <- numeric(k)

# Perform k-fold cross-validation
for (i in 1:(k-1)) {  
  # Split the data into training and testing sets for the current fold
  set.seed (2316529)
  indices <- sample(1:nrow(reduced_data_T_Test), replace = FALSE)
  split <- cut(indices, breaks = k, labels = FALSE)
  
  # Train the LDA model on the training set
  LDA_MODEL_KFOLD <- lda(Class ~ ., data = reduced_data_T_Test[split != i, ])
  
  # Predict on the testing set
  LDA_MODEL_PREDICT <- predict(LDA_MODEL_KFOLD, newdata = reduced_data_T_Test[split == i, ])$class
  
  # Create confusion matrix
  LDA_confusion_Matrix_summary <- confusionMatrix(as.factor(LDA_MODEL_PREDICT), as.factor(reduced_data_T_Test$Class[split == i]))
  
  # Calculate misclassification error for the current fold
  LDA_misclassification_error <- 1 - LDA_confusion_Matrix_summary$overall["Accuracy"]
  LDA_misclassification_errors[i] <- LDA_misclassification_error
}

# Calculate the average misclassification error across all folds
average_misclassification_error_LDA <- mean(LDA_misclassification_errors)
average_misclassification_error_LDA

# Re sampling Technique of k-fold cross-validation Random Forest MODEL #  

library(randomForest)

# Set the number of folds for cross-validation
k <- 5
min_errors_rf <- 1
error_rf=list()

# Create an empty vector to store misclassification errors for each fold
RF_misclassification_errors <- numeric(k)

# Perform k-fold cross-validation
for (j in 1:50){
  for (i in 1 :(k-1)){ 
    # Split the data into training and testing sets for the current fold
    set.seed (2316529)
    indices <- sample(1:nrow(reduced_data_T_Test), replace = FALSE)
    
    split <- cut(indices, breaks = k, labels = FALSE)
    
    # Train the Random Forest model on the training set
    RF_MODEL_KFOLD <- randomForest(Class ~ ., data = reduced_data_T_Test[split != i, ], ntree = j, set.seed (2316529))
    
    # Predict on the testing set
    RF_PREDICTIONS_KFOLD <- predict(RF_MODEL_KFOLD, newdata = reduced_data_T_Test[split == i, ], type = "response")
    RF_PREDICTED_CLASS_KFOLD <- ifelse(RF_PREDICTIONS_KFOLD > 0.5, 1, 0)
    
    # Create confusion matrix
    RF_confusion_Matrix_summary <- confusionMatrix(as.factor(RF_PREDICTED_CLASS_KFOLD), as.factor(reduced_data_T_Test$Class[split == i]))
    
    # Calculate misclassification error for the current fold
    RF_misclassification_error <- 1 - RF_confusion_Matrix_summary$overall["Accuracy"]
    
    #Append error to list
    error_rf<- append(error_rf,RF_misclassification_error)
  }
  
  error_rf
  aplha1<-unlist(error_rf)
  mean_aplha1<-mean(aplha1)
  mean_aplha1
  if (mean_aplha1<min_errors_rf){
    min_errors_rf<-mean_aplha1
    V<-j
  }
}
min_errors_rf
V

# Re sampling Technique of k-fold cross-validation for KNN MODEL #

# Set the number of folds for cross-validation
k <- 5
min_errorsknn <- 1
error_knn=list()

# Create an empty vector to store classification errors
KNN_misclassification_errors <- numeric(k)
for (j in 1:50){
  for (i in 1 :(k-1)){
    
    # Split the data into training and testing sets for the current fold
    set.seed(2316529)
    indices <- sample(1:nrow(reduced_data_T_Test), replace = FALSE)
    
    split <- cut(indices, breaks = k, labels = FALSE)
    
    # Train the k-NN model on the training set
    
    KNN_MODEL_KFOLD <- knn(train = reduced_data_T_Test[split != i, -which(names(reduced_data_T_Test) == "Class")], 
                           test = reduced_data_T_Test[split == i, -which(names(reduced_data_T_Test) == "Class")], 
                           cl = reduced_data_T_Test$Class[split != i],k=j)  
    
    # Predict on the testing set
    KNN_PREDICTIONS_KFOLD <- KNN_MODEL_KFOLD
    
    # Create confusion matrix
    KNN_confusion_Matrix_summary <- confusionMatrix(as.factor(KNN_PREDICTIONS_KFOLD), as.factor(reduced_data_T_Test$Class[split == i]))
    KNN_confusion_Matrix_summary
    
    # Calculate misclassification error for the current fold
    KNN_misclassification_error <- 1 - KNN_confusion_Matrix_summary$overall["Accuracy"]
    
    #Append error to list
    error_knn<- append(error_knn,KNN_misclassification_error)
  }
  
  error_knn
  aplha<-unlist(error_knn)
  mean_aplha<-mean(aplha)
  mean_aplha
  if (mean_aplha<min_errorsknn){
    min_errorsknn<-mean_aplha
    K<-j
  }
}
min_errorsknn
k

# Re sampling Technique of k-fold cross-validation for SVM MODEL #

# Set the number of folds for cross-validation
k <- 5

# Create an empty vector to store classifications errors
SVM_misclassification_errors <- numeric(k)

# Perform k-fold cross-validation
for (i in 1:(k-1)) {  # Iterate from 1 to 8 instead of 1 to 9
  # Split the data into training and testing sets for the current fold
  set.seed(2316529)
  indices <- sample(1:nrow(reduced_data_T_Test), replace = FALSE)
  split <- cut(indices, breaks = k, labels = FALSE)
  
  # Train the SVM model on the training set
  SVM_MODEL_KFOLD <- svm(Class ~ ., data = reduced_data_T_Test[split != i, ], kernel = "radial")
  
  # Predict on the testing set
  SVM_PREDICTIONS_KFOLD <- predict(SVM_MODEL_KFOLD, newdata = reduced_data_T_Test[split == i, ])
  SVM_PREDICTED_CLASS_KFOLD<-ifelse(SVM_PREDICTIONS_KFOLD > 0.5, 1, 0)
  SVM_confusion_Matrix_summary <- confusionMatrix(as.factor(SVM_PREDICTED_CLASS_KFOLD), as.factor(reduced_data_T_Test$Class[split == i]))
  
  # Calculate misclassification error for the current fold
  SVM_misclassification_error <- 1 - SVM_confusion_Matrix_summary$overall["Accuracy"]
  SVM_misclassification_errors[i]<-SVM_misclassification_error
}

# Calculate the average misclassification error across all folds
SVM_average_misclassification_error <- unlist(SVM_misclassification_errors)
SVM_average_misclassification_error
mean(SVM_average_misclassification_error)

#Split the data into training and testing sets using 2000 Data set(ALL 2000 DATSET)
library(caTools)
SPLIT_DATASET <- sample.split(New_Dataset$Class, SplitRatio = 0.7, set.seed (2316529))
TRAIN_DATASET_2000 <- subset(New_Dataset, SPLIT_DATASET == TRUE)
TEST_DATATSET_2000 <- subset(New_Dataset, SPLIT_DATASET == FALSE)

# FOR ALL 2000 SUBSET DATA LOGISTICS REGRISSION MODEL # 

# Fit a logistic regression model using training data
LOGISTICS_MODEL_2000 <- glm(Class ~ ., data = TRAIN_DATASET_2000, family = 'binomial')

# Predict probabilities on test data
LOGISTICS_MODEL_PREDICT_2000 <- predict(LOGISTICS_MODEL_2000, newdata = TEST_DATATSET_2000, type = "response")

# Convert predicted probabilities to binary predictions (0 or 1)
LOGISTICS_PREDICTED_CLASS_2000 <- ifelse(LOGISTICS_MODEL_PREDICT_2000 > 0.5, 1, 0)

# Create confusion matrix
LOGISTICS_confusion_Matrix_summary_2000 <- confusionMatrix(as.factor(LOGISTICS_PREDICTED_CLASS_2000), as.factor(TEST_DATATSET_2000$Class))
LOGISTICS_confusion_Matrix_summary_2000

# Or if you want to use table function
LOGISTICS_confusion_Matrix_2000 <- table(LOGISTICS_PREDICTED_CLASS_2000, TEST_DATATSET_2000$Class)
LOGISTICS_confusion_Matrix_2000

# Calculate misclassification error
LOGISTICS_misclassification_error_2000 <- 1 - LOGISTICS_confusion_Matrix_summary_2000$overall["Accuracy"]
LOGISTICS_misclassification_error_2000

# LDA MODEL for 2000 subset #
# Fit LDA model
LDA_MODEL_2000 <- lda(Class ~ ., data = TRAIN_DATASET_2000)

# Predict on test data
LDA_MODEL_PREDICT_2000 <- predict(LDA_MODEL_2000, newdata = TEST_DATATSET_2000)$class

# Create confusion matrix
LDA_confusion_Matrix_summary_2000 <- confusionMatrix(as.factor(LDA_MODEL_PREDICT_2000), as.factor(TEST_DATATSET_2000$Class))
LDA_confusion_Matrix_summary_2000

# Calculate misclassification error
misclassification_error_LDA <- 1 - LDA_confusion_Matrix_summary_2000$overall["Accuracy"]
misclassification_error_LDA

## RAMDOM FOREST MODEL FOR ALL 2000 SUBSET #

library(randomForest)


# Train the Random Forest model
RF_MODEL_2000 <- randomForest(Class ~ ., data = TRAIN_DATASET_2000,ntree  = 3,set.seed(2316529))

# Predict on test data
RF_PREDICTIONS_2000 <- predict(RF_MODEL_2000, newdata = TEST_DATATSET_2000,type= "response")
RF_PREDICTIONS_2000
RF_PREDICTED_CLASS_2000<-ifelse(RF_PREDICTIONS_2000 > 0.5, 1, 0)

# Create confusion matrix
RF_confusion_Matrix_summary_2000 <- confusionMatrix(as.factor(RF_PREDICTED_CLASS_2000), as.factor(TEST_DATATSET_2000$Class))
RF_confusion_Matrix_summary_2000

# Or if you want to use table function
RF_Conf_matrix_2000 <- table(RF_PREDICTED_CLASS_2000, TEST_DATATSET_2000$Class)
RF_Conf_matrix_2000

# Calculate misclassification error
RF_misclassification_error_2000 <- 1 - RF_confusion_Matrix_summary_2000$overall["Accuracy"]
RF_misclassification_error_2000

# KNN MODEL FOR 2000 SUBSET OF DATA #

library(class)

KNN_MODEL_2000 <- knn(train = TRAIN_DATASET_2000[, -which(names(TRAIN_DATASET_2000) == "Class")], 
                      test = TEST_DATATSET_2000[, -which(names(TEST_DATATSET_2000) == "Class")], 
                      cl = TRAIN_DATASET_2000$Class, 
                      k = 3)  

# Predict on test data
KNN_PREDICTIONS_2000 <- KNN_MODEL_2000

# Create confusion matrix
KNN_confusion_Matrix_summary_2000 <- confusionMatrix(as.factor(KNN_PREDICTIONS_2000), as.factor(TEST_DATATSET_2000$Class))
KNN_confusion_Matrix_summary_2000

# Create confusion matrix for table
KNN_confusion_Matrix_2000 <- table(KNN_PREDICTIONS_2000, TEST_DATATSET_2000$Class)
KNN_confusion_Matrix_2000

# Calculate misclassification error
KNN_misclassification_error_2000 <- 1 - KNN_confusion_Matrix_summary_2000$overall["Accuracy"]
KNN_misclassification_error_2000

# SVM MODEL FOR 2000 SUBSET OF DATA #

# Load the required package
library(e1071)

# Train the SVM model
SVM_MODEL_2000 <- svm(Class ~ ., data = TRAIN_DATASET_2000, kernel = "radial")

# Predict on test data
SVM_PREDICTIONS_2000 <- predict(SVM_MODEL_2000, newdata = TEST_DATATSET_2000)
SVM_PREDICTED_CLASS_2000<-ifelse(SVM_PREDICTIONS_2000 > 0.5, 1, 0)
# Create confusion matrix
SVM_confusion_Matrix_summary_2000 <- confusionMatrix(as.factor(SVM_PREDICTED_CLASS_2000), as.factor(TEST_DATATSET_2000$Class))
SVM_confusion_Matrix_summary_2000

# Create confusion matrix
SVM_confusion_Matrix_2000 <- table(SVM_PREDICTED_CLASS_2000, TEST_DATATSET_2000$Class)
SVM_confusion_Matrix_2000

# Calculate misclassification error
SVM_misclassification_error_2000 <- 1 - SVM_confusion_Matrix_summary_2000$overall["Accuracy"]
SVM_misclassification_error_2000

# LOGISTICS REGRISSION MODEL T-test ##

# Fit a logistic regression model using training data
LOGISTICS_MODEL <- glm(Class ~ ., data = TRAIN_TTest_REDUCED_DATASET, family = 'binomial')

# Predict probabilities on test data
LOGISTICS_MODEL_PREDICT <- predict(LOGISTICS_MODEL, newdata = TEST_TTest_REDUCED_DATASET, type = "response")

# Convert predicted probabilities to binary predictions (0 or 1)
LOGISTICS_PREDICTED_CLASS <- ifelse(LOGISTICS_MODEL_PREDICT > 0.5, 1, 0)

# Create confusion matrix
LOGISTICS_confusion_Matrix_summary <- confusionMatrix(as.factor(LOGISTICS_PREDICTED_CLASS), as.factor(TEST_TTest_REDUCED_DATASET$Class))
LOGISTICS_confusion_Matrix_summary

# Calculate misclassification error
LOGISTICS_misclassification_error <- 1 - LOGISTICS_confusion_Matrix_summary$overall["Accuracy"]
LOGISTICS_misclassification_error

# LDA MODEL #

# Fit LDA model
LDA_MODEL <- lda(Class ~ ., data = TRAIN_TTest_REDUCED_DATASET)

# Predict on test data
LDA_MODEL_PREDICT <- predict(LDA_MODEL, newdata = TEST_TTest_REDUCED_DATASET)$class

# Create confusion matrix
LDA_confusion_Matrix_summary <- confusionMatrix(as.factor(LDA_MODEL_PREDICT), as.factor(TEST_TTest_REDUCED_DATASET$Class))
LDA_confusion_Matrix_summary

# Calculate misclassification error
misclassification_error_LDA <- 1 - LDA_confusion_Matrix_summary$overall["Accuracy"]
misclassification_error_LDA

# QDA MODEL NOT WORKING #

# Fit a QDA model using training data
QDA_MODEL <- qda(Class ~ ., data = TRAIN_TTest_REDUCED_DATASET)

# Predict using the QDA model
QDA_PREDICTION <- predict(QDA_MODEL, newdata = TEST_TTest_REDUCED_DATASET)$class
QDA_PREDICTION

# Create confusion matrix
confusion_Matrix_QDA <- table(QDA_PREDICTION, TEST_TTest_REDUCED_DATASET$Class)
confusion_Matrix_QDA

# Calculate misclassification error for QDA
misclassification_error_QDA <- mean(QDA_PREDICTION != TEST_TTest_REDUCED_DATASET$Class)
misclassification_error_QDA

# RANDOM FOREST MODEL #

# Train the Random Forest model
RF_MODEL <- randomForest(Class ~ ., data = TRAIN_TTest_REDUCED_DATASET,ntree  = 3,set.seed(2316529))

# Predict on test data
RF_PREDICTIONS <- predict(RF_MODEL, newdata = TEST_TTest_REDUCED_DATASET,type= "response")
RF_PREDICTIONS
RF_PREDICTED_CLASS<-ifelse(RF_PREDICTIONS > 0.5, 1, 0)

# Create confusion matrix
RF_confusion_Matrix_summary <- confusionMatrix(as.factor(RF_PREDICTED_CLASS), as.factor(TEST_TTest_REDUCED_DATASET$Class))
RF_confusion_Matrix_summary

# Or if you want to use table function
RF_Conf_matrix <- table(RF_PREDICTED_CLASS, TEST_TTest_REDUCED_DATASET$Class)
RF_Conf_matrix

# Calculate misclassification error
RF_misclassification_error <- 1 - RF_confusion_Matrix_summary$overall["Accuracy"]
RF_misclassification_error

## KNN MODEL ##
# Train the k-NN model
KNN_MODEL <- knn(train = TRAIN_TTest_REDUCED_DATASET[, -which(names(TRAIN_TTest_REDUCED_DATASET) == "Class")], 
                 test = TEST_TTest_REDUCED_DATASET[, -which(names(TEST_TTest_REDUCED_DATASET) == "Class")], 
                 cl = TRAIN_TTest_REDUCED_DATASET$Class, 
                 k = 3)  

# Predict on test data
KNN_PREDICTIONS <- KNN_MODEL

# Create confusion matrix
KNN_confusion_Matrix_summary <- confusionMatrix(as.factor(KNN_PREDICTIONS), as.factor(TEST_TTest_REDUCED_DATASET$Class))
KNN_confusion_Matrix_summary

# Create confusion matrix
KNN_confusion_Matrix <- table(KNN_PREDICTIONS, TEST_TTest_REDUCED_DATASET$Class)
KNN_confusion_Matrix

# Calculate misclassification error
KNN_misclassification_error <- 1 - KNN_confusion_Matrix_summary$overall["Accuracy"]
KNN_misclassification_error


## SVM MODEL ##
# Train the SVM model
SVM_MODEL <- svm(Class ~ ., data = TRAIN_TTest_REDUCED_DATASET, kernel = "radial")

# Predict on test data
SVM_PREDICTIONS <- predict(SVM_MODEL, newdata = TEST_TTest_REDUCED_DATASET)
SVM_PREDICTED_CLASS<-ifelse(SVM_PREDICTIONS > 0.5, 1, 0)
# Create confusion matrix
SVM_confusion_Matrix_summary <- confusionMatrix(as.factor(SVM_PREDICTED_CLASS), as.factor(TEST_TTest_REDUCED_DATASET$Class))
SVM_confusion_Matrix_summary

# Create confusion matrix
SVM_confusion_Matrix <- table(SVM_PREDICTED_CLASS, TEST_TTest_REDUCED_DATASET$Class)
SVM_confusion_Matrix

# Calculate misclassification error
SVM_misclassification_error <- 1 - SVM_confusion_Matrix_summary$overall["Accuracy"]
SVM_misclassification_error

table(New_Dataset$Class)
New_Dataset$Class[New_Dataset$Class==1]<- 0
New_Dataset$Class[New_Dataset$Class==2]<- 1

## Re sampling Technique of k-fold cross-validation LOGISTICS REGRISSION MODEL for 2000 DATASET ### 
# Set the number of folds for cross-validation
k <- 5

# Create an empty vector to store classification errors
LOGISTICS_misclassification_error_KFOLD <- numeric(k)

# Perform k-fold cross-validation
for (i in 1:(k-1)) {  
  # Split the data into training and testing sets for the current fold
  set.seed (2316529)
  indices <- sample(1:nrow(New_Dataset), replace = FALSE)
  split <- cut(indices, breaks = k, labels = FALSE)
  
  # Fit a logistic regression model using training data
  LOGISTICS_MODEL_KFOLD <- glm(Class ~ ., data = New_Dataset[split != i, ], family = 'binomial')
  
  # Predict probabilities on the testing set
  LOGISTICS_MODEL_PREDICT_KFOLD <- predict(LOGISTICS_MODEL_KFOLD, newdata = New_Dataset[split == i, ], type = "response")
  
  # Convert predicted probabilities to binary predictions (0 or 1)
  LOGISTICS_PREDICTED_CLASS_KFOLD <- ifelse(LOGISTICS_MODEL_PREDICT_KFOLD > 0.5, 1, 0)
  
  # Create confusion matrix
  LOGISTICS_confusion_Matrix_summary_KFOLD <- confusionMatrix(as.factor(LOGISTICS_PREDICTED_CLASS_KFOLD), as.factor(New_Dataset$Class[split == i]))
 
  # Calculate misclassification error for the current fold
  misclassification_error_KFOLD <- 1 - LOGISTICS_confusion_Matrix_summary_KFOLD$overall["Accuracy"]
  LOGISTICS_misclassification_error_KFOLD[i] <- misclassification_error_KFOLD
}

# Calculate the average misclassification error across all folds
AVERAGE_misclassification_error_logit <- mean(LOGISTICS_misclassification_error_KFOLD)
AVERAGE_misclassification_error_logit

## Re sampling Technique of k-fold cross-validation LDA MODEL for 2000 dataset ##

library(MASS)
# Set the number of folds for cross-validation
k <- 5

# Create an empty vector to store classification errors
LDA_misclassification_errors <- numeric(k)

# Perform k-fold cross-validation
for (i in 1:(k-1)) {  
  # Split the data into training and testing sets for the current fold
  set.seed (2316529)
  indices <- sample(1:nrow(New_Dataset), replace = FALSE)
  split <- cut(indices, breaks = k, labels = FALSE)
  
  # Train the LDA model on the training set
  LDA_MODEL_KFOLD <- lda(Class ~ ., data = New_Dataset[split != i, ])
  
  # Predict on the testing set
  LDA_MODEL_PREDICT <- predict(LDA_MODEL_KFOLD, newdata = New_Dataset[split == i, ])$class
  
  # Create confusion matrix
  LDA_confusion_Matrix_summary <- confusionMatrix(as.factor(LDA_MODEL_PREDICT), as.factor(New_Dataset$Class[split == i]))
  
  # Calculate misclassification error for the current fold
  LDA_misclassification_error <- 1 - LDA_confusion_Matrix_summary$overall["Accuracy"]
  LDA_misclassification_errors[i] <- LDA_misclassification_error
}

# Calculate the average misclassification error across all folds
average_misclassification_error_LDA <- mean(LDA_misclassification_errors)
average_misclassification_error_LDA


# Re sampling Technique of k-fold cross-validation Random Forest MODEL for 2000 dataset ##  

library(randomForest)

# Set the number of folds for cross-validation
k <- 5
min_errors_rf <- 1
error_rf=list()

# Create an empty vector to store misclassification errors for each fold
RF_misclassification_errors <- numeric(k)

# Perform k-fold cross-validation
for (j in 1:50){
  for (i in 1 :(k-1)){ 
    # Split the data into training and testing sets for the current fold
    set.seed (2316529)
    indices <- sample(1:nrow(New_Dataset), replace = FALSE)
    
    split <- cut(indices, breaks = k, labels = FALSE)
    
    # Train the Random Forest model on the training set
    RF_MODEL_KFOLD <- randomForest(Class ~ ., data = New_Dataset[split != i, ], ntree = j, set.seed (2316529))
    
    # Predict on the testing set
    RF_PREDICTIONS_KFOLD <- predict(RF_MODEL_KFOLD, newdata = New_Dataset[split == i, ], type = "response")
    RF_PREDICTED_CLASS_KFOLD <- ifelse(RF_PREDICTIONS_KFOLD > 0.5, 1, 0)
    
    # Create confusion matrix
    RF_confusion_Matrix_summary <- confusionMatrix(as.factor(RF_PREDICTED_CLASS_KFOLD), as.factor(New_Dataset$Class[split == i]))
    
    # Calculate misclassification error for the current fold
    RF_misclassification_error <- 1 - RF_confusion_Matrix_summary$overall["Accuracy"]
    
    #Append error to list
    error_rf<- append(error_rf,RF_misclassification_error)
  }
  
  error_rf
  aplha1<-unlist(error_rf)
  mean_aplha1<-mean(aplha1)
  mean_aplha1
  if (mean_aplha1<min_errors_rf){
    min_errors_rf<-mean_aplha1
    V<-j
  }
}
min_errors_rf
V

# Re sampling Technique of k-fold cross-validation for KNN MODEL for 2000 dataset #

# Set the number of folds for cross-validation
k <- 5
min_errorsknn <- 1
error_knn=list()

# Create an empty vector to store classification errors
KNN_misclassification_errors <- numeric(k)
for (j in 1:50){
  for (i in 1 :(k-1)){
    
    # Split the data into training and testing sets for the current fold
    set.seed(2316529)
    indices <- sample(1:nrow(New_Dataset), replace = FALSE)
    
    split <- cut(indices, breaks = k, labels = FALSE)
    
    # Train the k-NN model on the training set
    
    KNN_MODEL_KFOLD <- knn(train = New_Dataset[split != i, -which(names(New_Dataset) == "Class")], 
                           test = New_Dataset[split == i, -which(names(New_Dataset) == "Class")], 
                           cl = New_Dataset$Class[split != i],k=j)  
    
    # Predict on the testing set
    KNN_PREDICTIONS_KFOLD <- KNN_MODEL_KFOLD
    
    # Create confusion matrix
    KNN_confusion_Matrix_summary <- confusionMatrix(as.factor(KNN_PREDICTIONS_KFOLD), as.factor(New_Dataset$Class[split == i]))
    KNN_confusion_Matrix_summary
    
    # Calculate misclassification error for the current fold
    KNN_misclassification_error <- 1 - KNN_confusion_Matrix_summary$overall["Accuracy"]
    
    #Append error to list
    error_knn<- append(error_knn,KNN_misclassification_error)
  }
  
  error_knn
  aplha<-unlist(error_knn)
  mean_aplha<-mean(aplha)
  mean_aplha
  if (mean_aplha<min_errorsknn){
    min_errorsknn<-mean_aplha
    K<-j
  }
}
min_errorsknn
K

# Re sampling Technique of k-fold cross-validation for SVM MODEL for 2000 dataset #

# Set the number of folds for cross-validation
k <- 5

# Create an empty vector to store classifications errors
SVM_misclassification_errors <- numeric(k)

# Perform k-fold cross-validation
for (i in 1:(k-1)) {  # Iterate from 1 to 8 instead of 1 to 9
  # Split the data into training and testing sets for the current fold
  set.seed(2316529)
  indices <- sample(1:nrow(New_Dataset), replace = FALSE)
  split <- cut(indices, breaks = k, labels = FALSE)
  
  # Train the SVM model on the training set
  SVM_MODEL_KFOLD <- svm(Class ~ ., data = New_Dataset[split != i, ], kernel = "radial")
  
  # Predict on the testing set
  SVM_PREDICTIONS_KFOLD <- predict(SVM_MODEL_KFOLD, newdata = New_Dataset[split == i, ])
  SVM_PREDICTED_CLASS_KFOLD<-ifelse(SVM_PREDICTIONS_KFOLD > 0.5, 1, 0)
  SVM_confusion_Matrix_summary <- confusionMatrix(as.factor(SVM_PREDICTED_CLASS_KFOLD), as.factor(New_Dataset$Class[split == i]))
  
  # Calculate misclassification error for the current fold
  SVM_misclassification_error <- 1 - SVM_confusion_Matrix_summary$overall["Accuracy"]
  SVM_misclassification_errors[i]<-SVM_misclassification_error
}

# Calculate the average misclassification error across all folds
SVM_average_misclassification_error <- unlist(SVM_misclassification_errors)
SVM_average_misclassification_error
mean(SVM_average_misclassification_error)

# Split the dataset into 0 and 1
table(reduced_data_T_Test$Class)
reduced_data_T_Test$Class[reduced_data_T_Test$Class==1]<- 0
reduced_data_T_Test$Class[reduced_data_T_Test$Class==2]<- 1

# Re sampling Technique of k-fold cross-validation LOGISTICS REGRISSION MODEL # 

#reduced_data_T_Test

# Set the number of folds for cross-validation
k <- 5

# Create an empty vector to store classification errors
LOGISTICS_misclassification_error_KFOLD <- numeric(k)

# Perform k-fold cross-validation
for (i in 1:(k-1)) {  
  # Split the data into training and testing sets for the current fold
  set.seed (2316529)
  indices <- sample(1:nrow(reduced_data_T_Test), replace = FALSE)
  split <- cut(indices, breaks = k, labels = FALSE)
  
  # Fit a logistic regression model using training data
  LOGISTICS_MODEL_KFOLD <- glm(Class ~ ., data = reduced_data_T_Test[split != i, ], family = 'binomial')
  
  # Predict probabilities on the testing set
  LOGISTICS_MODEL_PREDICT_KFOLD <- predict(LOGISTICS_MODEL_KFOLD, newdata = reduced_data_T_Test[split == i, ], type = "response")
  
  # Convert predicted probabilities to binary predictions (0 or 1)
  LOGISTICS_PREDICTED_CLASS_KFOLD <- ifelse(LOGISTICS_MODEL_PREDICT_KFOLD > 0.5, 1, 0)
  
  # Create confusion matrix
  LOGISTICS_confusion_Matrix_summary_KFOLD <- confusionMatrix(as.factor(LOGISTICS_PREDICTED_CLASS_KFOLD), as.factor(reduced_data_T_Test$Class[split == i]))
  
  # Calculate misclassification error for the current fold
  misclassification_error_KFOLD <- 1 - LOGISTICS_confusion_Matrix_summary_KFOLD$overall["Accuracy"]
  LOGISTICS_misclassification_error_KFOLD[i] <- misclassification_error_KFOLD
}

# Calculate the average misclassification error across all folds
AVERAGE_misclassification_error_logit <- mean(LOGISTICS_misclassification_error_KFOLD)
AVERAGE_misclassification_error_logit

# Re sampling Technique of k-fold cross-validation LDA MODEL #  

# Set the number of folds for cross-validation
k <- 5

# Create an empty vector to store classification errors
LDA_misclassification_errors <- numeric(k)

# Perform k-fold cross-validation
for (i in 1:(k-1)) {  
  # Split the data into training and testing sets for the current fold
  set.seed (2316529)
  indices <- sample(1:nrow(reduced_data_T_Test), replace = FALSE)
  split <- cut(indices, breaks = k, labels = FALSE)
  
  # Train the LDA model on the training set
  LDA_MODEL_KFOLD <- lda(Class ~ ., data = reduced_data_T_Test[split != i, ])
  
  # Predict on the testing set
  LDA_MODEL_PREDICT <- predict(LDA_MODEL_KFOLD, newdata = reduced_data_T_Test[split == i, ])$class
  
  # Create confusion matrix
  LDA_confusion_Matrix_summary <- confusionMatrix(as.factor(LDA_MODEL_PREDICT), as.factor(reduced_data_T_Test$Class[split == i]))
  
  # Calculate misclassification error for the current fold
  LDA_misclassification_error <- 1 - LDA_confusion_Matrix_summary$overall["Accuracy"]
  LDA_misclassification_errors[i] <- LDA_misclassification_error
}

# Calculate the average misclassification error across all folds
average_misclassification_error_LDA <- mean(LDA_misclassification_errors)
average_misclassification_error_LDA

##### Re sampling Technique of k-fold cross-validation Random Forest MODEL ###  

library(randomForest)

# Set the number of folds for cross-validation
k <- 5
min_errors_rf <- 1
error_rf=list()

# Create an empty vector to store misclassification errors for each fold
RF_misclassification_errors <- numeric(k)

# Perform k-fold cross-validation
for (j in 1:50){
  for (i in 1 :(k-1)){ 
    # Split the data into training and testing sets for the current fold
    set.seed (2316529)
    indices <- sample(1:nrow(reduced_data_T_Test), replace = FALSE)
    
    split <- cut(indices, breaks = k, labels = FALSE)
    
    # Train the Random Forest model on the training set
    RF_MODEL_KFOLD <- randomForest(Class ~ ., data = reduced_data_T_Test[split != i, ], ntree = j, set.seed (2316529))
    
    # Predict on the testing set
    RF_PREDICTIONS_KFOLD <- predict(RF_MODEL_KFOLD, newdata = reduced_data_T_Test[split == i, ], type = "response")
    RF_PREDICTED_CLASS_KFOLD <- ifelse(RF_PREDICTIONS_KFOLD > 0.5, 1, 0)
    
    # Create confusion matrix
    RF_confusion_Matrix_summary <- confusionMatrix(as.factor(RF_PREDICTED_CLASS_KFOLD), as.factor(reduced_data_T_Test$Class[split == i]))
    
    # Calculate misclassification error for the current fold
    RF_misclassification_error <- 1 - RF_confusion_Matrix_summary$overall["Accuracy"]
    
    #Append error to list
    error_rf<- append(error_rf,RF_misclassification_error)
  }
  
  error_rf
  aplha1<-unlist(error_rf)
  mean_aplha1<-mean(aplha1)
  mean_aplha1
  if (mean_aplha1<min_errors_rf){
    min_errors_rf<-mean_aplha1
    V<-j
  }
}
min_errors_rf
V

# Re sampling Technique of k-fold cross-validation for KNN MODEL ##

# Set the number of folds for cross-validation
k <- 5
min_errorsknn <- 1
error_knn=list()

# Create an empty vector to store classification errors
KNN_misclassification_errors <- numeric(k)
for (j in 1:50){
  for (i in 1 :(k-1)){
    
    # Split the data into training and testing sets for the current fold
    set.seed(2316529)
    indices <- sample(1:nrow(reduced_data_T_Test), replace = FALSE)
    
    split <- cut(indices, breaks = k, labels = FALSE)
    
    # Train the k-NN model on the training set
    
    KNN_MODEL_KFOLD <- knn(train = reduced_data_T_Test[split != i, -which(names(reduced_data_T_Test) == "Class")], 
                           test = reduced_data_T_Test[split == i, -which(names(reduced_data_T_Test) == "Class")], 
                           cl = reduced_data_T_Test$Class[split != i],k=j)  
    
    # Predict on the testing set
    KNN_PREDICTIONS_KFOLD <- KNN_MODEL_KFOLD
    
    # Create confusion matrix
    KNN_confusion_Matrix_summary <- confusionMatrix(as.factor(KNN_PREDICTIONS_KFOLD), as.factor(reduced_data_T_Test$Class[split == i]))
    KNN_confusion_Matrix_summary
    
    # Calculate misclassification error for the current fold
    KNN_misclassification_error <- 1 - KNN_confusion_Matrix_summary$overall["Accuracy"]
    
    #Append error to list
    error_knn<- append(error_knn,KNN_misclassification_error)
  }
  
  error_knn
  aplha<-unlist(error_knn)
  mean_aplha<-mean(aplha)
  mean_aplha
  if (mean_aplha<min_errorsknn){
    min_errorsknn<-mean_aplha
    K<-j
  }
}
min_errorsknn
k

# Re sampling Technique of k-fold cross-validation for SVM MODEL #

# Set the number of folds for cross-validation
k <- 5

# Create an empty vector to store classifications errors
SVM_misclassification_errors <- numeric(k)

# Perform k-fold cross-validation
for (i in 1:(k-1)) {  # Iterate from 1 to 8 instead of 1 to 9
  # Split the data into training and testing sets for the current fold
  set.seed(2316529)
  indices <- sample(1:nrow(reduced_data_T_Test), replace = FALSE)
  split <- cut(indices, breaks = k, labels = FALSE)
  
  # Train the SVM model on the training set
  SVM_MODEL_KFOLD <- svm(Class ~ ., data = reduced_data_T_Test[split != i, ], kernel = "radial")
  
  # Predict on the testing set
  SVM_PREDICTIONS_KFOLD <- predict(SVM_MODEL_KFOLD, newdata = reduced_data_T_Test[split == i, ])
  SVM_PREDICTED_CLASS_KFOLD<-ifelse(SVM_PREDICTIONS_KFOLD > 0.5, 1, 0)
  SVM_confusion_Matrix_summary <- confusionMatrix(as.factor(SVM_PREDICTED_CLASS_KFOLD), as.factor(reduced_data_T_Test$Class[split == i]))
  
  # Calculate misclassification error for the current fold
  SVM_misclassification_error <- 1 - SVM_confusion_Matrix_summary$overall["Accuracy"]
  SVM_misclassification_errors[i]<-SVM_misclassification_error
}


# Calculate the average misclassification error across all folds
SVM_average_misclassification_error <- unlist(SVM_misclassification_errors)
SVM_average_misclassification_error
mean(SVM_average_misclassification_error)

#t-distributed Stochastic Neighbor Embedding REDUCTION TECHQUE
# Load necessary libraries

library(Rtsne)

# Perform t-SNE dimensionality reduction
tsne <- Rtsne(Remove_Null_values,perplexity=10)
reduced_data_tsne<-data.frame(tsne[['Y']])
reduced_data_tsne$Class <- Remove_Null_values$Class


# Split the data into training and testing sets(Model Building)
library(caTools)
split_TTest_REDUCED_DATASET <- sample.split(reduced_data_tsne$Class, SplitRatio = 0.7, set.seed (2316529))
TRAIN_TTest_REDUCED_DATASET  <- subset(reduced_data_tsne, split_TTest_REDUCED_DATASET  == TRUE)
TEST_TTest_REDUCED_DATASET  <- subset(reduced_data_tsne, split_TTest_REDUCED_DATASET  == FALSE)

table(TRAIN_TTest_REDUCED_DATASET$Class)
TRAIN_TTest_REDUCED_DATASET$Class[TRAIN_TTest_REDUCED_DATASET$Class==1]<- 0
TRAIN_TTest_REDUCED_DATASET$Class[TRAIN_TTest_REDUCED_DATASET$Class==2]<- 1
TEST_TTest_REDUCED_DATASET$Class[TEST_TTest_REDUCED_DATASET$Class==1]<- 0
TEST_TTest_REDUCED_DATASET$Class[TEST_TTest_REDUCED_DATASET$Class==2]<-1

 #import libraries#

library(caret)
library(MASS)
library(randomForest)
library(class)
library(e1071)
# LOGISTICS REGRISSION MODEL T-test ##

# Fit a logistic regression model using training data
LOGISTICS_MODEL <- glm(Class ~ ., data = TRAIN_TTest_REDUCED_DATASET, family = 'binomial')

# Predict probabilities on test data
LOGISTICS_MODEL_PREDICT <- predict(LOGISTICS_MODEL, newdata = TEST_TTest_REDUCED_DATASET, type = "response")

# Convert predicted probabilities to binary predictions (0 or 1)
LOGISTICS_PREDICTED_CLASS <- ifelse(LOGISTICS_MODEL_PREDICT > 0.5, 1, 0)

# Create confusion matrix
LOGISTICS_confusion_Matrix_summary <- confusionMatrix(as.factor(LOGISTICS_PREDICTED_CLASS), as.factor(TEST_TTest_REDUCED_DATASET$Class))
LOGISTICS_confusion_Matrix_summary

# Calculate misclassification error
LOGISTICS_misclassification_error <- 1 - LOGISTICS_confusion_Matrix_summary$overall["Accuracy"]
LOGISTICS_misclassification_error

# LDA MODEL ##

# Fit LDA model
LDA_MODEL <- lda(Class ~ ., data = TRAIN_TTest_REDUCED_DATASET)

# Predict on test data
LDA_MODEL_PREDICT <- predict(LDA_MODEL, newdata = TEST_TTest_REDUCED_DATASET)$class

# Create confusion matrix
LDA_confusion_Matrix_summary <- confusionMatrix(as.factor(LDA_MODEL_PREDICT), as.factor(TEST_TTest_REDUCED_DATASET$Class))
LDA_confusion_Matrix_summary

# Calculate misclassification error
misclassification_error_LDA <- 1 - LDA_confusion_Matrix_summary$overall["Accuracy"]
misclassification_error_LDA

# QDA MODEL #

# Fit a QDA model using training data
QDA_MODEL <- qda(Class ~ ., data = TRAIN_TTest_REDUCED_DATASET)

# Predict using the QDA model
QDA_PREDICTION <- predict(QDA_MODEL, newdata = TEST_TTest_REDUCED_DATASET)$class
QDA_PREDICTION

# Create confusion matrix
confusion_Matrix_QDA <- table(QDA_PREDICTION, TEST_TTest_REDUCED_DATASET$Class)
confusion_Matrix_QDA

# Calculate misclassification error for QDA
misclassification_error_QDA <- mean(QDA_PREDICTION != TEST_TTest_REDUCED_DATASET$Class)
misclassification_error_QDA

# RANDOM FOREST MODEL #
# Train the Random Forest model
RF_MODEL <- randomForest(Class ~ ., data = TRAIN_TTest_REDUCED_DATASET,ntree  = 3,set.seed(2316529))

# Predict on test data
RF_PREDICTIONS <- predict(RF_MODEL, newdata = TEST_TTest_REDUCED_DATASET,type= "response")
RF_PREDICTIONS
RF_PREDICTED_CLASS<-ifelse(RF_PREDICTIONS > 0.5, 1, 0)

# Create confusion matrix
RF_confusion_Matrix_summary <- confusionMatrix(as.factor(RF_PREDICTED_CLASS), as.factor(TEST_TTest_REDUCED_DATASET$Class))
RF_confusion_Matrix_summary

# Or if you want to use table function
RF_Conf_matrix <- table(RF_PREDICTED_CLASS, TEST_TTest_REDUCED_DATASET$Class)
RF_Conf_matrix

# Calculate misclassification error
RF_misclassification_error <- 1 - RF_confusion_Matrix_summary$overall["Accuracy"]
RF_misclassification_error

# KNN MODEL #

# Train the k-NN model
KNN_MODEL <- knn(train = TRAIN_TTest_REDUCED_DATASET[, -which(names(TRAIN_TTest_REDUCED_DATASET) == "Class")], 
                 test = TEST_TTest_REDUCED_DATASET[, -which(names(TEST_TTest_REDUCED_DATASET) == "Class")], 
                 cl = TRAIN_TTest_REDUCED_DATASET$Class, 
                 k = 3)  

# Predict on test data
KNN_PREDICTIONS <- KNN_MODEL

# Create confusion matrix
KNN_confusion_Matrix_summary <- confusionMatrix(as.factor(KNN_PREDICTIONS), as.factor(TEST_TTest_REDUCED_DATASET$Class))
KNN_confusion_Matrix_summary

# Create confusion matrix
KNN_confusion_Matrix <- table(KNN_PREDICTIONS, TEST_TTest_REDUCED_DATASET$Class)
KNN_confusion_Matrix

# Calculate misclassification error
KNN_misclassification_error <- 1 - KNN_confusion_Matrix_summary$overall["Accuracy"]
KNN_misclassification_error


# SVM MODEL #

# Train the SVM model
SVM_MODEL <- svm(Class ~ ., data = TRAIN_TTest_REDUCED_DATASET, kernel = "radial")

# Predict on test data
SVM_PREDICTIONS <- predict(SVM_MODEL, newdata = TEST_TTest_REDUCED_DATASET)
SVM_PREDICTED_CLASS<-ifelse(SVM_PREDICTIONS > 0.5, 1, 0)

# Create confusion matrix
SVM_confusion_Matrix_summary <- confusionMatrix(as.factor(SVM_PREDICTED_CLASS), as.factor(TEST_TTest_REDUCED_DATASET$Class))
SVM_confusion_Matrix_summary

# Create confusion matrix
SVM_confusion_Matrix <- table(SVM_PREDICTED_CLASS, TEST_TTest_REDUCED_DATASET$Class)
SVM_confusion_Matrix

# Calculate misclassification error
SVM_misclassification_error <- 1 - SVM_confusion_Matrix_summary$overall["Accuracy"]
SVM_misclassification_error


# # Split the data set into 0 and 1
table(reduced_data_tsne$Class)
reduced_data_tsne$Class[reduced_data_tsne$Class==1]<- 0
reduced_data_tsne$Class[reduced_data_tsne$Class==2]<- 1

# Re sampling Technique of k-fold cross-validation LOGISTICS REGRISSION MODEL #  

#reduced_data_T_Test

# Set the number of folds for cross-validation
k <- 5

# Create an empty vector to store classification errors
LOGISTICS_misclassification_error_KFOLD <- numeric(k)

# Perform k-fold cross-validation
for (i in 1:(k-1)) {  
  # Split the data into training and testing sets for the current fold
  set.seed (2316529)
  indices <- sample(1:nrow(reduced_data_tsne), replace = FALSE)
  split <- cut(indices, breaks = k, labels = FALSE)
  
  # Fit a logistic regression model using training data
  LOGISTICS_MODEL_KFOLD <- glm(Class ~ ., data = reduced_data_tsne[split != i, ], family = 'binomial')
  
  # Predict probabilities on the testing set
  LOGISTICS_MODEL_PREDICT_KFOLD <- predict(LOGISTICS_MODEL_KFOLD, newdata = reduced_data_tsne[split == i, ], type = "response")
  
  # Convert predicted probabilities to binary predictions (0 or 1)
  LOGISTICS_PREDICTED_CLASS_KFOLD <- ifelse(LOGISTICS_MODEL_PREDICT_KFOLD > 0.5, 1, 0)
  
  # Create confusion matrix
  LOGISTICS_confusion_Matrix_summary_KFOLD <- confusionMatrix(as.factor(LOGISTICS_PREDICTED_CLASS_KFOLD), as.factor(reduced_data_tsne$Class[split == i]))
  
  # Calculate misclassification error for the current fold
  misclassification_error_KFOLD <- 1 - LOGISTICS_confusion_Matrix_summary_KFOLD$overall["Accuracy"]
  LOGISTICS_misclassification_error_KFOLD[i] <- misclassification_error_KFOLD
}

# Calculate the average misclassification error across all folds
AVERAGE_misclassification_error_logit <- mean(LOGISTICS_misclassification_error_KFOLD)
AVERAGE_misclassification_error_logit

# Re sampling Technique of k-fold cross-validation LDA MODEL #

# Set the number of folds for cross-validation
k <- 5

# Create an empty vector to store classification errors
LDA_misclassification_errors <- numeric(k)

# Perform k-fold cross-validation
for (i in 1:(k-1)) {  
  # Split the data into training and testing sets for the current fold
  set.seed (2316529)
  indices <- sample(1:nrow(reduced_data_tsne), replace = FALSE)
  split <- cut(indices, breaks = k, labels = FALSE)
  
  # Train the LDA model on the training set
  LDA_MODEL_KFOLD <- lda(Class ~ ., data = reduced_data_tsne[split != i, ])
  
  # Predict on the testing set
  LDA_MODEL_PREDICT <- predict(LDA_MODEL_KFOLD, newdata = reduced_data_tsne[split == i, ])$class
  
  # Create confusion matrix
  LDA_confusion_Matrix_summary <- confusionMatrix(as.factor(LDA_MODEL_PREDICT), as.factor(reduced_data_tsne$Class[split == i]))
  
  # Calculate misclassification error for the current fold
  LDA_misclassification_error <- 1 - LDA_confusion_Matrix_summary$overall["Accuracy"]
  LDA_misclassification_errors[i] <- LDA_misclassification_error
}

# Calculate the average misclassification error across all folds
average_misclassification_error_LDA <- mean(LDA_misclassification_errors)
average_misclassification_error_LDA

# Re sampling Technique of k-fold cross-validation Random Forest MODEL # 

library(randomForest)

# Set the number of folds for cross-validation
k <- 5
min_errors_rf <- 1
error_rf=list()

# Create an empty vector to store misclassification errors for each fold
RF_misclassification_errors <- numeric(k)

# Perform k-fold cross-validation
for (j in 1:50){
  for (i in 1 :(k-1)){ 
    # Split the data into training and testing sets for the current fold
    set.seed (2316529)
    indices <- sample(1:nrow(reduced_data_tsne), replace = FALSE)
    
    split <- cut(indices, breaks = k, labels = FALSE)
    
    # Train the Random Forest model on the training set
    RF_MODEL_KFOLD <- randomForest(Class ~ ., data = reduced_data_tsne[split != i, ], ntree = j, set.seed (2316529))
    
    # Predict on the testing set
    RF_PREDICTIONS_KFOLD <- predict(RF_MODEL_KFOLD, newdata = reduced_data_tsne[split == i, ], type = "response")
    RF_PREDICTED_CLASS_KFOLD <- ifelse(RF_PREDICTIONS_KFOLD > 0.5, 1, 0)
    
    # Create confusion matrix
    RF_confusion_Matrix_summary <- confusionMatrix(as.factor(RF_PREDICTED_CLASS_KFOLD), as.factor(reduced_data_tsne$Class[split == i]))
    
    # Calculate misclassification error for the current fold
    RF_misclassification_error <- 1 - RF_confusion_Matrix_summary$overall["Accuracy"]
    
    #Append error to list
    error_rf<- append(error_rf,RF_misclassification_error)
  }
  
  error_rf
  aplha1<-unlist(error_rf)
  mean_aplha1<-mean(aplha1)
  mean_aplha1
  if (mean_aplha1<min_errors_rf){
    min_errors_rf<-mean_aplha1
    V<-j
  }
}
min_errors_rf
V

# Re sampling Technique of k-fold cross-validation for KNN MODEL #

# Set the number of folds for cross-validation
k <- 5
min_errorsknn <- 1
error_knn=list()

# Create an empty vector to store classification errors
KNN_misclassification_errors <- numeric(k)
for (j in 1:50){
  for (i in 1 :(k-1)){
    
    # Split the data into training and testing sets for the current fold
    set.seed(2316529)
    indices <- sample(1:nrow(reduced_data_tsne), replace = FALSE)
    
    split <- cut(indices, breaks = k, labels = FALSE)
    
    # Train the k-NN model on the training set
    
    KNN_MODEL_KFOLD <- knn(train = reduced_data_tsne[split != i, -which(names(reduced_data_tsne) == "Class")], 
                           test = reduced_data_tsne[split == i, -which(names(reduced_data_tsne) == "Class")], 
                           cl = reduced_data_tsne$Class[split != i],k=j)  
    
    # Predict on the testing set
    KNN_PREDICTIONS_KFOLD <- KNN_MODEL_KFOLD
    
    # Create confusion matrix
    KNN_confusion_Matrix_summary <- confusionMatrix(as.factor(KNN_PREDICTIONS_KFOLD), as.factor(reduced_data_tsne$Class[split == i]))
    KNN_confusion_Matrix_summary
    
    # Calculate misclassification error for the current fold
    KNN_misclassification_error <- 1 - KNN_confusion_Matrix_summary$overall["Accuracy"]
    
    #Append error to list
    error_knn<- append(error_knn,KNN_misclassification_error)
  }
  
  error_knn
  aplha<-unlist(error_knn)
  mean_aplha<-mean(aplha)
  mean_aplha
  if (mean_aplha<min_errorsknn){
    min_errorsknn<-mean_aplha
    K<-j
  }
}
min_errorsknn
k

# Re sampling Technique of k-fold cross-validation for SVM MODEL #

# Set the number of folds for cross-validation
k <- 5

# Create an empty vector to store classifications errors
SVM_misclassification_errors <- numeric(k)

# Perform k-fold cross-validation
for (i in 1:(k-1)) {  # Iterate from 1 to 8 instead of 1 to 9
  # Split the data into training and testing sets for the current fold
  set.seed(2316529)
  indices <- sample(1:nrow(reduced_data_tsne), replace = FALSE)
  split <- cut(indices, breaks = k, labels = FALSE)
  
  # Train the SVM model on the training set
  SVM_MODEL_KFOLD <- svm(Class ~ ., data = reduced_data_tsne[split != i, ], kernel = "radial")
  
  # Predict on the testing set
  SVM_PREDICTIONS_KFOLD <- predict(SVM_MODEL_KFOLD, newdata = reduced_data_tsne[split == i, ])
  SVM_PREDICTED_CLASS_KFOLD<-ifelse(SVM_PREDICTIONS_KFOLD > 0.5, 1, 0)
  SVM_confusion_Matrix_summary <- confusionMatrix(as.factor(SVM_PREDICTED_CLASS_KFOLD), as.factor(reduced_data_tsne$Class[split == i]))
  
  # Calculate misclassification error for the current fold
  SVM_misclassification_error <- 1 - SVM_confusion_Matrix_summary$overall["Accuracy"]
  SVM_misclassification_errors[i]<-SVM_misclassification_error
}

# Calculate the average misclassification error across all folds
SVM_average_misclassification_error <- unlist(SVM_misclassification_errors)
SVM_average_misclassification_error
mean(SVM_average_misclassification_error)

# Re sampling Technique of k-fold cross-validation QDA MODEL # 

# Set the number of folds for cross-validation
k <- 5

# Create an empty vector to store classification errors
QDA_misclassification_errors <- numeric(k)

# Perform k-fold cross-validation
for (i in 1:(k-1)) {  
  # Split the data into training and testing sets for the current fold
  set.seed (2316529)
  indices <- sample(1:nrow(reduced_data_tsne), replace = FALSE)
  split <- cut(indices, breaks = k, labels = FALSE)
  
  # Train the LDA model on the training set
  QDA_MODEL_KFOLD <- qda(Class ~ ., data = reduced_data_tsne[split != i, ])
  
  # Predict on the testing set
  QDA_MODEL_PREDICT <- predict(QDA_MODEL_KFOLD, newdata = reduced_data_tsne[split == i, ])$class
  
  # Create confusion matrix
  QDA_confusion_Matrix_summary <- confusionMatrix(as.factor(QDA_MODEL_PREDICT), as.factor(reduced_data_tsne$Class[split == i]))
  
  # Calculate misclassification error for the current fold
  QDA_misclassification_error <- 1 - QDA_confusion_Matrix_summary$overall["Accuracy"]
  QDA_misclassification_errors[i] <- QDA_misclassification_error
}

# Calculate the average misclassification error across all folds
average_misclassification_error_QDA <- mean(QDA_misclassification_errors)
average_misclassification_error_QDA

