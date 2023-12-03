setwd('D:/MSc Statistics & Data Science/Second Year/Statistical Consulting/Client meetings/Final data analysis')

library(foreign)
library(lattice)
library(geepack)
library(aod)
library(GLMMadaptive)
library(nlme)
library(corrplot)
library(ggplot2)

# Read data
csv_data <- read.csv('sumscores_output_scaled_lesions_only.csv')
PCA_object = prcomp(csv_data[,11:34])
object_scores = PCA_object$x[,1:3]



meta_part = csv_data[,1:10]
PCA_part = object_scores
CATPCA_data = cbind(meta_part,PCA_part)


# Set up the layout for the plots
par(mfrow = c(1, ncol(PCA_part)))  # 1 row, 3 columns

# Plot each column separately with different colors for each category
for (col in names(PCA_part)) {
  plot(PCA_part[[col]], col = as.factor(meta_part$types ), 
       main = paste("Column", col), pch = 16)
#  legend("topright", legend = levels(as.factor(meta_part$types )), 
#         col = 1:length(levels(as.factor(meta_part$types ))), pch = 16, title = "Type")
}
# Change the components names according to interpretations
names(CATPCA_data)[11]='Membrane'
names(CATPCA_data)[12]='Storage'
names(CATPCA_data)[13]='Other'

#######################################################################
#................  Fitting multi GEE models ............................
#######################################################################
df = CATPCA_data
df=df[order(df$donor, df$X), ]

comparisons_9types=list()
for (i in 1:length(unique(df$types))){
  # prepare the datasets that correspond to each of the comparisons
  comparisons_9types[[i]]=df
  names(comparisons_9types)[[i]]=unique(df$types)[i]
  
  # setting the type columns to 1 for the studied type and 0 otherwise
  comparisons_9types[[i]][comparisons_9types[[i]]$types!=unique(df$types)[i],'types']=0
  comparisons_9types[[i]][comparisons_9types[[i]]$types==unique(df$types)[i],'types']=1
  
  # setting the factor variables
  comparisons_9types[[i]]$types=as.integer(comparisons_9types[[i]]$types)
  comparisons_9types[[i]]$sex = as.factor(comparisons_9types[[i]]$sex)
  
}

# Getting the prob for each sample on each class, use highest prob as prediction
probabilities=data.frame(index=1:nrow(df))


for (i in 1:length(comparisons_9types)) {
  gee <- geeglm(types ~ -1 +age +Membrane +Storage + Other ,
                data = comparisons_9types[[i]],
                id = donor, family = binomial("logit"),
                corstr = "exchangeable")
  
  
  # Predicted probabilities
  predictions <- predict(gee, newdata = df, type = "response")
  probabilities[,names(comparisons_9types)[i]]=predictions
  
}
# Getting rid of the index col
probabilities=probabilities[,-1]
for (i in 1:nrow(probabilities) ){
  idx = which.max(probabilities[i,])
  pred = names(idx)
  df[i,'prediction']=pred
}
accuracy = sum(df$types==df$prediction)/nrow(df)
accuracy # 0.54

#################################################################
#............... with cross validation ..........................
#################################################################
# To calculate the accuracy with tolerance; to see if the model predicted the right 
# type as the one with maximum probability or with second or third maximum probability.

# This function returns the three types with the highest probability
get_max_values <- function(row) {
  max_values <- sort(row, decreasing = TRUE)[1]
  max_values = rep(max_values,length(row))
  indx = which(max_values==row)
  max_column = names(probabilities)[indx]
  
  max2_values <- sort(row, decreasing = TRUE)[2]
  max2_values = rep(max2_values,length(row))
  indx = which(max2_values==row)
  max2_column = names(probabilities)[indx]
  
  max3_values <- sort(row, decreasing = TRUE)[3]
  max3_values = rep(max3_values,length(row))
  indx = which(max3_values==row)
  max3_column = names(probabilities)[indx]
  
  return(c(max_column,max2_column,max3_column))
}

# Creating a dataframe to save the prob of every sample
probabilities=data.frame(index=1:nrow(df))

# To extract the coefficients with their CI
coef_df = data.frame(type = rep(0,6),age=rep(0,6),Membrane=rep(0,6),Storage=rep(0,6))

# Applying leave-one-out cross validation
for (j in 1:nrow(df)){
  for (i in 1:length(comparisons_9types)) {
    gee <- geeglm(types ~ -1 + age +Membrane +Storage ,
                  data = comparisons_9types[[i]][-j,],
                  id = donor, family = binomial("logit"),
                  corstr = "exchangeable")
    
    # The extraction of the coefficients with CI
    estim = exp(summary(gee)$coefficients[,1])
    estim_min = round(estim - 1.96*exp(summary(gee)$coefficients[,2]),2)
    estim_max = round(estim + 1.96*exp(summary(gee)$coefficients[,2]),2)
    ci_coef = rep(0,4)
    ci_coef[1]=names(comparisons_9types)[[i]]
    for (co in 1:length(estim)){
      ci_coef[co+1] = paste('[',estim_min[co],',',estim_max[co],']',' Point Estimate: ',round(estim[co],2))
    } 
    coef_df[i,]=ci_coef
    
    
    
    # Predicted probabilities
    predictions <- predict(gee, newdata = comparisons_9types[[i]][j,], type = "response")
    probabilities[j,names(comparisons_9types)[i]]=predictions
    
  }}

# Deleting the index column
probabilities=probabilities[,-1]

# Getting the names of the types with the highest probs

# To get the one type with the highest prob
predictions = data.frame(index=1:nrow(df))
# To get the first three types with the highest prob
predictions_with_tolerance = data.frame(index=1:nrow(df),
                                        score_with_tolerance_two=rep(0,nrow(df)),
                                        score_with_tolerance_three=rep(0,nrow(df)))

for (i in 1:nrow(probabilities) ){
  # type with max prob
  idx = which.max(probabilities[i,])
  pred = names(idx)
  predictions[i,'prediction']=pred
  # to check if the true type is within the two types with highest prob
  row = probabilities[i,]
  max_two_col = get_max_values(row)
  predictions_with_tolerance[i,'first_prediction']=max_two_col[1]
  predictions_with_tolerance[i,'second_prediction']=max_two_col[2]
  predictions_with_tolerance[i,'third_prediction']=max_two_col[3]
  if(df[i,'types'] %in% max_two_col[1:2]){
    predictions_with_tolerance[i,'score_with_tolerance_two']=1
  }
  
  if(df[i,'types'] %in% max_two_col){
    predictions_with_tolerance[i,'score_with_tolerance_three']=1
  }
}
# Adding the true labels to the predictions so the view is complete
predictions_with_tolerance[,'True type']=df$types
predictions[,'True type']=df$types

# returning the accuracies
accuracy = sum(df$types==predictions$prediction)/nrow(df)
accuracy #  0.38
accuracy_with_tolerance_two = mean(predictions_with_tolerance$score_with_tolerance_two)
accuracy_with_tolerance_three = mean(predictions_with_tolerance$score_with_tolerance_three)
accuracy_with_tolerance_two # 62
accuracy_with_tolerance_three # 82

# the result is that lipids have the potential to classify, but more samples are needed.
# do contengency table with foamy ordering
# do table for e^(coef) with confidence intervals
# see how each true type was classified into what three types?

# contengency table with foamy ordering

foamy_types = unique(CATPCA_data[CATPCA_data$Morphology.microglia=='foamy','types'])
non_foamy_types = unique(CATPCA_data[CATPCA_data$Morphology.microglia=='non_foamy','types'])
inactive_types = unique(CATPCA_data[CATPCA_data$Morphology.microglia=='inactive','types'])

contengency_table = data.frame('3.3'=rep(0,6),'2.3'=rep(0,6),'2.1_2.2'=rep(0,6),'3.1_3.2'=rep(0,6),'6'=rep(0,6),'4'=rep(0,6))
names(contengency_table) = c('3.3','2.3','2.1_2.2','3.1_3.2','6','4')
rownames(contengency_table) = c('3.3','2.3','2.1_2.2','3.1_3.2','6','4')
for (i in 1:nrow(predictions)){
  pred_index = which(colnames(contengency_table)==predictions$prediction[i])
  true_index = which(colnames(contengency_table)==predictions$`True type`[i])
  contengency_table[pred_index,true_index] = contengency_table[pred_index,true_index] + 1
  contengency_table
}
contengency_table[,"predictions sum"]=rowSums(contengency_table)
contengency_table = rbind(contengency_table,colSums(contengency_table))
rownames(contengency_table)[7]='True sum'
contengency_table

