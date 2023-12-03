setwd('D:/MSc Statistics & Data Science/Second Year/Statistical Consulting/Client meetings/Final data analysis')
#loading packages
library(tidyverse)
library(data.table)
library(ggplot2)
library(cowplot)
library(patchwork)
library(readxl)
library(dplyr)
library(geepack)
Sensitivity_analysis_min_impute= function(){
  #loading the metadata + some cleaning up and refining
  metadata <- readxl::read_xlsx("Metadata table2.xlsx") %>% data.frame() 
  row.names(metadata) <- metadata$Unifying_code
  metadata$Morphology.microglia[is.na(metadata$Morphology.microglia)] <- "inactive"
  metadata$Morphology.microglia <- factor(metadata$Morphology.microglia, levels=c("inactive","non_foamy","foamy"))
  Lesion.types.6 <- c("CWM","NAWM","2","3","4","6") 
  metadata$Lesion_type_6 <- factor(metadata$Lesion_type_6, levels=Lesion.types.6)
  Lesion.types.9 <- c("CWM","NAWM","PLWM","2.1_2.2","2.3","3.1_3.2","3.3","4","6") 
  metadata$Lesion_type_9 <- factor(metadata$Lesion_type_9, levels=Lesion.types.9)
  metadata <- metadata %>% #making the pmd a numerical decimal instead of 08:30 character
    separate(pmd,sep=":",into = c("H", "M")) %>%
    mutate(H = as.numeric(H)) %>% 
    mutate(M = as.numeric(M)/60) %>% 
    mutate(pmd = H+M) %>% 
    select(-"H",-"M")
  
  #filtered metadata for lipidomics
  #not all samples have been analyzed by lipidomics, so we filter them 
  metadata.lipids <- metadata %>% 
    filter(sample_lipidomics != "") %>% 
    arrange(sample_lipidomics)
  metadata.ABPP <- metadata %>% 
    filter(sample_ABPP != "") %>% 
    arrange(sample_ABPP)
  
  #Color schemes for nice plots 
  Colors.lesion.6g <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FFF400","#FF7F00")
  names(Colors.lesion.6g) <- Lesion.types.6
  Colors.lesion.9g <- c("#E41A1C","#0072B2","#A6CEE3","#B2DF8A","green4", "#b391c4", "darkmagenta","#FFF400", "#D55E00")
  names(Colors.lesion.9g) <- Lesion.types.9
  Colors.morph <- c("#E1E1E1","#b58709","#0940b5")
  names(Colors.morph) <- c("inactive","non_foamy","foamy")
  
  # loading the categorization of lipid into lipid-classes 
  # only necessary for the heatmaps with lipids in there
  lipidclass <- readxl::read_xlsx("lipidclass MS.xlsx")
  lipidclass$class <- factor(lipidclass$class,levels=c("PC","PE","PS","PI","PG & BMP",
                                                       "LysoPC","LysoPE","LysoPS","LysoPI","LysoPG",
                                                       "NAPE","LysoNAPE","GPNAE","NAE",
                                                       "TAG","DAG","MAG","FFA","Oxylipins",
                                                       "SM","CER","Cerebroside","Sterols","CE"))
  lipidclass$class.broad <- factor(lipidclass$`general class`,
                                   levels=c("Phospholipids",
                                            "Lysophospholipids",
                                            "Sphingolipids",
                                            "NAPEs",
                                            "Sterols",
                                            "CEs",
                                            "TAGs",
                                            "DAGs",
                                            "eCBs",
                                            "FFAs",          
                                            "Oxylipins"))
  
  # Import data frame
  lipidomics.data.notimputed <- readxl::read_xlsx("20221115 Xinyu DATA.xlsx") %>%
    arrange(`Sample code`) %>%
    select(-"Sample code",-"type",-"Protein conc. (mg/mL)") %>% #remove non-lipid variables
    replace(.==0,NA) %>% #zero's are not really 0, they are not measured so NA 
    log2() %>% 
    scale(center=T,scale=F) %>% 
    data.frame(check.names = F)#read in the data file 
  
  #setting the right names, again this is because not all samples have been measured, 
  #so there are multiple names for the same samples. This sets the rownames in line with the metadata
  rownames(lipidomics.data.notimputed) <- metadata %>% 
    filter(sample_lipidomics != "") %>% 
    arrange(sample_lipidomics) %>% 
    row.names()
  
  # CLEANING UP DATA ##
  ## data imputation with sampling from a normal distribution around half of the minimum value of a certain lipid with an SD 1/3 of the original SD
  #this is based on the assumption that if a lipid was not quantified, that it was below the treshold of detection
  #We can not know the treshold, but at least it isLOWER than the lowest value in the dataset
  set.seed(123)
  imputation <- lipidomics.data.notimputed %>% as.matrix() #log2 transformation ensures normal distribution (assumption)
  for (i in 1:ncol(imputation)) { 
    imputation[,i] <- replace(imputation[,i], #select the a column in the data
                              is.na(imputation[,i]), #determine missing values
                              rnorm(   # replace these missing values with: 
                                sum(is.na(imputation[,i])),   #determining how many missing values are needed for a specific lipid (how many values have to be generated?)
                                mean=min(imputation[,i],na.rm=T), #mean is min
                                sd=sd(imputation[,i],na.rm=T)/3)) #SD is original SD divided by 3
  }
  lipidomics.data <- imputation %>%
    scale(center=T, scale=F) %>% #center the data, so the average per lipid is zero. NO SCALING for unit variance
    data.frame(check.names = FALSE)
  
  #write.csv(lipidomics.data,"lipidomics data.csv")
  
  #The following lines are to produce a dataframe that is easier to handle than the original dataset. 
  #To this end, all lipids that belong to the same class are summed. This dataset is quite representative of the real dataset. 
  #however, subtle changes might be missed here. Only for exploratory purposes.
  #it produces a data frame called "results.summed" with only 24 variables (the lipid classes) in the same 100 samples.
  sumdata <- matrix(ncol =length(levels(lipidclass$class)),nrow=100) %>% data.table()
  colnames(sumdata) <- levels(lipidclass$class) 
  
  for (i in levels(lipidclass$class)) { 
    sumdata[,i] <- imputation[,lipidclass$class == i] %>% 
      2^. %>% 
      apply(1,sum,na.rm=T) }
  lipidomics.data.summed <- sumdata %>% log2() %>% scale(center=T,scale=F) %>% data.frame()
  lipidomics.data.summed.scaled <- sumdata %>% log2() %>% scale(center=T,scale=T) %>% data.frame() # here we scale the sumscores
  
  rownames(lipidomics.data.summed) = rownames(lipidomics.data)
  rownames(lipidomics.data.summed.scaled) = rownames(lipidomics.data)
  
  #write.csv(lipidomics.data.summed,"sumdata.csv")
  #write.csv(lipidomics.data.summed.scaled,"sumdata_scaled.csv")
  
  max(lipidomics.data)
  max(sumdata)
  max(lipidomics.data.summed)
  max(lipidomics.data.summed.scaled)
  
  #######################
  ### donors metadata ###
  #######################
  # Load the donor ID
  meta_data <- read_excel("Metadata table2.xlsx")
  # matching each sample with its donor through the variable 'sample_lipidomics':
  #  and then retrieving the corresponding metadata
  donors=data.frame('sample_lipidomics'=c(),'NBB donor ID'=c(),'sex'=c(),'age'=c(),'brain weight'=c(),'pmd'=c(),"pH CSF"=c(),'Lesion_type_9'=c(),"Morphology microglia"=c())
  for (i in 1:nrow(lipidomics.data.summed)){
    sample_id=rownames(lipidomics.data.summed)[i]
    index=which(meta_data$Unifying_code==sample_id)
    current_donor=meta_data[index,c('NBB donor ID','sex','age','brain weight','pmd',"pH CSF",'Lesion_type_9',"Morphology microglia")]
    current_donor=cbind(sample_id=sample_id,current_donor)
    donors=rbind(donors,current_donor)
  }
  
  
  #######################
  ### some missingness ##
  #######################
  # Filling the missing values in the "Morphology microglia" with the 'inactive' label
  sum(is.na(donors$"Morphology microglia"))
  donors$"Morphology microglia" <- replace(donors$"Morphology microglia", is.na(donors$"Morphology microglia"), 'inactive')
  sum(is.na(donors$"Morphology microglia"))
  
  # Filling the 5 missing values in the 'pH CSF' with the mean value
  sum(is.na(donors$"pH CSF"))
  donors$`pH CSF` <- replace(donors$`pH CSF`, is.na(donors$`pH CSF`), mean(donors$`pH CSF`, na.rm = TRUE))
  sum(is.na(donors$"pH CSF"))
  
  # There is one donor missing                         ########################
  length(unique(donors$`NBB donor ID`))                ##### question here ####
  ########################   
  
  # Initiate the final dataframe that will have the components
  #  beginning with the metadata info
  PCA_per_class=donors
  # a variable we will need later on (ease of retrieval)
  types=PCA_per_class$Lesion_type_9
  
  
  
  ##########################################################
  #####  Processing of varibles: donorID,sampleID,pmd ######
  ##########################################################
  dim(PCA_per_class)
  # Changing the sample ID and the donor ID to numeric 
  #  so that the GEE model runs properly.
  # First step taking out the non-numeric characters
  for (i in 1:nrow(PCA_per_class)){
    # processing the donor id
    original_string <- PCA_per_class$`NBB donor ID`[i]
    # Create a new string with the third character removed
    modified_string <- paste0(substr(original_string, 1, 4), 
                              substr(original_string, 6, nchar(original_string)))
    # Print the modified string
    PCA_per_class$`NBB donor ID`[i] <- modified_string
    
    
    # processing the sample code
    original_string <- PCA_per_class$sample_id[i]
    # Create a new string with the third character removed
    modified_string <- paste0(substr(original_string, 6, 8))
    # Print the modified string
    PCA_per_class$sample_id[i] <- modified_string
    
    
    # processing the pmd
    if (nchar(PCA_per_class$pmd[i])==5) {
      original_string <- PCA_per_class$pmd[i]
      # Create a new string with the third character removed
      pmd_hours <- as.numeric(paste0(substr(original_string, 1, 2)))
      pmd_minutes <-as.numeric(paste0(substr(original_string, 4, 5)))
      modified_value <- (pmd_hours*60) + pmd_minutes
      # Print the modified string
      PCA_per_class$pmd[i] <- modified_value
    }
    else if (nchar(PCA_per_class$pmd[i])==4){
      original_string <- PCA_per_class$pmd[i]
      # Create a new string without the ':' character
      pmd_hours <- as.numeric(paste0(substr(original_string, 1,1)))
      pmd_minutes <-as.numeric(paste0(substr(original_string, 3, 4)))
      modified_value <- (pmd_hours*60) + pmd_minutes
      # Print the modified string
      PCA_per_class$pmd[i] <- modified_value
    }
  }
  # Second step converting to numeric
  PCA_per_class$`NBB donor ID`= as.numeric(PCA_per_class$`NBB donor ID`)
  PCA_per_class$sample_id = as.numeric(PCA_per_class$sample_id)
  PCA_per_class$pmd = as.numeric(PCA_per_class$pmd)
  PCA_per_class$sex = as.factor(PCA_per_class$sex)
  
  
  #######################
  ### final touches    ##
  #######################
  # changing some columns names:
  names(PCA_per_class)
  names(PCA_per_class)[c(1,2)]=c('sample','donor')
  names(PCA_per_class)[8]='types'
  names(PCA_per_class)
  #######################
  ###   final check    ##
  #######################
  # checking that there are no NAs
  sum(is.na(PCA_per_class))
  
  
  # The final output
  sumscores_output = cbind(PCA_per_class,lipidomics.data.summed)
  sumscores_output_scaled = cbind(PCA_per_class,lipidomics.data.summed.scaled)
  unique(sumscores_output_scaled$types)
  
  # Selecting only the lesion types
  sumscores_output_scaled_lesions_only = sumscores_output_scaled[which(! sumscores_output_scaled$types %in% c('CWM','PLWM','NAWM')) ,]
  
  #write.csv(sumscores_output,"sumscores_output.csv")
  #write.csv(sumscores_output_scaled_lesions_only,"sumscores_output_scaled_lesions_only.csv")
  
  ### PCA analysis ###
  pca <- prcomp(sumscores_output_scaled_lesions_only[,10:ncol(sumscores_output_scaled_lesions_only)],scale=F,center=F)
  pca.var <- pca$sdev^2
  pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
  pca.data <- data.frame(PC1=pca$x[,1],
                         PC2=pca$x[,2])
  
  spss_data <- pca.data
  csv_data <- sumscores_output_scaled_lesions_only
  meta_part = csv_data[,1:10]
  CATPCA_data = cbind(meta_part,spss_data)
  #######################################################################
  #................  Trying multi GEE models ............................
  #######################################################################
  df = CATPCA_data
  names(df)[1]='X'
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
  
  
  probabilities=data.frame(index=1:nrow(df))
  for (i in 1:length(comparisons_9types)) {
    gee <- geeglm(types ~ age +PC1 +PC2 ,
                  data = comparisons_9types[[i]],
                  id = donor, family = binomial("logit"),
                  corstr = "exchangeable")
    
    
    # Predicted probabilities
    predictions <- predict(gee, newdata = df, type = "response")
    probabilities[,names(comparisons_9types)[i]]=predictions
    
  }
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
  # type as the one with maximum probability or with second maximum probability.
  
  # This function returns the two types with the highest probability
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
  
  # Applying leave-one-out cross validation
  for (j in 1:nrow(df)){
    for (i in 1:length(comparisons_9types)) {
      gee <- geeglm(types ~ age +PC1 +PC2 ,
                    data = comparisons_9types[[i]][-j,],
                    id = donor, family = binomial("logit"),
                    corstr = "exchangeable")
      
      
      # Predicted probabilities
      predictions <- predict(gee, newdata = comparisons_9types[[i]][j,], type = "response")
      probabilities[j,names(comparisons_9types)[i]]=predictions
      
    }}
  
  # Deleting the index column
  probabilities=probabilities[,-1]
  
  # Getting the names of the types with the highest probs
  
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
  accuracy #  0.36
  accuracy_with_tolerance_two = mean(predictions_with_tolerance$score_with_tolerance_two)
  accuracy_with_tolerance_three = mean(predictions_with_tolerance$score_with_tolerance_three)
  accuracy_with_tolerance_two
  accuracy_with_tolerance_three
  # the result is that lipids have the potential to classify, but more samples are needed.
  return(c(accuracy,accuracy_with_tolerance_two,accuracy_with_tolerance_three))
}

################################################################################
################################################################################
################################################################################
# impute with the maximum value of each lipid
Sensitivity_analysis_max_impute= function(){
  #loading the metadata + some cleaning up and refining
  metadata <- readxl::read_xlsx("Metadata table2.xlsx") %>% data.frame() 
  row.names(metadata) <- metadata$Unifying_code
  metadata$Morphology.microglia[is.na(metadata$Morphology.microglia)] <- "inactive"
  metadata$Morphology.microglia <- factor(metadata$Morphology.microglia, levels=c("inactive","non_foamy","foamy"))
  Lesion.types.6 <- c("CWM","NAWM","2","3","4","6") 
  metadata$Lesion_type_6 <- factor(metadata$Lesion_type_6, levels=Lesion.types.6)
  Lesion.types.9 <- c("CWM","NAWM","PLWM","2.1_2.2","2.3","3.1_3.2","3.3","4","6") 
  metadata$Lesion_type_9 <- factor(metadata$Lesion_type_9, levels=Lesion.types.9)
  metadata <- metadata %>% #making the pmd a numerical decimal instead of 08:30 character
    separate(pmd,sep=":",into = c("H", "M")) %>%
    mutate(H = as.numeric(H)) %>% 
    mutate(M = as.numeric(M)/60) %>% 
    mutate(pmd = H+M) %>% 
    select(-"H",-"M")
  
  #filtered metadata for lipidomics
  #not all samples have been analyzed by lipidomics, so we filter them 
  metadata.lipids <- metadata %>% 
    filter(sample_lipidomics != "") %>% 
    arrange(sample_lipidomics)
  metadata.ABPP <- metadata %>% 
    filter(sample_ABPP != "") %>% 
    arrange(sample_ABPP)
  
  #Color schemes for nice plots 
  Colors.lesion.6g <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FFF400","#FF7F00")
  names(Colors.lesion.6g) <- Lesion.types.6
  Colors.lesion.9g <- c("#E41A1C","#0072B2","#A6CEE3","#B2DF8A","green4", "#b391c4", "darkmagenta","#FFF400", "#D55E00")
  names(Colors.lesion.9g) <- Lesion.types.9
  Colors.morph <- c("#E1E1E1","#b58709","#0940b5")
  names(Colors.morph) <- c("inactive","non_foamy","foamy")
  
  # loading the categorization of lipid into lipid-classes 
  # only necessary for the heatmaps with lipids in there
  lipidclass <- readxl::read_xlsx("lipidclass MS.xlsx")
  lipidclass$class <- factor(lipidclass$class,levels=c("PC","PE","PS","PI","PG & BMP",
                                                       "LysoPC","LysoPE","LysoPS","LysoPI","LysoPG",
                                                       "NAPE","LysoNAPE","GPNAE","NAE",
                                                       "TAG","DAG","MAG","FFA","Oxylipins",
                                                       "SM","CER","Cerebroside","Sterols","CE"))
  lipidclass$class.broad <- factor(lipidclass$`general class`,
                                   levels=c("Phospholipids",
                                            "Lysophospholipids",
                                            "Sphingolipids",
                                            "NAPEs",
                                            "Sterols",
                                            "CEs",
                                            "TAGs",
                                            "DAGs",
                                            "eCBs",
                                            "FFAs",          
                                            "Oxylipins"))
  
  # Import data frame
  lipidomics.data.notimputed <- readxl::read_xlsx("20221115 Xinyu DATA.xlsx") %>%
    arrange(`Sample code`) %>%
    select(-"Sample code",-"type",-"Protein conc. (mg/mL)") %>% #remove non-lipid variables
    replace(.==0,NA) %>% #zero's are not really 0, they are not measured so NA 
    log2() %>% 
    scale(center=T,scale=F) %>% 
    data.frame(check.names = F)#read in the data file 
  
  #setting the right names, again this is because not all samples have been measured, 
  #so there are multiple names for the same samples. This sets the rownames in line with the metadata
  rownames(lipidomics.data.notimputed) <- metadata %>% 
    filter(sample_lipidomics != "") %>% 
    arrange(sample_lipidomics) %>% 
    row.names()
  
  # CLEANING UP DATA ##
  ## data imputation with sampling from a normal distribution around half of the minimum value of a certain lipid with an SD 1/3 of the original SD
  #this is based on the assumption that if a lipid was not quantified, that it was below the treshold of detection
  #We can not know the treshold, but at least it isLOWER than the lowest value in the dataset
  set.seed(123)
  imputation <- lipidomics.data.notimputed %>% as.matrix() #log2 transformation ensures normal distribution (assumption)
  for (i in 1:ncol(imputation)) { 
    imputation[,i] <- replace(imputation[,i], #select the a column in the data
                              is.na(imputation[,i]), #determine missing values
                              rnorm(   # replace these missing values with: 
                                sum(is.na(imputation[,i])),   #determining how many missing values are needed for a specific lipid (how many values have to be generated?)
                                mean=max(imputation[,i],na.rm=T), #mean is max
                                sd=sd(imputation[,i],na.rm=T)/3)) #SD is original SD divided by 3
  }
  lipidomics.data <- imputation %>%
    scale(center=T, scale=F) %>% #center the data, so the average per lipid is zero. NO SCALING for unit variance
    data.frame(check.names = FALSE)
  
  #write.csv(lipidomics.data,"lipidomics data.csv")
  
  #The following lines are to produce a dataframe that is easier to handle than the original dataset. 
  #To this end, all lipids that belong to the same class are summed. This dataset is quite representative of the real dataset. 
  #however, subtle changes might be missed here. Only for exploratory purposes.
  #it produces a data frame called "results.summed" with only 24 variables (the lipid classes) in the same 100 samples.
  sumdata <- matrix(ncol =length(levels(lipidclass$class)),nrow=100) %>% data.table()
  colnames(sumdata) <- levels(lipidclass$class) 
  
  for (i in levels(lipidclass$class)) { 
    sumdata[,i] <- imputation[,lipidclass$class == i] %>% 
      2^. %>% 
      apply(1,sum,na.rm=T) }
  lipidomics.data.summed <- sumdata %>% log2() %>% scale(center=T,scale=F) %>% data.frame()
  lipidomics.data.summed.scaled <- sumdata %>% log2() %>% scale(center=T,scale=T) %>% data.frame() # here we scale the sumscores
  
  rownames(lipidomics.data.summed) = rownames(lipidomics.data)
  rownames(lipidomics.data.summed.scaled) = rownames(lipidomics.data)
  
  #write.csv(lipidomics.data.summed,"sumdata.csv")
  #write.csv(lipidomics.data.summed.scaled,"sumdata_scaled.csv")
  
  max(lipidomics.data)
  max(sumdata)
  max(lipidomics.data.summed)
  max(lipidomics.data.summed.scaled)
  
  #######################
  ### donors metadata ###
  #######################
  # Load the donor ID
  meta_data <- read_excel("Metadata table2.xlsx")
  # matching each sample with its donor through the variable 'sample_lipidomics':
  #  and then retrieving the corresponding metadata
  donors=data.frame('sample_lipidomics'=c(),'NBB donor ID'=c(),'sex'=c(),'age'=c(),'brain weight'=c(),'pmd'=c(),"pH CSF"=c(),'Lesion_type_9'=c(),"Morphology microglia"=c())
  for (i in 1:nrow(lipidomics.data.summed)){
    sample_id=rownames(lipidomics.data.summed)[i]
    index=which(meta_data$Unifying_code==sample_id)
    current_donor=meta_data[index,c('NBB donor ID','sex','age','brain weight','pmd',"pH CSF",'Lesion_type_9',"Morphology microglia")]
    current_donor=cbind(sample_id=sample_id,current_donor)
    donors=rbind(donors,current_donor)
  }
  
  
  #######################
  ### some missingness ##
  #######################
  # Filling the missing values in the "Morphology microglia" with the 'inactive' label
  sum(is.na(donors$"Morphology microglia"))
  donors$"Morphology microglia" <- replace(donors$"Morphology microglia", is.na(donors$"Morphology microglia"), 'inactive')
  sum(is.na(donors$"Morphology microglia"))
  
  # Filling the 5 missing values in the 'pH CSF' with the mean value
  sum(is.na(donors$"pH CSF"))
  donors$`pH CSF` <- replace(donors$`pH CSF`, is.na(donors$`pH CSF`), mean(donors$`pH CSF`, na.rm = TRUE))
  sum(is.na(donors$"pH CSF"))
  
  # There is one donor missing                         ########################
  length(unique(donors$`NBB donor ID`))                ##### question here ####
  ########################   
  
  # Initiate the final dataframe that will have the components
  #  beginning with the metadata info
  PCA_per_class=donors
  # a variable we will need later on (ease of retrieval)
  types=PCA_per_class$Lesion_type_9
  
  
  
  ##########################################################
  #####  Processing of varibles: donorID,sampleID,pmd ######
  ##########################################################
  dim(PCA_per_class)
  # Changing the sample ID and the donor ID to numeric 
  #  so that the GEE model runs properly.
  # First step taking out the non-numeric characters
  for (i in 1:nrow(PCA_per_class)){
    # processing the donor id
    original_string <- PCA_per_class$`NBB donor ID`[i]
    # Create a new string with the third character removed
    modified_string <- paste0(substr(original_string, 1, 4), 
                              substr(original_string, 6, nchar(original_string)))
    # Print the modified string
    PCA_per_class$`NBB donor ID`[i] <- modified_string
    
    
    # processing the sample code
    original_string <- PCA_per_class$sample_id[i]
    # Create a new string with the third character removed
    modified_string <- paste0(substr(original_string, 6, 8))
    # Print the modified string
    PCA_per_class$sample_id[i] <- modified_string
    
    
    # processing the pmd
    if (nchar(PCA_per_class$pmd[i])==5) {
      original_string <- PCA_per_class$pmd[i]
      # Create a new string with the third character removed
      pmd_hours <- as.numeric(paste0(substr(original_string, 1, 2)))
      pmd_minutes <-as.numeric(paste0(substr(original_string, 4, 5)))
      modified_value <- (pmd_hours*60) + pmd_minutes
      # Print the modified string
      PCA_per_class$pmd[i] <- modified_value
    }
    else if (nchar(PCA_per_class$pmd[i])==4){
      original_string <- PCA_per_class$pmd[i]
      # Create a new string without the ':' character
      pmd_hours <- as.numeric(paste0(substr(original_string, 1,1)))
      pmd_minutes <-as.numeric(paste0(substr(original_string, 3, 4)))
      modified_value <- (pmd_hours*60) + pmd_minutes
      # Print the modified string
      PCA_per_class$pmd[i] <- modified_value
    }
  }
  # Second step converting to numeric
  PCA_per_class$`NBB donor ID`= as.numeric(PCA_per_class$`NBB donor ID`)
  PCA_per_class$sample_id = as.numeric(PCA_per_class$sample_id)
  PCA_per_class$pmd = as.numeric(PCA_per_class$pmd)
  PCA_per_class$sex = as.factor(PCA_per_class$sex)
  
  
  #######################
  ### final touches    ##
  #######################
  # changing some columns names:
  names(PCA_per_class)
  names(PCA_per_class)[c(1,2)]=c('sample','donor')
  names(PCA_per_class)[8]='types'
  names(PCA_per_class)
  #######################
  ###   final check    ##
  #######################
  # checking that there are no NAs
  sum(is.na(PCA_per_class))
  
  
  # The final output
  sumscores_output = cbind(PCA_per_class,lipidomics.data.summed)
  sumscores_output_scaled = cbind(PCA_per_class,lipidomics.data.summed.scaled)
  unique(sumscores_output_scaled$types)
  
  # Selecting only the lesion types
  sumscores_output_scaled_lesions_only = sumscores_output_scaled[which(! sumscores_output_scaled$types %in% c('CWM','PLWM','NAWM')) ,]
  
  #write.csv(sumscores_output,"sumscores_output.csv")
  #write.csv(sumscores_output_scaled_lesions_only,"sumscores_output_scaled_lesions_only.csv")
  
  ### PCA analysis ###
  pca <- prcomp(sumscores_output_scaled_lesions_only[,10:ncol(sumscores_output_scaled_lesions_only)],scale=F,center=F)
  pca.var <- pca$sdev^2
  pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
  pca.data <- data.frame(PC1=pca$x[,1],
                         PC2=pca$x[,2])
  
  spss_data <- pca.data
  csv_data <- sumscores_output_scaled_lesions_only
  meta_part = csv_data[,1:10]
  CATPCA_data = cbind(meta_part,spss_data)
  #######################################################################
  #................  Trying multi GEE models ............................
  #######################################################################
  df = CATPCA_data
  names(df)[1]='X'
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
  
  
  probabilities=data.frame(index=1:nrow(df))
  for (i in 1:length(comparisons_9types)) {
    gee <- geeglm(types ~ age +PC1 +PC2 ,
                  data = comparisons_9types[[i]],
                  id = donor, family = binomial("logit"),
                  corstr = "exchangeable")
    
    
    # Predicted probabilities
    predictions <- predict(gee, newdata = df, type = "response")
    probabilities[,names(comparisons_9types)[i]]=predictions
    
  }
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
  # type as the one with maximum probability or with second maximum probability.
  
  # This function returns the two types with the highest probability
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
  
  # Applying leave-one-out cross validation
  for (j in 1:nrow(df)){
    for (i in 1:length(comparisons_9types)) {
      gee <- geeglm(types ~ age +PC1 +PC2 ,
                    data = comparisons_9types[[i]][-j,],
                    id = donor, family = binomial("logit"),
                    corstr = "exchangeable")
      
      
      # Predicted probabilities
      predictions <- predict(gee, newdata = comparisons_9types[[i]][j,], type = "response")
      probabilities[j,names(comparisons_9types)[i]]=predictions
      
    }}
  
  # Deleting the index column
  probabilities=probabilities[,-1]
  
  # Getting the names of the types with the highest probs
  
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
  accuracy #  0.36
  accuracy_with_tolerance_two = mean(predictions_with_tolerance$score_with_tolerance_two)
  accuracy_with_tolerance_three = mean(predictions_with_tolerance$score_with_tolerance_three)
  accuracy_with_tolerance_two
  accuracy_with_tolerance_three
  # the result is that lipids have the potential to classify, but more samples are needed.
  return(c(accuracy,accuracy_with_tolerance_two,accuracy_with_tolerance_three))
}



################################################################################
################################################################################
################################################################################
################################################################################

# now repeating the analysis for different values of imputations.

Sensitivity_analysis_varied= function(imput_mean){
  #loading the metadata + some cleaning up and refining
  metadata <- readxl::read_xlsx("Metadata table2.xlsx") %>% data.frame() 
  row.names(metadata) <- metadata$Unifying_code
  metadata$Morphology.microglia[is.na(metadata$Morphology.microglia)] <- "inactive"
  metadata$Morphology.microglia <- factor(metadata$Morphology.microglia, levels=c("inactive","non_foamy","foamy"))
  Lesion.types.6 <- c("CWM","NAWM","2","3","4","6") 
  metadata$Lesion_type_6 <- factor(metadata$Lesion_type_6, levels=Lesion.types.6)
  Lesion.types.9 <- c("CWM","NAWM","PLWM","2.1_2.2","2.3","3.1_3.2","3.3","4","6") 
  metadata$Lesion_type_9 <- factor(metadata$Lesion_type_9, levels=Lesion.types.9)
  metadata <- metadata %>% #making the pmd a numerical decimal instead of 08:30 character
    separate(pmd,sep=":",into = c("H", "M")) %>%
    mutate(H = as.numeric(H)) %>% 
    mutate(M = as.numeric(M)/60) %>% 
    mutate(pmd = H+M) %>% 
    select(-"H",-"M")
  
  #filtered metadata for lipidomics
  #not all samples have been analyzed by lipidomics, so we filter them 
  metadata.lipids <- metadata %>% 
    filter(sample_lipidomics != "") %>% 
    arrange(sample_lipidomics)
  metadata.ABPP <- metadata %>% 
    filter(sample_ABPP != "") %>% 
    arrange(sample_ABPP)
  
  #Color schemes for nice plots 
  Colors.lesion.6g <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FFF400","#FF7F00")
  names(Colors.lesion.6g) <- Lesion.types.6
  Colors.lesion.9g <- c("#E41A1C","#0072B2","#A6CEE3","#B2DF8A","green4", "#b391c4", "darkmagenta","#FFF400", "#D55E00")
  names(Colors.lesion.9g) <- Lesion.types.9
  Colors.morph <- c("#E1E1E1","#b58709","#0940b5")
  names(Colors.morph) <- c("inactive","non_foamy","foamy")
  
  # loading the categorization of lipid into lipid-classes 
  # only necessary for the heatmaps with lipids in there
  lipidclass <- readxl::read_xlsx("lipidclass MS.xlsx")
  lipidclass$class <- factor(lipidclass$class,levels=c("PC","PE","PS","PI","PG & BMP",
                                                       "LysoPC","LysoPE","LysoPS","LysoPI","LysoPG",
                                                       "NAPE","LysoNAPE","GPNAE","NAE",
                                                       "TAG","DAG","MAG","FFA","Oxylipins",
                                                       "SM","CER","Cerebroside","Sterols","CE"))
  lipidclass$class.broad <- factor(lipidclass$`general class`,
                                   levels=c("Phospholipids",
                                            "Lysophospholipids",
                                            "Sphingolipids",
                                            "NAPEs",
                                            "Sterols",
                                            "CEs",
                                            "TAGs",
                                            "DAGs",
                                            "eCBs",
                                            "FFAs",          
                                            "Oxylipins"))
  
  # Import data frame
  lipidomics.data.notimputed <- readxl::read_xlsx("20221115 Xinyu DATA.xlsx") %>%
    arrange(`Sample code`) %>%
    select(-"Sample code",-"type",-"Protein conc. (mg/mL)") %>% #remove non-lipid variables
    replace(.==0,NA) %>% #zero's are not really 0, they are not measured so NA 
    log2() %>% 
    scale(center=T,scale=F) %>% 
    data.frame(check.names = F)#read in the data file 
  
  #setting the right names, again this is because not all samples have been measured, 
  #so there are multiple names for the same samples. This sets the rownames in line with the metadata
  rownames(lipidomics.data.notimputed) <- metadata %>% 
    filter(sample_lipidomics != "") %>% 
    arrange(sample_lipidomics) %>% 
    row.names()
  
  # CLEANING UP DATA ##
  ## data imputation with sampling from a normal distribution around half of the minimum value of a certain lipid with an SD 1/3 of the original SD
  #this is based on the assumption that if a lipid was not quantified, that it was below the treshold of detection
  #We can not know the treshold, but at least it isLOWER than the lowest value in the dataset
  set.seed(123)
  imputation <- lipidomics.data.notimputed %>% as.matrix() #log2 transformation ensures normal distribution (assumption)
  for (i in 1:ncol(imputation)) { 
    imputation[,i] <- replace(imputation[,i], #select the a column in the data
                              is.na(imputation[,i]), #determine missing values
                              rnorm(   # replace these missing values with: 
                                sum(is.na(imputation[,i])),   #determining how many missing values are needed for a specific lipid (how many values have to be generated?)
                                mean=imput_mean, #mean of the normal dist is a function parameter
                                sd=sd(imputation[,i],na.rm=T)/3)) #SD is original SD divided by 3
  }
  lipidomics.data <- imputation %>%
    scale(center=T, scale=F) %>% #center the data, so the average per lipid is zero. NO SCALING for unit variance
    data.frame(check.names = FALSE)
  
  #write.csv(lipidomics.data,"lipidomics data.csv")
  
  #The following lines are to produce a dataframe that is easier to handle than the original dataset. 
  #To this end, all lipids that belong to the same class are summed. This dataset is quite representative of the real dataset. 
  #however, subtle changes might be missed here. Only for exploratory purposes.
  #it produces a data frame called "results.summed" with only 24 variables (the lipid classes) in the same 100 samples.
  sumdata <- matrix(ncol =length(levels(lipidclass$class)),nrow=100) %>% data.table()
  colnames(sumdata) <- levels(lipidclass$class) 
  
  for (i in levels(lipidclass$class)) { 
    sumdata[,i] <- imputation[,lipidclass$class == i] %>% 
      2^. %>% 
      apply(1,sum,na.rm=T) }
  lipidomics.data.summed <- sumdata %>% log2() %>% scale(center=T,scale=F) %>% data.frame()
  lipidomics.data.summed.scaled <- sumdata %>% log2() %>% scale(center=T,scale=T) %>% data.frame() # here we scale the sumscores
  
  rownames(lipidomics.data.summed) = rownames(lipidomics.data)
  rownames(lipidomics.data.summed.scaled) = rownames(lipidomics.data)
  
  #write.csv(lipidomics.data.summed,"sumdata.csv")
  #write.csv(lipidomics.data.summed.scaled,"sumdata_scaled.csv")
  
  max(lipidomics.data)
  max(sumdata)
  max(lipidomics.data.summed)
  max(lipidomics.data.summed.scaled)
  
  #######################
  ### donors metadata ###
  #######################
  # Load the donor ID
  meta_data <- read_excel("Metadata table2.xlsx")
  # matching each sample with its donor through the variable 'sample_lipidomics':
  #  and then retrieving the corresponding metadata
  donors=data.frame('sample_lipidomics'=c(),'NBB donor ID'=c(),'sex'=c(),'age'=c(),'brain weight'=c(),'pmd'=c(),"pH CSF"=c(),'Lesion_type_9'=c(),"Morphology microglia"=c())
  for (i in 1:nrow(lipidomics.data.summed)){
    sample_id=rownames(lipidomics.data.summed)[i]
    index=which(meta_data$Unifying_code==sample_id)
    current_donor=meta_data[index,c('NBB donor ID','sex','age','brain weight','pmd',"pH CSF",'Lesion_type_9',"Morphology microglia")]
    current_donor=cbind(sample_id=sample_id,current_donor)
    donors=rbind(donors,current_donor)
  }
  
  
  #######################
  ### some missingness ##
  #######################
  # Filling the missing values in the "Morphology microglia" with the 'inactive' label
  sum(is.na(donors$"Morphology microglia"))
  donors$"Morphology microglia" <- replace(donors$"Morphology microglia", is.na(donors$"Morphology microglia"), 'inactive')
  sum(is.na(donors$"Morphology microglia"))
  
  # Filling the 5 missing values in the 'pH CSF' with the mean value
  sum(is.na(donors$"pH CSF"))
  donors$`pH CSF` <- replace(donors$`pH CSF`, is.na(donors$`pH CSF`), mean(donors$`pH CSF`, na.rm = TRUE))
  sum(is.na(donors$"pH CSF"))
  
  # There is one donor missing                         ########################
  length(unique(donors$`NBB donor ID`))                ##### question here ####
  ########################   
  
  # Initiate the final dataframe that will have the components
  #  beginning with the metadata info
  PCA_per_class=donors
  # a variable we will need later on (ease of retrieval)
  types=PCA_per_class$Lesion_type_9
  
  
  
  ##########################################################
  #####  Processing of varibles: donorID,sampleID,pmd ######
  ##########################################################
  dim(PCA_per_class)
  # Changing the sample ID and the donor ID to numeric 
  #  so that the GEE model runs properly.
  # First step taking out the non-numeric characters
  for (i in 1:nrow(PCA_per_class)){
    # processing the donor id
    original_string <- PCA_per_class$`NBB donor ID`[i]
    # Create a new string with the third character removed
    modified_string <- paste0(substr(original_string, 1, 4), 
                              substr(original_string, 6, nchar(original_string)))
    # Print the modified string
    PCA_per_class$`NBB donor ID`[i] <- modified_string
    
    
    # processing the sample code
    original_string <- PCA_per_class$sample_id[i]
    # Create a new string with the third character removed
    modified_string <- paste0(substr(original_string, 6, 8))
    # Print the modified string
    PCA_per_class$sample_id[i] <- modified_string
    
    
    # processing the pmd
    if (nchar(PCA_per_class$pmd[i])==5) {
      original_string <- PCA_per_class$pmd[i]
      # Create a new string with the third character removed
      pmd_hours <- as.numeric(paste0(substr(original_string, 1, 2)))
      pmd_minutes <-as.numeric(paste0(substr(original_string, 4, 5)))
      modified_value <- (pmd_hours*60) + pmd_minutes
      # Print the modified string
      PCA_per_class$pmd[i] <- modified_value
    }
    else if (nchar(PCA_per_class$pmd[i])==4){
      original_string <- PCA_per_class$pmd[i]
      # Create a new string without the ':' character
      pmd_hours <- as.numeric(paste0(substr(original_string, 1,1)))
      pmd_minutes <-as.numeric(paste0(substr(original_string, 3, 4)))
      modified_value <- (pmd_hours*60) + pmd_minutes
      # Print the modified string
      PCA_per_class$pmd[i] <- modified_value
    }
  }
  # Second step converting to numeric
  PCA_per_class$`NBB donor ID`= as.numeric(PCA_per_class$`NBB donor ID`)
  PCA_per_class$sample_id = as.numeric(PCA_per_class$sample_id)
  PCA_per_class$pmd = as.numeric(PCA_per_class$pmd)
  PCA_per_class$sex = as.factor(PCA_per_class$sex)
  
  
  #######################
  ### final touches    ##
  #######################
  # changing some columns names:
  names(PCA_per_class)
  names(PCA_per_class)[c(1,2)]=c('sample','donor')
  names(PCA_per_class)[8]='types'
  names(PCA_per_class)
  #######################
  ###   final check    ##
  #######################
  # checking that there are no NAs
  sum(is.na(PCA_per_class))
  
  
  # The final output
  sumscores_output = cbind(PCA_per_class,lipidomics.data.summed)
  sumscores_output_scaled = cbind(PCA_per_class,lipidomics.data.summed.scaled)
  unique(sumscores_output_scaled$types)
  
  # Selecting only the lesion types
  sumscores_output_scaled_lesions_only = sumscores_output_scaled[which(! sumscores_output_scaled$types %in% c('CWM','PLWM','NAWM')) ,]
  
  #write.csv(sumscores_output,"sumscores_output.csv")
  #write.csv(sumscores_output_scaled_lesions_only,"sumscores_output_scaled_lesions_only.csv")
  
  ### PCA analysis ###
  pca <- prcomp(sumscores_output_scaled_lesions_only[,10:ncol(sumscores_output_scaled_lesions_only)],scale=F,center=F)
  pca.var <- pca$sdev^2
  pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
  pca.data <- data.frame(PC1=pca$x[,1],
                         PC2=pca$x[,2])
  
  spss_data <- pca.data
  csv_data <- sumscores_output_scaled_lesions_only
  meta_part = csv_data[,1:10]
  CATPCA_data = cbind(meta_part,spss_data)
  #######################################################################
  #................  Trying multi GEE models ............................
  #######################################################################
  df = CATPCA_data
  names(df)[1]='X'
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
  
  
  probabilities=data.frame(index=1:nrow(df))
  for (i in 1:length(comparisons_9types)) {
    gee <- geeglm(types ~ age +PC1 +PC2 ,
                  data = comparisons_9types[[i]],
                  id = donor, family = binomial("logit"),
                  corstr = "exchangeable")
    
    
    # Predicted probabilities
    predictions <- predict(gee, newdata = df, type = "response")
    probabilities[,names(comparisons_9types)[i]]=predictions
    
  }
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
  # type as the one with maximum probability or with second maximum probability.
  
  # This function returns the two types with the highest probability
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
  
  # Applying leave-one-out cross validation
  for (j in 1:nrow(df)){
    for (i in 1:length(comparisons_9types)) {
      gee <- geeglm(types ~ age +PC1 +PC2 ,
                    data = comparisons_9types[[i]][-j,],
                    id = donor, family = binomial("logit"),
                    corstr = "exchangeable")
      
      
      # Predicted probabilities
      predictions <- predict(gee, newdata = comparisons_9types[[i]][j,], type = "response")
      probabilities[j,names(comparisons_9types)[i]]=predictions
      
    }}
  
  # Deleting the index column
  probabilities=probabilities[,-1]
  
  # Getting the names of the types with the highest probs
  
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
  accuracy #  0.36
  accuracy_with_tolerance_two = mean(predictions_with_tolerance$score_with_tolerance_two)
  accuracy_with_tolerance_three = mean(predictions_with_tolerance$score_with_tolerance_three)
  accuracy_with_tolerance_two
  accuracy_with_tolerance_three
  # the result is that lipids have the potential to classify, but more samples are needed.
  return(c(accuracy,accuracy_with_tolerance_two,accuracy_with_tolerance_three))
}



# there is no change in accuracy
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
# with more components
Sensitivity_analysis_3comp= function(){
  #loading the metadata + some cleaning up and refining
  metadata <- readxl::read_xlsx("Metadata table2.xlsx") %>% data.frame() 
  row.names(metadata) <- metadata$Unifying_code
  metadata$Morphology.microglia[is.na(metadata$Morphology.microglia)] <- "inactive"
  metadata$Morphology.microglia <- factor(metadata$Morphology.microglia, levels=c("inactive","non_foamy","foamy"))
  Lesion.types.6 <- c("CWM","NAWM","2","3","4","6") 
  metadata$Lesion_type_6 <- factor(metadata$Lesion_type_6, levels=Lesion.types.6)
  Lesion.types.9 <- c("CWM","NAWM","PLWM","2.1_2.2","2.3","3.1_3.2","3.3","4","6") 
  metadata$Lesion_type_9 <- factor(metadata$Lesion_type_9, levels=Lesion.types.9)
  metadata <- metadata %>% #making the pmd a numerical decimal instead of 08:30 character
    separate(pmd,sep=":",into = c("H", "M")) %>%
    mutate(H = as.numeric(H)) %>% 
    mutate(M = as.numeric(M)/60) %>% 
    mutate(pmd = H+M) %>% 
    select(-"H",-"M")
  
  #filtered metadata for lipidomics
  #not all samples have been analyzed by lipidomics, so we filter them 
  metadata.lipids <- metadata %>% 
    filter(sample_lipidomics != "") %>% 
    arrange(sample_lipidomics)
  metadata.ABPP <- metadata %>% 
    filter(sample_ABPP != "") %>% 
    arrange(sample_ABPP)
  
  #Color schemes for nice plots 
  Colors.lesion.6g <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FFF400","#FF7F00")
  names(Colors.lesion.6g) <- Lesion.types.6
  Colors.lesion.9g <- c("#E41A1C","#0072B2","#A6CEE3","#B2DF8A","green4", "#b391c4", "darkmagenta","#FFF400", "#D55E00")
  names(Colors.lesion.9g) <- Lesion.types.9
  Colors.morph <- c("#E1E1E1","#b58709","#0940b5")
  names(Colors.morph) <- c("inactive","non_foamy","foamy")
  
  # loading the categorization of lipid into lipid-classes 
  # only necessary for the heatmaps with lipids in there
  lipidclass <- readxl::read_xlsx("lipidclass MS.xlsx")
  lipidclass$class <- factor(lipidclass$class,levels=c("PC","PE","PS","PI","PG & BMP",
                                                       "LysoPC","LysoPE","LysoPS","LysoPI","LysoPG",
                                                       "NAPE","LysoNAPE","GPNAE","NAE",
                                                       "TAG","DAG","MAG","FFA","Oxylipins",
                                                       "SM","CER","Cerebroside","Sterols","CE"))
  lipidclass$class.broad <- factor(lipidclass$`general class`,
                                   levels=c("Phospholipids",
                                            "Lysophospholipids",
                                            "Sphingolipids",
                                            "NAPEs",
                                            "Sterols",
                                            "CEs",
                                            "TAGs",
                                            "DAGs",
                                            "eCBs",
                                            "FFAs",          
                                            "Oxylipins"))
  
  # Import data frame
  lipidomics.data.notimputed <- readxl::read_xlsx("20221115 Xinyu DATA.xlsx") %>%
    arrange(`Sample code`) %>%
    select(-"Sample code",-"type",-"Protein conc. (mg/mL)") %>% #remove non-lipid variables
    replace(.==0,NA) %>% #zero's are not really 0, they are not measured so NA 
    log2() %>% 
    scale(center=T,scale=F) %>% 
    data.frame(check.names = F)#read in the data file 
  
  #setting the right names, again this is because not all samples have been measured, 
  #so there are multiple names for the same samples. This sets the rownames in line with the metadata
  rownames(lipidomics.data.notimputed) <- metadata %>% 
    filter(sample_lipidomics != "") %>% 
    arrange(sample_lipidomics) %>% 
    row.names()
  
  # CLEANING UP DATA ##
  ## data imputation with sampling from a normal distribution around half of the minimum value of a certain lipid with an SD 1/3 of the original SD
  #this is based on the assumption that if a lipid was not quantified, that it was below the treshold of detection
  #We can not know the treshold, but at least it isLOWER than the lowest value in the dataset
  set.seed(123)
  imputation <- lipidomics.data.notimputed %>% as.matrix() #log2 transformation ensures normal distribution (assumption)
  for (i in 1:ncol(imputation)) { 
    imputation[,i] <- replace(imputation[,i], #select the a column in the data
                              is.na(imputation[,i]), #determine missing values
                              rnorm(   # replace these missing values with: 
                                sum(is.na(imputation[,i])),   #determining how many missing values are needed for a specific lipid (how many values have to be generated?)
                                mean=min(imputation[,i],na.rm=T)-1, #mean of the normal dist is half of the minimum value of a certain lipid (log2 scale so minus 1)
                                sd=sd(imputation[,i],na.rm=T)/3)) #SD is original SD divided by 3
  }
  lipidomics.data <- imputation %>%
    scale(center=T, scale=F) %>% #center the data, so the average per lipid is zero. NO SCALING for unit variance
    data.frame(check.names = FALSE)
  
  #write.csv(lipidomics.data,"lipidomics data.csv")
  
  #The following lines are to produce a dataframe that is easier to handle than the original dataset. 
  #To this end, all lipids that belong to the same class are summed. This dataset is quite representative of the real dataset. 
  #however, subtle changes might be missed here. Only for exploratory purposes.
  #it produces a data frame called "results.summed" with only 24 variables (the lipid classes) in the same 100 samples.
  sumdata <- matrix(ncol =length(levels(lipidclass$class)),nrow=100) %>% data.table()
  colnames(sumdata) <- levels(lipidclass$class) 
  
  for (i in levels(lipidclass$class)) { 
    sumdata[,i] <- imputation[,lipidclass$class == i] %>% 
      2^. %>% 
      apply(1,sum,na.rm=T) }
  lipidomics.data.summed <- sumdata %>% log2() %>% scale(center=T,scale=F) %>% data.frame()
  lipidomics.data.summed.scaled <- sumdata %>% log2() %>% scale(center=T,scale=T) %>% data.frame() # here we scale the sumscores
  
  rownames(lipidomics.data.summed) = rownames(lipidomics.data)
  rownames(lipidomics.data.summed.scaled) = rownames(lipidomics.data)
  
  #write.csv(lipidomics.data.summed,"sumdata.csv")
  #write.csv(lipidomics.data.summed.scaled,"sumdata_scaled.csv")
  
  max(lipidomics.data)
  max(sumdata)
  max(lipidomics.data.summed)
  max(lipidomics.data.summed.scaled)
  
  #######################
  ### donors metadata ###
  #######################
  # Load the donor ID
  meta_data <- read_excel("Metadata table2.xlsx")
  # matching each sample with its donor through the variable 'sample_lipidomics':
  #  and then retrieving the corresponding metadata
  donors=data.frame('sample_lipidomics'=c(),'NBB donor ID'=c(),'sex'=c(),'age'=c(),'brain weight'=c(),'pmd'=c(),"pH CSF"=c(),'Lesion_type_9'=c(),"Morphology microglia"=c())
  for (i in 1:nrow(lipidomics.data.summed)){
    sample_id=rownames(lipidomics.data.summed)[i]
    index=which(meta_data$Unifying_code==sample_id)
    current_donor=meta_data[index,c('NBB donor ID','sex','age','brain weight','pmd',"pH CSF",'Lesion_type_9',"Morphology microglia")]
    current_donor=cbind(sample_id=sample_id,current_donor)
    donors=rbind(donors,current_donor)
  }
  
  
  #######################
  ### some missingness ##
  #######################
  # Filling the missing values in the "Morphology microglia" with the 'inactive' label
  sum(is.na(donors$"Morphology microglia"))
  donors$"Morphology microglia" <- replace(donors$"Morphology microglia", is.na(donors$"Morphology microglia"), 'inactive')
  sum(is.na(donors$"Morphology microglia"))
  
  # Filling the 5 missing values in the 'pH CSF' with the mean value
  sum(is.na(donors$"pH CSF"))
  donors$`pH CSF` <- replace(donors$`pH CSF`, is.na(donors$`pH CSF`), mean(donors$`pH CSF`, na.rm = TRUE))
  sum(is.na(donors$"pH CSF"))
  
  # There is one donor missing                         ########################
  length(unique(donors$`NBB donor ID`))                ##### question here ####
  ########################   
  
  # Initiate the final dataframe that will have the components
  #  beginning with the metadata info
  PCA_per_class=donors
  # a variable we will need later on (ease of retrieval)
  types=PCA_per_class$Lesion_type_9
  
  
  
  ##########################################################
  #####  Processing of varibles: donorID,sampleID,pmd ######
  ##########################################################
  dim(PCA_per_class)
  # Changing the sample ID and the donor ID to numeric 
  #  so that the GEE model runs properly.
  # First step taking out the non-numeric characters
  for (i in 1:nrow(PCA_per_class)){
    # processing the donor id
    original_string <- PCA_per_class$`NBB donor ID`[i]
    # Create a new string with the third character removed
    modified_string <- paste0(substr(original_string, 1, 4), 
                              substr(original_string, 6, nchar(original_string)))
    # Print the modified string
    PCA_per_class$`NBB donor ID`[i] <- modified_string
    
    
    # processing the sample code
    original_string <- PCA_per_class$sample_id[i]
    # Create a new string with the third character removed
    modified_string <- paste0(substr(original_string, 6, 8))
    # Print the modified string
    PCA_per_class$sample_id[i] <- modified_string
    
    
    # processing the pmd
    if (nchar(PCA_per_class$pmd[i])==5) {
      original_string <- PCA_per_class$pmd[i]
      # Create a new string with the third character removed
      pmd_hours <- as.numeric(paste0(substr(original_string, 1, 2)))
      pmd_minutes <-as.numeric(paste0(substr(original_string, 4, 5)))
      modified_value <- (pmd_hours*60) + pmd_minutes
      # Print the modified string
      PCA_per_class$pmd[i] <- modified_value
    }
    else if (nchar(PCA_per_class$pmd[i])==4){
      original_string <- PCA_per_class$pmd[i]
      # Create a new string without the ':' character
      pmd_hours <- as.numeric(paste0(substr(original_string, 1,1)))
      pmd_minutes <-as.numeric(paste0(substr(original_string, 3, 4)))
      modified_value <- (pmd_hours*60) + pmd_minutes
      # Print the modified string
      PCA_per_class$pmd[i] <- modified_value
    }
  }
  # Second step converting to numeric
  PCA_per_class$`NBB donor ID`= as.numeric(PCA_per_class$`NBB donor ID`)
  PCA_per_class$sample_id = as.numeric(PCA_per_class$sample_id)
  PCA_per_class$pmd = as.numeric(PCA_per_class$pmd)
  PCA_per_class$sex = as.factor(PCA_per_class$sex)
  
  
  #######################
  ### final touches    ##
  #######################
  # changing some columns names:
  names(PCA_per_class)
  names(PCA_per_class)[c(1,2)]=c('sample','donor')
  names(PCA_per_class)[8]='types'
  names(PCA_per_class)
  #######################
  ###   final check    ##
  #######################
  # checking that there are no NAs
  sum(is.na(PCA_per_class))
  
  
  # The final output
  sumscores_output = cbind(PCA_per_class,lipidomics.data.summed)
  sumscores_output_scaled = cbind(PCA_per_class,lipidomics.data.summed.scaled)
  unique(sumscores_output_scaled$types)
  
  # Selecting only the lesion types
  sumscores_output_scaled_lesions_only = sumscores_output_scaled[which(! sumscores_output_scaled$types %in% c('CWM','PLWM','NAWM')) ,]
  
  #write.csv(sumscores_output,"sumscores_output.csv")
  #write.csv(sumscores_output_scaled_lesions_only,"sumscores_output_scaled_lesions_only.csv")
  
  ### PCA analysis ###
  pca <- prcomp(sumscores_output_scaled_lesions_only[,10:ncol(sumscores_output_scaled_lesions_only)],scale=F,center=F)
  pca.var <- pca$sdev^2
  pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
  pca.data <- data.frame(PC1=pca$x[,1],
                         PC2=pca$x[,2],
                         PC3=pca$x[,3])
  
  spss_data <- pca.data
  csv_data <- sumscores_output_scaled_lesions_only
  meta_part = csv_data[,1:10]
  CATPCA_data = cbind(meta_part,spss_data)
  #######################################################################
  #................  Trying multi GEE models ............................
  #######################################################################
  df = CATPCA_data
  names(df)[1]='X'
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
  
  
  probabilities=data.frame(index=1:nrow(df))
  for (i in 1:length(comparisons_9types)) {
    gee <- geeglm(types ~ age +PC1 +PC2+PC3 ,
                  data = comparisons_9types[[i]],
                  id = donor, family = binomial("logit"),
                  corstr = "exchangeable")
    
    
    # Predicted probabilities
    predictions <- predict(gee, newdata = df, type = "response")
    probabilities[,names(comparisons_9types)[i]]=predictions
    
  }
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
  # type as the one with maximum probability or with second maximum probability.
  
  # This function returns the two types with the highest probability
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
  
  # Applying leave-one-out cross validation
  for (j in 1:nrow(df)){
    for (i in 1:length(comparisons_9types)) {
      gee <- geeglm(types ~ age +PC1 +PC2 +PC3,
                    data = comparisons_9types[[i]][-j,],
                    id = donor, family = binomial("logit"),
                    corstr = "exchangeable")
      
      
      # Predicted probabilities
      predictions <- predict(gee, newdata = comparisons_9types[[i]][j,], type = "response")
      probabilities[j,names(comparisons_9types)[i]]=predictions
      
    }}
  
  # Deleting the index column
  probabilities=probabilities[,-1]
  
  # Getting the names of the types with the highest probs
  
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
  accuracy #  0.36
  accuracy_with_tolerance_two = mean(predictions_with_tolerance$score_with_tolerance_two)
  accuracy_with_tolerance_three = mean(predictions_with_tolerance$score_with_tolerance_three)
  accuracy_with_tolerance_two
  accuracy_with_tolerance_three
  # the result is that lipids have the potential to classify, but more samples are needed.
  return(c(accuracy,accuracy_with_tolerance_two,accuracy_with_tolerance_three))
}


# the predictibility got worse with 3 components.
# let's try it now with one component

########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
# with one component
Sensitivity_analysis_1comp= function(){
  #loading the metadata + some cleaning up and refining
  metadata <- readxl::read_xlsx("Metadata table2.xlsx") %>% data.frame() 
  row.names(metadata) <- metadata$Unifying_code
  metadata$Morphology.microglia[is.na(metadata$Morphology.microglia)] <- "inactive"
  metadata$Morphology.microglia <- factor(metadata$Morphology.microglia, levels=c("inactive","non_foamy","foamy"))
  Lesion.types.6 <- c("CWM","NAWM","2","3","4","6") 
  metadata$Lesion_type_6 <- factor(metadata$Lesion_type_6, levels=Lesion.types.6)
  Lesion.types.9 <- c("CWM","NAWM","PLWM","2.1_2.2","2.3","3.1_3.2","3.3","4","6") 
  metadata$Lesion_type_9 <- factor(metadata$Lesion_type_9, levels=Lesion.types.9)
  metadata <- metadata %>% #making the pmd a numerical decimal instead of 08:30 character
    separate(pmd,sep=":",into = c("H", "M")) %>%
    mutate(H = as.numeric(H)) %>% 
    mutate(M = as.numeric(M)/60) %>% 
    mutate(pmd = H+M) %>% 
    select(-"H",-"M")
  
  #filtered metadata for lipidomics
  #not all samples have been analyzed by lipidomics, so we filter them 
  metadata.lipids <- metadata %>% 
    filter(sample_lipidomics != "") %>% 
    arrange(sample_lipidomics)
  metadata.ABPP <- metadata %>% 
    filter(sample_ABPP != "") %>% 
    arrange(sample_ABPP)
  
  #Color schemes for nice plots 
  Colors.lesion.6g <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FFF400","#FF7F00")
  names(Colors.lesion.6g) <- Lesion.types.6
  Colors.lesion.9g <- c("#E41A1C","#0072B2","#A6CEE3","#B2DF8A","green4", "#b391c4", "darkmagenta","#FFF400", "#D55E00")
  names(Colors.lesion.9g) <- Lesion.types.9
  Colors.morph <- c("#E1E1E1","#b58709","#0940b5")
  names(Colors.morph) <- c("inactive","non_foamy","foamy")
  
  # loading the categorization of lipid into lipid-classes 
  # only necessary for the heatmaps with lipids in there
  lipidclass <- readxl::read_xlsx("lipidclass MS.xlsx")
  lipidclass$class <- factor(lipidclass$class,levels=c("PC","PE","PS","PI","PG & BMP",
                                                       "LysoPC","LysoPE","LysoPS","LysoPI","LysoPG",
                                                       "NAPE","LysoNAPE","GPNAE","NAE",
                                                       "TAG","DAG","MAG","FFA","Oxylipins",
                                                       "SM","CER","Cerebroside","Sterols","CE"))
  lipidclass$class.broad <- factor(lipidclass$`general class`,
                                   levels=c("Phospholipids",
                                            "Lysophospholipids",
                                            "Sphingolipids",
                                            "NAPEs",
                                            "Sterols",
                                            "CEs",
                                            "TAGs",
                                            "DAGs",
                                            "eCBs",
                                            "FFAs",          
                                            "Oxylipins"))
  
  # Import data frame
  lipidomics.data.notimputed <- readxl::read_xlsx("20221115 Xinyu DATA.xlsx") %>%
    arrange(`Sample code`) %>%
    select(-"Sample code",-"type",-"Protein conc. (mg/mL)") %>% #remove non-lipid variables
    replace(.==0,NA) %>% #zero's are not really 0, they are not measured so NA 
    log2() %>% 
    scale(center=T,scale=F) %>% 
    data.frame(check.names = F)#read in the data file 
  
  #setting the right names, again this is because not all samples have been measured, 
  #so there are multiple names for the same samples. This sets the rownames in line with the metadata
  rownames(lipidomics.data.notimputed) <- metadata %>% 
    filter(sample_lipidomics != "") %>% 
    arrange(sample_lipidomics) %>% 
    row.names()
  
  # CLEANING UP DATA ##
  ## data imputation with sampling from a normal distribution around half of the minimum value of a certain lipid with an SD 1/3 of the original SD
  #this is based on the assumption that if a lipid was not quantified, that it was below the treshold of detection
  #We can not know the treshold, but at least it isLOWER than the lowest value in the dataset
  set.seed(123)
  imputation <- lipidomics.data.notimputed %>% as.matrix() #log2 transformation ensures normal distribution (assumption)
  for (i in 1:ncol(imputation)) { 
    imputation[,i] <- replace(imputation[,i], #select the a column in the data
                              is.na(imputation[,i]), #determine missing values
                              rnorm(   # replace these missing values with: 
                                sum(is.na(imputation[,i])),   #determining how many missing values are needed for a specific lipid (how many values have to be generated?)
                                mean=min(imputation[,i],na.rm=T)-1, #mean of the normal dist is half of the minimum value of a certain lipid (log2 scale so minus 1)
                                sd=sd(imputation[,i],na.rm=T)/3)) #SD is original SD divided by 3
  }
  lipidomics.data <- imputation %>%
    scale(center=T, scale=F) %>% #center the data, so the average per lipid is zero. NO SCALING for unit variance
    data.frame(check.names = FALSE)
  
  #write.csv(lipidomics.data,"lipidomics data.csv")
  
  #The following lines are to produce a dataframe that is easier to handle than the original dataset. 
  #To this end, all lipids that belong to the same class are summed. This dataset is quite representative of the real dataset. 
  #however, subtle changes might be missed here. Only for exploratory purposes.
  #it produces a data frame called "results.summed" with only 24 variables (the lipid classes) in the same 100 samples.
  sumdata <- matrix(ncol =length(levels(lipidclass$class)),nrow=100) %>% data.table()
  colnames(sumdata) <- levels(lipidclass$class) 
  
  for (i in levels(lipidclass$class)) { 
    sumdata[,i] <- imputation[,lipidclass$class == i] %>% 
      2^. %>% 
      apply(1,sum,na.rm=T) }
  lipidomics.data.summed <- sumdata %>% log2() %>% scale(center=T,scale=F) %>% data.frame()
  lipidomics.data.summed.scaled <- sumdata %>% log2() %>% scale(center=T,scale=T) %>% data.frame() # here we scale the sumscores
  
  rownames(lipidomics.data.summed) = rownames(lipidomics.data)
  rownames(lipidomics.data.summed.scaled) = rownames(lipidomics.data)
  
  #write.csv(lipidomics.data.summed,"sumdata.csv")
  #write.csv(lipidomics.data.summed.scaled,"sumdata_scaled.csv")
  
  max(lipidomics.data)
  max(sumdata)
  max(lipidomics.data.summed)
  max(lipidomics.data.summed.scaled)
  
  #######################
  ### donors metadata ###
  #######################
  # Load the donor ID
  meta_data <- read_excel("Metadata table2.xlsx")
  # matching each sample with its donor through the variable 'sample_lipidomics':
  #  and then retrieving the corresponding metadata
  donors=data.frame('sample_lipidomics'=c(),'NBB donor ID'=c(),'sex'=c(),'age'=c(),'brain weight'=c(),'pmd'=c(),"pH CSF"=c(),'Lesion_type_9'=c(),"Morphology microglia"=c())
  for (i in 1:nrow(lipidomics.data.summed)){
    sample_id=rownames(lipidomics.data.summed)[i]
    index=which(meta_data$Unifying_code==sample_id)
    current_donor=meta_data[index,c('NBB donor ID','sex','age','brain weight','pmd',"pH CSF",'Lesion_type_9',"Morphology microglia")]
    current_donor=cbind(sample_id=sample_id,current_donor)
    donors=rbind(donors,current_donor)
  }
  
  
  #######################
  ### some missingness ##
  #######################
  # Filling the missing values in the "Morphology microglia" with the 'inactive' label
  sum(is.na(donors$"Morphology microglia"))
  donors$"Morphology microglia" <- replace(donors$"Morphology microglia", is.na(donors$"Morphology microglia"), 'inactive')
  sum(is.na(donors$"Morphology microglia"))
  
  # Filling the 5 missing values in the 'pH CSF' with the mean value
  sum(is.na(donors$"pH CSF"))
  donors$`pH CSF` <- replace(donors$`pH CSF`, is.na(donors$`pH CSF`), mean(donors$`pH CSF`, na.rm = TRUE))
  sum(is.na(donors$"pH CSF"))
  
  # There is one donor missing                         ########################
  length(unique(donors$`NBB donor ID`))                ##### question here ####
  ########################   
  
  # Initiate the final dataframe that will have the components
  #  beginning with the metadata info
  PCA_per_class=donors
  # a variable we will need later on (ease of retrieval)
  types=PCA_per_class$Lesion_type_9
  
  
  
  ##########################################################
  #####  Processing of varibles: donorID,sampleID,pmd ######
  ##########################################################
  dim(PCA_per_class)
  # Changing the sample ID and the donor ID to numeric 
  #  so that the GEE model runs properly.
  # First step taking out the non-numeric characters
  for (i in 1:nrow(PCA_per_class)){
    # processing the donor id
    original_string <- PCA_per_class$`NBB donor ID`[i]
    # Create a new string with the third character removed
    modified_string <- paste0(substr(original_string, 1, 4), 
                              substr(original_string, 6, nchar(original_string)))
    # Print the modified string
    PCA_per_class$`NBB donor ID`[i] <- modified_string
    
    
    # processing the sample code
    original_string <- PCA_per_class$sample_id[i]
    # Create a new string with the third character removed
    modified_string <- paste0(substr(original_string, 6, 8))
    # Print the modified string
    PCA_per_class$sample_id[i] <- modified_string
    
    
    # processing the pmd
    if (nchar(PCA_per_class$pmd[i])==5) {
      original_string <- PCA_per_class$pmd[i]
      # Create a new string with the third character removed
      pmd_hours <- as.numeric(paste0(substr(original_string, 1, 2)))
      pmd_minutes <-as.numeric(paste0(substr(original_string, 4, 5)))
      modified_value <- (pmd_hours*60) + pmd_minutes
      # Print the modified string
      PCA_per_class$pmd[i] <- modified_value
    }
    else if (nchar(PCA_per_class$pmd[i])==4){
      original_string <- PCA_per_class$pmd[i]
      # Create a new string without the ':' character
      pmd_hours <- as.numeric(paste0(substr(original_string, 1,1)))
      pmd_minutes <-as.numeric(paste0(substr(original_string, 3, 4)))
      modified_value <- (pmd_hours*60) + pmd_minutes
      # Print the modified string
      PCA_per_class$pmd[i] <- modified_value
    }
  }
  # Second step converting to numeric
  PCA_per_class$`NBB donor ID`= as.numeric(PCA_per_class$`NBB donor ID`)
  PCA_per_class$sample_id = as.numeric(PCA_per_class$sample_id)
  PCA_per_class$pmd = as.numeric(PCA_per_class$pmd)
  PCA_per_class$sex = as.factor(PCA_per_class$sex)
  
  
  #######################
  ### final touches    ##
  #######################
  # changing some columns names:
  names(PCA_per_class)
  names(PCA_per_class)[c(1,2)]=c('sample','donor')
  names(PCA_per_class)[8]='types'
  names(PCA_per_class)
  #######################
  ###   final check    ##
  #######################
  # checking that there are no NAs
  sum(is.na(PCA_per_class))
  
  
  # The final output
  sumscores_output = cbind(PCA_per_class,lipidomics.data.summed)
  sumscores_output_scaled = cbind(PCA_per_class,lipidomics.data.summed.scaled)
  unique(sumscores_output_scaled$types)
  
  # Selecting only the lesion types
  sumscores_output_scaled_lesions_only = sumscores_output_scaled[which(! sumscores_output_scaled$types %in% c('CWM','PLWM','NAWM')) ,]
  
  #write.csv(sumscores_output,"sumscores_output.csv")
  #write.csv(sumscores_output_scaled_lesions_only,"sumscores_output_scaled_lesions_only.csv")
  
  ### PCA analysis ###
  pca <- prcomp(sumscores_output_scaled_lesions_only[,10:ncol(sumscores_output_scaled_lesions_only)],scale=F,center=F)
  pca.var <- pca$sdev^2
  pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
  pca.data <- data.frame(PC1=pca$x[,1])
  
  spss_data <- pca.data
  csv_data <- sumscores_output_scaled_lesions_only
  meta_part = csv_data[,1:10]
  CATPCA_data = cbind(meta_part,spss_data)
  #######################################################################
  #................  Trying multi GEE models ............................
  #######################################################################
  df = CATPCA_data
  names(df)[1]='X'
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
  
  
  probabilities=data.frame(index=1:nrow(df))
  for (i in 1:length(comparisons_9types)) {
    gee <- geeglm(types ~ age +PC1,
                  data = comparisons_9types[[i]],
                  id = donor, family = binomial("logit"),
                  corstr = "exchangeable")
    
    
    # Predicted probabilities
    predictions <- predict(gee, newdata = df, type = "response")
    probabilities[,names(comparisons_9types)[i]]=predictions
    
  }
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
  # type as the one with maximum probability or with second maximum probability.
  
  # This function returns the two types with the highest probability
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
  
  # Applying leave-one-out cross validation
  for (j in 1:nrow(df)){
    for (i in 1:length(comparisons_9types)) {
      gee <- geeglm(types ~ age +PC1 ,
                    data = comparisons_9types[[i]][-j,],
                    id = donor, family = binomial("logit"),
                    corstr = "exchangeable")
      
      
      # Predicted probabilities
      predictions <- predict(gee, newdata = comparisons_9types[[i]][j,], type = "response")
      probabilities[j,names(comparisons_9types)[i]]=predictions
      
    }}
  
  # Deleting the index column
  probabilities=probabilities[,-1]
  
  # Getting the names of the types with the highest probs
  
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
  accuracy #  0.36
  accuracy_with_tolerance_two = mean(predictions_with_tolerance$score_with_tolerance_two)
  accuracy_with_tolerance_three = mean(predictions_with_tolerance$score_with_tolerance_three)
  accuracy_with_tolerance_two
  accuracy_with_tolerance_three
  # the result is that lipids have the potential to classify, but more samples are needed.
  return(c(accuracy,accuracy_with_tolerance_two,accuracy_with_tolerance_three))
}
Sensitivity_analysis_1comp()
Sensitivity_analysis_3comp()
Sensitivity_analysis_min_impute()
Sensitivity_analysis_max_impute()
# so two components are the best
result = data.frame(mean_imput=1:20)
for (i in 1:20){
  round = Sensitivity_analysis_varied(i)
  result[,'without tolerance'] = round[1]
  result[,'with tolerance two'] = round[2]
  result[,'with tolerance three'] = round[3]
}
result = rbind(result,c('min of each lipid',0.36,0.64,0.88),
                 c('min of each lipid',0.38,0.68,0.84))
