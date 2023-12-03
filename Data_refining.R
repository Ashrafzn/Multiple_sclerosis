
#loading packages
library(tidyverse)
library(data.table)
library(ggplot2)
library(cowplot)
library(patchwork)
library(readxl)
library(dplyr)


setwd('D:/MSc Statistics & Data Science/Second Year/Statistical Consulting/Client meetings/Final data analysis')
################################################################
#........ Loading data, metadata, lipid classes, types .........
################################################################


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
###########################################################
#........ transformation , imputation, sumscores ..........
###########################################################

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


#The following lines are to produce a dataframe that is easier to handle than the original dataset. 
#To this end, all lipids that belong to the same class are summed. This dataset is quite representative of the real dataset. 
#it produces a data frame called "results.summed" with only 24 variables (the lipid classes) in the same 100 samples.
sumdata <- matrix(ncol =length(levels(lipidclass$class)),nrow=100) %>% data.table()
colnames(sumdata) <- levels(lipidclass$class) 

for (i in levels(lipidclass$class)) { 
  sumdata[,i] <- imputation[,lipidclass$class == i] %>% 
    2^. %>% 
    apply(1,sum,na.rm=T) }
lipidomics.data.summed.scaled <- sumdata %>% log2() %>% scale(center=T,scale=T) %>% data.frame() # here we scale the sumscores
rownames(lipidomics.data.summed.scaled) = rownames(lipidomics.data)


#######################
### donors metadata ###
#######################
# Load the donor ID
meta_data <- read_excel("Metadata table2.xlsx")
# matching each sample with its donor through the variable 'sample_lipidomics':
#  and then retrieving the corresponding metadata
donors=data.frame('sample_lipidomics'=c(),'NBB donor ID'=c(),'sex'=c(),'age'=c(),'brain weight'=c(),'pmd'=c(),"pH CSF"=c(),'Lesion_type_9'=c(),"Morphology microglia"=c())
for (i in 1:nrow(lipidomics.data.summed.scaled)){
  sample_id=rownames(lipidomics.data.summed.scaled)[i]
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
sumscores_output_scaled = cbind(PCA_per_class,lipidomics.data.summed.scaled)
unique(sumscores_output_scaled$types)

# Selecting only the lesion types
sumscores_output_scaled_lesions_only = sumscores_output_scaled[which(! sumscores_output_scaled$types %in% c('CWM','PLWM','NAWM')) ,]

write.csv(sumscores_output_scaled_lesions_only,"sumscores_output_scaled_lesions_only.csv")

# Types distribution
table(sumscores_output_scaled_lesions_only$types)

table_donors = table( sumscores_output_scaled_lesions_only$types)

table_donors = as.data.frame(table_donors)
names(table_donors)=c('Lesion type','Frequency')
table_donors

# Age distribution
c(min(sumscores_output_scaled_lesions_only$age),max(sumscores_output_scaled_lesions_only$age))
mean(sumscores_output_scaled_lesions_only$age)
sd(sumscores_output_scaled_lesions_only$age)

# donors and samples
donor_sample = as.data.frame(table(table(sumscores_output_scaled_lesions_only$donor)))
names(donor_sample)= c('Samples taken per donor','number of donor')
donor_sample

# Plot each class separately with different colors for each lesion type
for (col in names(sumscores_output_scaled_lesions_only)[10:33]) {
  plot(sumscores_output_scaled_lesions_only[[col]], col = as.factor(sumscores_output_scaled_lesions_only$types ), 
       main = paste("Class: ", col), pch = 16)
}
