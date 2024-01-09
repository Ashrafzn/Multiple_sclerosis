# Multiple_sclerosis
Decoding Multiple Sclerosis Lesions: Exploring Lipid Composition as a Classifier in Lesion Typing through Post-Mortem Analysis.
Due to privacy issues, the data is not shared. But all the codes that are used to import, process and analyze the data are available here.
Files:
Data_refining.R : The code that imports the data and metadata, performs the imputations, structure the resulting dataframe as the sumscores of the classes of the lipids.
PCAsum_KNN.Rmd : The code that imports the output sumscores of the previous file (Data_refining.R), and performs PCA on them.
Analysis_multgee.R :  The code that imports the output sumscores of the previous file (Data_refining.R), and performs multinomial GEE analysis on them, with the leave one out cross validation.
Sensitivity_analysis.R: The code that re-does the whole analysis for variety of imputation distributions and monitor the corresponding change in the results.
