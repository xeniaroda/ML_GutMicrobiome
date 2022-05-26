# ML_GutMicrobiome
This repository contains the code for a taxonomic, pathway abundance and gene family abundance analysis from different metagenomic programs of controls and Parkinson patients and a main script that trains different models from different ML approaches to discriminate between control and diseased patients. It also contains a script that takes as input the microbiome profile of a sample and the generated models.

# Name
ML scrip classification from the gut microbiome to discriminate between control and Parkinson's patients.

# Description
Parkinson's disease (PD) is a progressive multisystem neurodegenerative disease affecting people mainly in later years of life. The aetiology of PD is probably multifactorial and there is no available treatment that will stop the progression of the disease. For that reason, discovering important biomarkers for the disease could bring us some ways to mitigate or even halt PD.


** my **
The script ML_product1_TFG.py is based on data mining and classification techniques to detect Parkinson's disease (PD). We developed an automated machine learning solution to detect PD by exploring 10 taxonomic variables selected by the Boruta algorithm (which is a feature selection algorithm). The dataset used for evaluation consist of 40 individuals: 20 stool samples from control subjects and 20 from PD patients.
We use these data to build two different classification algorithms: Support Vector Machine (SVM) and Random Forest (RF). By using this models we will predict the PD and non PD status of the patients.
Accuracy results of different models are compared and from the comparison the best model with the relevant algorithm is shown, thus you will be able to select the best model to predict the PD status of patients. In addition, the program outputs the prediction for the given samples to be PD or control for both algorithms trained, being 0 control and 1 PD.
