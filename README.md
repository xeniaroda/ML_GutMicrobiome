# ML_GutMicrobiome
This repository contains the code for a taxonomic, pathway abundance and gene family abundance analysis from different metagenomic programs of controls and Parkinson patients and a main script that trains different models from different ML approaches to discriminate between control and diseased patients. It also contains a script that takes as input the microbiome profile of samples and the generated models.

## Name
Analysis of the gut microbiome with Artificial Intelligence techniques for Parkinson's disease.

## Description
Parkinson's disease (PD) is a progressive multisystem neurodegenerative disease affecting people mainly in later years of life. The aetiology of PD is probably multifactorial and there is no available treatment that will stop the progression of the disease. For that reason, discovering important biomarkers for the disease could bring us some ways to mitigate or even halt PD.


**Methods:** We have made a taxonomic analysis with the program Centrifuge and a Pathway and Gene family abundance analysis with HUMAnN3 program. Subsequetly, the three matrices generated for each analysis were trained with different ML algorithms to evaluate the classification accuracy to discrimnate between control or Parkinson patients.

The script ML_product1_TFG.py is based on data mining and classification techniques to detect Parkinson's disease (PD). We developed an automated machine learning solution to detect PD by exploring 10 taxonomic variables selected by the Boruta algorithm (which is a feature selection algorithm). The dataset used for evaluation consist of 40 individuals: 20 stool samples from control subjects and 20 from PD patients.
We use these data to build two different classification algorithms: Support Vector Machine (SVM) and Random Forest (RF). By using this models we will predict the PD and non PD status of the patients.
Accuracy results of different models are compared and from the comparison the best model with the relevant algorithm is shown, thus you will be able to select the best model to predict the PD status of patients. In addition, the program outputs the prediction for the given samples to be PD or control for both algorithms trained, being 0 control and 1 PD.

## Scripts generated

- Taxonomic analysis: centrifuge_TFG.R
- Patway and Gene abundance analysis: humann3_TFG.R
- Principal Coordinate Analysis: PCoA_TFG.R
- Evaluation of Machine Learning Algorithms with LOOCV_1: LOOCV_TFG.R
- Evaluation of Machine Learning Algorithms with LOOCV_2: LOOCV_python_TFG.ipynb
- ML-Discrimination script: ML_product1_TFG.py

## Input/Output files
In the folder named: Input Example, there is an example of input data for our ML-Discrimination script.
These are the Input and the Output files generated from our main script:
- Input 1: A taxonomic matrix in csv format having samples as columns and taxons as rows.
- Input 2: A metadata in csv format having 2 columns corresponding to their samples name and the condition (1 for PD and 0 for control).
- Output: Two models in Pickle format. One corresponding to SVM and the other from RF.

## Visuals
ML-Discrimination script output:

[![asciicast](https://asciinema.org/a/xrGI4AFP4QYUv8FzaFY2quMDN.svg)](https://asciinema.org/a/xrGI4AFP4QYUv8FzaFY2quMDN)

A workflow of the program can be seen in the folder named: workflow.

Workflow of Machine Learning Classifier for PD disease: splitted in two main paths: (a) having as input a taxonomic matrix, (a1) extracting rows containing the important features selected by Boruta algorithm. (a2) A preprocessing step taking as input the matrix with the extracted features and the metadata. (a3) Finally making the prediction step with Random Forest (RF) and Support Vector Machine (SVM) taking as input the taxonomic matrix already preprocessed and (b1) the model in pickle format generated for RF and SVM. If we wanted to classify one sample, we would skip (b) path and perform only (a) path without inputting metadata as we will already have our model saved in Pickle format.

## Version requirements

- R version (4.1.2)
- Python (3.6.9)
- Jupyter Notebook (3.7.33)
- Centrifuge (v1.0.4_beta)
- HUMAnN (3.0.1)

