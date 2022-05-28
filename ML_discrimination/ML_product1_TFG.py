#!/usr/bin/python


import argparse
import shutil
import glob
import json
import sys
import re
import os
import pandas as pd
import numpy as np
import pickle
from sklearn.model_selection import LeaveOneOut
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import cross_val_score
from sklearn.svm import SVC
from sklearn import svm
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
from numpy import mean
from numpy import absolute
from numpy import sqrt
from numpy import std


def InputPar():
	#parsing arguments
	parser = argparse.ArgumentParser(description="Trimming script")
	parser.add_argument("-m", "--matrix", metavar="PLATFORM", type=str,
		help="Sequencing platform, it can be: 'illumina', 'roche454', 'iontorrent' or 'oxfordnanopore'.", required=True)
	parser.add_argument("-e", "--metadata", metavar="DIRECTORY", type=str, required=True,
		help="Directory where the input .fq/.fq.gz/.fastq/fastq.gz are in.")
	return parser.parse_args()


# Random Forest training

def trainRF(matrix, metadata): # matrix with features selected and preprocessing done; metadada with binary column
  

	# taxons as columns (in order to train the ML model)
	matrixT = matrix.transpose()
	matrixT.rename(columns=matrixT.iloc[0], inplace = True) # convert first row in columns name
	matrixFinal = matrixT.iloc[1: , :] # remove first row (it is a duplicate)
	
	# in metadata select only the binary column
	metadata = metadata.loc[:, metadata.columns.str.startswith('condition1')]

	# train Random Forest by LOOCV
	X = matrixFinal # matrix
	y = metadata # labels
	y = np.ravel(y) # change the labels from a column format to be in a row


	cv = LeaveOneOut()
	model = RandomForestClassifier()

	scores_rf = cross_val_score(model,  X , y , cv = cv, scoring='neg_mean_squared_error',n_jobs = -1) 
	scores_acc_rf = cross_val_score(model,  X , y , cv = cv, scoring='accuracy',n_jobs = -1) 
	
	#view RMSE
	mean_squared_error = sqrt(mean(absolute(scores_rf))) # mean and standard deviation
	# accuracy
	accuracy = (mean(scores_acc_rf))


	# save model in pickle format
	
	# To fit our stimator (call fit in it explicity with provided dataset)
	model.fit(X, np.ravel(y, order = 'C'))

	# To serialize (store the model)
	with open('modelRF.pkl', 'wb') as fid:
			pickle.dump(model, fid, protocol = 2)
			
	
	return(mean_squared_error,accuracy)
	


  

# SVM training

def trainSVM(matrix, metadata):

	# taxons as columns (in order to train the ML model)
	matrixT = matrix.transpose()
	matrixT.rename(columns=matrixT.iloc[0], inplace = True) # convert first row in columns name
	matrixFinal = matrixT.iloc[1: , :] # remove first row (it is a duplicate)
	
	# in metadata select only the binary column
	metadata = metadata.loc[:, metadata.columns.str.startswith('condition1')]

	# train SVM by LOOCV
	X = matrixFinal # matrix
	y = metadata # labels
	y = np.ravel(y) # change the labels from a column format to be in a row

	cv = LeaveOneOut()
	modelSVM = SVC()

	scores_svm = cross_val_score(modelSVM,  X , y , cv = cv, scoring='neg_mean_squared_error',n_jobs = -1) 
	scores_svmAcc = cross_val_score(modelSVM,  X , y , cv = cv, scoring='accuracy',n_jobs = -1) 
	
	# RMSE
	mean_squared_errorSVM = sqrt(mean(absolute(scores_svm))) # mean and standard deviation
	# accuracy
	accuracySVM = (mean(scores_svmAcc))



	# save model in pickle format

	# fit our stimator
	modelSVM = modelSVM.fit(X, np.ravel(y, order = 'C'))

	# To serialize (store the model)
	with open('modelSVM.pkl', 'wb') as fid:
		pickle.dump(modelSVM, fid, protocol = 2)

	 
	
	return(mean_squared_errorSVM,accuracySVM) # return the MSE and accuracy



# Extraction of the variables that we consider important (made with Boruta algorithm)
def featureExtraction(initialMatrix): # initial matrix that we have
	dataSelected = initialMatrix[initialMatrix['Unnamed: 0'].str.contains('Arthrobacter phage KellEzio|Escherichia phage PBECO4|Flexistipes sinusarabici|Lactococcus phage 50101|Salmonella phage vB_SosS_Oslo|Cellulophaga phage phi38:1|Escherichia phage vB_EcoM_Alf5|Geobacillus sp. JS12|Mycobacterium phage Charlie|Yellowstone lake mimivirus')] 

	return(dataSelected)


def preprocessing(data, meta): # prepare input data for the model
	# taxons as columns (in order to train the ML model)
	data = data.transpose()
	data.rename(columns=data.iloc[0], inplace = True) # convert first row in columns name
	dataFinal = data.iloc[1: , :] # remove first row (it is a duplicate)
	
	# in metadata select only the binary column
	meta = meta.loc[:, meta.columns.str.startswith('condition1')]
	return(dataFinal, meta)

def predictionRF(inputMatrix):
	# load the RF model
	with open("modelRF.pkl", "rb") as fid:
		modelRF = pickle.load(fid)

	y_predict = modelRF.predict(inputMatrix)
	return(y_predict)

def predictionSVM(inputMatrix):
	# load the SVM model
	with open("modelSVM.pkl", "rb") as fid:
		modelSVM = pickle.load(fid)

	y_predict = modelSVM.predict(inputMatrix)
	return(y_predict)





def main():
	#take arguments
	
	
	#matrixHuge.csv: samples and taxa
	#metadataHuge.csv: correspondence control or parkinson

	data = pd.read_csv(sys.argv[1]) # "matrixHuge.csv"

	#same for conditions
	meta = pd.read_csv(sys.argv[2]) # "metadataHuge.csv"
	
	# get only the important features in our input matrix
	dataExtracted = featureExtraction(data)

	#merge and transpose
	dataInput, meta = preprocessing(dataExtracted, meta) # this is the correct format to train

	# we have already the correct format to train

	

	# results
	
	resultSVM = trainSVM(dataExtracted, meta)
	MSE_svm = resultSVM[0] # get MSE result from the tuple
	acc_svm = resultSVM[1] # get accuracy result
	

	MSE_svm.astype(np.int64) # transform into a int
	acc_svm.astype(np.int64)
	

	MSE_svm = "{0:.3f}".format(MSE_svm) # show only 3 decimals
	acc_svm = "{0:.3f}".format(acc_svm)

	# train and load the model that is stored in pickle format
	classificationSVM = predictionSVM(dataInput)

	print("Prediction by SVM: " + str(classificationSVM)) # classification prediction

	print('-----------------------------------------')

	#print("Metrics RF: "+ str(trainRF("matrix", "metadata")))
	resultRF = trainRF(dataExtracted, meta) #
	MSE_rf = resultRF[0]
	acc_rf = resultRF[1]
	
	MSE_rf.astype(np.int64) # convert numpy.float64 into int
	acc_rf.astype(np.int64)
	
	MSE_rf = "{0:.3f}".format(MSE_rf) # show only 3 decimals
	acc_rf = "{0:.3f}".format(acc_rf)
	
	# train and load the model that is stored in pickle format
	classificationRF = predictionRF(dataInput)

	print("Prediction by RF: " + str(classificationRF))
	print("-----------------------------------------")

	metaNew = np.ravel(meta)
	print("Correct classification: " + str(metaNew))
	print("-----------------------------------------")



	# Calculate ROC_AUC

	y_true = metaNew
	y_pred_SVM = classificationSVM
	y_pred_RF = classificationRF

	roc_auc_svm = "{0:.3f}".format(roc_auc_score(y_true, y_pred_SVM))
	roc_auc_rf = "{0:.3f}".format(roc_auc_score(y_true, y_pred_RF))

	# Calculate precision
	precision_svm = "{0:.3f}".format(precision_score(y_true, y_pred_SVM))
	precision_rf = "{0:.3f}".format(precision_score(y_true, y_pred_RF))

	# Calculate recall
	recall_svm = "{0:.3f}".format(recall_score(y_true, y_pred_SVM))
	recall_rf = "{0:.3f}".format(recall_score(y_true, y_pred_RF))

	# Calculate f1
	f1_svm = "{0:.3f}".format(f1_score(y_true, y_pred_SVM))
	f1_rf = "{0:.3f}".format(f1_score(y_true, y_pred_RF))

	# print table results
	d = {1: ["Random Forest", MSE_rf, acc_rf, roc_auc_rf, precision_rf, recall_rf, f1_rf],
	2: ["SVM", MSE_svm, acc_svm, roc_auc_svm, precision_svm, recall_svm, f1_svm]}
	print(("{:<14} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10}".format('Algorithm', 'MSE', 'Accuracy', 'ROC_AUC', 'Precision', 'Recall', 'f1')))
	for k, v in d.items():
		Algorithm, MSE, Accuracy, ROC_AUC, Precision, Recall, f1 = v
		print("{:<14} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10}".format(Algorithm, MSE, Accuracy, ROC_AUC, Precision, Recall, f1))

	# print the best algorithm based on the metric ROC_AUC
	if roc_auc_rf >= roc_auc_svm:
		print("The best model for your dataset is: Random Forest, with a ROC_AUC of " + str(roc_auc_rf))
	else:
		print("The best model for your data is: Super Vector Machine, with a ROC_AUC of " + str(roc_auc_svm))


print(main())
