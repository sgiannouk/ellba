### ANALYSING LIQUID BIOPSY BASED RNASEQ DATA ###
## ESR11 - Stavros Giannoukakos 
## Script: classification.py

#Version of the program
__version__ = "0.3.8"

# Importing general libraries
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import os
import glob
import gtfast
import argparse
import itertools
import subprocess
import numpy as np
import pandas as pd
from math import isclose
from gtfparse import read_gtf
from datetime import datetime
from itertools import compress
from collections import Counter, OrderedDict
# Importing Machine Learning libraries
from sklearn.pipeline import make_pipeline
from genetic_selection import GeneticSelectionCV
from fast_ml.feature_selection import get_constant_features
from sklearn.preprocessing import RobustScaler, Binarizer, MinMaxScaler, StandardScaler
from sklearn.model_selection import train_test_split, StratifiedKFold, KFold
from sklearn.metrics import balanced_accuracy_score, f1_score, jaccard_score, precision_recall_fscore_support
from sklearn.metrics import RocCurveDisplay, roc_auc_score, accuracy_score, confusion_matrix, auc
from sklearn.ensemble import ExtraTreesClassifier, GradientBoostingClassifier, RandomForestClassifier, AdaBoostClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn import linear_model
# import imblearn
from sklearn.svm import SVC
from numpy import mean, std
import pickle
# Using Matplotlib and Seaborn libraries
import matplotlib
matplotlib.use("Agg") #Needed to save figures
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import seaborn as sns
random_s = 12345


startTime = datetime.now()


### ADDITIONAL FILES
ellba_dir			   = os.path.dirname(os.path.realpath(__file__))
parent_dir			   = os.path.dirname(ellba_dir)
references_dir		   = os.path.join(parent_dir, "references")
rscripts 			   = os.path.join(ellba_dir, "additional_files")
### REFERENCE GENOME and ANNOTATION
referenceannot 	   	   = os.path.join(references_dir, "gencode.v35.primary_assembly.annotation.gtf")





usage       = "classification.py [options]"
epilog      = " -- October 2020 | Stavros Giannoukakos -- "
description = "DESCRIPTION"

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, usage=usage, description=description, epilog=epilog)
# Directory with the all training set matrices
parser.add_argument('-td', '--trainingdir', dest='trainingdir',
                	help="Path of the directory where the\nfiltered training matrices are stored\n(e.g /dataset/data_analysis/filtered_matrices)")
# Directory with the all external validation set matrices
parser.add_argument('-vd', '--validationdir', dest='validationdir',
                	help="Path of the directory where the filtered\nexternal validation matrices are stored\n(e.g /dataset/validationdata_analysis/filtered_matrices)")
# Control group
parser.add_argument('-ctrl', '--control', dest='control', default="NonCancer", 
                	help="In the two-group comparison, which\nlabel should be considered as CONTROL\n(e.g. Healthy)")
# Condition group
parser.add_argument('-cond', '--condition', dest='condition', default="Cancer",
                	help="In the two-group comparison, which\nlabel should be considered as CONDITION\n(e.g. Cancer)")
# Project type
parser.add_argument('-pj', '--project', dest='project', default="MyProject",
                	help="Indicate your project with one word. For example,\nin a cancer diagnosis project, you can use the\ncancer type e.g. NSCLC")
# Reference annotation
parser.add_argument('-ra', '--refannot', dest='refannot', default=referenceannot,
                	help="Path of the reference annotation in gtf format\n(also preferably from GENCODE)")
# Outlier samples
parser.add_argument('-ro', '--removeoutliers', dest='removeoutliers', action="store_true",
                	help="Based on the initial exploratory analysis on the data\nremove detected outlier samples (default %(default)s)")
# Mean training set cutoff
parser.add_argument('-c', '--cutoff', dest='cutoff', default=0.65, type=float,
                	help="Minimum cutoff auc value in order for the\nfeature matrix to be utilised in the majority\nvoting classifier (default %(default)s)")
# Mean training set cutoff
parser.add_argument('-cc', '--crosscut', dest='crosscut', default=0.5, type=float,
                	help="Crossover cutoff threshold for the ensemble classifier (default %(default)s)")
# Number of threads/CPUs to be used
parser.add_argument('-t', '--threads', dest='threads', default=5, type=int, metavar='', 
                	help="Number of threads to be used in the analysis (default %(default)s)")
# Display the version of the pipeline
parser.add_argument('-v', '--version', action='version', version='%(prog)s {0}'.format(__version__))
# Get the options and return them
args = parser.parse_args()





# Main folder hosting the analysis
analysis_dir 	 = os.path.dirname(os.path.realpath(args.trainingdir))
# Classification analysis
ml_dir 			 = os.path.join(analysis_dir, "classification")
ind_dir			 = os.path.join(ml_dir, "individual_features")
voting_dir 		 = os.path.join(ml_dir, "voting_classification")
biorel_dir 	 	 = os.path.join(ml_dir, "biological_relevance")
conditions	 	 = [args.control, args.condition]

estimators = []
XX_train = []
yy_train = []
XX_test = []
yy_test = []

def extract_sets():
	# Obtaining the clinical data and splitting training/test set
	initial_labels = pd.read_csv(os.path.join(args.trainingdir, "clinical_data.filtered.tsv"), index_col=0, sep='\t', usecols=["SampleName","Classification"])  # Read the clinical data and randomly shuffling the dataframe
	
	# Removing potential outlier samples
	if args.removeoutliers:
		initial_labels = remove_outlier_samples(initial_labels)
	initial_labels = initial_labels.sample(frac=1, random_state=random_s).fillna(-1)  # Shuffle the labels
	
	# Checking and report class balance
	check_class_balance(initial_labels['Classification'].to_numpy())


	trainingset_labels = None
	testset_labels = None
	# If external validation, the whole dataset will be used for 
	# analysis and the external validation for independent testing
	if args.validationdir:
		print("\nATTENTION: External validation set was detected and will be used as test set!\n")
		trainingset_labels = initial_labels
		testset_labels = pd.read_csv(os.path.join(args.validationdir, "clinical_data.filtered.tsv"), index_col=0, sep='\t', usecols=["SampleName","Classification"])  # Read the clinical data into a dataframe
	# Else, the input set will be split into two partitions for training and test sets
	else:
		print("\nATTENTION: No external validation set was given, so we will randomly split\nthe data into training and test set in a 70%-30% fashion!\n")
		# Splitting the input data to training and testing in a 70%-30% fashion
		trainingset_labels, testset_labels = train_test_split(initial_labels, test_size=0.30, random_state=random_s, stratify=initial_labels["Classification"])
		# testset_labels.pop('Condition')
	
	# print("Training set:",trainingset_labels.shape, "\nTest set:",testset_labels.shape)
	# print(trainingset_labels, testset_labels)
	return trainingset_labels, testset_labels

def remove_outlier_samples(data):
	# When performing the gene and isoform expression transformation, we also perform an IQR analysis on the 
	# samples. The removal of the  outlier samples in based on that analysis and the samples can be found 
	# in the plots directory
	print('Removing outlier samples:')

	outlier_samples = []
	gene_outlier = os.path.join(args.trainingdir, 'plots', 'potential_outlier_samples.gene.tsv')
	isoform_outlier = os.path.join(args.trainingdir, 'plots', 'potential_outlier_samples.isoform.tsv')
	
	if os.path.exists(gene_outlier):
		with open(gene_outlier) as geoin:
			for line in geoin:
				if not line.startswith("x"):
					outlier_samples.append(line.strip())

	if os.path.exists(isoform_outlier):
		with open(isoform_outlier) as geoin:
			for line in geoin:
				if not line.startswith("x"):
					outlier_samples.append(line.strip())
	
	outlier_samples = list(set(outlier_samples))
	if outlier_samples:
		print(f"WARNING: The following {len(outlier_samples)} samples, will be excluded from\nanalysis since they got characterised as outliers!")
		# print(outlier_samples)
		data = data.drop(outlier_samples)
	else:
		print("WARNING: NO outlier samples were detected!")
	print()
	return data

def check_class_balance(y_train):
	print('Identifying class balance in the initial dataset:')

	global class_balance

	vals = []
	for i, c in enumerate(conditions, 1):
		n_examples = len(y_train[y_train==c])
		vals.append(n_examples)
		percent = n_examples / len(y_train) * 100
		print(f"condition{i} {c}: {n_examples}/{len(y_train)} ({percent:.1f}%)")

	if len(vals) == 2:
		 class_balance = isclose(vals[0], vals[1], abs_tol=50)
	else:
		print(f"ATTENTION: It seems that {len(vals)} classes were imported ({','.join(vals)}),\tbut the algorithm was disigned for only two..")
		sys.exit()
	return

def get_matrix_ordered():
	# Obtaining all expression matrices, but with 
	# specific order 
	training_matrices = [None, None, None, None, None, None]
	for matrix in os.listdir(args.trainingdir):
		if matrix.startswith("gene_expression"):
			training_matrices[0] = os.path.join(args.trainingdir, matrix)
		elif matrix.startswith("isoform_expression"):
			training_matrices[1] = os.path.join(args.trainingdir, matrix)
		elif matrix.startswith("alternative_isoform_expression"):
			training_matrices[2] = os.path.join(args.trainingdir, matrix)
		elif matrix.startswith("gene_fusion_expression"):
			training_matrices[3] = os.path.join(args.trainingdir, matrix)
		elif matrix.startswith("RNAediting_expression"):
			training_matrices[4] = os.path.join(args.trainingdir, matrix)
		elif matrix.startswith("snp_expression"):
			training_matrices[5] = os.path.join(args.trainingdir, matrix)

	validation_matrices	= [None, None, None, None, None, None]
	if args.validationdir:
		for valmatrix in os.listdir(args.validationdir):
			if valmatrix.startswith("gene_expression"):
				validation_matrices[0] = os.path.join(args.validationdir, valmatrix)
			elif valmatrix.startswith("isoform_expression"):
				validation_matrices[1] = os.path.join(args.validationdir, valmatrix)
			elif valmatrix.startswith("alternative_isoform_expression"):
				validation_matrices[2] = os.path.join(args.validationdir, valmatrix)
			elif valmatrix.startswith("gene_fusion_expression"):
				validation_matrices[3] = os.path.join(args.validationdir, valmatrix)
			elif valmatrix.startswith("RNAediting_expression"):
				validation_matrices[4] = os.path.join(args.validationdir, valmatrix)
			elif valmatrix.startswith("snp_expression"):
				validation_matrices[5] = os.path.join(args.validationdir, valmatrix)
	return zip(training_matrices, validation_matrices)

def read_expression_matrices(trainingMat, validationMat, trainingset_labels, testset_labels):
	method = None
	expression_matrix = None
	validation_expression_matrix = None
	### Gene expression
	if 'gene_expression_matrix.filtered.scaled' in trainingMat:
		# Reading the input expression file
		method = 'gene_expression'
		expression_matrix = pd.read_csv(trainingMat, index_col=0, sep='\t') 
		if validationMat != None:
			validation_expression_matrix = pd.read_csv(validationMat, index_col=0, sep='\t')
	### Isoform expression
	elif 'isoform_expression_matrix.filtered.scaled' in trainingMat:
		# Reading the input expression file
		method = 'isoform_expression'
		expression_matrix = pd.read_csv(trainingMat, index_col=0, sep='\t')
		if validationMat != None:
			validation_expression_matrix = pd.read_csv(validationMat, index_col=0, sep='\t')
	### Alternative isoform expression
	elif 'alternative_isoform' in trainingMat:
		# Reading the input expression file
		method = 'alternative_isoform'
		expression_matrix = pd.read_csv(trainingMat, index_col=0, sep='\t')  
		if validationMat != None:
			validation_expression_matrix = pd.read_csv(validationMat, index_col=0, sep='\t')
	### Gene fusion expression
	elif 'gene_fusion' in trainingMat:
		# Reading the input expression file
		method = 'gene_fusion'
		expression_matrix = pd.read_csv(trainingMat, index_col=0, sep='\t')
		if validationMat != None:
			validation_expression_matrix = pd.read_csv(validationMat, index_col=0, sep='\t')
	### RNA editing events expression
	elif "rnaediting" in trainingMat.lower():
		# Reading the input expression file
		method = 'rnaediting'
		expression_matrix = pd.read_csv(trainingMat, index_col=False, sep='\t')
		expression_matrix.insert(0, 'mutationID', expression_matrix[expression_matrix.columns[0:3]].apply(lambda x: '.'.join(x.dropna().astype(str)),axis=1))
		expression_matrix = expression_matrix.drop(['Chromosome', 'Position', 'Substitution'], axis=1)  # Remove the positional columns 
		expression_matrix = expression_matrix.set_index('mutationID')  # Setting mutationID as index
		if validationMat != None:
			validation_expression_matrix = pd.read_csv(validationMat, index_col=False, sep='\t')
			validation_expression_matrix.insert(0, 'mutationID', validation_expression_matrix[validation_expression_matrix.columns[0:3]].apply(lambda x: '.'.join(x.dropna().astype(str)),axis=1))
			validation_expression_matrix = validation_expression_matrix.drop(['Chromosome', 'Position', 'Substitution'], axis=1)  # Remove the positional columns 
			validation_expression_matrix = validation_expression_matrix.set_index('mutationID')  # Setting mutationID as index
	### Single nucleotide variance expression
	elif "snp_expression" in trainingMat.lower():
		# Reading the input expression file
		method = 'snp'
		expression_matrix = pd.read_csv(trainingMat, index_col=False, sep='\t')
		expression_matrix.insert(0, 'mutationID', expression_matrix[expression_matrix.columns[0:4]].apply(lambda x: '.'.join(x.dropna().astype(str)),axis=1))
		expression_matrix = expression_matrix.drop(['Chromosome', 'Position', 'Ref', 'Alt'], axis=1)  # Remove the positional columns 
		expression_matrix = expression_matrix.set_index('mutationID')  # Setting mutationID as index
		if validationMat != None:
			validation_expression_matrix = pd.read_csv(validationMat, index_col=False, sep='\t')
			validation_expression_matrix.insert(0, 'mutationID', validation_expression_matrix[validation_expression_matrix.columns[0:4]].apply(lambda x: '.'.join(x.dropna().astype(str)),axis=1))
			validation_expression_matrix = validation_expression_matrix.drop(['Chromosome', 'Position', 'Ref', 'Alt'], axis=1)  # Remove the positional columns 
			validation_expression_matrix = validation_expression_matrix.set_index('mutationID')  # Setting mutationID as index


	# Keep only features common between the two independent datasets
	if validationMat != None:
		common_features = list(set(expression_matrix.index.values.tolist()) & set(validation_expression_matrix.index.values.tolist()))
		expression_matrix = expression_matrix.loc[common_features]
		validation_expression_matrix = validation_expression_matrix.loc[common_features]
	

	return training_testing_validation_sets(expression_matrix.transpose(), validation_expression_matrix, trainingset_labels, testset_labels, method)

def training_testing_validation_sets(expression_matrix, validation_expression_matrix, trainingset_labels, testset_labels, method):

	### Training Set
	# Remove the rows that exist in the training set
	training_mat = expression_matrix[expression_matrix.index.isin(trainingset_labels.index)]
	# Incorporating the labels at the end of the matrix
	training_mat = training_mat.join(trainingset_labels, how='left')
	training_mat = training_mat.sort_index()
	# Obtaining the final data and their labels
	X_train = rename_ids(training_mat.drop('Classification', axis=1), method)
	y_train = training_mat['Classification'].values.ravel()


	X_test = None
	y_test = None
	### Testing Set
	if args.validationdir:
		validation_mat = validation_expression_matrix.transpose().join(testset_labels, how='left')
		validation_mat = validation_mat.sort_index()
		# Obtaining the final data and their labels
		X_test = rename_ids(validation_mat.drop(['Classification'], axis=1), method)
		y_test = validation_mat[['Classification']].values.ravel()
	else:
		# Remove the rows that exist in the testing set
		testing_mat = expression_matrix[expression_matrix.index.isin(testset_labels.index)]
		# Incorporating the labels at the end of the matrix
		testing_mat = testing_mat.join(testset_labels, how='left')
		testing_mat = testing_mat.sort_index()
		# Obtaining the final data and their labels
		X_test = rename_ids(testing_mat.drop(['Classification'], axis=1), method)
		y_test = testing_mat[['Classification']].values.ravel()
	return X_train, X_test, y_train, y_test

def rename_ids(X, method):
	# Renaming all matrix feature IDs so we can identify 
	# which feature comes from which matrix
	if method == 'gene_expression':
		X = X.add_prefix('g')
	if method == 'isoform_expression':
		X = X.add_prefix('t')
	if method == 'alternative_isoform':
		X = X.add_prefix('a')
	if method == 'gene_fusion':
		X = X.add_prefix('f')
	if method == 'rnaediting':
		X = X.add_prefix('r')
	if method == 'snp':
		X = X.add_prefix('s')
	return X

def remove_correlated_features(X):
	### Removing highly correlated, outliers or quasi-constant features

	# Create correlation matrix
	print("Removing correlated and quasi-constant features..")
	corr_matrix = X.corr().abs()
	# Select upper triangle of correlation matrix
	upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
	# Find index of feature columns with correlation greater than 0.80
	to_drop = [column for column in upper.columns if any(upper[column] >= 0.80)]
	Xuncorr = X.drop(X[to_drop], axis=1)

	# Also remove quasi-constant features
	Xuncorr.drop(columns = get_constant_features(Xuncorr, threshold=0.99)['Var'].to_list(), inplace=True)
	return Xuncorr

def feature_selection(X_train, y_train, X_test, y_test, method, iteration):
	from catboost import CatBoostClassifier
	from sklearn.ensemble import BaggingClassifier
	from sklearn.preprocessing import PowerTransformer, QuantileTransformer, RobustScaler, StandardScaler, Normalizer, MaxAbsScaler


	method_id = method.lower().replace(" ","_")

	if not os.path.exists(ind_dir): os.makedirs(ind_dir)
	if not os.path.exists(voting_dir): os.makedirs(voting_dir)

	# Keep only features common between the two independent datasets
	print(f"Initial # of features in the training set: {X_train.shape[1]}")
	common_features = X_train[X_train.columns.intersection(X_test.columns)].columns.values.tolist()
	X_train = X_train[common_features]
	print(f"Remained common (with val set) # of features in the training set: {X_train.shape[1]}")
	
	X_train = remove_correlated_features(X_train)
	print(f'After correlation removal: {X_train.shape[1]}\n')

	if X_train.shape[1] < 5:
		# Safekeeping the selected features
		with open(os.path.join(ml_dir, "selected_features.tsv"), 'a') as fsout:  # Saving the features in a txt file
			fsout.write('{0}\tNone\n'.format(method.lower().replace(" ","_")))
		return

	
		
	fspipeline = None
	clpipeline = None
	if method == "Gene Expression":
		fspipeline = make_pipeline(MinMaxScaler(), RandomForestClassifier(random_state=random_s, class_weight='balanced'))
		clpipeline = make_pipeline(MinMaxScaler(), AdaBoostClassifier(ExtraTreesClassifier(random_state=random_s, class_weight='balanced'), random_state=random_s))
	
	elif method == "Isoform Expression": 
		fspipeline = make_pipeline(MinMaxScaler(), SVC(random_state=random_s, C=10, kernel='linear', class_weight='balanced', probability=True))
		clpipeline = make_pipeline(MinMaxScaler(), AdaBoostClassifier(ExtraTreesClassifier(random_state=random_s, class_weight='balanced'), random_state=random_s))

	elif method == "Alternative Isoform Expression":
		fspipeline = clpipeline = make_pipeline(StandardScaler(), LogisticRegression(solver='liblinear', C=1, class_weight='balanced', random_state=random_s))

	elif method == "Gene Fusion Expression":
		fspipeline = clpipeline =  LogisticRegression(solver='liblinear', class_weight='balanced', C=0.1, random_state=random_s)		
	
	elif method == "RNA-editing Expression":
		fspipeline = clpipeline = LogisticRegression(solver='liblinear', class_weight='balanced', C=1, random_state=random_s)

	elif method == "Mutation Profiling Expression":
		fspipeline = clpipeline = LogisticRegression(solver='liblinear', class_weight='balanced', C=5, random_state=random_s)
		


	skf = StratifiedKFold(n_splits=20, shuffle=True ,random_state=random_s)
	model = GeneticSelectionCV(
	    fspipeline, cv=skf, verbose=0,
	    scoring="accuracy", max_features=50,
	    n_population=130, n_generations=130, mutation_proba=0.2,
		n_jobs=15)
	
	model.fit(X_train, y_train)

	selected_features = X_train.columns[model.support_].tolist()
	print(f'{len(selected_features)} Features:', selected_features)

	# Safekeeping the selected features
	with open(os.path.join(ml_dir, "selected_features.tsv"), 'a') as fsout:  # Saving the features in a txt file
		fsout.write('{0}\t{1}\n'.format(method.lower().replace(" ","_"), ",".join(selected_features)))
	return

def model_training(X_train, y_train, X_test, y_test, method, iteration):
	
	method_id = method.lower().replace(" ","_")

	if not os.path.exists(ind_dir): os.makedirs(ind_dir)
	if not os.path.exists(voting_dir): os.makedirs(voting_dir)

	
		
	clpipeline = None
	if method == "Gene Expression":
		clpipeline = make_pipeline(MinMaxScaler(), AdaBoostClassifier(ExtraTreesClassifier(random_state=random_s, class_weight='balanced'), random_state=random_s))
	
	elif method == "Isoform Expression": 
		clpipeline = make_pipeline(MinMaxScaler(), AdaBoostClassifier(ExtraTreesClassifier(random_state=random_s, class_weight='balanced'), random_state=random_s))

	elif method == "Alternative Isoform Expression":
		clpipeline = make_pipeline(StandardScaler(), LogisticRegression(solver='liblinear', C=1, class_weight='balanced', random_state=random_s))

	elif method == "Gene Fusion Expression":
		clpipeline =  make_pipeline(Binarizer(threshold=0.4), LogisticRegression(solver='liblinear', class_weight='balanced', C=0.1, random_state=random_s))			
	
	elif method == "RNA-editing Expression":
		clpipeline = make_pipeline(Binarizer(threshold=0.4), LogisticRegression(solver='liblinear', C=5, class_weight='balanced', random_state=random_s))

	elif method == "Mutation Profiling Expression":
		clpipeline = make_pipeline(Binarizer(threshold=0.6), LogisticRegression(solver='liblinear', C=10, class_weight='balanced', random_state=random_s))
		
	
	selected_features = None
	# Retrieving the selected features
	with open(os.path.join(ml_dir, "selected_features.tsv")) as fin:
		for line in fin:
			if line.strip().split("\t")[0] == method_id:
				try:
					if method_id == "mutation_profiling_expression":
						selected_features = list(line.strip().split("\t")[1].split(",s"))
						selected_features = ['s'+item if not item.startswith('s') else item for item in selected_features]
								
					else:
						selected_features = list(line.strip().split("\t")[1].split(","))
				except IndexError:
					print(f'{method_id} did not output any features..')
					return

	# Choosing the selected features on each matrix
	X_train_selected = X_train[np.intersect1d(X_train.columns, selected_features)]
	X_test_selected = X_test[np.intersect1d(X_test.columns, selected_features)]


	tprs = []
	aucs = []
	accs = []
	sens = []
	spec = []
	f1 = []
	misclassifications = {}
	mean_fpr = np.linspace(0, 1, 100)
	fig, ax = plt.subplots(figsize=(8,6))
	plt.rcParams.update({'font.size': 10})
	skf = StratifiedKFold(n_splits=5, shuffle=True ,random_state=random_s)
	for i, (train, test) in enumerate(skf.split(X_train_selected, y_train), 1):
		
		# Obtaining the training, test and labels for each fold
		X1_train, X1_test, y1_train, y1_test = X_train_selected.iloc[train], X_train_selected.iloc[test], y_train[train], y_train[test]

		# Fitting the data to the classifier
		clpipeline.fit(X1_train, y1_train)

		# Evaluating our model with the test set
		training_acc_score, training_auc_score = calculate_scores(X1_test, y1_test, clpipeline)

		# # Compute ROC curve and area under the curve metrics
		viz = RocCurveDisplay.from_estimator(clpipeline, X1_test, y1_test, label=f"Fold.{i} | AUC = {training_auc_score:.2f}", alpha=0.2, lw=1, ax=ax, pos_label=conditions[1])
		interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
		interp_tpr[0] = 0.0
		tprs.append(interp_tpr)
		aucs.append(training_auc_score)
		accs.append(training_acc_score)

		# Calculating sensitivity and specificity per fold
		sensitivity, specificity, *_ = calculate_basic_metrics(X1_test, y1_test, clpipeline)
		f1_sc = f1_score(y1_test, clpipeline.predict(X1_test), pos_label=args.condition)
		f1.append(f1_sc)
		sens.append(sensitivity)
		spec.append(specificity)
		# Output accuracy and AUC scores
		# Output AUC, accuracy, sensitivity and specificity scores
		print(f'Fold.{i} - AUC:{training_auc_score:.2f}  Acc.:{training_acc_score:.2f}  Sen.:{sensitivity:.2f}  Spe.:{specificity:.2f}  F1:{f1_sc:.2f}') 
			
	# Obtaining the overall mean scores
	print(f'\n(Avg) Training set:\tAUC: {mean(aucs):.2f}  Acc.: {mean(accs):.2f}  Sen.: {mean(sens):.2f}  Spe.: {mean(spec):.2f}  F1:{mean(f1):.2f}')
	# Output stats
	with open(os.path.join(ml_dir,"metrics_overview.txt"), 'a') as oout:
		if iteration == 0:
			oout.write("InputMatrix\tUsedSet\tAUC\tAccuracy\tSensitivity\tSpecificity\tF1score\n")
		oout.write(f"{method_id}\tTrainingSet\t{mean(aucs):.2f}\t{mean(accs):.2f}\t{mean(sens):.2f}\t{mean(spec):.2f}\t{mean(f1):.2f}\n")


	# Train the model with the whole training set
	clpipeline.fit(X_train_selected, y_train)
	# Saving the trained model
	pickle.dump(clpipeline, open(os.path.join(ml_dir, f'trained_model.{method_id}.sav'), 'wb'))

	# Saving necessary info for the meta-classifier
	estimators.append((method_id, clpipeline))
	XX_train.append((method_id, X_train_selected))
	yy_train.append((method_id, y_train))
	XX_test.append((method_id, X_test_selected))
	yy_test.append((method_id, y_test))


	# Predicting the labels of the testing set!
	testing_acc_score, testing_auc_score = calculate_scores(X_test_selected, y_test, clpipeline)
	testing_sensitivity, testing_specificity, *_ = calculate_basic_metrics(X_test_selected, y_test, clpipeline)
	f1score = f1_score(y_test, clpipeline.predict(X_test_selected), pos_label=args.condition)


	print(f'Testing set:\t\tAUC: {testing_auc_score:.2f}  Acc.: {testing_acc_score:.2f}  Sen.: {testing_sensitivity:.2f}  Spe.: {testing_specificity:.2f}  F1:{f1score:.2f}')
	output_classification_scores(X_test_selected, y_test, clpipeline, voting_dir, method_id, 'Testing Set', iteration)
	confusion_matrix_plot(X_test_selected, y_test, clpipeline, ind_dir, method_id, 'testingset')
	with open(os.path.join(ml_dir,"metrics_overview.txt"), 'a') as oout:
		oout.write(f"{method_id}\tTestingSet\t{testing_auc_score:.2f}\t{testing_acc_score:.2f}\t{testing_sensitivity:.2f}\t{testing_specificity:.2f}\t{f1score:.2f}\n")

	# ROC plot
	# Obtaining the ROC plot where each fold is represented by a ROC curve
	ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='#7D1935', label='Chance', alpha=.5)
	mean_tpr = np.mean(tprs, axis=0)
	mean_tpr[-1] = 1.0
	mean_auc = auc(mean_fpr, mean_tpr)
	std_auc = np.std(aucs)
	ax.plot(mean_fpr, mean_tpr, color='#22577A', label=r'Training Set | AUC = %0.2f $\pm$ %0.2f' % (mean_auc, std_auc), lw=2, alpha=.6)
	std_tpr = np.std(tprs, axis=0)
	tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
	tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
	ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='#C9CCD5', alpha=.2, label=r'$\pm$ 1 std. dev.')
	ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05], title=method)
	ax.legend(loc="lower right")
	RocCurveDisplay.from_estimator(clpipeline, X_test_selected, y_test, label=f"Testing Set | AUC = {testing_auc_score:.2f}", color='#C36839', alpha=0.6, lw=2, ax=ax, pos_label=conditions[1])
	
	# Plotting the ROC curve
	ax.legend(loc='lower right', prop={'size': 8})
	fig.savefig(os.path.join(ind_dir,f'{method_id}.ROCplot.png'), format='png', dpi=600)
	plt.close(fig)
	return

def calculate_scores(X_data, y_data, classifier):
	acc_score = accuracy_score(y_data, classifier.predict(X_data), normalize=True)  # Accuracy score
	auc_score = roc_auc_score(y_data, classifier.predict_proba(X_data)[:,1])  # AUC score
	return acc_score, auc_score

def calculate_addional_scores(X_data, y_data, classifier):
	balanced_acc = balanced_accuracy_score(y_data, classifier.predict(X_data))
	f1 = f1_score(y_data, classifier.predict(X_data), pos_label=conditions[1])
	jac_score = jaccard_score(y_data, classifier.predict(X_data), pos_label=conditions[1])
	return balanced_acc, f1, jac_score

def calculate_basic_metrics(X_data, y_data, classifier):
	tn, fp, fn, tp = confusion_matrix(y_data, classifier.predict(X_data), labels=conditions).ravel()
	sensitivity = (tp / (tp + fn))
	specificity = (tn / (tn + fp))
	return sensitivity, specificity, tn, fp, fn, tp

def confusion_matrix_plot(X_data, y_data, classifier, wheretosave, methodname, mattype):
	# Confusion matrix
	# Obtaining all necessary metrics for the plot
	acc_score, auc_score = calculate_scores(X_data, y_data, classifier)
	_, _, tn, fp, fn, tp = calculate_basic_metrics(X_data, y_data, classifier)
	misclassified = fp + fn
	total = tn + fp + fn + tp
	# Creating custom fully informative labels for the plot
	fig2, ax2 = plt.subplots(figsize=(8,6))
	group_names = ['True Positive','False Negative','False Positive','True Negative']
	group_counts = [tp, fn, fp, tn]
	group_percentages = ["{0:.2%}".format(value) for value in group_counts/np.sum(group_counts)]
	group_counts = [f"{tp}/{np.sum(group_counts)}", f"{fn}/{np.sum(group_counts)}", 
					f"{fp}/{np.sum(group_counts)}", f"{tn}/{np.sum(group_counts)}"]
	labels = [f"{v1}\n{v2}\n{v3}" for v1, v2, v3 in zip(group_names, group_counts, group_percentages)]
	labels = np.asarray(labels).reshape(2,2)
	# Using the seaborn library to make the confusion matrix plot
	ax2 = sns.heatmap([[tp, fn],[fp, tn]], annot=labels, fmt='', cmap='BuPu')
	ax2.set_title(f'{methodname.replace("_"," ").title()} ({mattype}) Confusion Matrix\nAUC: {auc_score}  Acc.: {acc_score}  Miscls: {misclassified}/{total} ({(float(misclassified)/float(total)*100):.1f}%)')
	ax2.set_xlabel('\nPredicted Class')
	ax2.set_ylabel('Actual Class')

	## Ticket labels
	ax2.xaxis.set_ticklabels(reversed(conditions))
	ax2.yaxis.set_ticklabels(reversed(conditions))
	fig2.savefig(os.path.join(wheretosave, f'{methodname}.ConfusionMatrixPlot.{mattype}.png'), format='png', dpi=600)

	### Confusion Matrix ###
	print(f"Confusion matrix ({mattype}):")
	print(f"TP: {tp}\tFN: {fn}\nFP: {fp}\tTN: {tn}")
	return

def output_classification_scores(X_data, y_data, classifier, wheretosave, methodname, mattype, iteration):
	### Here we will save all kinds of metrics and scores that aim to 
	### better understand the classification output and will be the 
	### input for the ensemble voting classification
	# Obtain the misclassified labels and output them in a file

	misclassified = X_data.index.values[np.where(classifier.predict(X_data) != y_data)]
	print(f"{mattype}:\t\t{len(misclassified)}/{len(y_data)} ({(float(len(misclassified))/float(len(y_data))*100):.1f}%) samples were misclassified")
	with open(os.path.join(wheretosave, 'misclassifications.txt'), 'a') as fout: 
		fout.write("{0},{1}\n".format(methodname, ",".join(misclassified)))

	### Obtaining the statistics from the set
	pred = classifier.predict(X_data)
	class_labels = classifier.classes_.tolist()
	stats = zip(X_data.index.values.tolist(), 
		        y_data.tolist(),
		        pred.tolist(), 
		        np.round(classifier.predict_proba(X_data)[:,0],2).tolist(),
		        np.round(classifier.predict_proba(X_data)[:,1],2).tolist())
	stats = pd.DataFrame(stats)
	stats.columns = ["Sample", "True_label", "Predicted_label", f"Probability{conditions[1]}", f"Probability{conditions[0]}"]
	stats['Results'] = stats['True_label']==stats['Predicted_label']
	stats.to_csv(os.path.join(wheretosave, f'{methodname}.PredResults.txt'), index=False, sep='\t')
	return

def features_to_use():
	### We want to verify which feature matrices are adequate
	### to be used in the ensemble classifier
	adequate_features = []
	with open(os.path.join(ml_dir,"metrics_overview.txt")) as fin:
		for line in fin:
			if not line.startswith("InputMatrix"):
				if line.strip().split("\t")[1] == "TrainingSet" and float(line.strip().split("\t")[2]) >= args.cutoff:
					adequate_features.append(line.strip().split("\t")[0])
	if not adequate_features:
		print(f"\n\nUnfortunately, no feature matrix was able to score above the chosen cutoff ({args.cutoff}).\nThe default gene and isoform expression will be chosen.")
		adequate_features = ['Gene Expression matrix', 'Isoform Expression matrix']
	else:
		print(f"\n\nThe following feature matrices passed the training cutoff value of {args.cutoff}\nand will be used in the ensembl classifier (soft/majority voting)")
		for i, mat in enumerate(adequate_features,1):
			print(f"{i}.", mat.replace("_"," ").title(), "matrix")
	return adequate_features

def hard_voting(featurestouse):
	### The hard voting classification is being done based on the 
	### majority class lables from all feature matrices that 
	### passed the minimum cutoff score during training
	print("\nMAJORITY VOTING")
	multidf = []
	for i, filename in enumerate(glob.glob(os.path.join(voting_dir, f'*.PredResults.txt'))):
		if os.path.basename(filename).split(".")[0] in featurestouse:
			df = pd.read_csv(filename, index_col=['Sample','True_label'], header=0, sep="\t")
			df = df.sort_index()
			df.drop([f'Probability{args.condition}',f'Probability{args.control}','Results'], inplace=True, axis=1)
			multidf.append(df)

	hardVoting_results = pd.concat(multidf, axis=1, ignore_index=False)
	hardVoting_output = hardVoting_results.mode(axis=1).iloc[:, 0].tolist()
	hardVoting_results = hardVoting_results.drop(['Predicted_label'], axis=1)
	hardVoting_results['Predicted_label'] = hardVoting_output
	hardVoting_results = hardVoting_results.reset_index()

	hardVoting_results['Results'] = hardVoting_results['True_label'] == hardVoting_results['Predicted_label']
	hardVoting_results.to_csv(os.path.join(voting_dir, "hardVoting.TestingSet.results.tsv"), index=False, sep='\t')

	class_labels = [args.control, args.condition]
	label_to_idx = {label: idx for idx, label in enumerate(class_labels)}
	true_labels_numeric = hardVoting_results['True_label'].map(label_to_idx)
	predicted_lables_numeric = hardVoting_results['Predicted_label'].map(label_to_idx)
	
	auc = roc_auc_score(true_labels_numeric, predicted_lables_numeric)
	accuracy = accuracy_score(true_labels_numeric, predicted_lables_numeric)
	f1score = f1_score(true_labels_numeric, predicted_lables_numeric)
	tn, fp, fn, tp = confusion_matrix(true_labels_numeric, predicted_lables_numeric).ravel()
	misclas = fp + fn
	total = tn + fp + fn + tp

	print(f"AUC: {auc:.2f}", f"Acc.: {accuracy:.2f}", f"Sen.:{(tp /(tp + fn)):.2f}", f"Spe.:{(tn /(tn + fp)):.2f}", f"F1:{f1score:.2f}")
	print(f"Misclassified samples: {misclas}/{total} ({(misclas/total)*100:.1f}%)")
	return

def soft_voting(featurestouse):
	### The voting classification is being done based on the 
	### average probabilities of each feature matrix that 
	### passed the minimum cutoff score during training
	print("\nSOFT VOTING")

	multidf = []
	for i, filename in enumerate(glob.glob(os.path.join(voting_dir, f'*.PredResults.txt'))):
		if os.path.basename(filename).split(".")[0] in featurestouse:
			df = pd.read_csv(filename, index_col=['Sample','True_label'], header=0, sep="\t")
			df.drop(['Predicted_label','Results'], inplace=True, axis=1)
			multidf.append(df)

	softVoting_results = pd.concat(multidf, axis=1, ignore_index=False)
	softVoting_results = softVoting_results.groupby(by=softVoting_results.columns, axis=1).mean()
	softVoting_results = softVoting_results.reset_index()
	
	softVoting_results['Predicted_label'] = np.where(softVoting_results[f'Probability{args.condition}'] >= args.crosscut, args.condition, args.control)
	softVoting_results['Results'] = softVoting_results['True_label']==softVoting_results['Predicted_label']
	softVoting_results = softVoting_results[['Sample', 'True_label', 'Predicted_label', f'Probability{args.condition}', f'Probability{args.control}', 'Results']]
	softVoting_results.to_csv(os.path.join(voting_dir, "softVoting.TestingSet.results.tsv"), index=False, sep='\t')

	auc, accuracy = cm_plot("Testing Set", "soft Voting", softVoting_results['True_label'].values, softVoting_results[f'Probability{conditions[0]}'].values, softVoting_results[f'Probability{conditions[1]}'].values, softVoting_results['Predicted_label'].values)
	
	# ROC plot
	fig, ax = plt.subplots(figsize=(8,6))
	RocCurveDisplay.from_predictions(softVoting_results['True_label'].values, softVoting_results[f'Probability{conditions[1]}'].values, label=f"Testing Set | AUC = {auc:.2f} (Accuracy = {accuracy:.2f})", color='#937DC2', alpha=0.9, lw=2, ax=ax, pos_label=conditions[1])
	ax.legend(loc='lower right', prop={'size': 8})
	plt.title(f'Testing set Soft Voting ensemle learning classification')
	plt.savefig(os.path.join(voting_dir, 'softVoting.ROCplot.new.png'), dpi=900)
	plt.clf()
	plt.cla()
	return

def stacking(featurestouse):
	### The meta-learning classification is being done  
	### based on the average probabilities of each feature matrix that 
	### passed the minimum cutoff score during training
	print("\nSTACKING")


	estimators_ade = choose_adequate_data(estimators, featurestouse, 'estimators')
	XX_train_ade = choose_adequate_data(XX_train, featurestouse, 'Xdata')
	yy_train_ade = choose_adequate_data(yy_train, featurestouse, 'ydata')
	XX_test_ade = choose_adequate_data(XX_test, featurestouse, 'Xdata')
	yy_test_ade = choose_adequate_data(yy_test, featurestouse, 'ydata')


	meta_classifier = RandomForestClassifier(n_estimators=100, random_state=random_s)

	stacking_classifier = StackingClassifier(estimators=estimators_ade, final_estimator=meta_classifier)

	# Train the StackingClassifier
	stacking_classifier.fit(XX_train_ade, yy_train_ade)

	# Make predictions using the trained StackingClassifier
	stacking_predictions = stacking_classifier.predict(XX_test_ade)

	# Evaluate the performance
	stacking_accuracy, stacking_auc = calculate_scores(XX_test_ade, yy_test_ade, stacking_classifier)
	stacking_sensitivity, stacking_specificity, *_ = calculate_basic_metrics(XX_test_ade, yy_test_ade, stacking_classifier)
	stacking_f1score = f1_score(yy_test_ade, stacking_classifier.predict(XX_test_ade), pos_label=args.condition)
	
	print(f'Stacking:\t\tAUC: {stacking_auc:.2f}  Acc.: {stacking_accuracy:.2f}  Sen.: {stacking_sensitivity:.2f}  Spe.: {stacking_specificity:.2f}  F1:{stacking_f1score:.2f}')
	
	misclassified = XX_test_ade.index.values[np.where(stacking_classifier.predict(XX_test_ade) != yy_test_ade)]
	print(f"\t\t\t{len(misclassified)}/{len(yy_test_ade)} ({(float(len(misclassified))/float(len(yy_test_ade))*100):.1f}) samples were misclassified\n")
	# print()
	# print(classification_report(yy_test_ade, stacking_predictions))
	return

def choose_adequate_data(data, features, typeofdata):

	if typeofdata == 'estimators':
		for elms in data:
			if not elms[0] in features:
				data.remove((elms[0],elms[1]))
	
	elif typeofdata == 'Xdata':
		data_new = []
		df_final = pd.DataFrame()
		for elms in data:
			if not elms[0] in features:
				data.remove((elms[0],elms[1]))
		for elms in data:
			data_new.append(elms[1])
		
		data = data_new
		data = pd.concat(data, axis=1)
	
	elif typeofdata == 'ydata':
		data = data[0][1]

	return data

def cm_plot(usedset, votingtype, true_lbls, control_pred_probs, condition_pred_probs, pred_lbls):
	# Output stats
	auc, accuracy, sensitivity, specificity, tn, fp, fn, tp = metrics_based_on_predictions(true_lbls, control_pred_probs, condition_pred_probs, pred_lbls)
	f1score = f1_score(true_lbls, pred_lbls, pos_label=args.condition)
	print(f"AUC: {auc:.2f}  Acc.: {accuracy:.2f}  Sen.: {sensitivity:.2f}  Spe.: {specificity:.2f}  F1:{f1score:.2f}")
	
	# ROC plot
	fig, ax = plt.subplots(figsize=(8,6))
	RocCurveDisplay.from_predictions(true_lbls, condition_pred_probs, label=f"{usedset} | AUC = {auc:.2f} (Accuracy = {accuracy:.2f})", color='#937DC2', alpha=0.9, lw=2, ax=ax, pos_label=conditions[1])
	plt.title(f'{usedset} {votingtype.capitalize()} ensemle learning classification')
	ax.legend(loc='lower right', prop={'size': 8})
	plt.savefig(os.path.join(voting_dir, f'{votingtype.replace(" ","")}.{usedset.replace(" ","")}.ROCplot.png'), dpi=900)
	plt.clf()
	plt.cla()

	misclassified = fp + fn
	total = tn + fp + fn + tp
	print(f"Misclassified samples: {misclassified}/{total} ({(misclassified/total)*100:.1f}%)")
	# Creating custom fully informative labels for the plot
	fig2, ax2 = plt.subplots(figsize=(8,6))
	group_names = ['True Positive','False Negative','False Positive','True Negative']
	group_counts = [tp, fn, fp, tn]
	group_percentages = ["{0:.2%}".format(value) for value in group_counts/np.sum(group_counts)]
	group_counts = [f"{tp}/{np.sum(group_counts)}", f"{fn}/{np.sum(group_counts)}", 
					f"{fp}/{np.sum(group_counts)}", f"{tn}/{np.sum(group_counts)}"]
	labels = [f"{v1}\n{v2}\n{v3}" for v1, v2, v3 in zip(group_names, group_counts, group_percentages)]
	labels = np.asarray(labels).reshape(2,2)
	# Using the seaborn library to make the confusion matrix plot
	ax2 = sns.heatmap([[tp, fn],[fp, tn]], annot=labels, fmt='', cmap='BuPu')
	ax2.set_title(f'{usedset} {votingtype} Confusion Matrix\nAUC: {auc:.2f}  Acc.: {accuracy:.2f}  Miscls: {misclassified}/{total} ({(float(misclassified)/float(total)*100):.1f}%)')
	ax2.set_xlabel('\nPredicted Class')
	ax2.set_ylabel('Actual Class')

	## Ticket labels
	ax2.xaxis.set_ticklabels(reversed(conditions))
	ax2.yaxis.set_ticklabels(reversed(conditions))
	fig2.savefig(os.path.join(voting_dir, f'{votingtype.replace(" ","")}.{usedset.replace(" ","")}.ConfusionMatrixPlot.png'), format='png', dpi=600)
	
	with open(os.path.join(ml_dir,"metrics_overview.txt"), 'a') as oout:
		oout.write(f"{votingtype}\t{usedset}\t{auc:.2f}\t{accuracy:.2f}\t{sensitivity:.2f}\t{specificity:.2f}\n")
	return auc, accuracy

def metrics_based_on_predictions(y_data, control_pred_probs, condition_pred_probs, predicted_lbls):
	auc_score = roc_auc_score(y_data, control_pred_probs)  # AUC score
	acc_score = accuracy_score(y_data, predicted_lbls, normalize=True)  # Accuracy score
	tn, fp, fn, tp = confusion_matrix(y_data.ravel(), predicted_lbls, labels=conditions).ravel()
	sensitivity = tp / (tp + fn)
	specificity = tn / (tn + fp)
	return auc_score, acc_score, sensitivity, specificity, tn, fp, fn, tp

class biological_relevance:

	def __init__(self, featurestouse):

		if not os.path.exists(biorel_dir): os.makedirs(biorel_dir)

		gene_db, transcript_db, reference, findgeneid = self.retrieve_reference_annot()

		connect_features = []
		for features in featurestouse:
			if features == 'gene_expression':
				connect_features = self.retrieve_gene_details('gene_expression', gene_db, connect_features)
				selected_features = ','.join([lst[2] for lst in connect_features if lst[0] == 'gene'])
				self.run_gprofiler('GeneExpression', selected_features)

			elif features == 'isoform_expression':
				connect_features = self.retrieve_transcript_details('isoform_expression', transcript_db, connect_features)
				selected_features = ','.join([lst[2] for lst in connect_features if lst[0] == 'isoform'])
				self.run_gprofiler('IsoformExpression', selected_features)

			elif features == 'alternative_isoform_expression':
				connect_features = self.retrieve_gene_details('alternative_isoform_expression', gene_db, connect_features)
				selected_features = ','.join([lst[2] for lst in connect_features if lst[0] == 'alternative_isoform'])
				self.run_gprofiler('AlternativeIsoformExpression', selected_features)

			elif features == 'gene_fusion_expression':
				connect_features = self.retrieve_gene_details('gene_fusion_expression', gene_db, connect_features)
				selected_features = ','.join([lst[2] for lst in connect_features if lst[0] == 'gene_fusion'])
				self.run_gprofiler('GeneFusionExpression', selected_features)

			elif features == 'rna-editing_expression':
				connect_features = self.retrieve_enaediting_details('rna-editing_expression', reference, connect_features)
				selected_features = ','.join([lst[2] for lst in connect_features if lst[0] == 'rna-editing'])
				self.run_gprofiler('RNAediting', selected_features)

			elif features == 'mutation_profiling_expression':
				connect_features = self.retrieve_variant_details('mutation_profiling_expression', reference, connect_features)
				selected_features = ','.join([lst[2] for lst in connect_features if lst[0] == 'snv'])
				self.run_gprofiler('SNV', selected_features)

		self.extract_summary(connect_features)
		return


	def retrieve_reference_annot(self):

		reference = []
		findgeneid = {}
		gene_db = {}
		transcript_db = {}
		with open(args.refannot) as refin:
			for line in refin:
				if not line.startswith("#") and line.strip().split("\t")[2] == 'transcript':
					gene_name = line.strip().split("gene_name \"")[1].split("\";")[0]
					gene_id = line.strip().split("gene_id \"")[1].split("\";")[0]
					transcript_id = line.strip().split("transcript_id \"")[1].split("\";")[0]
					gene_type = line.strip().split("gene_type \"")[1].split("\";")[0]
					chromosome = line.strip().split("\t")[0]
					start = line.strip().split("\t")[3]
					end = line.strip().split("\t")[4]
					strand = line.strip().split("\t")[6]
					findgeneid[gene_id.split(".")[0]] = gene_id

					gene_db[gene_id] = [gene_name, gene_id, transcript_id, gene_type, chromosome, start, end, strand]
					transcript_db[transcript_id] = [gene_name, gene_id, transcript_id, gene_type, chromosome, start, end, strand]
					reference.append([gene_name, gene_id, transcript_id, gene_type, chromosome, start, end, strand])
					# if transcript_id == "ENST0000666013.1":
					# 	print([gene_name, gene_id, transcript_id, gene_type, chromosome, start, end, strand])
		return gene_db, transcript_db, reference, findgeneid

	def retrieve_gene_details(self, feature, gene_db, connect_features):

		with open(os.path.join(ml_dir, 'selected_features.tsv')) as fin:
			for line in fin:
				feature_name = line.strip().split("\t")[0]
				if feature_name in feature:
					for gene in line.strip().split("\t")[1].split(","):
						geneID = gene[1:]
						if geneID in gene_db:
							add_elm = [feature_name.replace("_expression",""), gene_db[geneID][0],gene_db[geneID][1],'None',gene_db[geneID][3],gene_db[geneID][4],'None','None',gene_db[geneID][7]]
							if not add_elm in connect_features:
								connect_features.append(add_elm)
					
		return connect_features

	def retrieve_transcript_details(self, feature, transcript_db, connect_features):
		
		not_found = {}
		with open(os.path.join(ml_dir, 'selected_features.tsv')) as fin:
			for line in fin:
				feature_name = line.strip().split("\t")[0]
				if feature_name in feature:
					for transcript in line.strip().split("\t")[1].split(","):
						transcriptID = transcript[1:]
						if transcriptID in transcript_db:
							add_elm = [feature_name.replace("_expression","")] + transcript_db[transcriptID]
							if not add_elm in connect_features:
								connect_features.append(add_elm)
							# else:
							# 	not_found[transcriptID] = None

		# print(not_found)
		return connect_features

	def retrieve_enaediting_details(self, feature, reference, connect_features):
		
		gene_expression = {}
		with open(os.path.join(args.trainingdir, 'gene_expression_matrix.filtered.scaled.tsv')) as gin:
			for line in gin:
				if not line.startswith("featureID"):
					gene_expression[line.strip().split("\t")[0]] = sum([float(i) for i in line.strip().split("\t")[1:]])	

		with open(os.path.join(ml_dir, 'selected_features.tsv')) as fin:
			for line in fin:
				feature_name = line.strip().split("\t")[0]
				if feature_name in feature:
					for elms in line.strip().split("\t")[1].split(","):
						pot_genes = []
						chrom = elms[1:].split(".")[0]
						pos = elms[1:].split(".")[1]
						for sublst in reference:
							if chrom == sublst[4]:
								if int(sublst[5]) <= int(pos) <= int(sublst[6]):
									pot_genes.append(f'{sublst[1]}_{elms}')
						
						if len(list(set(pot_genes))) == 1:
							geneid = list(set(pot_genes))[0].split("_")[0]
							edit = list(set(pot_genes))[0].split("_")[1][1:]
							edit = edit[:-1] + '>' + edit[-1:]
							for sublst in reference:
								if geneid in sublst:
									add_elm = ([feature_name.replace("_expression",""),sublst[0],sublst[1],edit,sublst[3],sublst[4],'None','None',sublst[7]])
									if add_elm not in connect_features:
										connect_features.append(add_elm)
						
						elif len(list(set(pot_genes))) > 1:
							geneid = None
							geneid1 = list(set(pot_genes))[0].split("_")[0]
							geneid2 = list(set(pot_genes))[1].split("_")[0]
							edit = list(set(pot_genes))[0].split("_")[1][1:]
							edit = edit[:-1] + '>' + edit[-1:]
							
							if not geneid1 in gene_expression:
								gene_expression[geneid1] = 0
							if not geneid2 in gene_expression:
								gene_expression[geneid2] = 0

							if gene_expression[geneid1] > gene_expression[geneid2]:
								geneid = geneid1
							elif gene_expression[geneid1] <= gene_expression[geneid2]:
								geneid = geneid2
							for sublst in reference:
								if geneid in sublst:
									add_elm = ([feature_name.replace("_expression",""),sublst[0],sublst[1],edit,sublst[3],sublst[4],'None','None',sublst[7]])
									if add_elm not in connect_features:
										connect_features.append(add_elm)
		return connect_features

	def retrieve_variant_details(self, feature, reference, connect_features):
		
		gene_expression = {}
		with open(os.path.join(args.trainingdir, 'gene_expression_matrix.filtered.scaled.tsv')) as gin:
			for line in gin:
				if not line.startswith("featureID"):
					gene_expression[line.strip().split("\t")[0]] = sum([float(i) for i in line.strip().split("\t")[1:]])	

		with open(os.path.join(ml_dir, 'selected_features.tsv')) as fin:
			for line in fin:
				feature_name = line.strip().split("\t")[0]
				if feature_name in feature:
					for elms in line.strip().split("\t")[1].split(",s"):
						elms = 's'+elms if not elms.startswith('s') else elms
						pot_genes = []
						chrom = elms[1:].split(".")[0]
						pos = elms[1:].split(".")[1]

						for sublst in reference:
							if chrom == sublst[4]:
								if int(sublst[5]) <= int(pos) <= int(sublst[6]):
									pot_genes.append(f'{sublst[1]}_{elms}')
						
						if len(list(set(pot_genes))) == 1:
							geneid = list(set(pot_genes))[0].split("_")[0]
							edit = list(set(pot_genes))[0].split("_")[1][1:]
							chromosome = edit.strip().split(".")[0]
							position = edit.strip().split(".")[1]
							fromvar = edit.strip().split(".")[2]
							tovar = edit.strip().split(".")[3]
							edit = f'{chromosome}.{position}.{fromvar}>{tovar}'
							for sublst in reference:
								if geneid in sublst:
									add_elm = (["snv",sublst[0],sublst[1],edit,sublst[3],sublst[4],'None','None',sublst[7]])
									if add_elm not in connect_features:
										connect_features.append(add_elm)
						
						elif len(list(set(pot_genes))) > 1:
							pot_genes = [*set(pot_genes)]
							geneid = None
							geneid1 = list(set(pot_genes))[0].split("_")[0]
							geneid2 = list(set(pot_genes))[1].split("_")[0]

							# print(geneid1, geneid2)
							edit = list(set(pot_genes))[0].split("_")[1][1:]
							chromosome = edit.strip().split(".")[0]
							position = edit.strip().split(".")[1]
							fromvar = edit.strip().split(".")[2]
							tovar = edit.strip().split(".")[3]
							edit = f'{chromosome}.{position}.{fromvar}>{tovar}'

							if not geneid1 in gene_expression:
								gene_expression[geneid1] = 0
							if not geneid2 in gene_expression:
								gene_expression[geneid2] = 0

							if gene_expression[geneid1] > gene_expression[geneid2]:
								geneid = geneid1
							elif gene_expression[geneid1] <= gene_expression[geneid2]:
								geneid = geneid2
							
							for sublst in reference:
								if geneid in sublst:
									add_elm = (["snv",sublst[0],sublst[1],edit,sublst[3],sublst[4],'None','None',sublst[7]])
									if add_elm not in connect_features:
										connect_features.append(add_elm)
		return connect_features

	def run_gprofiler(self, feature_mat, selected_features):

		# Filtering the Gene expression matrix 
		geneprofiler = " ".join([
		"Rscript",  # Call Rscript
		os.path.join(rscripts, "enrichment_analysis.R"),  # Calling the R script
		selected_features,  # Input gene list
		biorel_dir,  # Output dir
		feature_mat,  # Feature matrix
		"2>>", os.path.join(biorel_dir, "feature_mat_gProfiler-report.txt")])  # Directory where all reports reside
		subprocess.run(geneprofiler, shell=True)
		return

	def extract_summary(self, connect_features):

		output_file = os.path.join(biorel_dir, "biological_relevance.tsv")

		all_features = [sublst[2] for sublst in connect_features]

		table = {"GeneID": [], "Gene Name": [], "Gene Type": [], "Occurrence": [], "Matrix of origin": [], "Feature of Origin": []}
		feature_stats = dict(OrderedDict(Counter(all_features).most_common()))
		for key, value in feature_stats.items():
			details = [sublts for sublts in connect_features if key in sublts]
			if len(details) == 1:
				table["GeneID"].extend([key])
				table["Occurrence"].extend([value])
				table["Gene Name"].extend([details[0][1]])
				table["Gene Type"].extend([details[0][4].replace("_"," ")])
				table["Matrix of origin"].extend([details[0][0].replace("_"," ")])
				if details[0][0] in ['isoform_expression','rna-editing','snv']:
					table["Feature of Origin"].extend([details[0][3]])
				else:
					table["Feature of Origin"].extend([details[0][2]])
			if len(details) > 1:
				matoforigin = []
				featoforigin = []
				table["GeneID"].extend([key])
				table["Occurrence"].extend([value])
				table["Gene Name"].extend([details[0][1]])
				table["Gene Type"].extend([details[0][4].replace("_"," ")])
				matoforigin.extend([sublst[0].replace("_"," ") for sublst in details])
				table["Matrix of origin"].extend(['/'.join(matoforigin)])
				for sublst in details:
					# print(sublst)
					if sublst[0] in ('gene','alternative_isoform','gene_fusion'):
						featoforigin.extend([sublst[2]])
					elif sublst[0] in ('isoform','rna-editing','snv'):
						featoforigin.extend([sublst[3]])
				table["Feature of Origin"].extend(['/'.join(featoforigin)])
		df = pd.DataFrame(table)
		df.to_csv(output_file, index=False, sep='\t')


		# Filtering the Gene expression matrix 
		# print(f'\n{datetime.now().strftime("%d.%m.%Y %H:%M")}  Preparing geneID list for online enrichment analysis using gProfiler...')
		geneprofiler = " ".join([
		"Rscript",  # Call Rscript
		os.path.join(rscripts, "enrichment_analysis.R"),  # Calling the R script
		','.join(df['GeneID'].tolist()),  # Input features
		biorel_dir,  # Output file
		"AllFeaturesAllMatrices",  # All features are given
		"2>>", os.path.join(biorel_dir, "biological_rel_gProfiler-report.txt")])  # Directory where all reports reside
		subprocess.run(geneprofiler, shell=True)

		print(f"Please open the following file to proceed with an online enrichment analysis of your data:\n{os.path.join(biorel_dir, 'gProfiler_click_link.txt')}")
		return



def main():
	print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} CLASSIFICATION ANALYSIS (MACHINE LEARNING APPROACH)')
	print(f'\t\t{project} PROJECT')

	if not os.path.exists(ml_dir): os.makedirs(ml_dir)
	

	# We will split our datasets to training and testing in an 70%-30% 
	# fashion, as well as getting an external validation set, if available
	trainingset_labels, testset_labels = extract_sets()


	for iteration, (trainingMat, validationMat) in enumerate(get_matrix_ordered()):
		# Gene expression
		if os.path.basename(trainingMat) == 'gene_expression_matrix.filtered.scaled.tsv':
			print(f'\n{datetime.now().strftime("%d.%m.%Y %H:%M")} GENE EXPRESSION')
			X_train, X_test, y_train, y_test = read_expression_matrices(trainingMat, validationMat, trainingset_labels, testset_labels)
			feature_selection(X_train, y_train, X_test, y_test, 'Gene Expression', iteration)
			model_training(X_train, y_train, X_test, y_test, 'Gene Expression', iteration)

		# Transcript expression
		elif os.path.basename(trainingMat) == 'isoform_expression_matrix.filtered.scaled.tsv':
			print(f'\n{datetime.now().strftime("%d.%m.%Y %H:%M")} ISOFORM EXPRESSION')
			X_train, X_test, y_train, y_test = read_expression_matrices(trainingMat, validationMat, trainingset_labels, testset_labels)
			feature_selection(X_train, y_train, X_test, y_test, 'Isoform Expression', iteration)
			model_training(X_train, y_train, X_test, y_test, 'Isoform Expression', iteration)

		# Alternative Isoform expression
		elif os.path.basename(trainingMat) == 'alternative_isoform_expression_matrix.filtered.tsv':
			print(f'\n{datetime.now().strftime("%d.%m.%Y %H:%M")} ALTERNATIVE ISOFORM EXPRESSION')
			X_train, X_test, y_train, y_test = read_expression_matrices(trainingMat, validationMat, trainingset_labels, testset_labels)
			feature_selection(X_train, y_train, X_test, y_test, 'Alternative Isoform Expression', iteration)
			model_training(X_train, y_train, X_test, y_test, 'Alternative Isoform Expression', iteration)

		# Gene fusion
		elif os.path.basename(trainingMat) == 'gene_fusion_expression_matrix.filtered.tsv':
			print(f'\n{datetime.now().strftime("%d.%m.%Y %H:%M")} GENE FUSION EXPRESSION')
			X_train, X_test, y_train, y_test = read_expression_matrices(trainingMat, validationMat, trainingset_labels, testset_labels)
			feature_selection(X_train, y_train, X_test, y_test, 'Gene Fusion Expression', iteration)
			model_training(X_train, y_train, X_test, y_test, 'Gene Fusion Expression', iteration)

		# RNA editing
		elif os.path.basename(trainingMat) == 'RNAediting_expression_matrix.filtered.tsv':
			print(f'\n{datetime.now().strftime("%d.%m.%Y %H:%M")} RNA-EDITING EXPRESSION')
			X_train, X_test, y_train, y_test = read_expression_matrices(trainingMat, validationMat, trainingset_labels, testset_labels)
			feature_selection(X_train, y_train, X_test, y_test, 'RNA-editing Expression', iteration)
			model_training(X_train, y_train, X_test, y_test, 'RNA-editing Expression', iteration)	

		# SNP matrix
		elif os.path.basename(trainingMat) == 'snp_expression_matrix.genotype.filtered.tsv':
			print(f'\n{datetime.now().strftime("%d.%m.%Y %H:%M")} MUTATION PROFILING EXPRESSION')
			X_train, X_test, y_train, y_test = read_expression_matrices(trainingMat, validationMat, trainingset_labels, testset_labels)
			feature_selection(X_train, y_train, X_test, y_test, 'Mutation Profiling Expression', iteration)
			model_training(X_train, y_train, X_test, y_test, 'Mutation Profiling Expression', iteration)


	adequate_features = features_to_use()
	# hard_voting(adequate_features)
	# stacking(adequate_features) 
	soft_voting(adequate_features)
	biological_relevance(adequate_features)
	

	print(f'The pipeline finished after {datetime.now() - startTime}')
	
if __name__ == "__main__": main()
