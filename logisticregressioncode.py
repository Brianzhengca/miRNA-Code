# Code for Logistic Regression Classifier
import sklearn 
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_validate
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import precision_recall_fscore_support
from sklearn.metrics import accuracy_score
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import mean_squared_error
from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import log_loss
from sklearn import metrics
from sklearn import svm
from sklearn.naive_bayes import GaussianNB
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import MultinomialNB
from sklearn.svm import LinearSVC

from sklearn.preprocessing import MinMaxScaler

import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
import sklearn

from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from math import sqrt

import keras
from keras.models import Sequential
from keras.layers import Dense
from keras.utils import to_categorical 
from keras.optimizers import SGD
from imblearn.over_sampling import SMOTE

design_matrix = "colorectaltobreastdesign"
exp_matrix = "colorectaltobreastexp"
cv = KFold(n_splits=10, random_state=1, shuffle=True)
skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=1)

def pred():
    design = pd.read_csv(design_matrix, sep="\t", index_col=0)
    exp = pd.read_csv(exp_matrix, sep="\t", index_col=0)

    matrix = []
    gene_list = [line.rstrip() for line in open('sig_genes.txt')]

    matrix.append(design["colorectal cancer"].tolist())

    for gene in gene_list:
        l = []
        for index in design.index:
            l.append(exp.at[gene, index])
        matrix.append(l)

    matrix = pd.DataFrame(matrix)
    matrix = matrix.T
    matrix.columns = ["state"] + gene_list
    matrix = matrix.astype({"state":int})

    Y = matrix[["state"]]
    X = matrix.loc[:, matrix.columns!="state"]

    x_train, x_test, y_train, y_test = sklearn.model_selection.train_test_split(X, Y, test_size=0.3, random_state=1)

    logit = LogisticRegression(solver="lbfgs", max_iter=1000)
    lst_accu_stratified = []
    for train_index, test_index in skf.split(X, Y): 
        X_train_fold, X_test_fold = X.iloc[train_index,:], X.iloc[test_index,:] 
        y_train_fold, y_test_fold = Y.iloc[train_index,:], Y.iloc[test_index,:]
        sm = SMOTE()
        X_train_oversampled, y_train_oversampled = sm.fit_resample(X_train_fold, y_train_fold)
        logit.fit(X_train_oversampled, y_train_oversampled.values.ravel())
        lst_accu_stratified.append(logit.score(X_test_fold, y_test_fold))
    logit.fit(x_train.values, y_train.values.reshape(-1,).ravel())
    y_pred = logit.predict(x_test.values)
    pred_prob1 = logit.predict_proba(x_test.values)
    fpr1, tpr1, thresh1 = roc_curve(y_test, pred_prob1[:,1], pos_label=1)   
    return(balanced_accuracy_score(y_test, y_pred)), fpr1, tpr1, roc_auc_score(y_test, pred_prob1[:, 1]), sum(lst_accu_stratified)/len(lst_accu_stratified), confusion_matrix(y_test, y_pred)

def pred_unique():
    design = pd.read_csv(design_matrix, sep="\t", index_col=0)
    exp = pd.read_csv(exp_matrix, sep="\t", index_col=0)

    matrix = []
    gene_list = [line.rstrip() for line in open('sig_genes_extended')]

    matrix.append(design["colorectal cancer"].tolist())

    for gene in gene_list:
        l = []
        for index in design.index:
            l.append(exp.at[gene, index])
        matrix.append(l)

    matrix = pd.DataFrame(matrix)
    matrix = matrix.T
    matrix.columns = ["state"] + gene_list
    matrix = matrix.astype({"state":int})
    #print(matrix.head())

    Y = matrix[["state"]]
    X = matrix.loc[:, matrix.columns!="state"]

    x_train, x_test, y_train, y_test = sklearn.model_selection.train_test_split(X, Y, test_size=0.3, random_state=1)

    logit = LogisticRegression(solver="lbfgs", max_iter=1000)
    lst_accu_stratified = []
    for train_index, test_index in skf.split(X, Y): 
        X_train_fold, X_test_fold = X.iloc[train_index,:], X.iloc[test_index,:] 
        y_train_fold, y_test_fold = Y.iloc[train_index,:], Y.iloc[test_index,:]
        sm = SMOTE()
        X_train_oversampled, y_train_oversampled = sm.fit_resample(X_train_fold, y_train_fold)
        logit.fit(X_train_oversampled, y_train_oversampled.values.ravel())
        lst_accu_stratified.append(logit.score(X_test_fold, y_test_fold))
    logit.fit(x_train.values, y_train.values.reshape(-1,).ravel())
    y_pred = logit.predict(x_test.values)
    pred_prob1 = logit.predict_proba(x_test.values)
    fpr1, tpr1, thresh1 = roc_curve(y_test, pred_prob1[:,1], pos_label=1)
    return(balanced_accuracy_score(y_test, y_pred)), fpr1, tpr1, roc_auc_score(y_test, pred_prob1[:, 1]), sum(lst_accu_stratified)/len(lst_accu_stratified), confusion_matrix(y_test, y_pred)

def bio_pred():
    design = pd.read_csv(design_matrix, sep="\t", index_col=0)

    matrix = []

    matrix.append(design["colorectal cancer"].tolist())
    matrix.append(design["age"].tolist())
    matrix.append(design["sex"].tolist())
    matrix = pd.DataFrame(matrix)
    matrix = matrix.T
    matrix.columns = ["state", "age", "sex"]
    matrix = matrix.astype({"state":int})
    #print(matrix.head())

    Y = matrix[["state"]]
    X = matrix.loc[:, matrix.columns!="state"]

    x_train, x_test, y_train, y_test = sklearn.model_selection.train_test_split(X, Y, test_size=0.3, random_state=1)

    logit = LogisticRegression(solver="lbfgs", max_iter=1000)
    lst_accu_stratified = []
    for train_index, test_index in skf.split(X, Y): 
        X_train_fold, X_test_fold = X.iloc[train_index,:], X.iloc[test_index,:] 
        y_train_fold, y_test_fold = Y.iloc[train_index,:], Y.iloc[test_index,:]
        sm = SMOTE()
        X_train_oversampled, y_train_oversampled = sm.fit_resample(X_train_fold, y_train_fold)
        logit.fit(X_train_oversampled, y_train_oversampled.values.ravel())
        lst_accu_stratified.append(logit.score(X_test_fold, y_test_fold))
    logit.fit(x_train.values, y_train.values.reshape(-1,).ravel())
    pred_prob1 = logit.predict_proba(x_test.values)
    y_pred = logit.predict(x_test.values)
    fpr1, tpr1, thresh1 = roc_curve(y_test, pred_prob1[:,1], pos_label=1)
    return balanced_accuracy_score(y_test, y_pred), fpr1, tpr1, roc_auc_score(y_test, pred_prob1[:, 1]), sum(lst_accu_stratified)/len(lst_accu_stratified), confusion_matrix(y_test, y_pred)

normal, fpr1, tpr1, auc1, cross1, conf1 = pred()
unique, fpr2, tpr2, auc2, cross2, conf2 = pred_unique()
baseline, fpr3, tpr3, auc3, cross3, conf3 = bio_pred()

# plot roc curves
plt.plot(fpr1, tpr1, linestyle='--',color='orange', label='Classifier w/ Non-unique Genes')
plt.plot(fpr2, tpr2, linestyle='--',color='green', label='Classifier w/ Unique Genes')
plt.plot(fpr3, tpr3, linestyle='--', color="red", label="Baseline w/ Age and Sex")
# title
plt.title('ROC curve')
# x label
plt.xlabel('False Positive Rate')
# y label
plt.ylabel('True Positive rate')

plt.legend(loc='best')
plt.savefig('ROC',dpi=300)
print("Classifier w/ Non-unique Genes Accuracy:", normal, "AUC Score:", auc1, "Cross-Validation with SMOTE:", cross1)
print("Classifier w/ Unique Genes Accuracy:", unique, "AUC Score:", auc2, "Cross-Validation with SMOTE", cross2)
print("Baseline:", baseline, "AUC Score:", auc3, "Cross-Validation with SMOTE:", cross3)
print("Classifier w/ Non-unique Genes")
print(conf1)
print("Classifier w/ Unique Genes")
print(conf2)
print("Baseline")
print(conf3)
plt.show()