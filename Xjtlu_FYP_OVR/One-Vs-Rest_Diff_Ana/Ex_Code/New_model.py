
import pandas as pd
from sklearn.model_selection import train_test_split, cross_val_score, KFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix, recall_score, precision_score, classification_report
import math
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, cross_val_score, KFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report
from xgboost import XGBClassifier
import xgboost as xgb
from catboost import CatBoostClassifier, Pool, cv
# Load your data into df here
df=pd.read_csv('H:/FYP/New_Data/Chayi_gene/Gene_matrix_200', index_col=[0], sep='\t')
df=df.T
x=df.drop(["Label"],axis=1)
for k in x.columns:
    a = k.split('|')
    x.rename(columns={k:a[0]},inplace=True)

df1=np.zeros((x.shape[0],x.shape[1]),dtype=np.float32)
for a in range(0,x.shape[0]-1):
    for b in range(0,x.shape[1]):
        if float(x.iloc[a,b])<=1:
            x.iloc[a,b]=1
        df1[a,b]=math.log2(float(x.iloc[a,b]))

y=df['Label']
df2=pd.DataFrame(df1,columns=x.columns,index=x.index)

# X_train,X_test,y_train,y_test=train_test_split(df2,y,test_size=0.2,random_state=42)
# Split into features and target

# Split the data into train and test sets
X_train, X_test, y_train, y_test = train_test_split(df2, y, test_size=0.2, random_state=42)



# Initialize the CatBoost model
catboost_model = CatBoostClassifier(loss_function='MultiClass', classes_count=33, random_seed=42)

# Perform 10-fold cross-validation on the training set
cv_data = Pool(X_train, y_train)
cv_scores = cv(cv_data, catboost_model.get_params(), fold_count=10, plot=False)

print("Cross-validation scores: ", cv_scores)
print("Average cross-validation score: ", np.mean(cv_scores['test-MultiClass-mean']))

# Train the model
catboost_model.fit(X_train, y_train, verbose=False)

# Predict labels for the training set and calculate accuracy
train_predictions = catboost_model.predict(X_train)
train_accuracy = accuracy_score(y_train, train_predictions)
print("Training set accuracy: ", train_accuracy)

# Predict labels for the test set and calculate accuracy
test_predictions = catboost_model.predict(X_test)
test_accuracy = accuracy_score(y_test, test_predictions)
print("Test set accuracy: ", test_accuracy)

# Calculate and print the confusion matrix and classification report for the test set
test_confusion_matrix = confusion_matrix(y_test, test_predictions)
test_classification_report = classification_report(y_test, test_predictions)

print("Test set confusion matrix: \n", test_confusion_matrix)
print("Test set classification report: \n", test_classification_report)

# Calculate False Negative Rate for each class
FN = test_confusion_matrix.sum(axis=1) - np.diag(test_confusion_matrix)
TP = np.diag(test_confusion_matrix)
FNR = FN / (TP + FN)

print("False Negative Rate : ", FNR)

# Calculate the total False Negative Rate
total_FNR = FN.sum() / (TP.sum() + FN.sum())
print("Total False Negative Rate : ", total_FNR)

# Calculate the mean accuracy of all classes (average sensitivity)
mean_accuracy = np.diag(test_confusion_matrix).sum() / test_confusion_matrix.sum()
print("Mean accuracy of all classes : ", mean_accuracy)

# Initialize the logistic regression model
# logistic_regression = LogisticRegression(C=10000, multi_class='ovr')
# xgb_model = xgb.XGBClassifier(objective='multi:softmax', num_class=33, seed=42)
# catboost_model = CatBoostClassifier(loss_function='MultiClass', classes_count=33, random_seed=42)
#
#
# # Perform 10-fold cross-validation on the training set
# kfold = KFold(n_splits=10, random_state=42, shuffle=True)
# cross_val_scores = cross_val_score(catboost_model, X_train, y_train, cv=kfold)
#
# print("Cross-validation scores: ", cross_val_scores)
# print("Average cross-validation score: ", cross_val_scores.mean())
#
# # Train the model
# # logistic_regression.fit(X_train, y_train)
# catboost_model.fit(X_train, y_train)
# # Predict labels for the training set and calculate accuracy
# # train_predictions = logistic_regression.predict(X_train)
# train_predictions = catboost_model.predict(X_train)
# train_accuracy = accuracy_score(y_train, train_predictions)
# print("Training set accuracy: ", train_accuracy)
#
# # Predict labels for the test set and calculate accuracy
# # test_predictions = logistic_regression.predict(X_test)
# test_predictions = catboost_model.predict(X_test)
# test_accuracy = accuracy_score(y_test, test_predictions)
# print("Test set accuracy: ", test_accuracy)
#
# # Calculate and print the confusion matrix and classification report for the test set
# test_confusion_matrix = confusion_matrix(y_test, test_predictions)
# test_classification_report = classification_report(y_test, test_predictions)
#
# print("Test set confusion matrix: \n", test_confusion_matrix)
# print("Test set classification report: \n", test_classification_report)
#
#
# nm=["ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH",
#      "KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
#      "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"]
#
#
# test_conf_matrix=confusion_matrix(y_test,test_predictions,labels=nm)
# fn=[]
# for i in range(len(nm)):
#     fn.append(sum(test_conf_matrix[i,:]) - test_conf_matrix[i,i])
#
# fnr = []
# for i in range(len(nm)):
#     fnr.append(fn[i] / sum(test_conf_matrix[i,:]))
# print(fnr)
# overall_fnr=sum(fn)/sum(sum(test_conf_matrix))
# print(overall_fnr)