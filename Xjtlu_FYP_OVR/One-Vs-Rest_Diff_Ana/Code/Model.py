import pandas as pd

import numpy as np

from sklearn import linear_model

import math
from sklearn.model_selection import train_test_split, cross_val_score, cross_val_predict
from sklearn.metrics import accuracy_score, confusion_matrix
import pickle
np.set_printoptions(threshold=np.inf)
from sklearn.metrics import roc_auc_score, recall_score, precision_score, f1_score, confusion_matrix, precision_recall_curve, auc
from sklearn.metrics import precision_recall_curve


rank=['50','100','150','200','250','300']

for q in rank:
    df=pd.read_csv('H:/FYP/summer/Data/Chayi/Gene_matrix_'+q, index_col=[0], sep='\t')
    df=df.T
    x=df.drop(["Label"],axis=1)
    #这里可以更改列名只保留symbol
    # for k in x.columns:
    #     a = k.split('|')
    #     x.rename(columns={k:a[0]},inplace=True)

    df1=np.zeros((x.shape[0],x.shape[1]),dtype=np.float32)
    for a in range(0,x.shape[0]-1):
        for b in range(0,x.shape[1]):
            if float(x.iloc[a,b])<=1:
                x.iloc[a,b]=1
            df1[a,b]=math.log2(float(x.iloc[a,b]))

    y=df['Label']
    df2=pd.DataFrame(df1,columns=x.columns,index=x.index)

    #部分情况可能要对标签进行编码才能用
    X_train,X_test,y_train,y_test=train_test_split(df2,y,test_size=0.2,random_state=42)
    # lr_model=cb.CatBoostClassifier()
    # lr_model=linear_model.LogisticRegression(C=10000, multi_class='ovr',max_iter=1000)
    lr_model=linear_model.LogisticRegression(C=10000, multi_class='ovr', solver='sag',max_iter=2000)#另一种lrmodel的参数

    scores=cross_val_score(lr_model,X_train,y_train,cv=10)
    lr_model.fit(X_train,y_train)
    y_pred=lr_model.predict(X_test)
    y_pred_proba = lr_model.predict_proba(X_test)




    nm=["ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH",
         "KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
         "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"]


    test_conf_matrix=confusion_matrix(y_test,y_pred,labels=nm)

    #计算FNR值
    fn=[]
    for i in range(len(nm)):
        fn.append(sum(test_conf_matrix[i,:]) - test_conf_matrix[i,i])

    fnr = []
    for i in range(len(nm)):
        fnr.append(fn[i] / sum(test_conf_matrix[i,:]))

    overall_fnr=sum(fn)/sum(sum(test_conf_matrix))


    train_accuracy = accuracy_score(y_train, lr_model.predict(X_train))
    test_accuracy = accuracy_score(y_test,lr_model.predict(X_test))



    # Convert labels to integers if they are not already
    auc_score = roc_auc_score(y_test, y_pred_proba,
                              multi_class='ovr')  # Change to 'multiclass' for multinomial logistic regression

    # Calculate Sensitivity (Recall)
    recall1 = recall_score(y_test, y_pred, average='weighted')  # Choose the appropriate average method

    # Calculate F1-Score
    f1 = f1_score(y_test, y_pred, average='weighted')

    precision2 = precision_score(y_test, y_pred, average='weighted')

    #Assuming y_pred_proba contains the predicted probabilities for each class
    precision = dict()
    recall = dict()
    prc_auc = dict()

    for i, class_name in enumerate(nm):
        precision[class_name], recall[class_name], _ = precision_recall_curve(
            y_test == class_name, y_pred_proba[:, i]
        )
        prc_auc[class_name] = auc(recall[class_name], precision[class_name])

    # Calculate average PRC AUC
    average_prc_auc = sum(prc_auc.values()) / len(prc_auc)



    specificities = {}

    for i, class_name in enumerate(nm):
        tn = test_conf_matrix[i, :i].sum() + test_conf_matrix[i, i + 1:].sum()  # Calculate True Negatives
        fp = test_conf_matrix[:, i].sum() - test_conf_matrix[i, i]  # Calculate False Positives
        specificity = tn / (tn + fp)
        specificities[class_name] = specificity

    # Calculate average specificity
    valid_specificities = [value for value in specificities.values() if not np.isnan(value)]
    average_specificity = sum(valid_specificities) / 33


    print("Rank:",q)
    print("train_accuracy:",train_accuracy)
    print("test_accuracy", test_accuracy)
    print("fnr:",fnr)
    print("overall_fnr:",overall_fnr)
    print("AUC:", auc_score)
    print("Sensitivity:", recall1)
    print("F1-Score:", f1)
    print("Precision:", precision2)
    print(f"Average PRC AUC: {average_prc_auc}")
    print(f"Average Specificity: {average_specificity}")
    sensitivity_scores = recall_score(y_test, y_pred, average=None)
    print("sensitivity_scores:",sensitivity_scores)
    precision1 = precision_score(y_test, y_pred, average=None)
    print("recision:",precision1)

    # # # #模型保存
    # pickle.dump(lr_model,open('/home/haochun/OVR/Model/cn/model_'+q,'wb'))




#绘制Curve ROC
# from sklearn.metrics import roc_curve, auc
# import matplotlib.pyplot as plt
# from sklearn.metrics import roc_curve, plot_roc_curve
# import numpy as np
# from sklearn.preprocessing import label_binarize
# #Calculate ROC curve and AUC for each class
# y_test_binarized = label_binarize(y_test, classes=np.unique(y))
# fpr = dict()
# tpr = dict()
# roc_auc = dict()
# n_classes = y_test_binarized.shape[1]
# for i in range(n_classes):
#     fpr[i], tpr[i], _ = roc_curve(y_test_binarized[:, i], y_pred_proba[:, i])
#     roc_auc[i] = auc(fpr[i], tpr[i])
#
# # Compute micro-average ROC curve and AUC
# fpr["micro"], tpr["micro"], _ = roc_curve(y_test_binarized.ravel(), y_pred_proba.ravel())
# roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

#Plot ROC curve for each class and micro-average ROC curve
# import matplotlib.pyplot as plt
# plt.figure(dpi=350,figsize=(20,14))
# plt.plot(fpr["micro"], tpr["micro"],
#          label='micro-average ROC curve (area = {0:0.2f})'
#                ''.format(roc_auc["micro"]),
#          color='deeppink', linestyle=':', linewidth=4)
#
# colors=['yellowgreen','yellow','peru','dodgerblue','gold','hotpink','firebrick','cyan','slateblue','forestgreen','chocolate','lawngreen','sandybrown','lightpink','lightskyblue','red','turquoise','darkcyan','darkgoldenrod','moccasin','darkkhaki','seagreen','mediumslateblue','khaki','wheat','crimson','olive','darkolivegreen','mediumvioletred','aqua','slateblue','papayawhip','sienna']
# #colors = ['blue', 'red', 'green', 'orange', 'purple','blue', 'red', 'green', 'orange', 'purple','blue', 'red', 'green', 'orange', 'purple','blue', 'red', 'green', 'orange', 'purple','blue', 'red', 'green', 'orange', 'purple','blue', 'red', 'green', 'orange', 'purple','blue','red', 'green']
# for i, color in zip(range(n_classes), colors):
#     plt.plot(fpr[i], tpr[i], color=color, lw=2,
#              label='ROC curve of class {0} (area = {1:0.4f})'
#              ''.format(i, roc_auc[i]))
#
# plt.plot([0, 1], [0, 1], 'k--', lw=2)
# plt.xlim([-0.05, 1.0])
# plt.ylim([0.0, 1.05])
# plt.xlabel('False Positive Rate')
# plt.ylabel('True Positive Rate')
# plt.title('Receiver Operating Characteristic')
# plt.legend(loc="lower right")
# plt.savefig('H:/FYP/New_Data/200_ROC_DPI_600_4')
# #plt.show()
#
#





#SHAP
# explainer = shap.Explainer(lr_model.predict_proba, X_train)
# shap_values = explainer(X_test,max_evals=1500)
# shap.summary_plot(shap_values[:, :, 1], X_test,max_display=20)
#
#
# explainer = shap.Explainer(lr_model.predict_proba, X_train)
# shap_values = explainer(X_test,max_evals=1500)
# shap.summary_plot(shap_values[:, :, 0], X_test,max_display=20)
#
# explainer = shap.Explainer(lr_model.predict_proba, X_train)
# shap_values = explainer(X_test,max_evals=1500)
# shap.summary_plot(shap_values[:, :, 2], X_test,max_display=20)










#看数据及分布
# dic1={"ACC":0,"BLCA":0,"BRCA":0,"CESC":0,"CHOL":0,"COAD":0,"DLBC":0,"ESCA":0,"GBM":0,"GBMLGG":0,"HNSC":0,"KICH":0,
#     "KIRC":0,"KIRP":0,"LAML":0,"LGG":0,"LIHC":0,"LUAD":0,"LUSC":0,"MESO":0,"OV":0,"PAAD":0,"PCPG":0,"KIPAN":0,
#     "PRAD":0,"READ":0,"SARC":0,"SKCM":0,"STAD":0,"STES":0,"TGCT":0,"THCA":0,"THYM":0,"UCEC":0,"UCS":0,"UVM":0}
# for k in y:
#     dic1[k]+=1
# print(dic1)
#
# dic2={"ACC":0,"BLCA":0,"BRCA":0,"CESC":0,"CHOL":0,"COAD":0,"DLBC":0,"ESCA":0,"GBM":0,"GBMLGG":0,"HNSC":0,"KICH":0,
#     "KIRC":0,"KIRP":0,"LAML":0,"LGG":0,"LIHC":0,"LUAD":0,"LUSC":0,"MESO":0,"OV":0,"PAAD":0,"PCPG":0,"KIPAN":0,
#     "PRAD":0,"READ":0,"SARC":0,"SKCM":0,"STAD":0,"STES":0,"TGCT":0,"THCA":0,"THYM":0,"UCEC":0,"UCS":0,"UVM":0}
#
# for s in Y_train:
#     dic2[s]+=1
# print(dic2)