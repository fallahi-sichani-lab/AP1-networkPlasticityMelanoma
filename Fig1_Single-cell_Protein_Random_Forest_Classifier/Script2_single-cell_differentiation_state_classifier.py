import pandas as pd
import numpy as np

from sklearn import datasets
from sklearn import model_selection
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn import preprocessing
from sklearn.pipeline import make_pipeline
from sklearn.metrics import plot_confusion_matrix
from sklearn.preprocessing import StandardScaler
from sklearn import metrics
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn import metrics
from sklearn.metrics import accuracy_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import label_binarize
from sklearn.metrics import roc_auc_score

import shap
import random
from time import time

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D
from plotnine import *
import colorcet as cc
%matplotlib inline
from IPython.display import set_matplotlib_formats
set_matplotlib_formats('retina')

file_path = 'Sampled_Gated_Cells_by_Diff_States.csv'
df = pd.read_csv(file_path)
df.shape
df

AP1s = ["cFos","Phospho_cFos","Fra1","Phospho_Fra1","Fra2",
         "cJun","Phospho_cJun","JunB", "JunD",
         "Phospho_ATF1","ATF2","Phospho_ATF2","ATF3","ATF4","Phospho_ATF4","ATF5","ATF6"]
X = df[AP1s]
y = pd.Series(df["Diff.class"], dtype="category")
y = y.cat.reorder_categories(['M','T','N','U.NGFR_Low'])#Melanocytic, Transitory, Neural crest-like, Undifferentiated
y


for diff_class in y.cat.categories:
    
    dff = df[df['Diff.class']==diff_class]
    
    dff_test = dff.groupby('cellline').sample(frac = 0.20,replace=False, random_state = 100) #Taking 20% of the cells as test set for each cell line
    dff_train=dff.drop(dff_test.index)
    
    if diff_class==y.cat.categories[0]:
        df_test = dff_test
        df_train = dff_train
    else:
        df_test = pd.concat([df_test,dff_test])
        df_train = pd.concat([df_train,dff_train])
        
df_train.reset_index(drop=True, inplace=True)
df_test.reset_index(drop=True, inplace=True)

#Get count
count_train = df_train.groupby(['Diff.class','cellline']).size().reset_index(name='counts')
count_train['Diff.class'].astype("category").cat.reorder_categories(['M','T','N','U.NGFR_Low'])
count_test = df_test.groupby(['Diff.class','cellline']).size().reset_index(name='counts')
count_test['Diff.class'].astype("category").cat.reorder_categories(['M','T','N','U.NGFR_Low'])

#Plot cell line representations of training and testing data across differentiation states
plt.subplots(1,2,figsize=(22,7),dpi = 200)

cmap = cc.glasbey_category10
res = {df['cellline'].unique()[i]: cmap[i] for i in range(len(df['cellline'].unique()))}

plt.subplot(1,2,1)
ax = sns.histplot(count_train, x='Diff.class', hue='cellline', weights='counts',
             multiple='stack', palette=res, shrink=0.8,linewidth=0.5)
ax.set_ylabel('Count',fontsize = 20)
ax.set_xlabel('Differentiation State',fontsize = 20)
ax.set_title('Training Data',fontsize = 20)
ax.legend([],[], frameon=False)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)

plt.subplot(1,2,2)
ax = sns.histplot(count_test, x='Diff.class', hue='cellline', weights='counts',
             multiple='stack', palette=res, shrink=0.8,linewidth=0.5)
ax.set_ylabel('Count',fontsize = 20)
ax.set_xlabel('Differentiation State',fontsize = 20)
ax.set_title('Test Data',fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)

# Fix the legend so it's not on top of the bars.
legend = ax.get_legend()
legend.set_bbox_to_anchor((1,1))
legend.set_title('Cell Line',prop={'size':20})
ltext  = legend.get_texts()
plt.setp(ltext, fontsize='xx-large')

plt.savefig('Train_Test_Data_Cell_Line_Counts.pdf', format='pdf')

X_train = df_train[AP1s]
X_train.reset_index(drop=True, inplace=True)
y_train = pd.Series(df_train["Diff.class"], dtype="category")
y_train.reset_index(drop=True, inplace=True)
y_train = y_train.cat.reorder_categories(['M','T','N','U.NGFR_Low'])


X_test = df_test[AP1s]
X_test.reset_index(drop=True, inplace=True)
y_test = pd.Series(df_test["Diff.class"], dtype="category")
y_test.reset_index(drop=True, inplace=True)
y_test = y_test.cat.reorder_categories(['M','T','N','U.NGFR_Low'])



seed_id = 42
random.seed(seed_id)
cv = StratifiedShuffleSplit(n_splits = 10,test_size = 0.2,random_state=seed_id)

scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

max_depths = np.linspace(1, 50, 50, endpoint=True)
train_acc_results = []
test_acc_results = []
train_roc_results = []
test_roc_results = []

for max_depth in max_depths:
    rf = RandomForestClassifier(max_depth=max_depth, random_state = seed_id, n_jobs=-1)
    rf.fit(X_train_scaled, y_train)
    
    #accuracy
    train_acc = cross_val_score(rf, X_train_scaled, y_train, cv=cv,scoring= "accuracy")
    train_acc_results.append( train_acc.mean())
    
    test_acc = accuracy_score(y_test, rf.predict(X_test_scaled))
    test_acc_results.append(test_acc)
    
    #roc auc
    train_roc = cross_val_score(rf, X_train_scaled, y_train, cv=cv,scoring= "roc_auc_ovr")
    train_roc_results.append(train_roc.mean())
    
    y_pred_mat = label_binarize(rf.predict(X_test_scaled),classes = ["M","T","N","U.NGFR_Low"])
    y_test_mat = label_binarize(y_test,classes = ["M","T","N","U.NGFR_Low"])
    test_roc = roc_auc_score(y_test_mat, y_pred_mat,multi_class='ovr')
    test_roc_results.append(test_roc)
    



null = 1/len(y_train.unique())
line1, = plt.plot(max_depths, train_acc_results, '.-',color = 'xkcd:light blue', label="Training Accuracy")
line2, = plt.plot(max_depths, test_acc_results, '.-',color = 'xkcd:pink', label="Test Accuracy")
line3, = plt.plot(max_depths, train_roc_results, '.-',color = 'xkcd:azure', label="Training ROC AUC")
line4, = plt.plot(max_depths, test_roc_results, '.-',color = 'xkcd:tomato', label="Test ROC AUC")
line5, = plt.plot(max_depths, np.ones(len(max_depths),)*null,'--',color = "xkcd:grey", label = "Null")
#plt.legend(handler_map={line1: HandlerLine2D(numpoints=2,)})
plt.legend(loc='best', bbox_to_anchor=(1,1))
plt.ylabel('Model Performance')
plt.xlabel('Tree depth')
plt.ylim((0,1))

plt.savefig('RF_Optimize_max_depth.pdf', format='pdf')

max_depth_opt = 14 #optimized max_depth

n_estimators = [10,100,1000,5000,10000]
train_acc_results = []
test_acc_results = []
train_roc_results = []
test_roc_results = []
for n_estimator in n_estimators:
    rf = RandomForestClassifier(max_depth=max_depth_opt,
                                n_estimators = n_estimator,
                                random_state = seed_id,
                                n_jobs=-1)
    rf.fit(X_train_scaled, y_train)
    
    #accuracy
    train_acc = cross_val_score(rf, X_train_scaled, y_train, cv=cv,scoring= "accuracy")
    train_acc_results.append( train_acc.mean())
    
    test_acc = accuracy_score(y_test, rf.predict(X_test_scaled))
    test_acc_results.append(test_acc)
    
    #roc auc
    train_roc = cross_val_score(rf, X_train_scaled, y_train, cv=cv,scoring= "roc_auc_ovr")
    train_roc_results.append(train_roc.mean())
    
    y_pred_mat = label_binarize(rf.predict(X_test_scaled),classes = ["M","T","N","U.NGFR_Low"])
    y_test_mat = label_binarize(y_test,classes = ["M","T","N","U.NGFR_Low"])
    test_roc = roc_auc_score(y_test_mat, y_pred_mat,multi_class='ovr')
    test_roc_results.append(test_roc)


null = 1/len(y_train.unique())
line1, = plt.plot(np.log10(n_estimators), train_acc_results, '.-',color = 'xkcd:light blue', label="Training Accuracy")
line2, = plt.plot(np.log10(n_estimators), test_acc_results, '.-',color = 'xkcd:pink', label="Test Accuracy")
line3, = plt.plot(np.log10(n_estimators), train_roc_results,  '.-',color = 'xkcd:azure',  label="Training ROC AUC")
line4, = plt.plot(np.log10(n_estimators), test_roc_results, '.-',color = 'xkcd:tomato',  label="Test ROC AUC")
line5, = plt.plot(np.log10(n_estimators), np.ones(len(n_estimators),)*null,'--',color = "xkcd:grey",label = "Null")
#plt.legend(handler_map={line1: HandlerLine2D(numpoints=2,)})
plt.legend(loc='best', bbox_to_anchor=(1,1))
plt.ylabel('Model Performance')
plt.xlabel('n_estimators')
plt.ylim((0,1))

plt.savefig('RF_Optimize_n_estimators.pdf', format='pdf')

n_estimators_opt = 100 #optimized n_estimators

# z-score data
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# model training
rf = RandomForestClassifier(max_depth=max_depth_opt, #optimaized n_estimators (=100)
                            n_estimators = n_estimators_opt, #optimized max_depth (=14)
                            random_state = seed_id,
                            n_jobs=-1)

# cross validation performance on training data
train_auc = cross_val_score(rf, X_train_scaled, y_train, cv=cv,scoring= "roc_auc_ovr")
print("ROC AUC: %0.3f" % train_auc.mean())
train_accuracy = cross_val_score(rf, X_train_scaled, y_train, cv=cv,scoring= "accuracy")
print("Accuracy: %0.3f" %train_accuracy.mean())



rf.fit(X_train_scaled,y_train)

#Calculate performance
print("Accuracy: %0.3f" % accuracy_score(y_test, rf.predict(X_test_scaled)))

y_pred_mat = label_binarize(rf.predict(X_test_scaled),classes = ["M","T","N","U.NGFR_Low"])
y_test_mat = label_binarize(y_test,classes = ["M","T","N","U.NGFR_Low"])
test_scores = roc_auc_score(y_test_mat, y_pred_mat,multi_class='ovr')


print("ROC AUC: %0.3f" % test_scores)



#confusion matrix of test data prediction
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
metrics.plot_confusion_matrix(rf, 
                        X_test_scaled, 
                        y_test, 
                        ax=axes, 
                        labels=["M", "T","N",'U.NGFR_Low'])
    
plt.sca(axes)
plt.xticks(rotation=50,ha = "right")
plt.rc('font', size=14)    
plt.savefig('Test_Confusion_Mat_M-T-N-U.pdf', format='pdf')

df_X_test_scaled = pd.DataFrame(X_test_scaled,columns=AP1s)

# Create Tree Explainer object that can calculate shap values
explainer = shap.TreeExplainer(rf)
shap_values = explainer.shap_values(df_X_test_scaled)

plt.subplots(1,1,figsize=(6,6),dpi = 200)
shap.summary_plot(shap_values, df_X_test_scaled, plot_type="bar",class_names=rf.classes_,plot_size = None,show = False)

plt.savefig('RF_SHAP_Summary_Plot.pdf', format='pdf')

y_test.unique()


rf.classes_

#Save SHAP values for each differentiation state

df_shap_val = pd.DataFrame(shap_values[0],columns=df_X_test_scaled.columns) 

df_shap_val['SHAP_Score_Type'] = [rf.classes_[0]]*2000 #shap values for Melanocytic class

for j in np.array([2,1,3]): #shap values for Transitory, Neural crest-like and Undifferentiated class
    dff = pd.DataFrame(shap_values[j],columns=df_X_test_scaled.columns)
    
    dff['SHAP_Score_Type'] = [rf.classes_[j]]*2000
    df_shap_val = pd.concat([df_shap_val,dff])
    
df_shap_val.to_csv('RF_SHAP_Value_Output.csv')

df_shap_val

#Save test data
df_test.to_csv('RF_Test_Data.csv')

# Get ordered features
vals= np.abs(shap_values).mean(0)

feature_importance = pd.DataFrame(list(zip(df_X_test_scaled.columns, sum(vals))), columns=['col_name','feature_importance_vals'])
feature_importance.sort_values(by=['feature_importance_vals'], ascending=False,inplace=True)
feature_importance.col_name



# Local (single-cell) explanations by AP-1(ordered by mean absolute value of the SHAP values for each feature)
plt.subplots(1,4,figsize=(18,6),dpi = 100)
c=1
for i in np.array([0,2,1,3]):
    plt.subplot(1,4,c)
    plt.title(rf.classes_[i])
    shap.summary_plot(shap_values[i], df_X_test_scaled,show=False,plot_size = None)
    c+=1
plt.tight_layout()     
# plt.show()

plt.savefig('SHAP-local-importance_M-T-N-U.pdf', format='pdf')

#Get sorted list (by mean absolute shap values) of AP-1s for each class
nclass = 4
df_s= {}

c=0
for i in [0,2,1,3]:
    mean_abs_shap = np.abs(shap_values[i]).mean(0)
    
    df = pd.DataFrame(data = list(zip(AP1s, mean_abs_shap)), 
                    columns=['AP-1','mean abs shap'])
    df.sort_values(by= 'mean abs shap', ascending=False,inplace=True)
    df.reset_index(drop=False, inplace=True)
    df_s[c] = df
    c+=1
df_mas= pd.concat(df_s,axis = 1)
df_mas

scaler = StandardScaler()

test_scores = np.zeros([len(AP1s),3])
feature_list = []
for r in range(len(AP1s)):
    ml = df_mas[0,'AP-1'][0:r+1]
    tl = df_mas[1,'AP-1'][0:r+1]
    nl = df_mas[2,'AP-1'][0:r+1]
    ul = df_mas[3,'AP-1'][0:r+1]
    
    ap1_list = list(set().union(ml,tl,nl,ul))
    feature_list.append(ap1_list)
    #feature_list = [feature_list,ap1_list]
    print(ap1_list)
    X_train_scaled = scaler.fit_transform(X_train[ap1_list])
    X_test_scaled = scaler.transform(X_test[ap1_list])

    rf = RandomForestClassifier(max_depth=max_depth_opt,
                                n_estimators = n_estimators_opt,
                                random_state = seed_id,
                                n_jobs=-1)

    rf.fit(X_train_scaled,y_train)


    #Calculate performance
    test_accuracy = accuracy_score(y_test, rf.predict(X_test_scaled))
    
    print("%d features - Accuracy: %0.3f" % (r+1, test_accuracy))

    y_pred_mat = label_binarize(rf.predict(X_test_scaled),classes = ["M","T","N","U.NGFR_Low"])
    y_test_mat = label_binarize(y_test,classes = ["M","T","N","U.NGFR_Low"])
    test_auc = roc_auc_score(y_test_mat, y_pred_mat,multi_class='ovr')

    
    print("%d features - ROC AUC: %0.3f" % (r+1, test_auc))
    print("-----------------------")
    
    test_scores[r,0] = test_accuracy
    test_scores[r,1] = test_auc
    test_scores[r,2] = len(ap1_list)
df_test_scores = pd.DataFrame(data = test_scores, columns = ["Accuracy", "ROC AUC","Number of Features"])
df_test_scores["Features"] = feature_list

#Save performance
df_test_scores["Features"] = feature_list
df_test_scores.to_csv('RF_Test_Data_Performance.csv')
df_test_scores

g = sns.scatterplot(data = df_test_scores, x = df_test_scores.index+1, y = df_test_scores["Accuracy"],s = 50)
sns.lineplot(data = df_test_scores, x = df_test_scores.index+1, y = df_test_scores["Accuracy"],label = "Accuracy")
sns.scatterplot(data = df_test_scores, x = df_test_scores.index+1, y = df_test_scores["ROC AUC"],s = 50)
sns.lineplot(data = df_test_scores, x = df_test_scores.index+1, y = df_test_scores["ROC AUC"],label = "ROC AUC")
g.set(ylim=(0, 1),xlabel='Number of Top AP-1s', ylabel='Performance')
sns.set(font_scale=2)
sns.set_style("ticks")


plt.savefig('Model_Performance_Select_Top_AP-1.pdf', format='pdf')
