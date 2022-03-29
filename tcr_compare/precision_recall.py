import pandas as pd
import numpy as np
import sys
sys.path.append('/Users/danhudson/Documents/Academic/Oxford/Oxford_DPhil/Rotation_Projects/Hashem/clusTCR')
from util_functions import write_lines
from sklearn.metrics import precision_score,recall_score,confusion_matrix

from collections import Counter

def get_targets(epitopes):
    return [x[0] for x in Counter(epitopes['Epitope'].values.tolist()).most_common(3)]

def relabel_test(epitopes,testsplit=0.2):
    '''Extract one donor and set as the test set
        for later retrieval
    :param ClusterObject: Cluster object with preloaded chains
    :type ClusterObject: Cluster class object
    :param targets: target epitopes, others will be grouped
    :type targets: list'''
    

    epitopes=epitopes.copy().sample(frac=1).reset_index(drop=True)
    start=int((1-testsplit)*len(epitopes))
    train=epitopes.iloc[:start]
    test=epitopes.iloc[start:]

    test['Epitope']=['Test_%s'%(x) for x in test['Epitope'].values]
    epitopes=pd.concat([train,test])

    return epitopes

def get_confmat(clusters,epitopes):

    gt = epitopes[epitopes["CDR3"].isin(clusters["CDR3"])]

    gt = gt.drop_duplicates()
    gt = pd.merge(left=gt, right=clusters, on="CDR3")
    gt["count"] = 1
    conf_mat = pd.pivot_table(gt, values='count',
                                index=gt["Epitope"],
                                columns=gt["cluster"],
                                aggfunc=np.sum,
                                fill_value=0)
    
    return conf_mat

def group_other(df,targets):

    other=df[~df.index.isin(targets)]
    other=pd.DataFrame(np.sum(other,axis=0)).T
    df=df[df.index.isin(targets)]
    other.index=['Other']
    df=pd.concat([df,other],axis=0)
    
    return df

def train_test(clusters, epitopes,targets):
    confmat = get_confmat(clusters,epitopes)
    subset = [c for c in confmat.index if 'Test' in c]
    test=confmat[confmat.index.isin(subset)]
    train = group_other(confmat[~confmat.index.isin(test.index)],targets)
    test.index = [t[5:] for t in test.index]
    test = group_other(test,targets)

    return train, test

def get_pred(train,test):
    out=[]
    for clust in train.columns:
        c = train[clust].values.tolist()
        TP = np.max(c)
        label=train.index[c.index(TP)]
        out.append(train[clust].index==label)

    train_pred=pd.DataFrame(out).T
    train_pred.index=train.index
    out=[]

    for clust in test.columns:
        c = train[clust].values.tolist()  # Use train as a reference

        TP = np.max(c)
        label=train.index[c.index(TP)]
        out.append(test[clust].index==label)

    test_pred=pd.DataFrame(out).T
    test_pred.index=test.index

    return train_pred, test_pred

def get_actual(df):

    out=[]
    for clust in df.columns:
        out.append(df[clust]!=0)

    actual=pd.DataFrame(out).T
    actual.index=df.index
    
    return actual

def get_metrics(df,actual,pred):
    outdict={}

    for t in df.index:
        y =actual.loc[t].values
        y_pred=pred.loc[t].values      
        precision = precision_score(y,y_pred)
        precision_w = precision_score(y,y_pred,average='weighted')
        recall = recall_score(y,y_pred)
        recall_w = recall_score(y,y_pred, average='weighted')
        
        try:
            tn, fp, fn, tp= confusion_matrix(y,y_pred,normalize=None).ravel()

        except ValueError:
            tn, fp, fn, tp = [0,0,0,0]

        outdict[t] = {x:y for (x,y) in [('Precision', precision),
                                        ('Precision_weighted', precision_w),
                                        ('Recall', recall),
                                        ('Recall_weighted', recall_w),
                                        ('TN',tn),
                                        ('FP',fp),
                                        ('FN',fn),
                                        ('TP',tp),
                                        ('Support',tn+tp+fn+fp),
                                        ('TPR_norm',0 if (tp+fn==0) | (tp ==0) else tp/(tp+fn)),
                                        ('FPR_norm',0 if (tn+fp==0) | (tn ==0) else tn/(tn+fp)),
                                        ('Accuracy',0 if (tp+tn+fn+fp==0) | (tp + tn ==0) else (tp+tn)/(tp+tn+fp+fn)),
                                        ('F1',(2*precision*recall)/(precision + recall)),
                                        ('F1_weighted',(2*precision_w*recall_w)/(precision_w + recall_w)),
                                        ]}
    
                                    
    return outdict

def analyse(clusters,epitopes,targets,name):
    
    train, test = train_test(clusters,epitopes,targets)
    train_pred, test_pred = get_pred(train, test)
    train_actual, test_actual = [get_actual(df) for df in [train,test]]
    train_metrics= get_metrics(train,train_actual,train_pred)
    test_metrics= get_metrics(test,test_actual,test_pred)
    
    return {'Model': name, 'Train_metrics': train_metrics, 'Test_metrics': test_metrics}

def write_pr(ClusterObject,name,results_file,result):

    m = ['Accuracy','Precision','Recall', 'F1', 
         'Precision_weighted','Recall_weighted','F1_weighted',
            'TPR_norm','FPR_norm','Support']
    out=[]

    header=[ClusterObject.params['experiment'],name]
    for sub in ['Train_metrics','Test_metrics']:
        for epitope in result[sub].keys():
            out.append(header+[sub]+[epitope]+[result[sub][epitope][metric] for metric in m])
    write_lines(results_file, out)
    print('Complete')

def get_precision_recall(result):
    m = ['Accuracy','Precision','Recall', 'F1', 
         'Precision_weighted','Recall_weighted','F1_weighted',
            'TPR_norm','FPR_norm','Support']
    out=[]
    for subset in ['Train_metrics','Test_metrics']:
        eps = list(result[subset].keys())
        out.extend([np.mean([result[subset][epitope][metric] for epitope in eps]) for metric in m])

    return out