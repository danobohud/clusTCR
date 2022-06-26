import pandas as pd
import numpy as np
import sys
sys.path.append('/Users/danhudson/Documents/Academic/Oxford/Oxford_DPhil/Rotation_Projects/Hashem/clusTCR')
from util_functions import write_lines
from sklearn.metrics import precision_score, recall_score, confusion_matrix
from collections import Counter


def get_targets(epitopes):
    """Find the most commonly occurring epitopes in the dataset"""
    return [x[0] for x in Counter(epitopes['Epitope'].values.tolist()).most_common(3)]


def relabel_test(epitopes, testsplit=0.2):
    """Generate random train and test splits
    :param epitopes: input dataset
    :type epitopes: DataFrame
    :param testsplit: percentage to retain as test
    :type testplit: float
    :return epitopes: relabelled dataset
    :rtype epitopes: DataFrame"""

    epitopes = epitopes.copy().sample(frac=1).reset_index(drop=True)    # Shuffle the dataframe
    start = int((1-testsplit) * len(epitopes))                          # Set train and test split points
    train = epitopes.iloc[:start]
    test = epitopes.iloc[start:]

    test['Epitope'] = ['Test_%s' % (x) for x in test['Epitope'].values]
    epitopes = pd.concat([train, test])

    return epitopes


def get_confmat(clusters, epitopes):
    """Produce a confusion matrix from clustered data and original inputs
    :param clusters: cluster outputs
    :type clusters: DataFrame
    :param epitopes: source dataset
    :type epitopes: DataFrame
    :return conf_mat: confusion matrix of epitope counts per cluster"""

    gt = epitopes[epitopes["CDR3"].isin(clusters["CDR3"])]
    gt = gt.drop_duplicates()
    gt = pd.merge(left=gt, right=clusters, on="CDR3")
    gt["count"] = 1
    conf_mat = pd.pivot_table(gt,
                              values='count',
                              index=gt["Epitope"],
                              columns=gt["cluster"],
                              aggfunc=np.sum,
                              fill_value=0)
    
    return conf_mat


def group_other(df, targets):
    """Group non-target epitopes under 'Other'
    :param df: input dataset
    :type df: DataFrame
    :param targets: top n epitopes
    :type targets: list
    :return df: reformatted dataframe"""

    other=df[~df.index.isin(targets)]
    other=pd.DataFrame(np.sum(other,axis=0)).T
    df=df[df.index.isin(targets)]
    other.index=['Other']
    df=pd.concat([df,other],axis=0)
    
    return df


def train_test(clusters, epitopes, targets):
    """Split cluster into train and test for performance analysis
    :param clusters: cluster outputs
    :type clusters: DataFrame
    :param epitopes: source dataset
    :type epitopes: DataFrame
    :param targets: top n epitopes
    :type targets: list
    """
    confmat = get_confmat(clusters,epitopes)
    subset = [c for c in confmat.index if 'Test' in c]
    test=confmat[confmat.index.isin(subset)]
    train = group_other(confmat[~confmat.index.isin(test.index)],targets)
    test.index = [t[5:] for t in test.index]
    test = group_other(test,targets)

    return train, test


def get_pred(train, test):
    """Get predicted value per cluster
    :param train: train dataset
    :param test: test dataset"""

    out = []
    for clust in train.columns:
        c = train[clust].values.tolist()        # Get counts per epitope
        TP = np.max(c)                          # Find the maximum count
        label=train.index[c.index(TP)]          # Find the corresponding epitope
        out.append(train[clust].index==label)   # Apply boolean mask if epitope==predicted max

    train_pred = pd.DataFrame(out).T
    train_pred.index = train.index

    out=[]

    for clust in test.columns:
        c = train[clust].values.tolist()            # Use train epitopes as a reference
        TP = np.max(c)                              # Find the maximum count
        label = train.index[c.index(TP)]              # Find the corresponding epitope
        out.append(test[clust].index==label)        # Apply boolean mask if epitope==predicted max

    test_pred = pd.DataFrame(out).T
    test_pred.index = test.index

    return train_pred, test_pred

def get_actual(df):
    """"Find nonzero count values
    :param df: confusion matrix
    :type df: DataFrame
    :return actual: Nonzero values"""

    out=[]
    for clust in df.columns:
        out.append(df[clust] != 0)

    actual = pd.DataFrame(out).T
    actual.index = df.index
    
    return actual


def get_metrics(df, actual, predicted):
    """"Compute performance metrics from actual vs. predicted score per epitope
    :param df: confusion matrix
    :type df: DataFrame
    :param actual: Nonzero values
    :type actual: DataFrame
    :param predicted: Predicted max values
    :type predicted: DataFrame
    :return outdict: scores per epitope
    :rtype outdict: dictionary"""

    outdict = {}
    for t in df.index:
        y = actual.loc[t].values
        y_pred = predicted.loc[t].values
        precision = precision_score(y, y_pred)
        precision_w = precision_score(y, y_pred,average='weighted')
        recall = recall_score(y, y_pred)
        recall_w = recall_score(y, y_pred, average='weighted')
        
        try:
            tn, fp, fn, tp = confusion_matrix(y, y_pred, normalize=None).ravel()

        except ValueError:
            tn, fp, fn, tp = [0, 0, 0, 0]

        outdict[t] = {x: y for (x,y) in [('Precision', precision),
                                         ('Precision_weighted', precision_w),
                                         ('Recall', recall),
                                         ('Recall_weighted', recall_w),
                                         ('TN', tn),
                                         ('FP', fp),
                                         ('FN', fn),
                                         ('TP', tp),
                                         ('Support', tn+tp+fn+fp),
                                         ('TPR_norm', 0 if (tp+fn == 0) | (tp == 0) else tp/(tp+fn)),
                                         ('FPR_norm', 0 if (tn+fp == 0) | (tn == 0) else tn/(tn+fp)),
                                         ('Accuracy', 0 if (tp+tn+fn+fp == 0)|(tp + tn == 0) else (tp+tn)/(tp+tn+fp+fn)),
                                         ('F1', (2*precision*recall)/(precision + recall)),
                                         ('F1_weighted', (2*precision_w*recall_w)/(precision_w + recall_w)),
                                         ]
                      }
    return outdict


def analyse(clusters, epitopes, targets, name):
    """Get performance metrics for train and test sets
    :param clusters: cluster outputs
    :type clusters: DataFrame
    :param epitopes: source dataset
    :type epitopes: DataFrame
    :param targets: top n epitopes
    :type targets: list
    :param name: model name
    :type name: str
    :return: dictionary of train and test results
    """
    train, test = train_test(clusters, epitopes, targets)
    train_pred, test_pred = get_pred(train, test)
    train_actual, test_actual = [get_actual(df) for df in [train, test]]
    train_metrics= get_metrics(train, train_actual, train_pred)
    test_metrics= get_metrics(test, test_actual, test_pred)
    
    return {'Model': name, 'Train_metrics': train_metrics, 'Test_metrics': test_metrics}


def write_pr(cluster_object, name, results_file, result):
    """Write performance by epitope to csv
    :param cluster_object: Cluster object
    :type cluster_object: Object of class Cluster
    :param name: model name
    :type name: str
    :param results_file: output csv file
    :type results_file: str"""

    m = ['Accuracy','Precision',
         'Recall', 'F1',
         'Precision_weighted','Recall_weighted',
         'F1_weighted','TPR_norm',
         'FPR_norm','Support']

    out=[]

    header=[cluster_object.params['experiment'],name]
    for sub in ['Train_metrics','Test_metrics']:
        for epitope in result[sub].keys():
            out.append(header+[sub]+[epitope]+[result[sub][epitope][metric] for metric in m])
    write_lines(results_file, out)
    print('Complete')


def get_precision_recall(result):
    """Get mean performance across epitopes

    :param result: performance results
    :type result: dict
    :return out: reformatted results
    :rtype out: array"""

    m = ['Accuracy', 'Precision',
         'Recall', 'F1',
         'Precision_weighted', 'Recall_weighted',
         'F1_weighted', 'TPR_norm',
         'FPR_norm', 'Support']
    out = []
    for subset in ['Train_metrics','Test_metrics']:
        eps = list(result[subset].keys())
        out.extend([np.mean([result[subset][epitope][metric] for epitope in eps]) for metric in m])

    return out
