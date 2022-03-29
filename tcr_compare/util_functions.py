
import pandas as pd
import csv, datetime, os
import scipy.stats as st
import numpy as np
import pickle

def write_lines(csv_file,listoflists,header=False):
    '''Write rows to a csv file
    :param csvfile: output csv filename
    :type csvfile: str:
    :param listoflists: data to be written to csv
    :type listoflists: list of lists
    :param header: boolean to enable writing of a single header line to a new csv
    :type header: bool'''

    
    with open(csv_file,'a') as f:
        writer=csv.writer(f)
        if header:
            writer.writerow(listoflists)
        else:
            print('Writing results to %s\n'%(csv_file))
            for row in listoflists:
                writer.writerow(row)

def get_confidence(data,return_mean=False):
    '''Computes mean and confidence intervals on a given array'''

    #Calculate the sample parameters
    conf = 0.95   # 95% CI given
    df = len(data)-1  # degree of freedom = sample size-1
    mu = np.mean(data)    #sample mean
    se = st.sem(data)  #sample standard error

    #create 95% confidence interval for the population mean
    cI = st.t.interval(alpha=conf, df=df, loc=mu, scale=se)


    if return_mean:
        return mu, cI
    else:
        return cI

def save_pickle(file,destination):
    print('Saving pickle to ', destination)
    with open(destination, 'wb') as f:
        pickle.dump(file,f,protocol=pickle.HIGHEST_PROTOCOL)
    print('Complete')


def load_pickle(loadfile):
    print('Loading pickle from', loadfile)
    with open(loadfile, 'rb') as f:
        obj = pickle.load(f)
    print('Complete')
    return obj

def get_time():
    return datetime.datetime.now().strftime('%Y%m%d_%H%M%S')

def make_resultsfile(path,parameters,pr=False):
    if not os.path.exists(path):
        if not pr:
            header=['Experiment',
                            'Model',
                            'N_chains',
                            'retention: actual',
                            'retention: baseline',
                            'retention: dif',
                            'purity: actual',
                            'purity: baseline',
                            'purity: dif',
                            'purity_90: actual',
                            'purity_90: baseline',
                            'purity_90: dif',
                            'consistency: actual',
                            'consistency: baseline',
                            'consistency: dif',
                            'runtime',
                            'Accuracy_train',
                            'Precision_train',
                            'Recall_train',
                            'F1_train',
                            'Precision_weighted_train',
                            'Recall_weighted_train',
                            'F1_weighted_train',
                            'TPR_norm_train',
                            'FPR_norm_train',
                            'Support_train',
                            'Accuracy_test',
                            'Precision_test',
                            'Recall_test',
                            'F1_test',
                            'Precision_weighted_test',
                            'Recall_weighted_test',
                            'F1_weighted_test',
                            'TPR_norm_test',
                            'FPR_norm_test',
                            'Support_test']
                            
            header = header+[x for x in parameters.keys() if x not in ['input_file',
                                                                        'results_file',
                                                                        'chain_selection',
                                                                        'model_selection',
                                                                        'muscle_path',
                                                                        'graphs',
                                                                        'save',
                                                                        'root',
                                                                        'wdir',
                                                                        'name',
                                                                        'experiment']]

        else:
            header= ['Experiment', 'Model',
                        'Train/Test', 'Epitope',
                        'Accuracy','Precision',
                        'Recall', 'F1', 
                        'Precision_weighted','Recall_weighted',
                        'F1_weighted',
                        'TPR_norm','FPR_norm','Support']

        
        os.system('touch {}'.format(path))
        write_lines(path,header,header=True)

def mu_conf(name,df_1d):
    mu = np.mean(df_1d)
    conf=st.t.interval(0.95, len(df_1d)-1, loc=mu, scale=st.sem(df_1d))
    out = '{}: {:.0f} +- {:.2f}'.format(name,mu,(conf[1]-conf[0])/2)
    print(out)
    return out, mu
    