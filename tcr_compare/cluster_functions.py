from asyncio import open_unix_connection
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import numpy as np
import os, csv, time,sys, multiprocessing, subprocess
from collections import Counter
sys.path.append('/Users/danhudson/Documents/Academic/Oxford/Oxford_DPhil/Rotation_Projects/Hashem/clusTCR')
from clustcr import Clustering
from clustcr.modules.ismart.ismart import iSMART
from clustcr.clustering.clustering import ClusteringResult
from modules.gliph2.gliph2 import GLIPH2
from os.path import join, dirname, abspath
import matplotlib.pyplot as plt 
from clustcr.input.vdjdb import parse_vdjdb
from modules.tcrdist3.pw_tcrdist import cluster_TCRDist_matrix # Rerouted to local copy of pw_tcrdist for return of clusters
from tcrdist.repertoire import TCRrep
from sklearn.cluster import DBSCAN, KMeans
from util_functions import load_pickle
from sys import platform

def load_vdjdb(vdjdb_path):
    '''Parse and load VDJDB data
    :param vdjdb_path: path to vdjdb dataset
    :type vdjdb_path: str'''

    out=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]

    with open(vdjdb_path, 'r') as f:
        for line in f.readlines():
            l=line.split('\t')
            for i in range(14):
                out[i].append(l[i])
            out[-2].append(l[19])
            out[-1].append(l[23])

    df2=pd.DataFrame(out).T
    cols = df2.iloc[0].tolist()

    df2=df2.iloc[1:]
    df2.columns=cols
    df2=df2[df2['species']=='HomoSapiens']
    for col in ['cdr3.alpha','v.alpha','j.alpha',
            'cdr3.beta','v.beta','d.beta','j.beta',
            'mhc.a','mhc.b','mhc.class']:
            df2[col].replace('',np.nan,inplace=True)
    df2=df2.rename(columns={'meta.subject.id': 'subject',
                            'method.verification': 'condition'})
    df2['subject:condition']=[df2.iloc[i]['subject']+':'+'vdjdb' for i in range(len(df2))]
    df2['antigen.epitope']=[df2.iloc[i]['mhc.a']+'_'+df2.iloc[i]['antigen.epitope']+'_'+df2.iloc[i]['antigen.gene']+'_'+df2.iloc[i]['antigen.species'] for i in range(len(df2))]
    return df2


def prepare_chains(ClusterObject):
    '''Prepares chain data for clustering by dropping duplicates 
        and reformatting column labels
    :param chain_selection: select chain input from ['alpha','beta','paired']
    :type chain_selection: str'''
    
    epitopes=ClusterObject.data

    chain_selection = ClusterObject.chain_selection

    # Handle vdjdb labels to align entries with GLIPH2 requirements
    if ('subject' in list(epitopes.columns)) & ('subject:condition' not in list(epitopes.columns)): 
        subj = ['%s:unknown'%(x) for x in epitopes['subject'].values.tolist()]
        epitopes['subject:condition']=subj
    if 'antigen.epitope' in epitopes.columns:
        epitopes=epitopes.rename(columns={'antigen.epitope':'Epitope'})
    if 'count' in epitopes.columns:
        epitopes=epitopes.rename(columns={'count':'Count'})
    if not 'Epitope' in epitopes.columns:
        epitopes['Epitope']=['None']*len(epitopes)
        
    epitopes['Epitope']=epitopes['Epitope'].replace(np.nan,'None')
    
    if ClusterObject.params['single_epitopes']==True:

            # Drop TCRs associated with more than one epitope
            epitopes=ClusterObject.check_epitopes(epitopes)

    # Reformat columns according to chain selection

    if chain_selection == 'alpha':
        if 'count_alpha' in epitopes.columns:
            epitopes = epitopes.rename(columns={'count_alpha':'Count'})
        epitopes = epitopes.drop(labels=[x for x in epitopes.columns if x not in ['cdr3.alpha', 
                                                                                'v.alpha',
                                                                                'j.alpha',
                                                                                'subject:condition',
                                                                                'Count',
                                                                                'Epitope']],axis=1)
        epitopes = epitopes.rename(columns={'cdr3.alpha':'CDR3',
                                            'v.alpha':'V',
                                            'j.alpha': 'J',
                                            'antigen.epitope':'Epitope'})
                            


    elif chain_selection == 'beta':
        if 'count_beta' in epitopes.columns:
            epitopes = epitopes.rename(columns={'count_beta':'Count'})

        epitopes = epitopes.drop(labels=[x for x in epitopes.columns if x not in ['cdr3.beta',
                                                                                    'v.beta',
                                                                                    'j.beta',
                                                                                    'subject:condition',
                                                                                    'Count',
                                                                                    'Epitope']],axis=1)
        epitopes = epitopes.rename(columns={'cdr3.beta':'CDR3',
                                            'v.beta':'V',
                                            'j.beta': 'J'})


    else:
        cnt=[1]*len(epitopes)
        epitopes['Count']=cnt

        epitopes = epitopes.rename(columns={'cdr3.beta':'CDR3',
                                            'v.beta':'V',
                                            'j.beta': 'J',
                                            'cdr3.alpha':'CDR3_alpha',
                                            'v.alpha':'V_alpha',
                                            'j.alpha': 'J_alpha'})




    chain = epitopes.drop(labels='Epitope',axis=1)

    return chain, epitopes

def get_epitope(cdr3,ref):
    eps = ref[ref['CDR3']==cdr3]['Epitope'].values.tolist()
    return Counter(eps).most_common(1)[0][0]

def trim_clusters(clusters,cutoff):

        out=pd.DataFrame()
        for c in clusters['cluster'].unique():
            sub=clusters[clusters['cluster']==c]
            if len(sub)>=cutoff:
                out=pd.concat([out,sub])
        out.columns=clusters.columns
        return out

def clusTCR(ClusterObject):
    
    clusTCR_data=ClusterObject.chain.CDR3
    chain_selection=ClusterObject.chain_selection

    # Run clusTCR
    print('\n*** Clustering %s %s chains with clusTCR **'%(len(clusTCR_data),chain_selection))
    t0 = time.time() 
    model = Clustering(n_cpus=ClusterObject.cpus).fit(clusTCR_data) 
    t1 = time.time()    
    t = t1 - t0

    # Record results

    clusters = model.clusters_df.copy()

    if ClusterObject.params['min_clustsize']:
        cutoff = ClusterObject.params['min_clustsize']
        clusters = trim_clusters(clusters,cutoff)

    epitopes=ClusterObject.epitopes[['CDR3','Epitope']]

    output = ClusteringResult(clusters).metrics(epitopes).summary()

    return output, t, model.clusters_df

def gliph2(ClusterObject):

    chain_data=ClusterObject.chain
    chain_selection=ClusterObject.chain_selection
    # Reformat input for GLIPH2
    cols=['CDR3','V','J','CDR3_alpha','subject:condition','Count']
    for x in cols[2:-1]:
        if x not in chain_data.columns:
            chain_data[x]=['NA']*len(chain_data)
        if 'V' not in chain_data.columns:
            raise KeyError(("Variable gene information is required to run gliph2"))
        if 'subject:condition' not in chain_data.columns:
            chain_data['subject:condition']=['vdjdb:none']*len(chain_data)
        if 'Count' not in chain_data.columns:
            chain_data['Count']=np.ones(len(chain_data)).astype(int)
   
    chain_data=chain_data[cols]
    chain_data=chain_data.drop(labels=[x for x in chain_data.columns if x not in cols],axis=1)

    # Run GLIPH2

    # NB: .centos executable required for Linux Centos, .OSX for Mac
    print('\n*** Clustering %s %s chains with GLIPH2 **'%(len(chain_data),chain_selection))
    GLIPH2_PATH = os.path.join(ClusterObject.params['wdir'],'modules/gliph2/lib')
    os.chdir(GLIPH2_PATH)
    chain_data.to_csv(os.path.join(GLIPH2_PATH,'metarepertoire.txt'), index=False, header=False, sep='\t')
    t0 = time.time()

    if sys.platform.lower() == 'darwin':
        os.system('./irtools.osx -c parameters.txt')
    elif sys.platform.lower() == 'linux':
        os.system('./irtools.centos -c parameters.txt')
    else:
        raise SystemError('GLIPH2 can run on Mac (Darwin) or Linux systems. Windows users are directed to the GLIPH2 webtool, a link to which is provided in the readme, or to select a different model')
    t1 = time.time()
    t = t1 - t0

    # Reformat gliph2 clustering results
    
    clusters = pd.DataFrame()
    with open('metarepertoire_cluster.txt', 'r') as f:
        results = f.read().splitlines()
    c = 0
    for line in results:
        columns = line.split(' ')
        motif = columns[3]
        cluster = columns[4:]
        if len(cluster) >= 2:
            nodes = pd.DataFrame({'CDR3':cluster})
            nodes['cluster'] = c
            nodes['motif'] = motif
            clusters = clusters.append(nodes)
            c += 1
    
    # Record results
    epitopes=ClusterObject.epitopes

    if ClusterObject.params['min_clustsize']:
        cutoff = ClusterObject.params['min_clustsize']
        clusters = trim_clusters(clusters,cutoff)

    output = ClusteringResult(clusters).metrics(epitopes).summary()
    os.chdir(ClusterObject.params['wdir'])

    return output, t, clusters

def join_cdr3_v(df, data_type='clusters_df'):
    '''Combine CDR3 and V information for iSMART performance analysis
    :param data_type: data type description, default clusters_df
    :type data_type: str'''

    joint_id = df['CDR3'].astype(str) + '_' + df['V'].astype(str)

    df['TCR'] = joint_id
    df = df.drop(labels=['CDR3', 'V'],axis=1)
    df = df.rename(columns={'TCR':'CDR3'})
    if data_type == 'clusters_df':
        df = df[['CDR3', 'cluster']]
    else:
        df = df[['CDR3', 'Epitope']]

    return df

def GIANA(ClusterObject):

    '''Cluster chains with GIANA'''

    chain_selection=ClusterObject.chain_selection
    # Reformat input for GIANA
    GIANA_data = ClusterObject.chain.drop(labels=[x for x in ClusterObject.chain.columns if x not in ['CDR3','V']],axis=1)

    # Run GIANA
    print('\n*** Clustering %s %s chains with GIANA **'%(len(GIANA_data),chain_selection))
    
    GIANA_PATH=os.path.join(ClusterObject.params['wdir'],'modules/GIANA/')

    os.chdir(GIANA_PATH)
    GIANA_data.to_csv('input.txt', index=False, header=False, sep='\t')

    print('Clustering {} sequences with GIANA.'.format(len(GIANA_data)))

    # Perform GIANA algorithm on test sequences
    t0 = time.time()
    os.system('python GIANA4.1.py -f input.txt -O input_clustered.txt -v True -N {}'.format(ClusterObject.cpus))       # To Do: test inclusion of v genes, referencing TRVB or TRAV files
    t1 = time.time()
    t = t1 - t0

    print('Elapsed time: {} seconds.'.format(t))

    with open(os.path.join(GIANA_PATH,'input_clustered.txt'), 'r') as f:
        clusters = f.read().splitlines()[3:]
        clusters = pd.DataFrame([x.split('\t') for x in clusters], columns=['CDR3', 'cluster','V'])
    
    epitopes=ClusterObject.epitopes[['CDR3','Epitope']]

    if ClusterObject.params['min_clustsize']:
        cutoff = ClusterObject.params['min_clustsize']
        clusters = trim_clusters(clusters,cutoff)

    output = ClusteringResult(join_cdr3_v(clusters)).metrics(join_cdr3_v(ClusterObject.epitopes, 
                    data_type='epitope_data')).summary()

    os.chdir(ClusterObject.params['wdir'])

    return output, t, clusters

def hamming_hash(cdr3):
    print('Calculating Hamming distances')
    # Adapted from clusTCR edge list production function

    cdr3hash = dict()
    for cdr in cdr3:
        for hash in (cdr[::2], cdr[1::2]):
            if hash not in cdr3hash:
                cdr3hash[hash] = set()
            cdr3hash[hash].add(cdr)
    
    edgedict={}
    for i,hash in enumerate(cdr3hash):
        if len(cdr3hash[hash]) >= 1:
            for cdr1 in cdr3hash[hash]:
                if cdr1 not in edgedict.keys():
                    edgedict[cdr1]=i
                    for cdr2 in cdr3hash[hash]:
                        if cdr2 not in edgedict.keys():
                            if cdr1 != cdr2:
                                if cdr1 <= cdr2:
                                    if sum(ch1 != ch2 for ch1, ch2 in zip(cdr1, cdr2)) == 1:
                                        edgedict[cdr2]=i
    
    return edgedict

def hamming(ClusterObject):

    eps=ClusterObject.epitopes
    chain=ClusterObject.chain
    chain_selection=ClusterObject.chain_selection
        
    # Reformat input 
    if chain_selection=='paired':
        eps=eps[['CDR3','CDR3_alpha','Epitope']]
        alpha=eps[['CDR3_alpha','Epitope']].rename(columns={'CDR3_alpha':'CDR3'})
        beta=eps[['CDR3','Epitope']]
        eps=pd.concat([beta,alpha])
        chain=eps.drop(labels='Epitope',axis=1)
        # Run comparison

    print('\n*** Clustering %s %s chains on Hamming distance **'%(len(chain),chain_selection))
    t0 = time.time() # Initialise runtime

    clusters=chain.copy()
    cdr3=clusters['CDR3'].unique()
    hamdict = hamming_hash(cdr3)
    clusters=clusters[clusters['CDR3'].isin(list(hamdict.keys()))]
    clusters['cluster']=[hamdict[seq] for seq in clusters['CDR3'].values]
    
    if ClusterObject.params['min_clustsize']:
        cutoff = ClusterObject.params['min_clustsize']
        clusters = trim_clusters(clusters,cutoff)

    epitopes=ClusterObject.epitopes[['CDR3','Epitope']]
    output = ClusteringResult(clusters).metrics(epitopes).summary()
    t1 = time.time()
    t = t1 - t0
    
    return output, t, clusters

def ismart(ClusterObject):

    '''Cluster chains with iSMART
    :param chain_selection: select chain input from ['alpha','beta','paired']
    :type chain_selection: str'''

    chain_selection=ClusterObject.chain_selection
    # Reformat input for iSMART
    ismart_data = ClusterObject.chain.drop(labels=[x for x in ClusterObject.chain.columns if x not in ['CDR3','V']],axis=1)

    # Run iSMART
    print('\n*** Clustering %s %s chains with iSMART **'%(len(ismart_data),chain_selection))
    
    ISMART_PATH=os.path.join(ClusterObject.params['root'],'clustcr/modules/ismart/lib/')

    os.chdir(ISMART_PATH)
    ismart_data.to_csv('input.txt', index=False, header=False, sep='\t')

    print('Clustering {} sequences with iSMART.'.format(len(ismart_data)))

    # Perform iSMART algorithm on test sequences
    t0 = time.time()
    os.system('python iSMARTf3.py -f input.txt -v True -N {}'.format(ClusterObject.cpus))       
    t1 = time.time()
    t = t1 - t0

    print('Elapsed time: {} seconds.'.format(t))

    with open(ISMART_PATH+'/input_clustered_v3.txt', 'r') as f:
        clusters = f.read().splitlines()[3:]
        clusters = pd.DataFrame([x.split('\t') for x in clusters], columns=['CDR3', 'V','cluster'])
    
    epitopes=ClusterObject.epitopes[['CDR3','V','Epitope']]
    
    if ClusterObject.params['min_clustsize']:
        cutoff = ClusterObject.params['min_clustsize']
        clusters = trim_clusters(clusters,cutoff)

    output = ClusteringResult(join_cdr3_v(clusters)).metrics(join_cdr3_v(epitopes, 
                    data_type='epitope_data')).summary()

    os.chdir(ClusterObject.params['wdir'])
    return output, t, clusters


def length_cluster(ClusterObject):
    '''Cluster chains with length comparator'''

    eps=ClusterObject.epitopes
    chain=ClusterObject.chain
    chain_selection=ClusterObject.chain_selection
        
    # Reformat input for length comparison
    if chain_selection=='paired':
        eps=eps[['CDR3','CDR3_alpha','Epitope']]
        alpha=eps[['CDR3_alpha','Epitope']].rename(columns={'CDR3_alpha':'CDR3'})
        beta=eps[['CDR3','Epitope']]
        eps=pd.concat([beta,alpha])
        chain=eps.drop(labels='Epitope',axis=1)

    # Run comparison
    out=pd.DataFrame()
    print('\n*** Clustering %s %s chains on length comparison **'%(len(chain),chain_selection))
    t0 = time.time() # Initialise runtime
    cdr3=sorted(chain['CDR3'].unique())
    lengthdict = {'CDR3': cdr3, 'Length': [len(c) for c in cdr3]}
    # Trim small clusters
    clusters = pd.DataFrame.from_dict(lengthdict)
    cid = {l:idx for idx,l in enumerate(sorted(clusters['Length'].unique()))}
    clusters['cluster']=[cid[l] for l in clusters['Length'].values.tolist()]

    
    if ClusterObject.params['min_clustsize']:
        cutoff = ClusterObject.params['min_clustsize']
        clusters = trim_clusters(clusters,cutoff)

    epitopes=ClusterObject.epitopes[['CDR3','Epitope']]
    output = ClusteringResult(clusters).metrics(epitopes).summary()
    t1 = time.time()
    t = t1 - t0
    
    return output, t, clusters

def vgene_map(sequences,vgenedict,chain_selection):
        vls=[]

        for v in sequences:
            if v not in vgenedict.keys():
                if chain_selection=='alpha':
                    vls.append('TRAV12-2*01')
                else:
                    vls.append('TRBV19*01')
            else:
                vls.append(vgenedict[v])
        return vls
        
def tcrdist(ClusterObject):
    '''Cluster chains with tcrdist3
    :param chain_selection: select chain input from ['alpha','beta','paired']
    :type chain_selection: str
    :param sparse: analyse using tcrdist3 sparse matrix function
    :type sparse: bool'''

    tcrd_data=ClusterObject.epitopes
    chain_selection=ClusterObject.chain_selection

    if len(tcrd_data)>10000:
        sparse=True
    else:
        sparse=False

    # Reformat input for tcrdist3
    cdr3a = 'cdr3_a_aa'
    va= 'v_a_gene'
    ja= 'j_a_gene'
    cdr3b = 'cdr3_b_aa'
    vb= 'v_b_gene'
    ja= 'j_a_gene'

    if not 'count' in tcrd_data.columns:
        cnt = [1]*len(tcrd_data)
        tcrd_data['count']=cnt
        
    if 'antigen.epitope' in tcrd_data.columns:
        tcrd_data=tcrd_data.rename({'antigen.epitope': 'Epitope'})
    
    if ('*' not in tcrd_data['V'].values.tolist()) & (ClusterObject.params['input_file']!='vdjdb'):
        vdict=load_pickle(os.path.join(ClusterObject.params['wdir'],'data/vgene_dict.pkl'))

        if chain_selection in ['alpha','beta']:
            tcrd_data['V']=vgene_map(tcrd_data['V'].values, vdict, chain_selection)
        else:
            tcrd_data['V']=vgene_map(tcrd_data['V'].values, vdict, 'beta')
            tcrd_data['V_alpha']=vgene_map(tcrd_data['V'].values, vdict, 'alpha')
    
    if chain_selection =='alpha':
        tcrd_data=tcrd_data.rename(columns={'CDR3': cdr3a,
                                            'V': va})
        df_epi = tcrd_data[[cdr3a,va,'Epitope','count']]
        gt = df_epi.rename(columns = {cdr3a:'CDR3',
                                va:'V'})
        

    elif chain_selection =='beta':

            tcrd_data=tcrd_data.rename(columns={'CDR3': cdr3b,
                                                'V': vb})
            df_epi = tcrd_data[[cdr3b,vb,'Epitope','count']]
            gt = df_epi.rename(columns = {cdr3b:'CDR3',
                                    vb:'V'})
    else:
        tcrd_data=tcrd_data.rename(columns={'CDR3': cdr3b,
                                            'V': vb,
                                            'CDR3_alpha': cdr3a,
                                            'V_alpha': va,
                                            })
        
        
        df_epi = tcrd_data[[va,cdr3a,vb,cdr3b,'Epitope','count']]
        gt = df_epi.rename(columns = {cdr3b:'CDR3',
                                    vb:'V',
                                    'antigen.epitope':'Epitope'})


    seq = df_epi.drop(columns = ['Epitope'],axis=1).reset_index(drop=True)

    if chain_selection in ['alpha','beta']:
        chain=[chain_selection]
    else:
        chain=['alpha','beta']


    # Run tcrdist3

    print('\n*** Clustering %s %s chains with tcrdist3'%(len(seq),chain_selection))

    t0 = time.time()

    if sparse:
        print('Implementing sparse matrix distance computation')
        tr = TCRrep(
                    cell_df = seq,
                    # seqs1=df[cdr3].values,
                    organism = 'human',
                    chains = chain,
                    db_file = 'alphabeta_gammadelta_db.tsv',
                    compute_distances = False)
    
        tr.cpus = ClusterObject.cpus
        tr.compute_sparse_rect_distances(radius = 50, chunk_size = 500)
        if chain_selection =='alpha':
            S = tr.rw_alpha
            seq=seq.rename(columns={cdr3a:'CDR3',
                        va:'V'})
            
        elif chain_selection =='beta':
            S = tr.rw_beta
            seq=seq.rename(columns={cdr3b:'CDR3',
                                vb:'V'})
        
        else:
            S = tr.rw_beta
            seq=seq.rename(columns={cdr3b:'CDR3',
                                vb:'V'})

    else:
        tr = TCRrep(cell_df = seq,
                    organism = 'human',
                    chains = chain,
                    db_file = 'alphabeta_gammadelta_db.tsv',
                    compute_distances = True)

        tr.cpus = ClusterObject.cpus
        
        if chain_selection=='alpha':
            S = tr.pw_cdr3_a_aa
            seq=seq.rename(columns={cdr3a:'CDR3',
                        va:'V'})
        elif chain_selection in ['beta','paired']:
            S=tr.pw_cdr3_b_aa
            seq=seq.rename(columns={cdr3b:'CDR3',
                                vb:'V'})

    
    # Record results
    t1=time.time()
    print('Distance matrix calculated in %s seconds, clustering'%(t1-t0))

    clusters, output = cluster_TCRDist_matrix(S, seq, gt,method='KMeans')

    if ClusterObject.params['min_clustsize']:
        cutoff = ClusterObject.params['min_clustsize']
        clusters = trim_clusters(clusters,cutoff)

    t2 = time.time()
    t = t2 - t0
    
    return output, t, clusters


def tcrai_clusts(ClusterObject):
    path = os.path.join(ClusterObject.params['wdir'],'modules/TCRAI/')
    
    data=ClusterObject.epitopes.drop(labels='Epitope',axis=1)

    if ClusterObject.chain_selection=='alpha':
        data=data.rename(columns={'CDR3': 'CDR3_alpha'})
    else:
        data=data.rename(columns={'CDR3': 'CDR3_beta'})
    input_file=os.path.join(path,'TCRAI.csv')
    output_file=os.path.join(path,'TCRAI_result.csv')
    print('Writing data to: ',input_file)
    data.to_csv(input_file)
    # To Do: define automatic lookup of current conda envt, upload tcrai .yml into modules folder
    cmd = '. /Users/danhudson/opt/miniconda3/etc/profile.d/conda.sh && conda activate tcrai && python {}tcrai_implement.py -r {} -i {} -o {}'.format(path, path,
                                                                                                                                                input_file,
                                                                                                                                                output_file)
    print('Running TCRAI')
    subprocess.call(cmd, shell=True, executable='/bin/zsh')
    print('Reading results')
    tcrai_clusts= pd.read_csv(output_file)
    tcrai_clusts=tcrai_clusts.rename(columns={'TRB_cdr3':'CDR3',
                            'TRB_v_gene':'V',
                            'TRB_j_gene':'J',
                            'TRA_cdr3': 'CDR3_alpha',
                            'TRA_v_gene': 'V_alpha',
                            'TRA_j_gene':'J_alpha',
                            'pmhc_code':'Epitope'}).drop(labels='Unnamed: 0',axis=1)

    t=tcrai_clusts['Time'].unique()[0]
    clusters=tcrai_clusts[['CDR3','cluster']]
    
    if ClusterObject.params['min_clustsize']:
        cutoff = ClusterObject.params['min_clustsize']
        clusters = trim_clusters(clusters,cutoff)

    epitopes=ClusterObject.epitopes[['CDR3','Epitope']]
    output = ClusteringResult(clusters).metrics(epitopes).summary()
    
    return output, t, clusters