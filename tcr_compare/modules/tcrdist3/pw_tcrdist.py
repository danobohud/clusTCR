import pandas as pd
import numpy as np
import networkx as nx
import os

from tcrdist.repertoire import TCRrep
from sklearn.cluster import DBSCAN, KMeans

from clustcr.clustering.clustering import ClusteringResult
from clustcr.input.vdjdb import parse_vdjdb

def normalize(edge):
    '''
    Introduce normalization property on an edge that is represented as a tuple.
    The normalization property will constraint the ordering of two nodes that
    make up an edge. This prevents duplicated edges.

    Parameters
    ----------
    edge : tuple
        Tuple of length two, containing two nodes, represented as integers.

    Returns
    -------
    (n1, n2) : tuple
        Normalized edge.
        
    '''
    n1, n2 = edge
    if n1 > n2:
        n1, n2 = n2, n1
    return (n1, n2)

def greedy_clustering(dm, threshold):
    '''
    Greedy graph clustering approach that uses a fixed distance-threshold to 
    assign nodes to cluster. The algorithm starts by computing all pairs 
    of sequences that satisfy the predefined distance threshold (edge list). 
    Next, it finds the sequence with the highest degree (i.e. the most neighbors), 
    assigns this as the cluster centre, and removes it and its neighbors 
    from the edge list. This process is repeated until all sequences are clustered.

    Parameters
    ----------
    dm : numpy.array
        Distance matrix.
    threshold : int
        Distance threshold for defining network edges.

    Returns
    -------
    res : pandas.DataFrame
        Dataframe containing clustering results.

    '''

    edges = np.argwhere(dm<=threshold)
    # print(len(edges))
    edges = set(map(normalize, edges)) # Remove duplicated edges
    edges = np.array(list(edges)) # Edgelist to array
    # print(len(edges))
    
    cid = 0
    res = pd.DataFrame()
    
    while len(edges) > 0:
        
        G = nx.from_edgelist(edges)
        degrees = pd.DataFrame(G.degree(), columns=['node', 'degree'])
        degrees = degrees.set_index('node')
        degrees = degrees.sort_values(by='degree', ascending=False)
        max_degree = degrees.idxmax().values
        
        cluster = edges[np.where(edges[:,0]==max_degree)]
        ids = np.unique(cluster)
        cids = [cid] * len(ids)

        if len(ids) <= 1:
            break
        
        cluster_iter = pd.DataFrame({'seq_id':ids,'cluster':cids})
        res = res.append(cluster_iter)
        
        for i in ids:
            edges = edges[np.where(edges[:,0]!=i)] # Remove from column 1
            edges = edges[np.where(edges[:,1]!=i)] # Remove from column 2
        
        cid += 1
            
    return res

def cluster_TCRDist_matrix(dm=None, cdr3=None, gt=None, method='KMeans'):
    '''
    Function for clustering of the TCRDist-calculated distance matrix.
    The function provides several methods for clustering, which include:
        - Greedy: fixed-distance threshold clustering method that groups
        of sequences that are connected in a graph.
        - DBSCAN: density-based clustering method that groups of densely
        packed points.
        - Kmeans: iterative clustering approach that partitions the data
        into k clusters.

    Parameters
    ----------
    dm : numpy.array, optional
        TCRDist-calculated distance matrix. The default is None.
    cdr3 : pandas.Series, optional
        pandas.Series containing the input CDR3 sequences. The default is None.
    gt : pandas.DataFrame, optional
        Ground truth data. The default is None.
    method : str, optional
        Clustering method used to cluster the TCRDist distance matrix. 
        Available methods include: greedy, DBSCAN and Kmeans.
        The default is 'DBSCAN'.

    Returns
    -------
    Clustering results

    '''
    # Available methods
    methods = ['GREEDY', 'DBSCAN', 'KMEANS']
    assert method.upper() in methods, r'Please choose one of the following: /n %s' % methods
    
    # If any of the parameters not specified, compute it using default settings
    if dm is None:
        dm, cdr3, gt = TCRDist(sparse=False)
    if gt is None:
        dm, cdr3, gt = TCRDist(sparse=False)
    if cdr3 is None:
        dm, cdr3, gt = TCRDist(sparse=False)
             
    if method.upper() == 'GREEDY':    
        # Greedy clustering
        clusters = greedy_clustering(dm, 12)
        clusters = clusters.rename(columns={'seq_id':'Index'})
        clusters = clusters.set_index('Index', drop=True)
        clusters = clusters.merge(right=cdr3, left_index=True, right_index=True)
        clusters = clusters.rename(columns={'cdr3_b_aa':'CDR3',
                                            'v_b_gene':'V'})
        metrics = ClusteringResult(clusters).metrics(gt)
        return clusters, metrics.summary()

    elif method.upper() == 'DBSCAN':
        # DBSCAN
        clustering = DBSCAN(eps=100, min_samples=2, n_jobs=-1).fit(dm)
        labels = clustering.labels_
        clusters = cdr3.rename(columns={'cdr3_b_aa':'CDR3',
                                        'v_b_gene':'V'})
        clusters['cluster'] = labels
        clusters = clusters[clusters['cluster']!=-1]
        metrics = ClusteringResult(clusters).metrics(gt)
        return clusters, metrics.summary()
    
    else:
        # K-Means

        kmeans = KMeans(n_clusters=500).fit(dm)
        labels = kmeans.labels_

        clusters = cdr3.rename(columns={'cdr3_b_aa':'CDR3',
                                        'v_b_gene':'V'})

        clusters['cluster'] = labels
        metrics = ClusteringResult(clusters).metrics(gt)
        return clusters, metrics.summary()