import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import os
import numpy as np

from cluster_functions import prepare_chains, load_vdjdb, get_epitope
from cluster_functions import clusTCR, GIANA, gliph2, ismart, hamming, length_cluster, tcrdist
from util_functions import write_lines, get_time, make_resultsfile
from precision_recall import relabel_test, get_targets, analyse, write_pr, get_precision_recall


def checkdir(parameters, directory):
    """Check directory structure to produce time-stamped results
    :param parameters: input parameters
    :type parameters: dictionary
    :param directory: name of directory:
    :type directory: str"""

    wdir = parameters['wdir']
    if directory not in os.listdir(wdir):
        # Create parent directory
        os.mkdir(os.path.join(wdir,directory))
    try:
        # Create a time and dataset-stamped folder for results
        os.mkdir(os.path.join(wdir, directory, parameters['experiment']))
    except FileExistsError("Folder already exists, renaming"):
        os.mkdir(os.path.join(wdir, directory, parameters['experiment'],'_1'))


def initialise(parameters):
    """Prepare a parameter dictionary for later clustering and annotation functionality
    :param parameters: input parameters
    :type parameters: dict
    :return: Dictionary of parameter values for Cluster object initialisation"""

    # Ensure directory structure is in place and user is in the correct working directory
    wdir, input_file = [parameters[x] for x in ['wdir', 'input_file']]
    print('Checking directories')
    
    if not input_file:
        input_file = 'vdjdb'

    if (input_file != 'vdjdb') & (not os.path.exists(input_file)):
        raise(FileNotFoundError("Input file not found in {}".format(input_file)))

    # Generate results folders
    if 'results' not in os.listdir(wdir):
        os.mkdir(os.path.join(wdir, 'results'))
    results_file = os.path.join(wdir, 'results/results.csv')
    results_file2 = os.path.join(wdir, 'results/results_pr.csv')
    make_resultsfile(results_file, parameters)
    make_resultsfile(results_file2, parameters, pr=True)
    
    for (key,value) in [('input_file', input_file),
                        ('name', input_file.split('/')[-1]),
                        ('results_file', results_file)]:
                        parameters[key] = value
    
    time = get_time()   # Generate timestamp
    parameters['experiment'] = time+'_%s' % (parameters['name'])

    for (name, direc) in [('n_logos', 'logos'),
                        ('graphs', 'cytoscape')]:

        if parameters[name]:
            checkdir(parameters, direc)

    return parameters


class Cluster:
    """Create a class object for running cluster analyses
    :param parameters: dictionary of input parameters defining experiments to be conducted
    :type parameters: dict"""
    
    def __init__(self, parameters):
        print('Initialising Cluster object\n')

        self.params = parameters                                # Read parameter dictionary to Cluster object
        self.chain_selection = self.params['chain_selection']
        self.clusters = dict()                                  # Initialise record of output clusters
        self.cpus = self.params['cpus']
        self.pr_results = dict()                                # Initalise accuracy results tracker
        self.results = pd.DataFrame()                           # Initialise cluster results tracker
        self.times = dict()                                     # Initialise runtime tracker
        self.chain = None
        self.epitopes = None
        self.targets = None
        self.N = None
        self.result_df = None
        self.alphabet = {v:k for k,v in enumerate(
            ['A', 'C', 'D','E','F',
            'G','H','I','K','L',
            'M','N','P','Q','R',
            'S','T','V','W','Y'])}

        print("Initialised with chain selection: ", self.chain_selection)

    def check_epitopes(self, data):
        """Drop TCRs associated with more than one Epitope in the database
            step 1: Find epitopes associated with a given CDR3
            step 2: Find CDR3s linked to more than one epitope
            step 3: Drop multi-epitope CDR3s
        :param data: TCR dataset
        :type data: DataFrame
        :return: cleaned dataset"""

        print('Dropping chains with more than one epitope assignment')
        l1 = len(data)

        if 'cdr3.alpha' in data.columns:
            ep_dict_alpha = {cdr3: data[data['cdr3.alpha'] == cdr3]['Epitope'].values for cdr3 in data['cdr3.alpha'].unique()}
            alpha = [k for k in ep_dict_alpha.keys() if len(ep_dict_alpha[k])>1]
            data=data[~data['cdr3.alpha'].isin(alpha)]

        if 'cdr3.beta' in data.columns:
            ep_dict_beta = {cdr3: data[data['cdr3.beta'] == cdr3]['Epitope'].values for cdr3 in data['cdr3.beta'].unique()}
            beta = [k for k in ep_dict_beta.keys() if len(ep_dict_beta[k]) > 1]
            
            data = data[~data['cdr3.beta'].isin(beta)]

        print('%s chains dropped' % (l1 - len(data)))

        return data

    def spike_in(self, dataset, n):
        """Spike in random sequences generated with OLGA
        :param dataset: input data
        :type dataset: pandas dataframe
        :param n: number of instances to add
        :type n: integer"""

        wdir = self.params['wdir']
        if (type(n) != int) | (n < 0) | (n > 1e7):
            raise TypeError("Please enter a spike-in value between 0 and 1000000")

        # Read in random TCRs
        TRA = pd.read_csv(os.path.join(wdir, 'modules/olga/TRA.csv'))
        TRA.columns = ['nt', 'cdr3.alpha', 'v.alpha', 'j.alpha']
        TRA = TRA.drop(columns=['nt']).iloc[:n]
        TRB = pd.read_csv(os.path.join(wdir, 'modules/olga/TRB.csv'))
        TRB.columns = ['nt', 'cdr3.beta', 'v.beta', 'j.beta']
        TRB = TRB.drop(columns=['nt']).iloc[:n]
        joint = pd.concat([TRA, TRB], axis=1)

        # Align formatting if allels are included in HLA descriptor
        if ('*' in dataset['v.alpha'].values.any()) | ('*' in dataset['v.beta'].values.any()):
            joint['v.alpha']=[x.split('/')[0]+'*01' for x in joint['v.alpha'].values]
            joint['v.beta']=[x.split('/')[0]+'*01' for x in joint['v.beta'].values]
        else:
            joint['v.alpha']=[x.split('/')[0] for x in joint['v.alpha'].values]
            joint['v.beta']=[x.split('/')[0] for x in joint['v.beta'].values]

        # Merge input and noise dataframes
        cols = [c for c in dataset.columns if c not in joint.columns]
        for c in cols:
            joint[c] = ['NA']*len(joint)
        joint = joint[dataset.columns]

        dataset = pd.concat([dataset, joint], axis=0).reset_index(drop=True)

        return dataset
    
    def annotate(self, epitopes):
        """Annotate CDR3s from a reference dataset based on exact match
        :param epitopes: input dataset
        :type epitopes: DataFrame"""

        ref = pd.read_csv('data/combined_cdr3_2.csv')
        print('Annotating missing epitopes from VDJDB, McPas, MIRA and GIANA reference databases')

        e = 0
        eps = []
        for i, cdr3 in enumerate(epitopes['CDR3'].values):
            if epitopes.iloc[i]['Epitope'] in ['None',np.nan]:
                if cdr3 in set(ref['CDR3'].values):
                    eps.append(get_epitope(cdr3, ref))
                    e += 1
                else:
                    eps.append('None')
            else:
                eps.append(epitopes.iloc[i]['Epitope'])
        
        epitopes['Epitope'] = eps
        print('Completed with %s sequences from the reference set' % (e))
        return epitopes

    def load_data(self):
        """Load database of TCR sequence and metadata"""

        vdjdb_on = False
        if self.params['input_file'].lower() == 'vdjdb':
            vdjdb_on = True
            # Parse VDJdb from clusTCR input files
            path = os.path.join(self.params['wdir'], 'data/vdjdb_trimmed_3.csv')
            data = pd.read_csv(path)
            print("VDJdb data loaded with {} instances".format(len(data)))

        else:
            print('Loading data from', self.params['input_file'])
            file_type = self.params['input_file'].split('.')[-1]

            if (file_type not in ['pkl','pickle','csv']) & (vdjdb_on == False):
                raise(ValueError("Input file must be .csv, .txt or pickled dataframe"))

            if file_type in ['pkl', 'pickle']:
                data=pd.read_pickle(self.params['input_file'])

            elif (file_type == 'csv')|(file_type == 'csv'):
                data=pd.read_csv(self.params['input_file'], index_col=0)

            else:
                raise(ValueError("File type %s is not supported, please use .csv, .txt or .pkl format" % (file_type)))
        
            print("Data loaded with {} instances".format(len(data)))

        if ('cdr3.beta' not in data.columns) & ('cdr3.alpha' not in data.columns):
            raise KeyError("Please ensure cdr3 columns are labelled as cdr3.alpha or cdr3.beta")
        
        if self.params['annotate']:
            print('Adding reference sequences for co-clustering analysis')
            # if self.params['chain_selection'] in ['alpha','paired']:
            reference = load_vdjdb(os.path.join(self.params['wdir'],'data/vdjdb_full.txt'))

            # else:
                # reference = pd.read_csv(os.path.join(self.params['wdir'],'data/vdjdb_mira.csv'))
            
            data = pd.concat([data, reference]).reset_index(drop=True)
            data['Epitope'].replace(np.nan,'Orphan TCR')
        
        original = len(data)
        
        if self.params['chain_selection']=='alpha':
            cols = ['cdr3.alpha', 'v.alpha','j.alpha']
        elif self.params['chain_selection']=='beta':
            cols = ['cdr3.beta', 'v.beta','j.beta']
        else:
            cols = ['cdr3.alpha', 'v.alpha','j.alpha',
                    'cdr3.beta', 'v.beta','j.beta']
        
        data = data.dropna(subset=cols)
        data=data.replace(np.nan,'NA')
        l2 = len(data)

        print('Dropping out-of-vocabulary values')

        cdr = [c for c in cols if 'cdr3' in c]
        for c in cdr:
            data=data[~data[c].apply(lambda x: len([aa for aa in x if aa not in self.alphabet])!=0)]
        l3 = len(data)


        data = data.drop_duplicates()

        print("Dropped {} NaNs, {} OOV and {} duplicates".format(original-l2, l2-l3, l3-len(data)))

        if self.params['spike_in']:
            # Add random instances produced with OLGA
            n = self.params['spike_in']
            print("Spiking in {} random TCRs".format(n))
            data = self.spike_in(data,n)

        # Add data to cluster object
        self.data = data

    def get_chain_data(self):
        """Prepares chain data for clustering by dropping duplicates
            and reformatting column labels"""

        chain, epitopes = prepare_chains(self)

        if self.params['annotate']:
            epitopes = self.annotate(epitopes)
            epitopes['Epitope'].replace(np.nan,'Orphan TCR',inplace=True)

        if bool(self.params['drop_nonbinders']):
            print('Dropping non-binders')
            orig = len(epitopes)    
            epitopes = epitopes[epitopes['Epitope'] != 'None'] # Include only TCRs with linked pMHC
            print('%s instances dropped with no epitope assignment' % (len(epitopes) - orig))

        self.chain = chain
        self.epitopes = relabel_test(epitopes)
        self.targets = get_targets(epitopes)
        self.N = len(epitopes)

    def run_model(self, modelfunc, name):
        """Run cluster models and record results
        :param modelfunc: model function
        :type modelfunc: function
        :param name: model name
        :type name: str"""

        if (('giana' in name)|('ismart' in name)) & ('paired' in name):
            raise KeyError("GIANA and iSMART do not currently accepted paired chain information")

        if ('tcrdist' in name) & (not self.params['single_epitopes']):
            raise KeyError('TCRDist3 should be run with TCRs binding only one epitope, set -se True in input command')

        output, t, clusters = modelfunc(self)
        output['method'] = name
        self.results = self.results.append(output)
        self.times[name] = t
        print('{} completed, elapsed time: {:.0f} seconds.'.format(name, t))
        return clusters

    def get_pr(self, cluster, name):
        """Get and record performance in epitope inference task
        :param cluster: cluster results
        :type cluster: DataFrame
        :param name: model name
        :type name: str"""

        self.clusters[name] = cluster
        pr_result = analyse(cluster,
                            self.epitopes[['CDR3', 'Epitope']],
                            self.targets,
                            name)
        write_pr(self, name, os.path.join(self.params['wdir'],
                                          'results/results_pr.csv'),
                 pr_result)
        self.pr_results[name] = get_precision_recall(pr_result)

    def cluster(self, model_selection):
        """Run cluster models from selection
        :param model_selection: models to evaluate
        :type model_selection: list"""

        for (model,func) in [('clustcr', clusTCR),
                             ('gliph2', gliph2),
                             ('hamming', hamming),
                             ('ismart', ismart),
                             ('giana', GIANA),
                             ('length', length_cluster),
                             ('tcrdist3', tcrdist)]:

            if model in model_selection:
                name = '_'.join([model, self.chain_selection])
                c=self.run_model(func, name)
                self.get_pr(c, name)

    def get_metrics(self):
        """Parse all results and parameters to read out to results file
        :return results: array of results for csv writer
        :rtype results: array"""

        expt = self.params['experiment']
        results=[]
        if 'method' not in self.results.columns:
            raise KeyError('No results detected. This may be because a model has been selected that cannot be used with the input dataset or chain selection. Please refer to the README for details')

        for method in self.results['method'].unique():
            # Extract results and pepare for csv write
            sub = self.results[self.results['method'] == method]
            out=[expt, method, self.N]
            for metric in ['retention', 'purity',
                           'purity_90', 'consistency']:
                sub2=sub[sub['metrics'] == metric]
                out.append(sub2['actual'].iloc[0])
                out.append(sub2['baseline'].iloc[0])
                out.append(sub2['actual'].iloc[0]-sub2['baseline'].iloc[0])
            out.append(self.times[method])
            out.extend(self.pr_results[method])

            # Record parameter values
            vals = [x for i,
                    x in enumerate(self.params.values()) if list(self.params.keys())[i] not in ['input_file',
                                                                                                'chain_selection',
                                                                                                'model_selection',
                                                                                                'muscle_path',
                                                                                                'graphs',
                                                                                                'root',
                                                                                                'wdir',
                                                                                                'name',
                                                                                                'experiment']]
            out.extend(vals)
            results.append(out)
        self.result_df = results
        return results
        
    def write_record(self, result_array):
        """Write results to csv
        :param result_array: results per model
        :type result_array: array"""

        print('Recording results and input parameters')
        write_lines(self.params['results_file'],
                    result_array,
                    header=False)