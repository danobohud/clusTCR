import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys,multiprocessing, os, argparse
from collections import Counter
global root, wdir

wdir=os.getcwd()
root='/'.join(wdir.split('/')[:-1])
sys.path.append(root)

from cluster import initialise, Cluster
from downstream import get_motifs, annotate_experimental

def arg_check(args):

    '''Check command line arguments for formatting or input errors
    :return: corrected arguments'''
    
    args.chain_selection=args.chain_selection.lower()
    args.model_selection=args.model_selection.lower()

    if args.chain_selection not in ['alpha','beta','paired','all']:
        raise KeyError("Please select a chain format from ['alpha','beta','paired']")
    
    elif args.chain_selection=='all':
        args.chain_selection=['alpha','beta','paired']

    else:
        args.chain_selection=[args.chain_selection]
    
    if str(args.model_selection) not in ['clustcr','giana','gliph2','ismart','hamming','length','tcrdist3','all']:
        raise KeyError("Please select a model from ['clustcr','giana','gliph2','hamming','ismart','length','tcrdist3','all']")
    
    elif args.model_selection=='all':
        args.model_selection=['clustcr','giana','gliph2','hamming','ismart','length','tcrdist3']
    

    return args

def cluster_analysis(parameters,model_selection):
    '''Runs clustering and records results for given 
        parameter set and model_selection
    :param parameter_dict: parameter dictionary
    :type parameter_dict: dictionary
    :param model_selection: list of models
    :type model_selection: list
    :return: Cluster object and clusters produced per model
    :rtype: Cluster object, dictionary'''

    clust = Cluster(parameters)         # Initialise
    clust.load_data()                   # Load in data
    clust.get_chain_data()              # Prepare chains for clustering
    clust.cluster(model_selection)      # Run clustering    
    results=clust.get_metrics()         # Format results
    clust.write_record(results)         # Write csv file with results of the analyis
    
    return clust

def run(parameters,n):
        
    '''Run clustering and analysis
    :param parameters: parameter dictionary:
    :type parameters: dict
    :param chain_selection: chains with which to run analysis
    :type chain_selection: list
    :param model_selection: models with which to run analysis
    :type model_selection: list
    :param graphs: include node and edge lists in output
    :type graphs: bool
    :param save: save cluster dictionary in .pkl
    :type save: bool

    :return: cluster dictionary, Cluster object
    '''
    chain_selection = parameters['chain_selection']
    model_selection= parameters['model_selection']

    clusterdict={}

    for chain in chain_selection:
        parameters['chain_selection']=chain
        clusterdict[chain]={}
        clust = cluster_analysis(parameters,model_selection) # Run clustering and read out metrics
        clusterdict[chain]['clusters']=clust.clusters # Record clusters

        if parameters['n_logos']:
            # Run mutliple sequence alignment on clusters and save weblogos
            get_motifs(clust.clusters,
                    outdir=os.path.join(wdir,'logos',parameters['experiment']),
                    n=parameters['n_logos'],
                    muscle_path=parameters['muscle_path'])

        if parameters['graphs']==True:
            # Generate features per cluster and save node/edgle list for Cytoscape
            clusterdict[chain]['features'] = annotate_experimental(clusters=clust.clusters,
                                                                epitopes=clust.epitopes,
                                                                parameters=parameters,
                                                                outdir=os.path.join(wdir,'cytoscape',parameters['experiment']),
                                                                pGen=parameters['pGen'])

if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Generate and annotate clusters from T Cell Receptor data')

    parser.add_argument('-i', '--input_file', type=str, required=False, default='vdjdb',
                        help='input file: .csv or .pkl')
    parser.add_argument('-m', '--model_selection', type=str,required=False, default = 'clustcr',
                        help='selection of models from ["clustcr","gliph2","ismart","length","tcrdist3","all"]')
    parser.add_argument('-c', '--chain_selection', required=False, type=str,default='beta',
                        help= 'chain selection from ["alpha","beta","paired","all"]')
    parser.add_argument('-sp', '--spike_in', type=int,required=False, default=None,
                        help='Add up to 1e7 random TCR sequences from OLGA as background noise')
    parser.add_argument('-dn', '--drop_nonbinders', type=bool,required=False, default=None,
                        help='Drop chains for which there is no Antigen or Epitope information')
    parser.add_argument('-a', '--annotate', type=bool,required=False, default=None,
                        help='Annotate and cocluster for epitope inference')                    
    parser.add_argument('-se', '--single_epitopes', type=bool,required=False, default=None,
                        help='Drop TCRs for which there is more than one epitope match')
    parser.add_argument('-re', '--repeats', type=int,required=False, default=1,
                        help='Repeat runs')
    parser.add_argument('-mc', '--min_clustsize', type=int,required=False, default=3,
                        help='Set the cluster size, default 3')
    parser.add_argument('-mp', '--muscle_path', type=str,required=False, default=None,
                        help='Add in path to muscle executable for WebLogo production')
    parser.add_argument('-n', '--n_logos', type=int,required=False, default=None,
                        help='Select number of logos to produce (top N clusters), or "all". NB: selecting all may generate thousands of logos, care is advised to avoid kernel timeout')
    parser.add_argument('-pg', '--pgen', type=bool,required=False, default=False,
                        help='Add pGen to nodelist features. NB: This may significantly increase runtime for node list production')
    parser.add_argument('-cp', '--cpus', type=int,required=False, default=1,
                        help='Set number of CPUs, default 1 for cluster deployment')
    parser.add_argument('-g', '--graphs', type=bool,required=False, default=None,
                        help='Output nodelist and edge list to cytoscape')
    parser.add_argument('-co', '--comment', type=str,required=False, default=None,
                        help='Annotate experiment with a decriptive comment')

    args = parser.parse_args()

    wdir =os.getcwd()
    root = '/'.join(os.getcwd().split('/')[:-1])
    print('Setting root at ',root)
    print('Modules loaded, working directory: ',wdir)

    args=arg_check(args)

    parameters = {'wdir': wdir,
            'root': root, 
        'chain_selection': args.chain_selection,
        'model_selection': args.model_selection,
        'input_file': args.input_file,
        'drop_nonbinders': args.drop_nonbinders,
        'single_epitopes': args.single_epitopes,
        'spike_in': args.spike_in,
        'annotate': args.annotate,
        'min_clustsize': args.min_clustsize,
        'muscle_path': args.muscle_path,
        'n_logos': args.n_logos,
        'pGen': args.pgen,
        'cpus': args.cpus,  
        'graphs': args.graphs,
        'Comment':args.comment}

    # Run
    for n in range(args.repeats):
        parameters['chain_selection']=args.chain_selection
        run(initialise(parameters),n)