import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys, os, argparse
global root, wdir
wdir=os.getcwd()
root='/'.join(wdir.split('/')[:-1])
sys.path.append(root)
from cluster import initialise, Cluster
from downstream import get_motifs, annotate_experimental



def arg_check(args):
    """Check command line arguments for formatting or input errors
    :return: corrected arguments"""
    
    args.chain_selection = args.chain_selection.lower()
    args.model_selection = args.model_selection.lower()

    if args.chain_selection not in ['alpha','beta','paired','all']:
        raise KeyError("Please select a chain format from ['alpha','beta','paired']")
    
    elif args.chain_selection == 'all':
        args.chain_selection = ['alpha','beta','paired']

    else:
        args.chain_selection = [args.chain_selection]
    
    if str(args.model_selection) not in ['clustcr','giana','gliph2',
                                         'ismart','hamming','length',
                                         'tcrdist3','all']:
        raise KeyError("Please select a model from ['clustcr','giana','gliph2'",
                       ",'hamming','ismart','length','tcrdist3','all']")
    
    elif args.model_selection == 'all':
        args.model_selection = ['clustcr','giana','gliph2',
                              'hamming','ismart','length',
                              'tcrdist3']

    return args


def cluster_analysis(params, model_selection):
    """Runs clustering and records results for given
        parameter set and model_selection
    :param params: parameter dictionary
    :type params: dictionary
    :param model_selection: list of models
    :type model_selection: list
    :return: Cluster object"""

    clust = Cluster(params)             # Initialise
    clust.load_data()                   # Load in data
    clust.get_chain_data()              # Prepare chains for clustering
    clust.cluster(model_selection)      # Run clustering    
    results = clust.get_metrics()       # Format results
    clust.write_record(results)         # Write csv file with results of the analysis
    
    return clust


def run(params):
    """Run clustering and analysis
    :param params: parameter dictionary
    :type params: dict"""

    chain_selection = params['chain_selection']
    model_selection = params['model_selection']
    clusterdict = {}
    for chain in chain_selection:
        params['chain_selection'] = chain
        clusterdict[chain] = {}
        clust = cluster_analysis(params, model_selection) # Run clustering and read out metrics
        clusterdict[chain]['clusters'] = clust.clusters # Record clusters

        if params['n_logos']:
            # Run multiple sequence alignment on clusters and save weblogos
            get_motifs(clust.clusters,
                    outdir = os.path.join(wdir,'logos', params['experiment']),
                    n = params['n_logos'],
                    muscle_path = params['muscle_path'])

        if parameters['graphs']:
            # Generate features per cluster and save node/edge list for Cytoscape
            clusterdict[chain]['features'] = annotate_experimental(clusters=clust.clusters,
                                                                   epitopes=clust.epitopes,
                                                                   parameters=params,
                                                                   outdir=os.path.join(wdir, 'cytoscape',
                                                                                      params['experiment']),
                                                                   pGen=params['pGen'])


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Generate and annotate clusters from T Cell Receptor data')

    parser.add_argument('-i', '--input_file', type=str, required=False, default = 'vdjdb',
                        help='input file: .csv or .pkl')
    parser.add_argument('-m', '--model_selection', type=str, required=False, default = 'clustcr',
                        help='selection of models from ["clustcr","gliph2","ismart","length","tcrdist3","all"]')
    parser.add_argument('-c', '--chain_selection', required=False, type=str,default='beta',
                        help='chain selection from ["alpha","beta","paired","all"]')
    parser.add_argument('-sp', '--spike_in', type=int, required=False, default=None,
                        help='Add up to 1e7 random TCR sequences from OLGA as background noise')
    parser.add_argument('-dn', '--drop_nonbinders', type=bool, required=False, default=None,
                        help='Drop chains for which there is no Epitope information')
    parser.add_argument('-a', '--annotate', type=bool, required=False, default=None,
                        help='Annotate and co-cluster for epitope inference')
    parser.add_argument('-se', '--single_epitopes', type=bool, required=False, default=None,
                        help='Drop TCRs for which there is more than one epitope match')
    parser.add_argument('-re', '--repeats', type=int, required=False, default=1,
                        help='Repeat runs')
    parser.add_argument('-mc', '--min_clustsize', type=int, required=False, default=2,
                        help='Set the cluster size, default 2')
    parser.add_argument('-mp', '--muscle_path', type=str, required=False, default=None,
                        help='Add in path to muscle executable for WebLogo production')
    parser.add_argument('-n', '--n_logos', type=int, required=False, default=None,
                        help='Select number of logos to produce (top N clusters)')
    parser.add_argument('-pg', '--pgen', type=bool, required=False, default=False,
                        help='Add pGen to nodelist features')
    parser.add_argument('-cp', '--cpus', type=int, required=False, default=1,
                        help='Set number of CPUs, default 1 for cluster deployment')
    parser.add_argument('-g', '--graphs', type=bool, required=False, default=None,
                        help='Output node and edge lists')
    parser.add_argument('-co', '--comment', type=str, required=False, default=None,
                        help='Annotate experiment with a descriptive comment')

    args = parser.parse_args()
    wdir = os.getcwd()
    root = '/'.join(os.getcwd().split('/')[:-1])
    print('Setting root at ',root)
    print('Modules loaded, working directory: ', wdir)

    args = arg_check(args)
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
                  'Comment': args.comment}

    # Run
    for n in range(args.repeats):
        parameters['chain_selection'] = args.chain_selection
        run(initialise(parameters))