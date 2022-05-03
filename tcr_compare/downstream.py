from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import subprocess, os, sys
from collections import Counter
from util_functions import write_lines
from clustcr.analysis.features import FeatureGenerator
import numpy as np

global root, wdir
wdir=os.getcwd()
root='/'.join(os.getcwd().split('/')[:-1])
sys.path.append(root)


def prepare_fasta(sequences, model, cluster, output_file):
    """Convert a list of sequences to a temporary fasta file
        in preparation for MSA
    :param sequences: list of input CDR3 sequences
    :type sequences: list
    :param model: model name
    :type model: str
    :param cluster: cluster number
    :type cluster: str
    :param output_file: location of temporary fasta"""

    print('Converting %s sequences to fasta' % (len(sequences)))
    out = []
    for i, sequence in enumerate(sequences):
        header = '{}_{:.0f}_{}'.format(i, float(cluster), model)
        seq = Seq(sequence)
        seq_rec = SeqRecord(seq,
                            id=header,
                            name='CDR3',
                            description='CDR3 sequence #%s of cluster %s from %s' % (i,
                                                                                     cluster,
                                                                                     model))
    
        out.append(seq_rec)
    SeqIO.write(out, output_file, "fasta")


def run_MSA(infile, output, muscle_path=None):
    """Run multiple sequence alignment with Muscle
    NB: Muscle needs to be downloaded and available in PATH
    Accessed from https://drive5.com/muscle5/
    :param infile: fasta input path
    :type infile: str
    :param output: msa output path:
    :type output: str
    :param muscle_path: path to muscle executable
    :type muscle_path: str"""

    print('Generating MSA')
    if not muscle_path:
        muscle_path='/opt/local/bin/muscle'

    cmd = '{} -align {} -output {}'.format(muscle_path, infile, output) # Muscle v3
    # cmd = '{} -in {} -out {}'.format(muscle_path,input,output)    # Muscle v5

    if sys.platform.lower() == 'darwin':
        subprocess.call(cmd, shell=True, executable='/bin/zsh')
    else:
        subprocess.call(cmd, shell=True, executable='/bin/bash')


def get_weblogo(msa, out, start=1, fmat='pdf'):
    """Produce Weblogo png image from MSA
    NB: weblogo .png production requires both Weblogo and Ghostscript: see README

    :param msa: input MSA file
    :type msa: str
    :param out: output png destination
    :type out: str
    :param start: Amino acid start position, default 2 to clip N terminal Cys
    :type start: int
    :param fmat: format
    :type fmat: str"""

    cmd = 'weblogo -f {} -l {} -o {} -F {}'.format(msa, start, out, fmat)
    os.system(cmd)


def get_motifs(clusters, outdir, n=5, muscle_path=None):
    """Generate weblogo motifs from input
    :param clusters: dictionary of cluster dataframes per method used
    :type clusters: dict
    :param outdir: directory for logos:
    type outdir: str
    :param n: number of logos to print
    :type n: int
    :param muscle_path: path to muscle executable
    :type muscle_path: str"""

    for method in clusters.keys():                      # Iterate over cluster results for each method
        df = clusters[method]
        if n > len(df['cluster'].unique()):
            raise KeyError('More logos selected than clusters')
        print('*** Generating motifs for %s ***' % method)
        input_fasta = os.path.join(wdir, 'logos/seqs.fa')
        output_fasta = os.path.join(wdir, 'logos/msa.fa')
        topn = df['cluster'].value_counts().index[:n].tolist()    # Select top N
        for cluster in topn:                                    # Iterate over clusters
            logo = os.path.join(outdir,'%s_c%s.pdf' % (method, cluster))
            seqs=df[df['cluster']== cluster]['CDR3'].values.tolist() # Extract CDR3s
            prepare_fasta(seqs, method, str(cluster), input_fasta)     # Generate combined fasta
            run_MSA(input_fasta, output_fasta, muscle_path)           # run MSA
            get_weblogo(output_fasta, logo)  # Create logo
    print('Complete')


def annotate(sequence, epitopes):
    """Annotate each CDR3 with the most common epitope observed in that dataset
    :param sequence: input CDR3
    :type sequence: str
    :param epitopes: reference database
    :return: most common epitope, or 'None'"""

    # to do: implement a database annotator e.g. TCRMatch
    eps = list(epitopes[epitopes['CDR3'] == sequence]['Epitope'].values)
    try:
        return Counter(eps).most_common(1)[0][0]
    except IndexError:
        return 'None'


def make_edgelist(nodes, output_file):
    """Create an edgelist for each node (sequence) in a graph
    :param nodes: database of sequences with cluster annotations
    :type nodes: DataFrame
    :param output_file: output text file:
    :type output_file: str"""

    out = []
    for cl in nodes['cluster'].unique():
        seqs = nodes[nodes['cluster'] == cl]['CDR3'].unique().tolist()
        edges = [[seqs[i], seqs[j]] for i in range(len(seqs)) for j in range(len(seqs)) if i != j]
        out.extend(edges)
    write_lines(output_file, out)


def annotate_experimental(clusters, epitopes, outdir, parameters, pGen=True):
    """Generate cytoscape graphs annotated with features and epitope binding data
    :param clusters: dictionary of cluster values
    :type clusters: dict
    :param epitopes: reference epitope database
    :type epitopes: DataFrame
    :param outdir: output file for edge and node lists
    :type outdir: str
    :param parameters: model parameters
    :type parameters: dict
    :param pGen: compute probability of generation
    :type pGen: bool
    :return: cluster features per model type
    :rtype: dict"""

    # Strip 'test' annotation in precision_recall analysis
    print(epitopes['Epitope'].unique())
    epitopes['Epitope'] = ['Orphan TCR' if type(x) ==float else x.strip('Test_') for x in epitopes['Epitope']]

    featuredict={}
    for key in clusters.keys(): # Iterate over models
        print('*** Generating edge list for %s ***' % (key))
        nodes = clusters[key]
        print('Creating node list with cluster features')

        if pGen:
            print('Annotating with pGen')
            print('NB: Annotating with pGen significantly increases the time to load')

        # Annotate clusters with features
        analysis = FeatureGenerator(nodes) 
        analysis.get_features(compute_pgen=pGen)
        motifs = analysis.clustermotif()
        nodes['motif'] = nodes['cluster'].map(motifs)
        featuredict[key] = nodes

        eps = epitopes.copy()
        mapdict = {}
        
        for col in eps.columns:
            if col not in nodes.columns:
                # Map features from input data to output clusters
                mapdict[col] = {cdr3: x for cdr3, x in zip(eps['CDR3'].values,eps[col].values)}
                nodes[col] = nodes['CDR3'].map(mapdict[col]).replace(np.nan,'NA').astype(str)
    
        nodes['Length']=[len(cdr3) for cdr3 in nodes['CDR3'].values]
        
        if ('_' in nodes['Epitope'].unique())&('HLA' in nodes['Epitope'].unique()):
            split = [nodes.iloc[i]['Epitope'].split('_') for i in range(len(nodes))]
            nodes['HLA'] = [x[0] for x in split]

        # Export networks
        print('Exporting network\n')
        make_edgelist(nodes, outdir+'/%s_%s_edgelist.txt' % (key, parameters['name']))
        nodes.to_csv(outdir+'/%s_%s_nodelist.txt' % (key, parameters['name']), index=False)
        print('Complete')
        
    return featuredict
