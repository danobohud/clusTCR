from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import subprocess, os, sys
from collections import Counter
from util_functions import get_time, write_lines
from clustcr.analysis.features import FeatureGenerator

global root, wdir

wdir=os.getcwd()
root='/'.join(os.getcwd().split('/')[:-1])
sys.path.append(root)


# from weblogo import *

def prepare_fasta(sequences,model,cluster,output_file):
    '''Convert a list of sequences to a temporary fasta file
        in preparation for MSA
        :param sequences: list of input CDR3 sequences
        :type sequences: list
        :param model: model name
        :type model: str
        :param cluster: cluster number
        :type cluster: str
        :param output_file: location of temporary fasta'''

    print('Converting %s sequences to fasta'%(len(sequences)))
    out=[]
    for i,sequence in enumerate(sequences):

        header = '{}_{:.0f}_{}'.format(i,float(cluster),model)
        seq=Seq(sequence)
        seq_rec = SeqRecord(seq, id = header , name='CDR3',
                        description='CDR3 sequence #%s of cluster %s from %s'%(i,cluster,model))
    
        out.append(seq_rec)
    SeqIO.write(out, output_file, "fasta")

def run_MSA(input,output,muscle_path=None):
    '''Run multiple sequence alignment with Muscle
    NB: Muscle needs to be downloaded and available in PATH
    Accessed from https://drive5.com/muscle5/
    :param input: fasta input path
    :type input: str
    :param output: msa output path:
    :type output: str
    :param muscle_path: path to muscle executable
    :type muscle_path: str
    '''

    print('Generating MSA')
    if not muscle_path:
        muscle_path='/opt/local/bin/muscle'
    # cmd = '/opt/local/bin/muscle -align {} -output {}'.format(input,output)
    cmd = '{} -in {} -out {}'.format(muscle_path,input,output)
    subprocess.call(cmd, shell=True, executable='/bin/zsh')

def get_weblogo(msa,out,start=2,format='png'):
    '''Produce Weblogo png image from MSA
    NB: weblogo .png production requires both Weblogo and Ghostscript

    :param msa: input MSA file
    :type msa: str
    :param out: output png destination
    :type out: str
    :param title: title for weblogo output
    :type title: str
    :param start: Amino acid start position, default 2 to clip N terminal Cys
    :type start: int'''

    cmd = 'weblogo -f {} -l {} -o {} -F {}'.format(msa,start,out,format)
    os.system(cmd)

def get_motifs(clusters,outdir,n=5,muscle_path=None):
    '''Generate weblogo motifs from input
    :param clusters: dictionary of cluster dataframes per method used
    :type clusters: dict
    :param outdir: directory for logos:
    type outdir: str
    :param n: number of logos to print
    :type n: int
    :param muscle_path: path to muscle executable
    :type muscle_path: str'''

    for method in clusters.keys():# Iterate over cluster results for each method
        df=clusters[method]
        print('*** Generating motifs for %s ***'%(method))
        input_fasta=os.path.join(wdir,'logos/seqs.fa')
        output_fasta = os.path.join(wdir,'logos/msa.fa')
        if n=='all':
            topn=df['cluster'].value_counts().index.tolist()        # Generate logos for all clusters
        else:
            topn=df['cluster'].value_counts().index[:n].tolist()    # Select top N
        for cluster in topn:                                    # Iterate over clusters
            logo = os.path.join(outdir,'%s_c%s.png'%(method,cluster))
            seqs=df[df['cluster']==cluster]['CDR3'].values.tolist() # Extract CDR3s
            prepare_fasta(seqs,method,str(cluster),input_fasta)     # Generate combiend fasta
            run_MSA(input_fasta,output_fasta,muscle_path)           # run MSA
            get_weblogo(output_fasta,logo)  # Create logo

    print('Complete')

def annotate(sequence,epitopes):
    '''Anotates each CDR3 with the most common epitope observed in that dataset
    :param sequence: input CDR3
    :type sequence: str
    :param epitopes: reference database
    :return: most common epitope, or "None"'''

    # to do: implement a database annotator e.g. TCRMatch
    eps=list(epitopes[epitopes['CDR3']==sequence]['Epitope'].values)
    try:
        return Counter(eps).most_common(1)[0][0]
    except IndexError:
        return 'None'

def make_edgelist(nodes,output_file):
    '''Create an edgelist for each node (sequence) in a graph
    :param nodes: database of sequences with cluster annotations
    :type nodes: DataFrame
    :param output_file: output text file:
    :type output_file: str'''

    out=[]
    for cl in nodes['cluster'].unique():
        seqs = nodes[nodes['cluster']==cl]['CDR3'].unique().tolist()
        edges = [[seqs[i],seqs[j]] for i in range(len(seqs)) for j in range(len(seqs)) if i!=j]
        out.extend(edges)
    write_lines(output_file,out)

def annotate_experimental(clusters,epitopes,outdir,parameters,pGen=True):
    '''Generate cytoscape graphs annotated with features and epitope binding data
    :param clusters: dictionary of cluster values
    :type clusters: dict
    :param epitopes: reference epitope database
    :type epitopes: DataFrame
    :param outdir: output file for edge and node lists
    :type outdir: str
    :param parameters: model parameters
    :type parameters: dict
    :return: cluster features per model type
    :rtype: dict'''

    epitopes['Epitope']=[x.strip('Test_') for x in epitopes['Epitope']] # Required to strip Test annotation in precision_recall analysis

    featuredict={}
    for key in clusters.keys(): # Iterate over models
        print('*** Generating edge list for %s ***'%(key))
        nodes=clusters[key]
        
        # To Do: record average length of cluster for every cluster
        print('Creating node list with cluster features')

        if pGen==True:
            print('NB: Annotating with pGen signficantly increases the time to load')

        # Annotate clusters with features
        analysis = FeatureGenerator(nodes) 
        analysis.get_features(compute_pgen=pGen)
        featuredict[key]=analysis
        analysis = FeatureGenerator(nodes) 
        analysis.get_features(compute_pgen=pGen)
        nodes=nodes.merge(epitopes[['CDR3','V','J','Epitope','subject:condition']])
        nodes['Length']=[len(cdr3) for cdr3 in nodes['CDR3'].values]

        split = [nodes.iloc[i]['Epitope'].split('_') for i in range(len(nodes))]
        nodes['HLA']=[x[0] for x in split]
        if '*' in nodes['HLA'].unique():
            nodes['HLA_super']=['None' if 'None' in x else x[0].split('*')[1].split(':')[0] for x in split]
        else:
            nodes['HLA_super'] = ['None' if 'None' in x else x[0][-2:] for x in split]

        # nums=Counter(nodes['Epitope'].values)
        # cnt = [c for c in nums.keys() if  nums[c]>parameters['min_clustsize']]
        # nodes=nodes[nodes['Epitope'].isin(cnt)]

        # Export networks
        print('Exporting network\n')
        make_edgelist(nodes, outdir+'/%s_%s_edgelist.txt'%(key,parameters['name']))
        nodes.to_csv(outdir+'/%s_%s_nodelist.txt'%(key,parameters['name']),index=False)
        print('Complete')
        
    return featuredict