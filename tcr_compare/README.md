*** Overview *** 

The TCR compare package has been developed as a sub module of the excellent ClusTCR platform developed by Sebastiaan Valkiers et al. (2021), to enable like-for-like comparison across a series of clustering algorithms accessible from the command line. 

The main amendments are as follows, all accessible via a command line interface (main.py):
- Integration of a suite of clustering models to allow simulataneous execution and comparison of clustering results across ClusTCR [1], GLIPH2[2], GIANA[3], iSMART[4] and TCRDist3[5]. Additional Hamming and CDR3 length approaches are also included for benchmarking purposes
- Annotation of orphan TCRs with exact matches from VDJdb
- Co-clustering analysis of orphan / user input data with reference VDJdb data
- Addition of measures of predictive performance in epitope specifcity inference  (Accuracy, Precision, Recall and F1-score) for datasets of experimentally detemined TCR-epitope pairs
- Automated recording of clustering results for ease of comparison
- Incorporation of alpha, beta or paired chain selections across models
- Export of Cytoscape network graphs including any feature from the input CSV, to permit customised feature analysis for private datasets

The core functionality of ClusTCR is unchanged from the original package, and so users are directed to the ClusTCR documentation for answers to questions on ClusTCR feature extraction and graph plotting, available at https://svalkiers.github.io/clusTCR/.

*** Installation ***

To install and run the tcr_compare module:

git clone https://github.com/danobohud/clusTCR

conda env create -f clusTCR_mod.yml

*** Running the package ***

Navigate to the ClusTCR package on your local machine
cd clusTCR/tcr_compare
conda activate clusTCR_mod

To view the command line syntax and defaults:
python main.py -h

To run:

i) an example clustering analysis on VDJdb with a single model and beta chain selection

python main.py -m clusTCR -c beta

ii) an example clustering analysis with example 10X data for all models and paired chain selections

python main.py -m all -c paired -i data/10X_dexnorm.csv

iii) with output of weblogo motifs for the top 5 largest clusters, and cytoscape network graphs (saved in logos and cytoscape folders respectively):

python main.py -g True -n 5

iv) with VDJdb annotation and co-clustering with VDJdb:

python main.py -i data/10X_dexnorm.csv -a True -g True

Results will be recorded in:
data/results.csv (clustering results and grouped accuracy scores over all epitopes) data/results_pr.csv (train and test accuracy scores for the top three epitopes)
logos (Weblogo motifs)
cytoscape/datetime (edge and node lists, see ClusTCR documentation for instructions on how to visualise these in Cytoscape)

To run tcr_compare on your own data, input .csv files should include as a minimum one column entitled 'cdr3.beta' or 'cdr3.alpha'. V, D and J gene usage may be included as columns entitled v.alpha, j.alpha, v.beta, d.beta or j.beta. Clonotype frequency information should be recorded under 'Count', subject information under 'subject:condition' and epitope information under 'Epitope' (all case-sensitive). Any other columns will be passed to the model as-is, and linked to input sequences via unique CDR3 keys post-clustering. 

*** Datasets provided ***

- VDJDb_trimmed_3.csv: A VDJDB reference [6]treated to remove any epitopes represented fewert han three times in the dataset
- 10X_dexnorm.csv: CD8+ multimodal single cell data from 4 healthy donors produced with 10X CellRanger and available at https://support.10xgenomics.com/single-cell-vdj/datasets/3.0.2/vdj_v1_hs_aggregated_donor[X]? where [X]=1,2,3 or 4. Ths data has been processed to improve dextramer UMI signal to noise ratios following the methodology of [7]
- combined_cdr3.csv: Combined inputs from VDJdb, McPas-TCR and GIANA for epitope annotation of orphan TCR datasets [3]

*** Dependencies ***

Software dependencies are listed in requirements.txt
Please note that MUSCLE v3.8.1551 is required for production of WebLogo motifs. MUSCLE should be accessible on the system path, or alternatively a path to the user's MUSCLE executable can be customised using the -mp flag from the main.py python executable.

GLIPH2 executables are provided for .osx and .centos machines in tcr_compare/modules/gliph2/lib/. A GLIPH2 webtool can also be accessed via http://50.255.35.37:8080/


*** Citation ***

For use of this package, please cite the original ClusTCR manuscript [1], and the following:

TCR Compare, co-developed by the Koohy Lab at the University of Oxford and the Basham group at the Rosalind Franklin Institute. Released April 8th 2022 at https://github.com/danobohud/clustcr.

Update 28/03/22: A detailed description of the methodology and the principle findings of this package will be released on Arxiv in due course.

*** References ***

1. Valkiers, S., van Houcke, M., Laukens, K., & Meysman, P. (2021). ClusTCR: a python interface for rapid clustering of large sets of CDR3 sequences with unknown antigen specificity. Bioinformatics, 37(24), 4865–4867. https://doi.org/10.1093/BIOINFORMATICS/BTAB446
2. Huang, H., & Wang, C. (n.d.). Analyzing the Mycobacterium tuberculosis immune response by T-cell receptor clustering with GLIPH2 and genome-wide antigen screening. Nature Biotechnology. https://doi.org/10.1038/s41587-020-0505-4
3. Zhang, H., Zhan, X., & Li, B. (2021). GIANA allows computationally-efficient TCR clustering and multi-disease repertoire classification by isometric transformation. Nature Communications 2021 12:1, 12(1), 1–11. https://doi.org/10.1038/s41467-021-25006-7 
4. Zhang, H., Liu, L., Zhang, J., Chen, J., Ye, J., Shukla, S., Qiao, J., Zhan, X., Chen, H., Wu, C. J., Fu, Y. X., & Li, B. (2020a). Investigation of Antigen-Specific T-Cell Receptor Clusters in Human Cancers. Clinical Cancer Research, 26(6), 1359–1371. https://doi.org/10.1158/1078-0432.CCR-19-3249
5. Mayer-Blackwell, K., Schattgen, S., Cohen-Lavi, L., Crawford, J. C., Souquette, A., Gaevert, J. A., Hertz, T., Thomas, P. G., Bradley, P. G., & Fiore-Gartland, A. (2021). Tcr meta-clonotypes for biomarker discovery with tcrdist3 enabled identification of public, hla-restricted clusters of sars-cov-2 tcrs. ELife, 10. https://doi.org/10.7554/ELIFE.68605
6. Bagaev, D. v., Vroomans, R. M. A., Samir, J., Stervbo, U., Rius, C., Dolton, G., Greenshields-Watson, A., Attaf, M., Egorov, E. S., Zvyagin, I. v., Babel, N., Cole, D. K., Godkin, A. J., Sewell, A. K., Kesmir, C., Chudakov, D. M., Luciani, F., & Shugay, M. (2020). VDJdb in 2019: database extension, new analysis infrastructure and a T-cell receptor motif compendium. Nucleic Acids Research, 48(D1), D1057–D1062. https://doi.org/10.1093/NAR/GKZ874
7. Zhang, W., Hawkins, P. G., He, J., Gupta, N. T., Liu, J., Choonoo, G., Jeong, S. W., Chen, C. R., Dhanik, A., Dillon, M., Deering, R., Macdonald, L. E., Thurston, G., & Atwal, G. S. (2021). A framework for highly multiplexed dextramer mapping and prediction of T cell receptor sequences to antigen specificity. Sci. Adv, 7, 5835–5849. https://www.science.org  