# DGPathinter
a novel model for identifying driver genes via knowledge-driven matrix factorization with prior knowledge from interactome and pathways
![image](https://github.com/JianingXi/DGPathinter/blob/master/image/splash.jpg)

Developer: Jianing Xi <xjn@mail.ustc.edu.cn> from Health Informatics Lab, School of Information Science and Technology, University of Science and Technology of China

## Instructions to DGPathinter (version 1.0.0)

Requirement
------------------------
* 4GB memory
* MATLAB R2015a or later

Prior knowledge
------------------------
1. Prior knowledge from pathways
The file `./PriorKnowledge/Pathway_map_bipartite.mat` contains the bipartite matrix of the mapping from genes to pathways, and the lists including genes and pathways. The pathway used here is [KEGG](https://academic.oup.com/nar/article/27/1/29/1238108/KEGG-Kyoto-Encyclopedia-of-Genes-and-Genomes) [1], [Reactome](https://link.springer.com/chapter/10.1007/11915034_95) [2] and [BioCarta](http://online.liebertpub.com/doi/pdf/10.1089/152791601750294344) [3], which are downloaded from a [previous study](https://academic.oup.com/bioinformatics/article/32/11/1643/1742725/An-integrative-somatic-mutation-analysis-to) [4].
2. Prior knowledge from interactome
The file `./PriorKnowledge/iRefIndex_adj_matrix.mat` contains the adjanceny matrix of interaction network [iRefIndex 9](http://irefindex.org) [5] and the list of genes in the network.


Input data descriptions
------------------------
We provide the input data files of TCGA somatic mutation data of three cancers from [cBioPortal](http://www.cbioportal.org/data_sets.jsp) [6] respectively.

        =================================================================================================
        | CANCER NAME                  | Abbreviation |  FILE DIRECTORY                                 |
        =================================================================================================
        |Breast invasive carcinoma     |    BRCA      |`./Input_data/BRCA.mat`                          |
        -------------------------------------------------------------------------------------------------
        |Glioblastoma multiforme       |     GBM      |`./Input_data/KIRC.mat`                          |
        -------------------------------------------------------------------------------------------------
        |Thyroid carcinoma             |     THCA     |`./Input_data/THCA.mat`                          |
        -------------------------------------------------------------------------------------------------

run DGPathinter
------------------------
To apply DGPathinter on the example input datasets with the default configurations, please run the script file `./example.m` and the result file will be automatically saved as `.mat` file in directory `./Output_data` when the program is finished.

Output data descriptions
------------------------
The descriptions of output variables of mCGfinder are provided below:

        =================================================================================================
        | VARIABLE NAME        | DESCRIPTION                                                            |
        =================================================================================================
        |GeneSelected          |The top 200 genes selected by DGPathinter as potential driver genes.      |
        -------------------------------------------------------------------------------------------------
        |U_sample_indicator    |The sample indicator matrix indicates the assignment of tumor samples to|
        |                      |each layer of the reconstructed matrix. The entry U_ij = 1 indicate that|
        |                      |the i-th samples are assigned to the j-th layer, and 0 otherwise.       |
        -------------------------------------------------------------------------------------------------
        |V_gene_score          |The gene score matrix for each layers of the reconstructed matrix and   |
        |                      |the investigated genes. Higher scores represent the larger possibility  |
        |                      |of the related genes to be potential driver genes.                      |
        -------------------------------------------------------------------------------------------------

Parameter configurations
------------------------
The configurations of DGPathinter can be set as Struct variables `Gene2Path` and `Net_Conf` in `./example.m`, with their descriptions provided below:

        =================================================================================================
        | PARAMETER NAME       | DESCRIPTION                                                            |
        =================================================================================================
        |Gene2Path.lambda_C    |The tuning parameter of prior knowledge from pathways.                  |
        -------------------------------------------------------------------------------------------------
        |Gene2Path.lambda_V    |The tuning parameter of pathway score vectors.                          |
        -------------------------------------------------------------------------------------------------
        |Net_Conf.lambda_L     |The tuning parameter of prior knowledge from interactome.               |
        -------------------------------------------------------------------------------------------------
        |MC                    |Maximum layer number or maximum rank of the reconstructed matrix.       |
        |                      |The default value is 5.                                                 |
        -------------------------------------------------------------------------------------------------
        |LLP                   |Least proportion of samples assigned to each layers. The default value  |
        |                      |is 15%.                                                                 |
        -------------------------------------------------------------------------------------------------

[1] Ogata H, Goto S, Sato K, et al. KEGG: Kyoto encyclopedia of genes and genomes[J]. Nucleic acids research, 1999, 27(1): 29-34.

[2] Schmidt E, Birney E, Croft D, et al. Reactome¨Ca knowledgebase of biological pathways[C]//OTM Confederated International Conferences" On the Move to Meaningful Internet Systems". Springer, Berlin, Heidelberg, 2006: 710-719.

[3] Nishimura D. BioCarta[J]. Biotech Software & Internet Report: The Computer Software Journal for Scient, 2001, 2(3): 117-120.

[4] Park S, Kim S J, Yu D, et al. An integrative somatic mutation analysis to identify pathways linked with survival outcomes across 19 cancer types[J]. Bioinformatics, 2015, 32(11): 1643-1651.

[5] Razick S, Magklaras G, Donaldson I M. iRefIndex: a consolidated protein interaction database with provenance[J]. BMC bioinformatics, 2008, 9(1): 405.

[6] Gao J, Aksoy B A, Dogrusoz U, et al. Integrative analysis of complex cancer genomics and clinical profiles using the cBioPortal[J]. Science signaling, 2013, 6(269): pl1.