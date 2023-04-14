This folder contains all data and code to perform the polyG analysis. 

Genotyping of polyGs was performed in two batches, and the raw data from these batches can be found in the "data" directory.
The code to generate all results can be found in the "src" directory and is written in R version 4.1.2

These R libraries are needed for the analysis and plotting: 
"BiocManager"
"tidyverse"
"data.table"
"ggtree"
"tidytree"
"ape"
"phangorn"
"Rphylip"
"ggbeeswarm"
"rstatix"
"ggpubr"

First, the R scripts need to be run in this order: "cor_batch1.R", and "cor_batch2.R", "t790m_tree.R", "cor_interpatient.R".

Then the jupyter notebook "PolyG Results.ipynb" can be run to generate all plots and results.
The jupyter notebook has been exported as an HTML file to view the results with accompanying code without having to rerun the analysis. 