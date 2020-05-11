# ElenMatchR v1.0.9003 comparative genomics tool for Eggerthella lenta and Coriobacteriia.

## Background
Welcome to ElenMatchR, a tool for matching phenotypes to genes and SNPs within the Coriobacteriia, and especially Eggerthella lenta. Two modes are available: genes and kmers. Genes describes clusters of orthologous genes across genomes and the BLASTP clustering thresholds can be modified as desired. An earlier version of this BLASTP method was used to identify the enzymes responsible for the metabolism of digoxin (Koppel et al., eLife 2018) and pinoresinol (Bess et al., Nat Microbiol 2019). The kmer-mode is an improvement to the method that is capable of both detecting the presence and absence of genes, and also nucleotide variation impacting phenotypes. Although compute time and interpretation is more involved, this method captures subtle variation including the SNP responsible for determining the ability of E. lenta to metabolize dopamine (Rekdal et al., Science 2019). While a web-accessable version is currently running, users are encouraged to download the ElenMatchR source and run locally on their own machine. Genomes and annotations can be extracted as below for direct access to the underlying information for the tool.

***

## Usage

### Options

The first important consideration is selecting the phenotype to study. Published datasets are included for demonstration purposes and are available under the **Preloaded phenotype to analyze:** drop down menu:

![inpheno](https://github.com/jbisanz/ElenMatchR/raw/master/images/inpheno.png)

To analyze a new phenotype, select **Download xlsx template**, modify the file as needed and reupload using the **Browse** button above.

![meta](https://github.com/jbisanz/ElenMatchR/raw/master/images/metadata.png)

Next, analysis parameters can be adjusted in the **Parameters** pane:

![mode](https://github.com/jbisanz/ElenMatchR/raw/master/images/mode.png)

A critical choice is whether the user wants to study the presence/absence of genes or kmers. Genes are easily interpreted, but kmers can also detect the assosciation of SNPs with a particular phenotype. *Note: Kmers will require considerably more computational power and may take longer to compute.*

![inpram](https://github.com/jbisanz/ElenMatchR/raw/master/images/in_pram.png)

Next the user can modify parameters for the random forest increasing the number of replications for the analysis and the number of trees built per analysis. More replications will stabilize the results but take longer to process. Default settings should provide robust and stable estimates of variable importance. **In Genes mode only** the user can also pick the clustering parameters used to define orthologs (the cluster of proteins across organisms taken to represent the same gene). When dealing with within species-analysis, a more stringent threshold is proposed, while across species a more lenient  threshold (i.e. 30% identity, 50% coverage) is suggested.

Finally, when all parameters have been adjusted, analysis is initiated via the **Execute Analysis** button and the status of the analysis is displayed in the bottom right hand corner of the screen.

![exec](https://github.com/jbisanz/ElenMatchR/raw/master/images/exec.png)

***

### Outputs
All outputs are available via tabs at the top of the screen. These tabs will be populated with after the analysis is complete and will be dynamically regenerated in response to changes within the tab.

#### Phenotypes
![phenos](https://github.com/jbisanz/ElenMatchR/raw/master/images/phenos.png)
This tab provides a table containing a copy of the input phenotypes. The table can be downloaded in csv or excel format to keep as a record.

#### Classifier Results
![confmat](https://github.com/jbisanz/ElenMatchR/raw/master/images/confmat.png)
This table contains the confusion matrix for each replication of the model building. In the provided example, of the 13 strains which did not display metabolism, 11 were correctly assigned to the no metabolism phenotype based on the generated model resulting in a classification error of 15.4% (1-(11/13)). These values are provided for reference; however, a poorly performing classifier does not necessarily indicate a failure to uncover features with predictive accuracy.

#### Importance Scatter Plot
![impscat](https://github.com/jbisanz/ElenMatchR/raw/master/images/impscat.png)
This plot ranks all features (genes/kmers) and shows their relative importance (mean Decrease GINI). The error bars represent the standard deviation of the estimate determined from running multiple replications of the classifier. In this example we see a set of 7 features which clearly have elevated importance compared to the next set of best predictors. The plot can be downloaded in .pdf format using the **Download Importance Scatter Plot** button and the scores themselves can be downloaded using **Download Importance Scores for all features** button.


#### Heatmap
![heatmap](https://github.com/jbisanz/ElenMatchR/raw/master/images/heatmap.png)
This tab may contain one of the most crucial visual tools for interpreting the results. In this example we see multiple genes which are uniquely present in only strains capable of metabolism. The default is to plot the 30 most important features (which can be adjusted using the **Number of features on plot** slider). Features on the Y-axis are ranked from highest to lowest score (top to bottom). A .pdf is available using the **Download Heatmap** button.


#### Manhattan Plot
![manhat](https://github.com/jbisanz/ElenMatchR/raw/master/images/manhat.png)
To understand the genomic context of hits and how they may be co-localized, a Manhattan plot is also generated. The default reference genome is the contiguous assembly of Eggerthella lenta DSM2243; however, it can be modified using the **Reference Genome** drop down box. The number of features plotted can also be modified using the **Number of features on plot** slider. A .pdf can be downloaded using the **Download Plot** button. In this case, we can see that the most predictive features are co-localized around ~2.9Mb in the genome.

#### Phylogenetic Tree
![phylotree](https://github.com/jbisanz/ElenMatchR/raw/master/images/phylotree.png)
To understand potential correlations between the phenotype and phylogeny, a tree is provided. Two trees are available from two tools: PhyloPhlAn (Segata et al 2013) and Roary (Page et al. 2015). The Roary tree is based on SNPs in Eggerthella-specific SNPs and as such will only plot Eggerthella genus members. The PhyloPhlAn tree is based on alignment of ~400 conserved marker genes conserved across all bacteria. The tree can also be plotted as linear or circular tree (using the **Orientation** drop down) and as a cladogram (using the **Branch Lengths** drop down). Finally, strains without phenotypes can be maintained in the tree using the **Prune Tree?** drop down.
The figure is then available for download as a .pdf (**Download Tree Plot**) or as a Newick format tree (**Download Newick Format Tree**).

#### Annotation Table
![annotations](https://github.com/jbisanz/ElenMatchR/raw/master/images/annotations.png)
Finally, information about every gene feature (Ortholog or SNP) is available in a table. The sequences of every gene can be downloaded using the **Download results in fasta (nucleotide) format** button. Columns in the table are as follows:
* Rank: The rank based on importance score.
* FeatureID: A unique identifier for the feature. OG_ denotes running in gene mode while ID60_COV80_ refers to the clustering thresholds for the run and 2997 is a the unique number of this cluster of genes.
* GenomeID: A unique identifier for the genome in which the feature is observed in. These match to the genome files which are available for extraction (see below).
* Phenotype: The user-assigned phenotype for the strain.
* Gene: The unique gene ID within the strain which is a member of the feature. If running in Kmer mode, not every Kmer will occur within a gene and this will be NA.
* Contig: The unique ID for the contig within the strain.
* start: the base position at which the gene/kmer starts.
* end: the base position at which the gene/kmer ends.
* strand: the orientation of the gene/kmer.
* attributes: information about annotations and other predictions for the feature.

***

## Installation
### Bypassing Installation
ElenMatchR is currently available via a [shinyapps instance](https://jbisanz.shinyapps.io/elenmatchr/) which does not require installation and runs on an external server. Certain jobs, including Kmer analysis, may exceed available memory and will require a local installation as below.
### Local Installation
ElenMatchR and dependencies can be installed via the devtools package as below:
```
library(devtools) # may need to install if not already present
install_github("jbisanz/ElenMatchR")
```
The application can then be initialized using the command:
```
ElenMatchR::run_ElenMatchR()
```

***
## Accessing Genomes
ElenMatchR contains 95 genomes and their accompanying annotations which are stored internally as RDS files (compressed serialized r objects). These can be exported to fasta and gff files using the an included utility script as below in a bash terminal:
```
cd LocationOfYourInstall
Rscript UnpackGenomes.R
```
The resulting files will be available in a new directory: `$PWD/genomes`. Genomes are also available from the NCBI under the following accessions: PRJNA412637, GCA_000403355.2, GCA_000478885.1, GCA_003725335.1, GCA_000422625.1, GCA_900169505.1, GCA_003726015.1, GCA_000169035.1, GCA_000195315.1, GCA_000023845.1, GCA_900110565.1, GCA_000763035.1, GCA_000185625.1, GCA_002148255.1, GCA_000191845.2, GCA_003438525.1, GCA_000270285.1, GCA_900184265.1, GCA_000311845.1, GCA_002899715.1, GCA_003788985.1, GCA_900170005.1, GCA_000210055.1, GCA_003788975.1 , GCA_000174015.1, GCA_900143685.1, GCA_900078545.1, GCA_000143845.1, GCA_003726035.1, GCA_003966955.1, GCA_900169485.1, GCA_900199545.1, GCA_900240215.1, GCA_002899695.1, GCA_000236865.1, GCA_003725995.1, GCA_000162875.1, GCA_003725295.1, GCA_000023885.1, GCA_003725955.1, GCA_000296445.1.

***

## Dependencies
* R >3.5.0
* tidyverse
* shiny
* randomForest
* ggtree
* ggdendro
* DT
* readxl
* Biostrings
* Matrix
* ape
* pryr
