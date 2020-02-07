# ElenMatchR v1.0dev: comparative genomics tool for Eggerthella lenta and Coriobacteriia.

Under Construction!!!!

## Background

## Usage



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

### Accessing Genomes
ElenMatchR contains 95 genomes and their accompanying annotations which are stored internally as RDS files (compressed serialized r objects). These can be exported to fasta and gff files using the an included utility script as below in a bash terminal:
```
cd LocationOfYourInstall
Rscript UnpackGenomes.R
```
The resulting files will be available in a new directory: `$PWD/genomes`.

### Dependencies
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
