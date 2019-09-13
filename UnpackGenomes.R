# UnpackGenomes.R
# script to get local assemblies and annotations for Coriobacteriia collection
# requires tidyverse and Biostrings

library(tidyverse)
library(Biostrings)

dir.create("genomes")

for(f in list.files("resources/fna/", pattern="RDS$", full.names=TRUE)){
  message(f)
  tmp<-readRDS(f)
  writeXStringSet(tmp, paste0("genomes/", gsub("RDS$","fna",basename(f))))
}

gff<-readRDS("resources/merged.gff")

gff<-split(gff, gff$GenomeID)

lapply(names(gff), function(x) gff[[x]] %>% filter(GenomeID==x) %>% select(-GenomeID) %>% write_tsv(., paste0("genomes/", x, ".gff3")) )

message(date(), " Complete, see $PWD/genomes for fasta files and tabular gff files with annotations.")