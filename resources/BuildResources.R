# BuildResources.R
# Generate the necessary databases for ElenMatchR
stop()

setwd("/Volumes/turnbaughlab/qb3share/jbisanz/ElGenomes2019/ElenMatchR_v1.0/ElenMatchR/")
library(tidyverse)
library(Biostrings)

dir.create("resources")
dir.create("resources/fna")

# Get tree
library(ape)
#tree<-read.tree("/Volumes/turnbaughlab/qb3share/jbisanz/ElGenomes2019/phylophlan/output/AllStrains/AllStrains.tree.nwk")
tree<-read.tree("/Volumes/turnbaughlab/qb3share/jbisanz/ElGenomes2019/phylophlan/output/AllStrains_nomags/AllStrains_nomags.tree.nwk")
tree<-root(tree, outgroup="GCA_000022965.1_ASM2296v1_protein", resolve.root=TRUE)
is.rooted(tree)
plot(tree)
saveRDS(tree,"resources/tree.RDS")

# Gene calls
gffs<-
list.files("/Volumes/turnbaughlab/qb3share/jbisanz/ElGenomes2019/annotation/gff/", pattern="gff", full.names = TRUE) %>%
  lapply(., function(x) {
    read_tsv(x, col_names=c("seqid","source","type","start","end","score","strand","phase","attributes")) %>%
      mutate(GenomeID=basename(x) %>% gsub("\\.gff","",.)) 
    }) %>%
  do.call(bind_rows, .) %>%
  select(GenomeID, everything())

saveRDS(gffs, "resources/merged.gff")

# contigs
for(f in list.files("/Volumes/turnbaughlab/qb3share/jbisanz/ElGenomes2019/annotation/scaffolds/", pattern="fna", full.names = TRUE)){
  message(f)
  readDNAStringSet(f) %>%
    saveRDS(., paste0("resources/fna/", gsub("\\.fna",".RDS", basename(f))))
}

# Gene orthologs
dir.create("resources/matrices")
for(f in list.files("/Volumes/turnbaughlab/qb3share/jbisanz/ElGenomes2019/proteinortho/outs/", pattern="\\.tsv$", full.names = TRUE)){
  message(f)
  id=basename(f) %>% str_split(., "_") %>% .[[1]] %>% .[1]
  cov=basename(f) %>% str_split(., "_") %>% .[[1]] %>% .[2]
  
  f<-
  read_tsv(f) %>%
    arrange(desc(`# Species`)) %>%
    select(-1,-2,-3) %>%
    mutate(OG=paste0("OG_",id,"_",cov,"_", 1:nrow(.))) %>%
    select(OG, everything())

  colnames(f)<-gsub("\\.faa","", colnames(f))
  
  f[f=="*"]<-NA
  
  saveRDS(f, paste0("resources/matrices/", id,cov,".RDS"))
}
  
#kmers

## build list of unique kmers across Eggerthella lenta, the kmer analysis will be limited to this species due to computational requirements, run on biggut
dir.create("resources/kmers")
library(tidyverse)
klist<-vector()
for (f in list.files("/labmainshare/qb3share/jbisanz/ElGenomes2019/kmers/", pattern = "31mers", full.names=TRUE) %>% grep("Eggerthella_lenta", ., value=T)){
  fk<-readLines(f)
  fk<-gsub(" [0-9]+$","", fk)
  fk<-fk[!fk %in% klist]
  klist<-c(klist, fk)
  message(f, ": ", length(fk),"kmers, total kmers: ", length(klist )," object size:", format(object.size(klist), "Gb"))
}
klist<-sort(klist)
saveRDS(klist, "resources/kmers/kmerlist.RDS")

## Now Generate Kmer Indicies on a per genome basis
for (f in list.files("/labmainshare/qb3share/jbisanz/ElGenomes2019/kmers/", pattern = "31mers", full.names=TRUE) %>% grep("Eggerthella_lenta", ., value=T)){
  message(f)
  fk<-readLines(f)
  fk<-gsub(" [0-9]+$","", fk)
  idx<-which(klist %in% fk)
  saveRDS(idx, paste0("resources/kmers/",gsub("\\.31mers\\.gz", ".idx", basename(f))))
}

# Now Generate Kmer sparse matrix
library(Matrix)
library(tidyverse)
klist<-readRDS("resources/kmers/kmerlist.RDS")
  idxs<-list.files("resources/kmers/", pattern="idx$", full.names = TRUE) # remove me
  ks<-lapply(idxs, function(x) readRDS(x))
  
  j=as.integer()
  for(i in 1:length(ks)){
    j<-c(j, rep(i, length(ks[[i]])))
  }
  
  kmat<-sparseMatrix(i = unlist(ks), j=j ,x=1)
  rm(j)

  colnames(kmat)<-sapply(idxs, function(x) basename(x) %>% gsub("\\.idx","", .))
  rownames(kmat)<-klist
  
  saveRDS(kmat, "../../kmers/kmat_fromelen_prefilt.RDS")
  
  #sum(rowSums(kmat)==48)
  #542496
  kmat<-kmat[rowSums(kmat)<48,]
  #remove kmers which are conserved across the species
  
  saveRDS(kmat, "resources/kmers/KmerMatrix.RDS")

  #Now need to collapse the co-occuring kmers which should hopefully help save considerably on memory, then store the dereplicated matrix and a kmer to id lookup table
  kmat<-readRDS("resources/kmers/KmerMatrix.RDS")
  
  
  
  kclusts<-tibble(kmer=rownames(kmat), cluster=apply(kmat, 1,  function(x) paste(which(x==1), collapse=","))) #store comma separated list of positive locations as this will prevent the storing of excessive 0s
  klookup<-kclusts %>% select(cluster) %>% distinct() %>% mutate(KClusterID=paste0("kc_", 1:nrow(.))) # actually only 96,825 unique combinations
  klookup<-kclusts %>% left_join(klookup) %>% select(KClusterID, kmer)
  saveRDS(klookup,"resources/kmers/kmer_lookup.RDS")
  
  keeps<-klookup %>% filter(!duplicated(KClusterID))
  
  fmat<-kmat[keeps$kmer,]
  rownames(fmat)<-keeps[match(keeps$kmer, rownames(fmat)),]$KClusterID
  
  saveRDS(fmat, "resources/kmers/kmer_derep.RDS")
#removed all the temporary .idx files and the kmer list as it is redundant with the table

#split kmers into two files to upload to biggut
kmers<-readRDS("/Volumes/turnbaughlab/qb3share/jbisanz/ElGenomes2019/ElenMatchR_v1.0/ElenMatchR/resources/kmers/kmer_lookup.RDS")
k1<-kmers[1:5319404,]
k2<-kmers[5319405:10638807,]

saveRDS(k1,"resources/kmers/kmer_lookup1.RDS")
saveRDS(k2,"resources/kmers/kmer_lookup2.RDS")


