library(shiny)
library(randomForest)
library(ggtree)
library(ggdendro)
library(DT)
library(tidyverse)
library(readxl)
library(Biostrings)
library(Matrix)
library(ape)
library(pryr)

#options(repos=BiocManager::repositories())
#library(rsconnect)
#deployApp()

#Static Vars on load
phenos<-read_excel("resources/ElenMatchR_phenotypes.xlsx", skip=1)
tree<-readRDS("resources/tree.RDS")
PWIDTH=5 #plot width
PHEIGHT=4 #plot height

# User interface ########## ########## ########## ########## ########## ########## ##########
ui <- navbarPage(HTML("ElenMatchR: Comparative Genomics Tool v1.0.9000"),
                 
      tabPanel("Input",
        h3("Instructions"),
          p("Welcome to ElenMatchR, a tool for matching phenotypes to genes and SNPs within the Coriobacteriia, and especially Eggerthella lenta. To conduct the primary analysis, select the appropriate parameters below and click execute analysis to begin. You will be notified of the progress, and upon completion, the results tabs will be populated with interactive figures."),
          p("Two modes are available: genes and kmers. Genes describes clusters of orthologous genes across genomes and the BLASTP clustering thresholds can be modified below. An earlier version of this BLASTP method was used to identify the enzymes responsible for the metabolism of digoxin (Koppel et al., eLife 2018) and pinoresinol (Bess et al., Nat Microbiol 2019). The kmer-mode is an improvement to the method that is capable of both detecting the presence and absence of genes, and also nucleotide variation impacting phenotypes. Although computing time and interpretation is more involved, this method captures subtle variation including the SNP responsible for determining the ability of E. lenta to metabolize dopamine (Rekdal et al., Science 2019)."),
          p("Note: Shinyapps.io sessions are limited to 1Gb of RAM and may fail due to a lack of memory for certain Kmer-based analyses. Users are encouraged to download the ElenMatchR source and run locally on their own machine (github.com/jbisanz/ElenMatchR). Genomes and annotations are also available for download within the source files."),
        hr(),
          actionButton(inputId="StartAnalysis", label="Execute Analysis"),
        hr(),
        h3("Phenotypes"),
          selectInput("UserPhenotype", "Preloaded phenotype to analyze:", selected="Tetracycline_Resistance", choices=c("Tetracycline_Resistance","Digoxin_Metabolism","Pinoresinol_Metabolism","Dopamine_Metabolism")),
          fileInput("phenofile","Or upload your own using the provided template:", multiple = FALSE),
          downloadButton("DLtemplate", label="Download xlsx template"),
        hr(),
        h3("Parameters"),
          selectInput("Mode", "Analysis mode:", selected="Genes", choices=c("Genes","Kmers")),
        h4("Random forest:"),
          sliderInput("Nreps", "Number of replications:", min=3, max=10, value=3),
          sliderInput("Ntrees", "Number of trees generated in model:", min=500, max=5000, value=1000, step=500),
        h4("BLAST (Genes mode only):"),
          selectInput("PID", "Minimum percent BLASTP identity:", selected=60, choices=c(seq(30,90,10), 95, 99)),
          selectInput("COV", "Minimum percent BLASTP coverage:", selected=80, choices=c(seq(50,90,10), 95, 99)),
        hr(),
        p("Citation: Bisanz et al., 2019 DOI:10.1101/304840")
      ),
      tabPanel("Phenotypes",
        dataTableOutput("PhenoTable"),
        p("Description: List of user-provided phenotypes. Please check for accuracy.")
      ),
      tabPanel("Classifier Results",
        dataTableOutput("ConfusionMatrices"),
        p("Description: This table represents the confusion matrix for most accurate classifier built during feature selection. Given that we are primarily interested in feature selection, these results are provided for informational purposes.")            
      ),
      tabPanel("Importance Scatter Plot",
        sliderInput("Scatter_Nfeats", "Number of features on plot:",  min = 1, max = 1000, value = 100),
        downloadButton("DLImportancePlot", label="Download Importance Scatter Plot"),
        downloadButton("DLImportanceRaw", label="Download Importance Scores for all features"),
        p("Description: This plot has ranked all features (gene or kmer clusters) by their estimated importance (Mean Decrease GINI)."),
        plotOutput("ImportancePlot")
      ),
      tabPanel("Heatmap",
        sliderInput("Heatmap_Nfeats", "Number of features on plot:",  min = 1, max = 100, value = 30),
        downloadButton("DLHeatmap", label="Download Heatmap"),
        p("Description: This plot has ranked all features (gene or kmer clusters) by their estimated importance (Mean Decrease GINI), and plotted their presense as a function of phenotype."),
        plotOutput("Heatmap")
      ),
      tabPanel("Manhattan Plot",
        selectInput("RefGenome", "Reference Genome", selected="Eggerthella lenta DSM 2243", choices=phenos$Strain_Name),
        sliderInput("MA_Nfeats", "Number of features on plot:",  min = 10, max = 100, value = 10),
        downloadButton("DLMAplot", label="Download Plot"),
        p("This plot shows the genomic location of important features within a selected reference genome. Where multiple scaffolds are present, each has been plotted individually. scaffolds without a hit are not plotted."),
        plotOutput("ManhattanPlot")
      ),
      tabPanel("Phylogenetic Tree",
        selectInput("TreeSource", "Tree Source", selected="Phylophlan", choices=c("Phylophlan","Eggerthella core genes")),
        selectInput("BranchLengths", "Branch Lengths", selected="Phylogenetic Tree", choices=c("Phylogenetic Tree","Cladogram")),
        selectInput("Orientation", "Orientation", selected="Linear", choices=c("Linear","Circular")),
        selectInput("Prune", "Prune Tree?", selected="Remove genomes without phenotype", choices=c("Plot all genomes in tree","Remove genomes without phenotype")),
        downloadButton("DLTreePlot", label="Download Tree Plot"),
        downloadButton("DLTree", label="Download Newick Format Tree"),
        p("Note: phylogenetic trees created using PhyloPhlAn or from alignment of core genes using Roary and FastTree. Citations: Segata N, et al. PhyloPhlAn is a new method for improved phylogenetic and taxonomic placement of microbes. Nat Commun. 2013;4:2304. doi:10.1038/ncomms3304; Page AJ et al. Roary: rapid large-scale prokaryote pan genome analysis. Bioinformatics 2-15. 31(22):3691-3693. DOI:10.1093/bioinformatics/btv421."),
        plotOutput("TreePlot")
      ),
      tabPanel("Annotation Table",
        sliderInput("AT_Nfeats", "Number of features recalled for table:",  min = 2, max = 100, value = 10),
        p("The table below contains the top ranked features, and where they occur in each reference genome, and available annotations at these loci. You can also download the sequences in fasta format."),
        downloadButton("DLfasta", label="Download results in fasta (nucleotide) format."),
        hr(),
        dataTableOutput("AnnotationTable")
      )
)




#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

server <- function(input, output, session) {
  options(shiny.maxRequestSize = 25*1024^2)
  tmps<<-list()
  
  #Functions ########## ########## ########## ########## ########## ########## ##########
  
  theme_ElenMatchR<- function () { 
    theme_classic(base_size=8, base_family="Helvetica") +
      theme(panel.border = element_rect(color="black", size=1, fill=NA)) +
      theme(axis.line = element_blank(), strip.background = element_blank())
  }
  
  
  DoRF<-function(phenofactor, pmatrix, nrep, ntrees){
    fits<-lapply(1:nrep, function(x){randomForest(x=t(as.matrix(pmatrix)), y=phenofactor, ntree = ntrees, importance=TRUE)})
    imps<-lapply(1:nrep, function(x){fits[[x]]$importance %>% as.data.frame() %>% rownames_to_column("FeatureID") %>% mutate(Iteration=x)})  %>% 
      do.call(bind_rows, .) %>%
      group_by(FeatureID) %>% 
      summarize(Importance=mean(MeanDecreaseGini), SD=sd((MeanDecreaseGini))) %>%
      arrange(desc(Importance)) %>%
      mutate(Rank=1:nrow(.))
    confs<-lapply(1:nrep, function(x){fits[[x]]$confusion %>% as.data.frame() %>% rownames_to_column("Contrast") %>% mutate(Iteration=x)}) %>%
      do.call(bind_rows, .)
    return(list(ImportanceTable=imps, ConfusionMatrices=confs))
  }
  
  DrawScatter<-function(importance, ntoplot, phenoname, mode){
    iplot<-
      importance %>%
      filter(Rank<=ntoplot) %>%
      ggplot(aes(x=Rank, y=Importance, ymin=Importance-SD, ymax=Importance+SD, label=FeatureID)) +
      geom_line() +
      geom_errorbar(width=0) +
      geom_point() +
      ylab("Importance (Mean Decrease GINI ± SD)") +
      theme_ElenMatchR() +
      ggtitle(paste("Importance scatter plot:", phenoname))
    
    if(ntoplot<40){iplot<-iplot+geom_text(angle=45, hjust=0)}
    
    return(list(Plot=iplot, Table=importance %>% select(Rank, everything())))
  }
  
  DrawHeatMap<-function(importance, phenotable, ntoplot, matrix, phenoname, mode){
    hm<-
      importance %>%
      filter(Rank<=ntoplot) %>%
      select(FeatureID) %>%
      left_join(
        matrix %>%
          as.matrix() %>%
          as.data.frame() %>%
          rownames_to_column("FeatureID") %>%
          gather(-FeatureID, key="GenomeID", value="Present")
      ) %>%
      mutate(FeatureID=factor(FeatureID, levels=rev(importance$FeatureID))) %>%
      left_join(phenotable) %>%
      mutate(Present=if_else(Present==1, "Present","Absent")) %>%
      ggplot(aes(x=Strain_Name, y=FeatureID, fill=Present)) +
      geom_tile() +
      facet_grid(~Phenotype, scales="free", space="free") +
      theme_ElenMatchR() +
      coord_cartesian(expand=F) +
      theme(axis.text.x=element_text(angle=45, hjust=1)) +
      ggtitle(paste("Feature heatmap:", phenoname)) +
      scale_fill_manual(values=c("grey90","cornflowerblue"), name=paste(mode, "presence")) +
      xlab("Strain ID")
    return(hm)
  }
  
  
  DrawMA<-function(refgenome, importance, ntoplot, mode, genelist, phenoname){
    genomeid<-read_excel("resources/ElenMatchR_phenotypes.xlsx", skip=1) %>% filter(Strain_Name==refgenome) %>% pull(GenomeID)
    fasta<-readRDS(paste0("resources/fna/",genomeid ,".RDS"))
    gff<-readRDS("resources/merged.gff") %>% filter(GenomeID==genomeid) %>% filter(type=="CDS") %>% mutate(Genes=gsub(";..+","", attributes) %>% gsub("ID=","", .))
    if(mode=="Genes"){
      maplot<-
        importance %>%
        filter(Rank<=ntoplot) %>%
        left_join(
          genelist %>% dplyr::rename(FeatureID=OG) %>% select(FeatureID, genomeid)
        ) %>%
        dplyr::rename(Genes=genomeid) %>%
        mutate(Genes=gsub(",..+","", Genes)) %>% #masking multiple copies
        left_join(gff) %>%
        filter(!is.na(start)) %>% #remove features not in the reference
        ggplot(aes(x=start, y=Importance, ymin=0, ymax=Importance, color=Importance)) +
        geom_errorbar(width=0) +
        geom_point() +
        theme_ElenMatchR() +
        geom_hline(yintercept = 0) +
        facet_grid(~seqid, space="free") +
        xlab("Base position") +
        ylab("Mean Decrease GINI") +
        scale_color_viridis_c() +
        theme(legend.position="none") +
        ggtitle(paste("Manhattan Plot:", phenoname))
    } else{
      #kmerlookup<-readRDS("resources/kmers/kmer_lookup.RDS")
      kmerlookup<-rbind(readRDS("resources/kmers/kmer_lookup1.RDS"),readRDS("resources/kmers/kmer_lookup2.RDS"))
      timp<-
        importance %>%
        filter(Rank<=ntoplot) %>%
        left_join(kmerlookup %>% dplyr::rename(FeatureID=KClusterID)) %>%
        select(FeatureID, Importance, kmer)
      
      hits<-sapply(timp$kmer, function(x) vmatchPattern(DNAString(x), fasta))
      locations<-
        lapply(names(hits), function(x){
          z=hits[[x]] %>% unlist()
          tibble(kmer=x, Contig=names(z), position=z@start)
        }) %>%
        do.call(bind_rows, .) %>%
        left_join(timp)
      maplot<-
        locations %>%
        ggplot(aes(x=position, y=Importance, ymin=0, ymax=Importance, color=Importance)) +
        geom_errorbar() +
        geom_point() +
        theme_ElenMatchR() +
        geom_hline(yintercept = 0) +
        facet_grid(~Contig, space="free") +
        xlab("Base position") +
        ylab("Mean Decrease GINI") +
        scale_color_viridis_c() +
        theme(legend.position="none") +
        ggtitle(paste("Manhattan Plot:", phenoname))
    }
    
    return(maplot)
  }
  
  DrawTree<-function(branches, ori, prune, phenotable, phenoname, tsource){
    phenotable<-phenotable %>% select(GenomeID, everything()) %>% right_join(phenos[,c("GenomeID","Strain_Name")])
    mtree<-tree
    if(tsource!="Phylophlan"){mtree<-readRDS("resources/eggtree.RDS")}
   if(prune=="Remove genomes without phenotype"){ttree<-drop.tip(mtree, mtree$tip.label[!mtree$tip.label %in% subset(phenotable, !is.na(Phenotype))$GenomeID])} else {ttree<-mtree}
  if(branches=="Phylogenetic Tree") {
    if(ori=="Linear"){
      tplot<-
      ggtree(ttree) %<+% phenotable +
       geom_tippoint(aes(fill=Phenotype), shape=21, size=3) +
       geom_tiplab(aes(label=Strain_Name), size=4) +
        scale_fill_discrete(name=phenoname) +
        theme(legend.position="right")
      tplot<-tplot + scale_x_continuous(limits=c(0, 1.5*max(tplot$data$x)))
     } else {
       tplot<-
       ggtree(ttree, layout="circular") %<+% phenotable +
         geom_tippoint(aes(fill=Phenotype), shape=21, size=3) +
         geom_tiplab2(aes(label=Strain_Name)) +
         scale_fill_discrete(name=phenoname) +
         theme(legend.position="right")
       tplot<-tplot + scale_x_continuous(limits=c(0, 1.5*max(tplot$data$x)))
     }
  } else {
    if(ori=="Linear"){
      tplot<-
      ggtree(ttree, branch.length = "none") %<+% phenotable +
        geom_tippoint(aes(fill=Phenotype), shape=21, size=3) +
        geom_tiplab(aes(label=Strain_Name), size=4) +
        scale_fill_discrete(name=phenoname) +
        theme(legend.position="right")
    } else {
      tplot<-
      ggtree(ttree, layout="circular", branch.length = "none") %<+% phenotable +
        geom_tippoint(aes(fill=Phenotype), shape=21, size=3) +
        geom_tiplab2(aes(label=Strain_Name)) +
        scale_fill_discrete(name=phenoname) +
        theme(legend.position="right")
    }
  }
   
  return(tplot)
}

MakeAnnoTable<-function(mode, importance, genelist, nfeats, phenotable){
  if(mode=="Genes"){
    gff<-readRDS("resources/merged.gff") %>% filter(type=="CDS") %>% mutate(Gene=gsub("ID=","", attributes) %>% gsub(";..+","", .))
    timp<-
    importance %>%
      filter(Rank<=nfeats) %>%
      left_join(genelist %>% dplyr::rename(FeatureID=OG)) %>%
      gather(-FeatureID, -Importance, -SD, -Rank, key="Genome", value="Gene") %>%
      filter(!is.na(Gene)) %>%
      left_join(gff) %>%
      left_join(phenotable[,c("GenomeID","Phenotype")]) %>%
      filter(!is.na(Phenotype)) %>%
      select(Rank, FeatureID, GenomeID, Phenotype, Gene, Contig=seqid, start, end, strand, attributes) %>%
      arrange(Rank)
    
    seqs<-
    timp %>%
      select(FeatureID, GenomeID, Gene, Contig, start, end, strand) %>%
      arrange(GenomeID)
    
    splitseqs<-split(seqs, seqs$GenomeID)
    splitseqs<-lapply(names(splitseqs), function(x){
      fa<-readRDS(paste0("resources/fna/", x, ".RDS"))
      tmp<-splitseqs[[x]]
      tmp$Seq<-apply(tmp, 1, function(y) subseq(fa[[as.character(y['Contig'])]], start=as.numeric(y['start']), end=as.numeric(y['end'])) %>% as.character())
      return(tmp)
    })
    seqs<-do.call(bind_rows, splitseqs)
    seqs$RC<-sapply(seqs$Seq, function(x) as.character(reverseComplement(DNAString(x))))
    seqs<-seqs %>% mutate(fasta=if_else(strand=="+", Seq, RC)) %>% select(-Seq, -RC)
    fa<-DNAStringSet(seqs$fasta)
    names(fa)<-paste0(seqs$FeatureID,"|", seqs$Gene)
    return(list(fasta=fa, annotable=timp))
  } else {
    #klookup<-readRDS("resources/kmers/kmer_lookup.RDS")
    klookup<-rbind(readRDS("resources/kmers/kmer_lookup1.RDS"),readRDS("resources/kmers/kmer_lookup2.RDS"))
    timp<-
      importance %>%
      filter(Rank<=nfeats) %>%
        left_join(
          klookup %>% dplyr::rename(FeatureID=KClusterID)
        ) %>%
      mutate(subid=1:nrow(.))
    fa<-DNAStringSet(timp$kmer)
    names(fa)<-paste0(timp$FeatureID,"|",timp$subid)
    rm(klookup)
    
    
    fasta<-sapply(phenotable$GenomeID, function(x) readRDS(paste0("resources/fna/",x,".RDS")) ) #make fasta, which is all the ref genomes of interest
    mfasta<-DNAStringSet()
    for(f in names(fasta)){mfasta<-c(mfasta, fasta[[f]])}
    rm(fasta)
    
    hits<-sapply(timp$kmer, function(x) vmatchPattern(DNAString(x), mfasta)) # get hits
    rm(mfasta)
    locations<-
      lapply(names(hits), function(x){
        z=hits[[x]] %>% unlist()
        tibble(kmer=x, Contig=names(z), position=z@start)
      }) %>%
      do.call(bind_rows, .)
    rm(hits)
    # find if start position is within a coding sequence
    gff<-readRDS("resources/merged.gff") %>% filter(GenomeID %in% phenotable$GenomeID) %>% filter(type=="CDS") %>% mutate(Gene=gsub("ID=","", attributes) %>% gsub(";..+","", .)) %>% filter(seqid %in% locations$Contig)
    
    hits<-tibble()
    for(i in 1:nrow(locations)){
      hits<-gff %>% filter(seqid==locations$Contig[i]) %>% filter((start<locations$position[i] & end>locations$position[i])) %>% mutate(kmer=locations$kmer[i]) %>% bind_rows(hits)
    }
    
    timp<-
    hits %>% full_join(timp) %>%
      select(FeatureID, kmer, Importance, SD, Rank, GenomeID, Contig=seqid, start, end, strand, attributes) %>%
      arrange(Rank)
    
    return(list(fasta=fa, annotable=timp))
    }
  
}
  
  ############# Main Analysis #############  #############  #############  #############
  observeEvent(input$StartAnalysis, {
    withProgress(message = 'Analyzing Data', value = 0, {    n=10
    message(date(),"-> Cleaning up from last run")
    message("Preclean memory: ", mem_used())
    tmps<<-NULL
    tmps<-NULL
    gc()
    message("Postclean memory: ", mem_used())
    
    
    message(date(),"-> Loading files for main analysis")

    incProgress(1/n, detail = "Loading data...")
    #Get name of phenotype and figure out if user provided their own, then subset table

    if(is.null(input$phenofile)){
      tmps$PhenoName<-input$UserPhenotype
      tmps$PhenoTable<-read_excel("resources/ElenMatchR_phenotypes.xlsx", skip=1) %>% select(Strain_Name, PublicationID, GenomeID, Phenotype=tmps$PhenoName)
      tmps$PhenoTable<-tmps$PhenoTable %>% filter(!is.na(Phenotype)) %>% arrange(Phenotype)
    } else {
      tmps$PhenoName<-read_excel(input$phenofile$datapath, skip=1) %>% colnames(.)
      tmps$PhenoName<-tmps$PhenoName[length(tmps$PhenoName)] #assume it must be last column
      tmps$PhenoTable<-read_excel(input$phenofile$datapath, skip=1) %>% select(Strain_Name, PublicationID, GenomeID, Phenotype=tmps$PhenoName)
      tmps$PhenoTable<-tmps$PhenoTable %>% filter(!is.na(Phenotype)) %>% arrange(Phenotype)
    }
    
    #load matrix and subset matrix
    if(input$Mode=="Genes"){
      tmps$GeneList<-readRDS(paste0("resources/matrices/ID",input$PID,"COV",input$COV,".RDS"))
      tmps$matrix<-tmps$GeneList %>% column_to_rownames("OG")
      tmps$matrix<-apply(tmps$matrix, 2, function(x){if_else(is.na(x), 0, 1)})
      rownames(tmps$matrix)<-tmps$GeneList$OG
      } else {
        tmps$GeneList<-""
        tmps$matrix<-readRDS("resources/kmers/kmer_derep.RDS")
        colnames(tmps$matrix)<-as.character(colnames(tmps$matrix))
        message("Warning-> Kmer mode detected non Eggerthella lenta genomes, these are being removed")
        tmps$PhenoTable<-tmps$PhenoTable %>% filter(grepl("Eggerthella lenta", Strain_Name))
      }
    tmps$matrix<-tmps$matrix[,tmps$PhenoTable$GenomeID] #remove the unneeded columns
    tmps$matrix<-tmps$matrix[rowSums(tmps$matrix)!=ncol(tmps$matrix),]# remove features that are all present
    tmps$matrix<-tmps$matrix[rowSums(tmps$matrix)!=0,]# remove features that are all absent
    
    message("Post data load memory: ", mem_used())
    
    
    ###
    message(date(), "-> Starting analysis of ",tmps$PhenoName, " on ", input$Mode)
    incProgress(2/n, detail = "Running random forests...")
    
    tmps$PhenoTable<-tmps$PhenoTable %>% mutate(Phenotype=factor(Phenotype, levels=unique(Phenotype)))
    tmps$RF<-DoRF(phenofactor=tmps$PhenoTable$Phenotype, pmatrix=tmps$matrix, nrep=input$Nreps, ntree=input$Ntrees)
    ###
    message(date(), "-> Filling in phenotypes")
    incProgress(3/n, detail = "Adding phenotypes...")
    
    output$PhenoTable<-renderDataTable({
      tmps$PhenoTable
      }, 
      extensions = 'Buttons', 
      options=list(pageLength=10,
                   "dom" = 'T<"clear">lBfrtip', 
                   buttons = list('copy', 'csv', 'excel')
                   ),
      rownames=FALSE
      )
    
    ###
    message(date(), "-> Filling in confusion matrices")
    incProgress(4/n, detail = "Adding confusion matrices")
    
    output$ConfusionMatrices<-renderDataTable({
      tmps$RF$ConfusionMatrices
    }, 
    extensions = 'Buttons', 
    options=list(pageLength=2,
                 "dom" = 'T<"clear">lBfrtip', 
                 buttons = list('copy', 'csv', 'excel'),
                 autoWidth=TRUE
                ),
                rownames=FALSE
    )
    
    
    ###
    message(date(), "-> Drawing scatter plots")
    incProgress(5/n, detail = "Drawing scatter plots...")
    
    tmps$DrawScatter<-DrawScatter(tmps$RF$ImportanceTable, input$Scatter_Nfeats, tmps$PhenoName)
    
    output$ImportancePlot <- renderPlot({tmps$DrawScatter$Plot}, height=700)
    output$DLImportancePlot <- downloadHandler(
      filename = function() { paste0("ScatterPlot_", tmps$PhenoName, "_cov", input$COV, "_pid", input$PID, '.pdf') },
      content = function(file) {
        ggsave(file, plot = tmps$DrawScatter$Plot, device = "pdf", width=PWIDTH, height=PHEIGHT, useDingbats=F)
      }
    )
    output$DLImportanceRaw <- downloadHandler(
      filename = function() { paste0("ImportanceScores_", tmps$PhenoName, "_cov", input$COV, "_pid", input$PID, '.tsv') },
      content = function(file) {
       write_tsv(file, x = tmps$DrawScatter$Table)
      }
    )
    
    ###
    message(date(), "-> Drawing Heatmap")
    incProgress(6/n, detail = "Drawing heatmap...")
    
    tmps$DrawHeatMap<-DrawHeatMap(tmps$RF$ImportanceTable, tmps$PhenoTable, input$Heatmap_Nfeats, tmps$matrix, tmps$PhenoName, input$Mode)
    output$Heatmap <- renderPlot({tmps$DrawHeatMap}, height=700)
    output$DLHeatmap <- downloadHandler(
      filename = function() { paste0("Heatmap_", tmps$PhenoName, "_cov", input$COV, "_pid", input$PID, '.pdf') },
      content = function(file) {
        ggsave(file, plot = tmps$DrawHeatMap, device = "pdf", width=PWIDTH, height=PHEIGHT, useDingbats=F)
      }
    )
    

    ###
    message(date(), "-> Drawing MA plot")
    incProgress(7/n, detail = "Drawing MA Plot...")
    
    tmps$DrawMA<-DrawMA(input$RefGenome, tmps$RF$ImportanceTable, input$MA_Nfeats, input$Mode, tmps$GeneList, tmps$PhenoName)
    output$ManhattanPlot <- renderPlot({tmps$DrawMA}, height=700)
    output$DLMAplot <- downloadHandler(
      filename = function() { paste0("MAPlot_", tmps$PhenoName, "_cov", input$COV, "_pid", input$PID, '.pdf') },
      content = function(file) {
        ggsave(file, plot = tmps$DrawMA, device = "pdf", width=PWIDTH, height=PHEIGHT, useDingbats=F)
      }
    )
 
    ###
    message(date(), "-> Drawing Tree")
    incProgress(8/n, detail = "Drawing tree...")
    
    tmps$DrawTree<-DrawTree(input$BranchLengths, input$Orientation, input$Prune, tmps$PhenoTable, tmps$PhenoName, input$TreeSource)
    
    output$TreePlot <- renderPlot({tmps$DrawTree}, height=700)
    output$DLTreePlot <- downloadHandler(
      filename = function() { paste0("PhyloTree_", tmps$PhenoName, '.pdf') },
      content = function(file) {
        ggsave(file, plot = tmps$DrawTree, device = "pdf", width=PWIDTH, height=PHEIGHT, useDingbats=F)
      }
    )   
    output$DLTree <- downloadHandler(
      filename = function() { paste0("PhyloTree_", tmps$PhenoName, '.nwk') },
      content = function(file) {
        write.tree(file, phy=tree)
      }
    )   
    
    ###
    message(date(), "-> Building Annotation Table")
    incProgress(9/n, detail = "Building annotation table and extracting features")
    
    tmps$MakeAnnoTable<-MakeAnnoTable(input$Mode, tmps$RF$ImportanceTable, tmps$GeneList, input$AT_Nfeats, tmps$PhenoTable)

    output$AnnotationTable<-renderDataTable({
      tmps$MakeAnnoTable$annotable
    }, 
    extensions = 'Buttons', 
    options=list(pageLength=10,
                 "dom" = 'T<"clear">lBfrtip', 
                 buttons = list('copy', 'csv', 'excel')
    ),
    rownames=FALSE
    )
  
    output$DLfasta <- downloadHandler(
      filename = function() { paste0("Sequences_", tmps$PhenoName, "_cov", input$COV, "_pid", input$PID, '.fasta') },
      content = function(file) {
        writeXStringSet(file, x=tmps$MakeAnnoTable$fasta)
      }
    )
    
    
    tmps<<-tmps # needed to make tmps available to event observations below
    incProgress(10/n, detail = "Finished!")
    })
    
    message(date(), "-> Main analysis done")
})


  
  
  
  
  
  
  
  
  ########## Update number of features on scatter plot
  observeEvent(input$Scatter_Nfeats, {
    message(date(), "-> Redrawing scatter plots")
    tmps$DrawScatter<-DrawScatter(tmps$RF$ImportanceTable, input$Scatter_Nfeats, tmps$PhenoName)
    
    output$ImportancePlot <- renderPlot({tmps$DrawScatter$Plot}, height=700)
    output$DLImportancePlot <- downloadHandler(
      filename = function() { paste0("ScatterPlot_", tmps$PhenoName, "_cov", input$COV, "_pid", input$PID, '.pdf') },
      content = function(file) {
        ggsave(file, plot = tmps$DrawScatter$Plot, device = "pdf", width=PWIDTH, height=PHEIGHT, useDingbats=F)
      }
    )
    output$DLImportanceRaw <- downloadHandler(
      filename = function() { paste0("ImportanceScores_", tmps$PhenoName, "_cov", input$COV, "_pid", input$PID, '.tsv') },
      content = function(file) {
        write_tsv(file, x = tmps$DrawScatter$Table)
      }
    )
  }, ignoreInit = TRUE)
  
  
  ########## Update number of features on heatmap
  observeEvent(input$Heatmap_Nfeats, {
    message(date(), "-> Redrawing Heatmap")
    tmps$DrawHeatMap<-DrawHeatMap(tmps$RF$ImportanceTable, tmps$PhenoTable, input$Heatmap_Nfeats, tmps$matrix, tmps$PhenoName, input$Mode)
    output$Heatmap <- renderPlot({tmps$DrawHeatMap}, height=700)
    output$DLHeatmap <- downloadHandler(
      filename = function() { paste0("Heatmap_", tmps$PhenoName, "_cov", input$COV, "_pid", input$PID, '.pdf') },
      content = function(file) {
        ggsave(file, plot = tmps$DrawHeatMap, device = "pdf", width=PWIDTH, height=PHEIGHT, useDingbats=F)
      }
    )
  }, ignoreInit = TRUE)
  
  ########## Update MA plot
  
  observeEvent(c(input$RefGenome,input$MA_Nfeats), {
    withProgress(message = 'Redrawing MA plot', value = 0, {    n=10
    message(date(), "-> Redrawing MA plot")
    incProgress(7/n, detail = "Extracting hits")
    
      tmps$DrawMA<-DrawMA(input$RefGenome, tmps$RF$ImportanceTable, input$MA_Nfeats, input$Mode, tmps$GeneList, tmps$PhenoName)
      output$ManhattanPlot <- renderPlot({tmps$DrawMA}, height=700)
      output$DLMAplot <- downloadHandler(
        filename = function() { paste0("MAPlot_", tmps$PhenoName, "_cov", input$COV, "_pid", input$PID, '.pdf') },
        content = function(file) {
          ggsave(file, plot = tmps$DrawMA, device = "pdf", width=PWIDTH, height=PHEIGHT, useDingbats=F)
        }
      )
    })
  }, ignoreInit = TRUE)
  
  
  ########## Update tree
  
  observeEvent(c(input$BranchLengths,input$Orientation, input$Prune, input$TreeSource), {
    message(date(), "-> Drawing Tree")
    tmps$DrawTree<-DrawTree(input$BranchLengths, input$Orientation, input$Prune, tmps$PhenoTable, tmps$PhenoName, input$TreeSource)
    
    output$TreePlot <- renderPlot({tmps$DrawTree}, height=700)
    output$DLTreePlot <- downloadHandler(
      filename = function() { paste0("PhyloTree_", tmps$PhenoName, '.pdf') },
      content = function(file) {
        ggsave(file, plot = tmps$DrawTree, device = "pdf", width=PWIDTH, height=PHEIGHT, useDingbats=F)
      }
    )   
  }, ignoreInit = TRUE)
  
  
  ########## Update table
observeEvent(input$AT_Nfeats, {
  message(date(), "-> Rebuilding Annotation Table")
  withProgress(message = 'Rebuilding Annotation Table', value = 0, {    n=10
  incProgress(7/n, detail = "Extracting hits")
  
  tmps$MakeAnnoTable<-MakeAnnoTable(input$Mode, tmps$RF$ImportanceTable, tmps$GeneList, input$AT_Nfeats, tmps$PhenoTable)
  
  output$AnnotationTable<-renderDataTable({
    tmps$MakeAnnoTable$annotable
  }, 
  extensions = 'Buttons', 
  options=list(pageLength=10,
               "dom" = 'T<"clear">lBfrtip', 
               buttons = list('copy', 'csv', 'excel')
  ),
  rownames=FALSE
  )
  
  output$DLfasta <- downloadHandler(
    filename = function() { paste0("Sequences_", tmps$PhenoName, "_cov", input$COV, "_pid", input$PID, '.fasta') },
    content = function(file) {
      writeXStringSet(file, x=tmps$MakeAnnoTable$fasta)
    }
  )
  })
}, ignoreInit=TRUE)
  

  output$DLtemplate <- downloadHandler(
    filename = function() {"ElenMatchR_tenmplate.xlsx"},
    content = function(file) {
      file.copy(to=file, from="resources/ElenMatchR_template.xlsx")
    }
  )


  
}


# Run the application 
shinyApp(ui = ui, server = server)
