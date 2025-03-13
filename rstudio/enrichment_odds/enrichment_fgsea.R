library(hypeR)
library(dplyr)
library(tibble)
require(pals)
suppressPackageStartupMessages(library(fgsea))
library(cowplot)
library(ggplot2)

##fgsea hack
.handle.genesets <- function(genesets) {
  if (is(genesets, "list")) {
    gsets.obj <- gsets$new(genesets, quiet=TRUE)
  }
  else if (is(genesets, "gsets") | is(genesets, "rgsets")) {
    gsets.obj <- genesets
  } 
  else {
    stop("Genesets must be gsets/rgsets object or named list of genesets")
  }
  return(gsets.obj)
}

fgsea.wrapper <- function(signature, genesets, sample.size=101, min.size=1, max.size=Inf, ...) {
  # Save original arguments
  args <- as.list(environment())
  
  # Save gsets object
  gsets.obj <- .handle.genesets(genesets)
  args$genesets <- gsets.obj
  
  # Run fgsea
  results <- fgsea::fgseaMultilevel(stats=signature, 
                                    pathways=gsets.obj$genesets, 
                                    sampleSize=sample.size, 
                                    minSize=min.size, 
                                    maxSize=max.size, 
                                    ...)
  
  data <- results %>%
    data.frame() %>%
    plyr::rename(c("pathway"="label", "padj"="fdr", "log2err"="lte", "size"="overlap", "leadingEdge"="le")) %>%
    dplyr::rename_with(tolower) %>%
    mutate(pval=signif(pval, 2)) %>%
    mutate(fdr=signif(fdr, 2)) %>%
    mutate(le=sapply(le, function(x) paste(x, collapse=','))) %>%
    mutate(signature=length(signature)) %>%
    mutate(geneset=sapply(label, function(x) length(gsets.obj$genesets[[x]]))) %>%
    dplyr::select(c("label", "pval", "fdr", "lte", "es", "nes", "signature", "geneset", "overlap", "le"))
  
  data.up <- data %>%
    dplyr::filter(es > 0) %>%
    dplyr::arrange(pval, es)
  
  data.dn <- data %>%
    dplyr::filter(es < 0) %>%
    dplyr::arrange(pval, es)    
  
  # Reproducibility information
  info <- list(fgsea=paste("v", packageVersion("fgsea"), sep=""),
               signature=length(signature), 
               genesets=args$genesets$info())
  
  info <- c(info, args[c("sample.size", "min.size", "max.size")])
  info <- lapply(info, as.character)
  
  # Wrap dataframe in hyp object
  hyp.up <- hyp$new(data=data.up, args=args, info=info)
  hyp.dn <- hyp$new(data=data.dn, args=args, info=info)
  mhyp <- multihyp$new(data=list("Meta"=hyp.up, "Isolate"=hyp.dn))
  return(mhyp)
}
#PFAM-GO
enrichment_table <- read.csv('pfam_odds.tsv', sep='\t')
pfam_vec <- as.vector(enrichment_table$PFAM)
signature <- enrichment_table %>% arrange(desc(log_odds)) %>% select(PFAM, log_odds) %>% deframe()
#signature <- as.vector(signature)


pfam2go_table <- read.csv('pfam2GO.tsv', sep='\t') 

grp <- pfam2go_table %>% filter(PFAM %in% pfam_vec) %>% group_by(GO_term) %>% summarise(PFAM_list = list(PFAM)) %>% deframe()

genesets <- gsets$new(grp, name='PFAM2GO', version='v1.0')

mhyp_obj <- fgsea.wrapper(signature, genesets)

a <-rbind(mhyp_obj$data$Meta$data, mhyp_obj$data$Isolate$data)
a <- a %>% mutate(Origin=if_else(es<0,'Isolate',"Meta"))
a <- a %>% filter(fdr <0.05)
a <- a %>% arrange(desc(nes), fdr)
b <- a %>% arrange(nes,fdr)
top_10 <- rbind(head(a,10), head(b,10))


p1 <- (ggplot(top_10, mapping=aes(Origin,label))+
        geom_point(aes(colour = nes, size=-log10(fdr)))+
        scale_colour_gradientn(colours=coolwarm(100), guide = "colourbar")+
        labs(title='PFAM-GO-FGSEA')
      +theme_minimal())



#KEGG Modules
enrichment_table <- read.csv('kofam_odds.tsv', sep='\t')
ko_vec <- as.vector(enrichment_table$KO)
signature <- enrichment_table %>% arrange(desc(log_odds)) %>% select(KO,log_odds) %>% deframe()
#signature <- as.vector(signature)

ko2go_table <- read.csv('KO2Module.tsv', sep='\t') 

grp <- ko2go_table %>% filter(KO %in% ko_vec) %>% group_by(Module_desc) %>% summarise(KO_list = list(KO)) %>% deframe()

genesets <- gsets$new(grp, name='KO2Module', version='v1.0')

mhyp_obj <- fgsea.wrapper(signature, genesets)

a <-rbind(mhyp_obj$data$Meta$data, mhyp_obj$data$Isolate$data)
a <- a %>% mutate(Origin=if_else(es<0,'Isolate',"Meta"))
a <- a %>% filter(fdr <0.05)
a <- a %>% arrange(desc(nes), fdr)
b <- a %>% arrange(nes,fdr)
top_10 <- rbind(head(a,10), head(b,10))


p2 <- (ggplot(top_10, mapping=aes(Origin,label))+
        geom_point(aes(colour = nes, size=-log10(fdr)))+
        scale_colour_gradientn(colours=coolwarm(100), guide = "colourbar")+
        labs(title='KEGG-Module-FGSEA')
      +theme_minimal())

#KEGG-GO
enrichment_table <- read.csv('kofam_odds.tsv', sep='\t')
ko_vec <- as.vector(enrichment_table$KO)
signature <- enrichment_table %>% arrange(desc(log_odds)) %>% select(KO, log_odds) %>% deframe()

ko2go_table <- read.csv('KO2GO.tsv', sep = '\t')

grp <- ko2go_table %>% filter(KO %in% ko_vec) %>% group_by(GO_term) %>% summarise(KO_list = list(KO)) %>% deframe()

genesets <- gsets$new(grp, name = 'KO2GO', version = 'v1.0')

mhyp_obj <- fgsea.wrapper(signature, genesets)

a <-rbind(mhyp_obj$data$Meta$data, mhyp_obj$data$Isolate$data)
a <- a %>% mutate(Origin=if_else(es<0,'Isolate',"Meta"))
a <- a %>% filter(fdr <0.05)
a <- a %>% arrange(desc(nes), fdr)
b <- a %>% arrange(nes,fdr)
top_10 <- rbind(head(a,10), head(b,10))


p3 <- (ggplot(top_10, mapping=aes(Origin,label))+
        geom_point(aes(colour = nes, size=-log10(fdr)))+
        scale_colour_gradientn(colours=coolwarm(100), guide = "colourbar")+
        labs(title='KEGG-GO-FGSEA')
      +theme_minimal())

#KEGG-map
enrichment_table <- read.csv('kofam_odds.tsv', sep='\t')
ko_vec <- as.vector(enrichment_table$KO)
signature <- enrichment_table %>% arrange(desc(log_odds)) %>% select(KO, log_odds) %>% deframe()

ko2go_table <- read.csv('KO2map.tsv', sep = '\t')

grp <- ko2go_table %>% filter(KO %in% ko_vec) %>% group_by(pathway_description) %>% summarise(KO_list = list(KO)) %>% deframe()

genesets <- gsets$new(grp, name = 'KO2map', version = 'v1.0')

mhyp_obj <- fgsea.wrapper(signature, genesets)

a <-rbind(mhyp_obj$data$Meta$data, mhyp_obj$data$Isolate$data)
a <- a %>% mutate(Origin=if_else(es<0,'Isolate',"Meta"))
a <- a %>% filter(fdr <0.05)
a <- a %>% arrange(desc(nes), fdr)
b <- a %>% arrange(nes,fdr)
top_10 <- rbind(head(a,10), head(b,10))


p4 <- (ggplot(top_10, mapping=aes(Origin,label))+
        geom_point(aes(colour = nes, size=-log10(fdr)))+
        scale_colour_gradientn(colours=coolwarm(100), guide = "colourbar")+
        labs(title='KO-map_FGSEA')
      +theme_minimal())


print(p1)
ggsave2("odds_pfam-go_fgsea.pdf", width=10, height = 10)
print(p2)
ggsave2("odds_kegg-module_fgsea.pdf", width=10, height = 10)
print(p3)
ggsave2("odds_kegg-go_fgsea.pdf", width=10, height = 10)
print(p4)
ggsave2("odds_kegg-map_fgsea.pdf", width=10, height = 10)
