library(hypeR)
library(dplyr)
library(tibble)

enrichment_table <- read.csv('pfam_t_test.tsv', sep = '\t')
pfam_vec <- as.vector(enrichment_table$PFAM)
signature_meta <- enrichment_table %>% filter(q < 0.001 &
                                                Cohen_d > 0) %>% select(PFAM) %>% deframe()
signature_isolate <- enrichment_table %>% filter(q < 0.001 &
                                                   Cohen_d < 0) %>% select(PFAM) %>% deframe()

pfam2go_table <- read.csv('pfam2GO.tsv', sep = '\t')

grp <- pfam2go_table %>% filter(PFAM %in% pfam_vec) %>% group_by(GO_term) %>% summarise(PFAM_list = list(PFAM)) %>% deframe()

genesets <- gsets$new(grp, name = 'PFAM2GO', version = 'v1.0')

hyp_obj_meta <- hypeR(signature_meta, genesets, test = "hypergeometric", fdr =
                        0.01)
hyp_obj_isolate <- hypeR(signature_isolate, genesets, test = "hypergeometric", fdr =
                           0.01)

meta_df <- hyp_obj_meta$data %>% mutate(
  origin = 'Meta',
  o_by_s = overlap / signature,
  geneset_total = geneset / length(pfam_vec),
  effect = o_by_s / geneset_total
) %>% arrange(desc(effect), fdr)
isolate_df <- hyp_obj_isolate$data %>% mutate(
  origin = 'Isolate',
  o_by_s = overlap / signature,
  geneset_total = geneset /
    length(pfam_vec),
  effect = o_by_s / geneset_total
) %>% arrange(desc(effect), fdr)

consolidated <- rbind(head(meta_df, 10), head(isolate_df, 10))


#print(ggplot(consolidated, aes(origin, label)) +
#        geom_point(aes(colour = fdr, size = effect)) +
#        ylab("GO"))

print(ggplot(consolidated, aes(fdr, label)) +
        geom_point(aes(colour = origin, size = effect)) +
        ylab("GO")
      +theme_minimal())

#print(hyp_dots(hyp_obj_meta))
#print(hyp_dots(hyp_obj_isolate))

#KEGG Modules
library(hypeR)
library(dplyr)
library(tibble)

enrichment_table <- read.csv('kofam_t_test.tsv', sep = '\t')
kofam_vec <- as.vector(enrichment_table$KO)
signature_meta <- enrichment_table %>% filter(q < 0.001 &
                                                Cohen_d > 0) %>% select(KO) %>% deframe()
signature_isolate <- enrichment_table %>% filter(q < 0.001 &
                                                   Cohen_d < 0) %>% select(KO) %>% deframe()

ko2module_table <- read.csv('KO2Module.tsv', sep = '\t')

grp <- ko2module_table %>% filter(KO %in% kofam_vec) %>% group_by(Module_desc) %>% summarise(KO_list = list(KO)) %>% deframe()

genesets <- gsets$new(grp, name = 'KO2Module', version = 'v1.0')

hyp_obj_meta <- hypeR(signature_meta, genesets, test = "hypergeometric", fdr =
                        0.01)
hyp_obj_isolate <- hypeR(signature_isolate, genesets, test = "hypergeometric", fdr =
                           0.01)

meta_df <- hyp_obj_meta$data %>% mutate(
  origin = 'Meta',
  o_by_s = overlap / signature,
  geneset_total = geneset / length(kofam_vec),
  effect = o_by_s / geneset_total
) %>% arrange(desc(effect), fdr)
isolate_df <- hyp_obj_isolate$data %>% mutate(
  origin = 'Isolate',
  o_by_s = overlap / signature,
  geneset_total = geneset /
    length(kofam_vec),
  effect = o_by_s / geneset_total
) %>% arrange(desc(effect), fdr)


consolidated <- rbind(head(meta_df, 10), head(isolate_df, 10))

#print(ggplot(consolidated, aes(origin, label)) +
#        geom_point(aes(colour = fdr, size = effect)) +
#        ylab("Module"))

print(ggplot(consolidated, aes(fdr, label)) +
        geom_point(aes(colour = origin, size = effect)) +
        ylab("Module")
      +theme_minimal())

#print(hyp_dots(hyp_obj_meta))
#print(hyp_dots(hyp_obj_isolate))


#KEGG GO
library(hypeR)
library(dplyr)
library(tibble)

enrichment_table <- read.csv('kofam_t_test.tsv', sep = '\t')
kofam_vec <- as.vector(enrichment_table$KO)
signature_meta <- enrichment_table %>% filter(q < 0.001 &
                                                Cohen_d > 0) %>% select(KO) %>% deframe()
signature_isolate <- enrichment_table %>% filter(q < 0.001 &
                                                   Cohen_d < 0) %>% select(KO) %>% deframe()

ko2module_table <- read.csv('KO2GO.tsv', sep = '\t')

grp <- ko2module_table %>% filter(KO %in% kofam_vec) %>% group_by(GO_term) %>% summarise(KO_list = list(KO)) %>% deframe()

genesets <- gsets$new(grp, name = 'KO2GO', version = 'v1.0')

hyp_obj_meta <- hypeR(signature_meta, genesets, test = "hypergeometric", fdr =
                        0.01)
hyp_obj_isolate <- hypeR(signature_isolate, genesets, test = "hypergeometric", fdr =
                           0.01)

meta_df <- hyp_obj_meta$data %>% mutate(
  origin = 'Meta',
  o_by_s = overlap / signature,
  geneset_total = geneset / length(kofam_vec),
  effect = o_by_s / geneset_total
) %>% arrange(desc(effect), fdr)
isolate_df <- hyp_obj_isolate$data %>% mutate(
  origin = 'Isolate',
  o_by_s = overlap / signature,
  geneset_total = geneset /
    length(kofam_vec),
  effect = o_by_s / geneset_total
) %>% arrange(desc(effect), fdr)


consolidated <- rbind(head(meta_df, 10), head(isolate_df, 10))

print(ggplot(consolidated, aes(fdr, label)) +
        geom_point(aes(colour = origin, size = effect)) +
        ylab("GO")
      +theme_minimal())
