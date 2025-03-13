library(dplyr)
library(tibble)
library(ggplot2)

df_descriptions = read.csv('PFAM_NCBIFAM_KOFAM_entries.tsv', sep = '\t')

df_pfam = read.csv('pfam_odds.tsv', sep = '\t')
df_kofam = read.csv('kofam_odds.tsv', sep = '\t')
df_ncbifam = read.csv('ncbifam_odds.tsv', sep = '\t')
df_pfam = df_pfam %>% rename(profile = PFAM) %>% left_join(df_descriptions, by = join_by(profile == Accession)) %>% mutate(full =paste(Name, '-', profile)) %>% arrange(desc(log_odds), padj)
df_kofam = df_kofam %>% rename(profile = KO) %>% left_join(df_descriptions, by = join_by(profile == Accession)) %>% mutate(full =paste(Name, '-', profile)) %>% arrange(desc(log_odds), padj)
df_ncbifam = df_ncbifam %>% rename(profile = NCBIFam) %>% left_join(df_descriptions, by = join_by(profile == Accession)) %>% mutate(full =paste(Name, '-', profile)) %>% arrange(desc(log_odds), padj)
df_all = rbind(df_pfam, df_kofam, df_ncbifam) 



df_all = df_all %>% arrange(desc(log_odds), padj)
write_csv(df_all, './all_ods.csv')
df_all_signif = df_all %>% filter(padj<0.05) %>% arrange(desc(log_odds), padj)
df_all_top = rbind(head(df_all_signif,30), tail(df_all_signif,30))

#df_exclusive_meta = rbind(df_pfam, df_kofam) %>% filter(Isolate == 0 & padj < 0.05) %>% left_join(df_descriptions, by = join_by(profile == Accession)) %>% mutate(full =
                                                                                                                                                                #paste(Name, '-', profile)) %>% arrange(desc(log_odds), padj)

p1 = (ggplot(df_all_top, mapping=aes(log_odds,reorder(full,log_odds), fill=log_odds >0))+
  geom_bar(stat='identity')+
    scale_x_continuous(breaks = seq(-10, 10, by = 2), limits = c(-10, 10))+
  theme_minimal())

#p2 = (ggplot(head(df_exclusive_meta,20), aes(log_odds, reorder(full, log_odds))) +
#  geom_bar(stat = 'identity')+
#theme_minimal())

#p3 = (ggplot(rbind(head(df_pfam,20), tail(df_pfam,20)), mapping=(aes(log_odds,reorder(full,log_odds), fill=log_odds >0)))+
#        geom_bar(stat='identity')+
#        scale_x_continuous(breaks = seq(-8, 8, by = 2), limits = c(-8, 8))+
#        theme_minimal())

#print(p1)
#print(p2)
print(p1)

