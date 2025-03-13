library(dplyr)
library(tibble)
library(ggplot2)

df_descriptions = read.csv('PFAM_NCBIFAM_KOFAM_entries.tsv', sep = '\t')

df_pfam = read.csv('pfam_t_test.tsv', sep = '\t')
df_kofam = read.csv('kofam_t_test.tsv', sep = '\t')
df_ncbifam = read.csv('ncbifam_t_test.tsv', sep = '\t')
df_pfam = df_pfam %>% rename(profile = PFAM)
df_kofam = df_kofam %>% rename(profile = KO)
df_ncbifam = df_ncbifam %>% rename(profile = NCBIFam)
df_all = rbind(df_pfam, df_kofam) %>% left_join(df_descriptions, by = join_by(profile == Accession)) %>% mutate(full =
                                                                                                                   paste(Name, '-', profile))

df_all = df_all %>% arrange(desc(Cohen_d), q)
df_all_top = rbind(head(df_all,20), tail(df_all,20))


df_exclusive_meta = rbind(df_pfam, df_kofam) %>% filter(Meta_n >=
                                                                      100 &
                                                                      Isolate_n == 0 &
                                                                      q < 0.05) %>% left_join(df_descriptions, by = join_by(profile == Accession)) %>% mutate(full =
                                                                                                                                                                paste(Name, '-', profile)) %>% arrange(desc(Cohen_d), q)

p1 = ggplot(df_all_top, aes(Cohen_d, reorder(full, Cohen_d), fill = Cohen_d >
                          0)) +
  geom_bar(stat = 'identity') +
  scale_x_continuous(breaks = seq(-2, 2, by = 1), limits = c(-2, 2))+
theme_minimal()

p2 = ggplot(head(df_exclusive_meta,20), aes(Cohen_d, reorder(full, Cohen_d))) +
  geom_bar(stat = 'identity')+
theme_minimal()

print(p1)
print(p2)
