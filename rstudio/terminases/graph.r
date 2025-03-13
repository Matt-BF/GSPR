library(ggplot2)
library(dplyr)
library(cowplot)
df = read.csv('report.tsv', sep = '\t', row.names = 1)

head(df)

lens = ggplot(df) +
  geom_histogram(aes(seq_len), fill='blue')

gap_rates = ggplot(df) +
  geom_histogram(aes(gap_rate))
median(df$gap_rate)
head(df)
plot_grid(lens,gap_rates)

#subsetted = df[,c(1,2)]

#subsetted = subsetted %>% mutate(gap_rate=gaps/seq_len)


