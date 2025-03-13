library(iNEXT)
library(ggplot2)
library(ggokabeito)
ptus_ab = read.csv('ptu_abundance_per_ecosystem.tsv', sep='\t', row.names = 1)
out = iNEXT(ptus_ab, endpoint=30000, nboot = 100)

ggiNEXT(out) +
  theme(legend.position = 'right') +
  scale_shape_manual(values=c(19,19,19,19,19,19,19,19,19))+
  scale_colour_okabe_ito()+
  scale_fill_okabe_ito() +
  theme_minimal() +
  labs(x="Number of plasmids", y='Number of PTUs')
