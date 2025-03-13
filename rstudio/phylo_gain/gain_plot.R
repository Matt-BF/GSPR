library(ggplot2)
library(tidyr)

df = read.csv('phylo_gain.csv')

df_long = df %>% pivot_longer(cols = c('soil_vs_all', 'soil_vs_soil'))


print(
  ggplot(df_long) +
    geom_linerange(
      aes(
        x = value,
        xmin = 0,
        xmax = value,
        y = Gene,
        color = name
      ),
      position = position_dodge(width = 0.3)
    ) +
    geom_point(aes(
      x = value, y = Gene, color = name
    ), position = position_dodge(width = 0.3))
  + ylab("") +
    xlab('Phylogenetic Gain')
  + theme_minimal()
)

print(
  ggplot(df_long, aes(
    y = Gene, x = value, fill = name
  )) +
    geom_bar(
      position = position_dodge(width = 0.2),
      stat = 'identity',
      width = 0.2
    )
  + ylab("") +
    xlab('Phylogenetic Gain')
  + theme_minimal()
)

