library(tidyverse)
library(ggbump)

df1 <- read_tsv('cog_categories.tsv') %>% mutate(Origin=if_else(Origin=='Rank_meta',2020,2021))


year <- rep(2019:2021, 4)
position <- c(4, 2, 2, 3, 1, 4, 2, 3, 1, 1, 4, 3)
player <- c("A", "A", "A",
            "B", "B", "B", 
            "C", "C", "C",
            "D", "D", "D")

df <- data.frame(x = year,
                 y = position,
                 group = player)


print(ggplot(df, aes(x = x, y = y, color = group)) +
  geom_bump(size = 1.5) +
  geom_point(size = 6) +
  geom_text(data = df %>% filter(x == min(x)),
            aes(x = x - 0.1, label = group),
            size = 5, hjust = 1) +
  geom_text(data = df %>% filter(x == max(x)),
            aes(x = x + 0.1, label = group),
            size = 5, hjust = 0) +
  scale_color_brewer(palette = "RdBu") +
  theme_void() +
  theme(legend.position = "none"))

print(ggplot(df1, aes(Origin, Rank, color=Category_description)) +
        geom_bump(size = 1.5) +
        geom_point(size = 6) +
        geom_text(data = df1 %>% filter(Origin == min(Origin)),
                  aes(x = Origin - 0.1, label = Category_description),
                  size = 5, hjust = 1) +
        geom_text(data = df1 %>% filter(Origin == max(Origin)),
                  aes(x = Origin + 0.1, label = Category_description),
                  size = 5, hjust = 0) +
        scale_y_reverse()+
        theme_void()+
        theme(legend.position = 'none'))
ggsave2("test.pdf")
