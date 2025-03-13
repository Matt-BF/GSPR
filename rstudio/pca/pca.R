library(PCAtools)
library(ggplot2)
library(ggokabeito)

df <- read.csv('pfam_pivot_for_pca.csv', row.names = 1)

metadata <- data.frame(row.names = colnames(df))

# Perform PCA using PCAtools
p <- pca(df, metadata=metadata, removeVar = 0.1)  # Adjust removeVar as needed


# Create the biplot
print(biplot(p, showLoadings = TRUE, pointSize = 5, legendPosition = 'right'))


df2 <- read.csv('kofam_pivot_for_pca.csv', row.names = 1)

metadata <- data.frame(row.names = colnames(df2))

# Perform PCA using PCAtools
p2 <- pca(df2, metadata=metadata, removeVar = 0.1)  # Adjust removeVar as needed


# Create the biplot
print(biplot(p2, showLoadings = TRUE, pointSize = 5, legendPosition = 'right'))

#print(screeplot(p))
