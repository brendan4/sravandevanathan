library(tidyverse)
theme_set(theme_bw(base_size=12)) # set default ggplot2 theme
library(egg) # required to arrange plots side-by-side
library(grid) # required to draw arrows
library(ggthemes) # for colorblind color scale

t(na.omit(expressed.genes)) %>%  #t() around this part only diff between transposed and untransposed
  scale() %>%            # scale to 0 mean and unit variance
  prcomp() ->            # do PCA
  PCA_expressed.genes2   # store result as `pca`

PCA_expressed.genes2_data <- data.frame(PCA_expressed.genes2$x)

ggplot(PCA_expressed.genes2_data, aes(x = PC1, y = PC2, label = rownames(PCA_expressed.genes2_data))) + 
  geom_text(size = 3.0) +
  scale_color_colorblind()+
  ggtitle("PCA- Genes: untransposed")

# capture the rotation matrix in a data frame
rotation_data_genes <- data.frame(
  PCA_expressed.genes2$rotation, 
  variable = row.names(PCA_expressed.genes2$rotation)
)

# define a pleasing arrow style
arrow_style <- arrow(
  length = unit(0.05, "inches"),
  type = "closed"
)

# now plot, using geom_segment() for arrows and geom_text() for labels
ggplot(rotation_data_genes) + 
  geom_segment(aes(xend = PC1, yend = PC2), x = 0, y = 0, arrow = arrow_style) + 
  geom_text(aes(x = PC1, y = PC2, label = variable), hjust = 0, size = 3, color = "red")+
  xlim(-.25, 1.) + 
  ylim(-.50, .50) +
  coord_fixed() +
  ggtitle("PCA rotation- Genes: untransposed")
