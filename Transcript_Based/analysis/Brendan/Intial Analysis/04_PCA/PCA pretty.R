library(tidyverse)
theme_set(theme_bw(base_size=12)) # set default ggplot2 theme
library(egg) # required to arrange plots side-by-side
library(grid) # required to draw arrows
library(ggthemes) # for colorblind color scale

t(na.omit(diff.genes)) %>% #t() around this part only diff between transposed and untransposed
  scale()%>% 
  prcomp() ->            # do PCA
  PCA    # store result as `pca`

PCA_plot <- data.frame(PCA$x)


ggplot(PCA_plot, aes(x = PC1, y = PC2, color = pheno$pheno ,label = rownames(PCA_plot))) + 
  geom_text(size = 3.0) +
  ggtitle("PCA- Genes: Transposed")

# capture the rotation matrix in a data frame
rotation_data_genes <- data.frame(
  PCA$rotation, 
  variable = rownames(PCA$rotation)
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
  ggtitle("PCA rotation- Genes: Transposed")

