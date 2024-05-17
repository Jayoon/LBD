source('functionsBonf.R')
## Real data analysis for paper
### DNA CNV -------------
library(cumSeg)
data(fibroblast)
df <- fibroblast[, c("Chromosome", "Genome.Order", "gm05296")]
df <- df[df$Chromosome <= 23, ]
df <- na.omit(df)
x <-df$gm05296
alpha <- 0.05

ci_results <- lbd(x, alpha = alpha)
ci_results$disjoint