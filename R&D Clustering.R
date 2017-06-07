library(ggplot2)
library(clusterSim)
library(data.table)
library(MixSim)
options(scipen = 100)


set.seed(1234)
Q <- MixSim(BarOmega = 0.01, K = 3, p = 2)
A <- simdataset(n = 5000, Pi = Q$Pi, Mu = Q$Mu, S = Q$S, n.out = 10)
colors <- c("red", "green", "blue", "brown", "magenta")

par(mar = c(0.1, 0.1, 0.1, 0.1))
plot(A$X, col = colors[A$id], pch = 19, cex = 0.8, xlab = "", ylab = "", axes = FALSE)
box()

