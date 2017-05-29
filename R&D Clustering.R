library(ggplot2)
library(clusterSim)
library(data.table)
options(scipen = 100)
grnd <- cluster.Gen(c(50,50), model=3, means=c(8,11,14),
                    cov=c(1,1),numNoisyVar=1)
clusterdata <- data.frame(grnd$data)
clusterdata$X3 <- rnorm(n = 100,mean = 4,sd = 1)
clusterdata$X4 <- abs(rexp(100))





set.seed(1234)
Q <- MixSim(BarOmega = 0.01, K = 4, p = 2)
A <- simdataset(n = 500, Pi = Q$Pi, Mu = Q$Mu, S = Q$S, n.out = 10)
colors <- c("red", "green", "blue", "brown", "magenta")

par(mar = c(0.1, 0.1, 0.1, 0.1))
plot(A$X, col = colors[A$id], pch = 19, cex = 0.8, xlab = "", ylab = "", axes = FALSE)
box()

library("MixSim")