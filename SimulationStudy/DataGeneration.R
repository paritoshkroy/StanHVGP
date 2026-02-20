obj_prev <- ls()
library(fields)
library(Matrix)
library(tidyverse)
library(magrittr)
library(tidybayes)
library(coda)
library(nleqslv)

coords <- unname(as.matrix(expand.grid(x = seq(-0.99, 0.99, length.out = 100), y = seq(-0.99, 0.99, length.out = 100))))
nsite <- nrow(coords); nsite
theta <- c(5,-1,1)
sigma <- 5
sigmasq <- sigma^2
lscale <- 0.5
tau <- 2
tausq <- tau^2

distMat <- fields::rdist(coords)

SigmaX <- 1*0.25^abs(outer(1:2,1:2,'-'))
#set.seed(1988)
X <- cbind(1,cbind(rnorm(n=nsite),rnorm(n=nsite)) %*% t(chol(SigmaX)))
muX <- drop(X %*% theta)
z <- drop(crossprod(chol(matern32(d = fields::rdist(coords), sigma = sigma, lscale = lscale) + diag(x=1e-9, nrow = nsite, ncol = nsite)), rnorm(nsite)))
linpred <- muX + z
y <- rnorm(n = nsite, mean = linpred, sd = tau)
nsize <- 5000
idSampled <- sample.int(n = nsite, size = nsize, replace = FALSE)
#set.seed(NULL)

obj_all <- ls()
obj_keep <- c("idSampled", "y", "z", "X", "theta", "sigma", "lscale", "tau", "coords", "nsize")
obj_drop <- obj_all[!obj_all %in% c(obj_prev, obj_keep)]
rm(list = obj_drop)
gc()
mean(y)
