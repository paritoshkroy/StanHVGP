rm(list=ls())
graphics.off()
library(fields)
library(Matrix)
library(tidyverse)
library(magrittr)
library(tidybayes)
library(coda)
library(nleqslv)
library(cmdstanr)

###########################################################################
# Local PC
###########################################################################
fpath <- "/home/ParitoshKRoy/git/ApproximateGP/HVGP/"
node <- 1
##########################################################################
# ARC Preparation
##########################################################################
fpath <- "/home/pkroy/ApproximateGP/HVGP/"
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
}
node <- as.numeric(args[1]) - 100
cat("Node is ", node, "\n")

######################################################################
# Generating the data
######################################################################
source(paste0(fpath,"UtilityFunctions.R"))
source(paste0(fpath,"DataGeneration.R"))

######################################################################
# partition as observed and predicted
######################################################################
obsCoords <- coords[idSampled,]
prdCoords <- coords[-idSampled,]
obsY <- y[idSampled]
prdY <- y[-idSampled]
obsX <- X[idSampled,]
prdX <- X[-idSampled,]
obsZ <- z[idSampled]
prdZ <- z[-idSampled]

obsDistMat <- fields::rdist(obsCoords)
str(obsDistMat)
obsDistVec <- obsDistMat[lower.tri(obsDistMat, diag = FALSE)]
obsMaxDist <- max(obsDistVec); obsMaxDist
obsMedDist <- median(obsDistVec); obsMedDist
obsMinDist <- min(obsDistVec); obsMinDist
rm(obsDistMat)

################################################################################
## HVGP Preparation
################################################################################
nNeighbors <- 10
mra.spec <- GPvecchia::vecchia_specify(locs = obsCoords, m = nNeighbors, conditioning = 'mra', ordering = "maxmin", verbose = TRUE)
obsY <- obsY[mra.spec$ord]
obsCoords <- mra.spec$locsord
obsX <- obsX[mra.spec$ord,]
NN_inds <- mra.spec$U.prep$revNNarray
NN_star <- sapply(1:nsize, function(i) min(which(!is.na(mra.spec$U.prep$revNNarray[i,]))))
Sraw <- Matrix::sparseMatrix(i = mra.spec$U.prep$rowpointers, j = mra.spec$U.prep$colindices, x = 1, dims = c(mra.spec$U.prep$size, mra.spec$U.prep$size))
dim(Sraw)
sparseS <- Sraw[-seq(2,mra.spec$U.prep$size,by=2),-seq(2,mra.spec$U.prep$size,by=2)]
str(sparseS)

NN_ind <- mra.spec$U.prep$revNNarray
dim(NN_ind)
NN_ind[is.na(NN_ind)] <- -1
startingID <- sapply(1:nsize, function(i) min(which(NN_ind[i,]!=-1)))
endingID <- sapply(1:nsize, function(i) max(which(NN_ind[i,]!=-1)))

## Position for tau
sparseDistMat <- rdist(obsCoords)*sparseS
diag(sparseDistMat) <- 111
obsDistVec <- rstan::extract_sparse_parts(sparseDistMat)$w
tauID <- which(obsDistVec == 111)
obsDistVec[tauID] <- 0
nnZero <- length(obsDistVec)

## Prior elicitation: Prior Set 2
lLimit <- quantile(obsDistVec[obsDistVec!=0], prob = 0.05)/2.75; lLimit
uLimit <- quantile(obsDistVec[obsDistVec!=0], prob = 0.99)/2.75; uLimit

## Exponential and PC prior
lambda_sigma <- -log(0.01)/sd(y); lambda_sigma
lambda_tau <- -log(0.01)/sd(y); lambda_tau

## Pr(tau > sd(y)) = 0.05 and Pr(tau > sd(y)) = 0.05
pexp(q = 1, rate = lambda_tau, lower.tail = TRUE)

## P(ell < lLimit) = 0.05
lambda_ell <- as.numeric(-log(0.01)*lLimit); lambda_ell
pfrechet(q = lLimit, alpha = 1, sigma = lambda_ell, lower.tail = TRUE)
quantile(rfrechet(n = 1000, alpha = 1, sigma = lambda_ell))

## Stan input
P <- 3
mu_theta <- c(mean(obsY),rep(0,P-1))
mu_theta
V_theta <- diag(c(100,rep(1,P-1)))

# Keep in mind that the data should be ordered following nearest neighbor settings
input <- list(N = nsize, K = ncol(NN_ind)-1, 
              P = P, y = obsY, X = obsX, 
              coords = obsCoords, 
              S = as.matrix(sparseS), 
              neiID = NN_ind, startingID = startingID, 
              endingID = endingID, tauID = tauID, nnZero = nnZero, 
              distVec = obsDistVec, mu_theta = mu_theta, 
              V_theta = V_theta, 
              #Prior Set 2: Exponential and Frechet
              rate_sigma = lambda_sigma, 
              rate_tau = lambda_tau, 
              scale_ell = lambda_ell)
str(input)
library(cmdstanr)
set_cmdstan_path(path = "~/cmdstan-2.36.0")
stan_file <- paste0(fpath,"HVGP_PriorSet2.stan")
mod <- cmdstan_model(stan_file, compile = TRUE)
mod$check_syntax(pedantic = TRUE)
mod$print()
cmdstan_fit <- mod$sample(data = input, 
                          chains = 4,
                          parallel_chains = 4,
                          iter_warmup = 1000,
                          iter_sampling = 1000)
elapsed_time <- cmdstan_fit$time()
elapsed_time
elapsed_time$total/60

cmdstan_fit$cmdstan_diagnose()
sampler_diag <- cmdstan_fit$sampler_diagnostics(format = "df")
str(sampler_diag)
table(sampler_diag$divergent__)
max(sampler_diag$treedepth__)

## Posterior summaries
pars <- c(paste0("theta[",1:P,"]"),"sigma","ell","tau")
pars_true_df <- tibble(variable = pars, true = c(theta,sigma,lscale,tau))
fit_summary <- cmdstan_fit$summary(NULL, c("mean","sd","quantile50","quantile2.5","quantile97.5","rhat","ess_bulk","ess_tail"))
fixed_summary <- inner_join(pars_true_df, fit_summary)
fixed_summary %>% print(digits = 3)

## Posterior draws
draws_df <- cmdstan_fit$draws(format = "df")
draws_df

library(bayesplot)
color_scheme_set("brewer-Spectral")
mcmc_trace(draws_df,  pars = pars, facet_args = list(ncol = 3)) + facet_text(size = 15)

save(fit_summary, file = paste0(fpath, "HVGPSimulationStudy10NeighPriorSet2", node, ".RData"))

