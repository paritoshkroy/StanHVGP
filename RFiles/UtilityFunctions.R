matern32 <- function(d, sigma, lscale){
  ds <- sqrt(3)*d/lscale
  (sigma^2) * (1 + ds) * exp(-ds)
}
matern12 <- function(d, sigma, lscale){
  ds <- d/lscale
  (sigma^2) * exp(-ds)
}

############################################################################
## Prior elicitation for scale parameter
############################################################################
# P(lower_pars < pars < upper_pars) = probability

getIGamma <- function(par, lRange, uRange, prob) {
  
  lProb <- (1 - prob)/2
  uProb <- 1 - lProb
  
  invlRange <- 1/lRange
  invuRange <- 1/uRange
  
  c(qgamma(p = lProb, shape = par[1], rate = par[2], lower.tail = TRUE) - invuRange, qgamma(p = uProb, shape = par[1], rate = par[2], lower.tail = TRUE) - invlRange)
}

############################################################################
## Frechet distribution
############################################################################

dfrechet <- function(x, alpha, sigma, log = FALSE){
  ldensity <- dweibull(1/x, shape = alpha, scale = 1/sigma, log = TRUE) - 2 * log(x)
  if (log){
    out <- ldensity
  } else{
    out <- exp(ldensity)
  }
  return(out)
}

qfrechet <- function(p, alpha, sigma, lower.tail = TRUE, log.p = FALSE){
  1/qweibull(1 - p, shape = alpha, scale = 1/sigma, lower.tail = lower.tail, log.p = log.p)
}

rfrechet <- function(n, alpha, sigma){
  1/rweibull(n, shape = alpha, scale = 1/sigma)
}

pfrechet <- function(q, alpha, sigma, lower.tail = TRUE, log.p = FALSE){
  pweibull(1/q, shape = alpha, scale = 1/sigma, lower.tail = !lower.tail, log.p = log.p)
}

############################################################################
## Inverse Gamma distribution
############################################################################
dinvgamma <- function (x, shape, scale, log = FALSE) {
  ldensity <- dgamma(1/x, shape, rate = scale, log = TRUE) - 2 * log(x)
  if (log){
    out <- ldensity
  } else{
    out <- exp(ldensity)
  }
  return(out)
}

qinvgamma <- function(p, shape, scale, lower.tail = TRUE, log.p = FALSE){
  1/qgamma(1 - p, shape, rate = scale, lower.tail = lower.tail, log.p = log.p)
}

rinvgamma <- function(n, shape, scale){
  1/rgamma(n, shape, rate = scale)
}

pinvgamma <- function(q, shape, scale, lower.tail = TRUE, log.p = FALSE){
  pgamma(1/q, shape, rate = scale, lower.tail = !lower.tail, log.p = log.p)
}


quantile2.5 <- function(x) quantile(x, prob = 0.025)
quantile50 <- function(x) quantile(x, prob = 0.50)
quantile97.5 <- function(x) quantile(x, prob = 0.975)


### GPvecchia function
## evaluate Vecchia likelihood based on U
vecchia_likelihood_U=function(z,U.obj) {
  
  ### output: loglikelihood (for z)
  U=U.obj$U
  latent=U.obj$latent
  zord=z[U.obj$ord.z]
  
  # constant
  const=sum(!latent)*log(2*pi)
  
  # numerator
  z1=Matrix::crossprod(U[!latent,],zord) # standardization using tau
  quadform.num=sum(z1^2)
  logdet.num=-2*sum(log(Matrix::diag(U))) # it for precision
  
  # denominator
  if(sum(latent)==0){ # no latents -> denominator not needed
    
    logdet.denom=quadform.denom=0
    
  } else {  # if latents, need denominator
    
    U.y=U[latent,]
    z2=as.numeric(U.y%*%z1)
    V.ord=U2V(U.obj)   # main factor when tau != 0
    z3=Matrix::solve(V.ord,rev(z2),system='L') # main factor when tau != 0
    quadform.denom=sum(z3^2)
    logdet.denom=-2*sum(log(Matrix::diag(V.ord)))
    
  }
  
  # putting everything together
  neg2loglik=logdet.num-logdet.denom+quadform.num-quadform.denom+const
  loglik=-neg2loglik/2
  return(loglik)
  
}

revMat <- function(mat) mat[nrow(mat):1,ncol(mat):1,drop=FALSE]

U2V <- function(U.obj){
  ### when changing this function make sure it returns a dtCMatrix!
  ### Otherwise solve in the parent function will be very slow
  
  U.y=U.obj$U[U.obj$latent,]
  
  if(U.obj$cond.yz=='zy') {
    
    V.ord=revMat(U.y[,U.obj$latent,drop=FALSE])
    
  } else if(U.obj$ord.pred!='obspred'){
    
    W=Matrix::tcrossprod(U.y)
    W.rev=revMat(W)
    
    if(U.obj$ic0){
      V.ord=Matrix::t(ichol(W.rev))
    } else {
      V.ord=Matrix::t(Matrix::chol(W.rev))
    }
    
    V.ord = methods::as(V.ord, 'dtCMatrix')
    
    
  } else {  # for obspred ordering
    
    last.obs=max(which(!U.obj$latent))
    latents.before=sum(U.obj$latent[1:last.obs])
    latents.after=sum(U.obj$latent[-(1:last.obs)])
    
    # pred columns are unchanged
    V.pr=revMat(U.y[,(last.obs+1):ncol(U.y),drop=FALSE])
    
    # have to compute cholesky for obs block
    U.oo=U.y[1:latents.before,1:last.obs]
    A=Matrix::tcrossprod(U.oo)
    A.rev=revMat(A)
    if(U.obj$ic0){ V.oor=Matrix::t(ichol(A.rev))
    }     else   V.oor=Matrix::t(Matrix::chol(A.rev))
    
    # combine the blocks into one matrix
    zeromat.sparse=Matrix::sparseMatrix(c(),c(),dims=c(latents.after,latents.before))
    V.or=rbind(zeromat.sparse,V.oor)
    
    V.ord=methods::as(cbind(V.pr,V.or),'dtCMatrix')
    
  }
  
  return(V.ord)
}
