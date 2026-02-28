functions {
  matrix preconditioning_crossprod_low_tri(matrix L, array[] int rid, array[] int cid){
  int n = rows(L);
  matrix[n,n] out = rep_matrix(0,n,n);
  int m = size(rid);
  for(i in 1:m){
    out[rid[i],cid[i]] = dot_product(L[1:n,rid[i]],L[1:n,cid[i]]);
     out[cid[i],rid[i]] = out[rid[i],cid[i]];
    }
  return out;
}

matrix hv_inverse_left_tri_low_to_upp(matrix L, array[] int rid, array[] int nid){
  int n = rows(L);
  matrix[n,n] X = rep_matrix(0,n,n);
  for(k in 1:n){
    X[k,k] = inv(L[k,k]);
    for(i in (nid[k]+2):nid[k+1]){
      X[k,rid[i]] = -L[rid[i],k:(rid[i]-1)] * to_vector(X[k,k:(rid[i]-1)]) * inv(L[rid[i],rid[i]]);
      }
    }
  return X;
}

  matrix my_multiply_lower_tri_self_transpose(matrix x){
    return multiply_lower_tri_self_transpose(x);
  }
  
  matrix my_csr_to_dense_matrix(int m, int n, vector w, array[] int v, array[] int u){
    return csr_to_dense_matrix(m, n, w, v, u);
  }

  matrix matern32_cross(array[] vector x, array[] vector y, real sigma, real lscale){
  return gp_matern32_cov(x, y, sigma, lscale);
  }
  
  matrix matern32_matrix(array[] vector coords, real sigma, real lscale){
    int N = dims(coords)[1];
    matrix[N,N] C = gp_matern32_cov(coords, sigma, lscale);
    return C;
  }
  
  // Forward substitution algorithm to solve Lx = b
  // Stan's equivalent is mdivide_left_tri_low
  // x = forward_sub_solve(L, b) is the solution to L x = b
  //     L must be a lower-triangular matrix
  //     b must be a vector of the same leading dimension as L
  vector forward_sub_solve(matrix L, vector b) {
    int n = rows(L);
    vector[n] x = rep_vector(0,n);
    for(i in 1:n){
      real tmp = 0;
      int j = 1;
      while(j<i){
        tmp = tmp + L[i,j] * x[j];
        j = j + 1;
        }
      x[i] = (b[i] - tmp)*inv(L[i,i]);
      }
    return x;
    }
    
  // Sparse forward substitution algorithm to solve Lx = b
  // Stan's equivalent is mdivide_left_tri_low
  // x = forward_sub(L, b) is the solution to L x = b
  //     L must be a lower-triangular sparse matrix
  //     b must be a vector of the same leading dimension as L
  vector nnforward_sub_solve(matrix L, vector b, array[,] int neiID, array[] int startingID, array[] int endingID) {
    
    int n = rows(L);
    vector[n] x = rep_vector(0,n);
    
    for(i in 1:n){
      
      real tmp = 0;
      if(i == 1){
        x[i] = (b[i] - tmp)*inv(L[i,i]);
      } else {
        int dim = endingID[i] - startingID[i] + 1;
        array[dim] int index = linspaced_int_array(dim, startingID[i], endingID[i]);
        array[dim] int thisID = neiID[i,index];
        for(j in 1:(dim-1)){
          tmp = tmp + L[i,thisID[j]] * x[thisID[j]];
          }
        x[i] = (b[i] - tmp)*inv(L[i,i]); // i and thisID[dim] are the same
      }
    }
  return x;
  }
  
  matrix inverse_left_tri_low(matrix L){
    int n = rows(L);
    matrix[n,n] X = rep_matrix(0,n,n);
    for(k in 1:n){
      X[k,k] = inv(L[k,k]);
      for(i in (k+1):n){
        X[i,k] = -L[i,k:(i-1)] * X[k:(i-1),k] * inv(L[i,i]);
        }
      }
    return X;
}

  matrix hv_inverse_left_tri_low(matrix L, array[] int rid, array[] int nid){
    int n = rows(L);
    matrix[n,n] X = rep_matrix(0,n,n);
    for(k in 1:n){
      X[k,k] = inv(L[k,k]);
      //for(i in (k+1):n){
    for(i in (nid[k]+2):nid[k+1]){
    //X[i,k] = -L[i,k:(i-1)] * X[k:(i-1),k] * inv(L[i,i]);
    X[rid[i],k] = -L[rid[i],k:(rid[i]-1)] * X[k:(rid[i]-1),k] * inv(L[rid[i],rid[i]]);
  }
  }
  return X;
}

//
  vector marginal_matern32_vec(vector d, real sigma, real lscale, real tau, array[] int tauID){
   vector[size(d)] out;
   vector[size(d)] ds = sqrt(3)*d*inv(lscale);
   out = square(sigma) * ((1 + ds) .* exp(-ds));
   out[tauID] = out[tauID] + square(tau);
   return out;
  }
  
  vector matern32_vec(vector d, real sigma, real lscale){
   vector[size(d)] out;
   vector[size(d)] ds = sqrt(3)*d*inv(lscale);
   out = square(sigma) * ((1 + ds) .* exp(-ds));
   return out;
  }
  
  vector mdivide_left_tri_upper(matrix U, vector b) { 
    int n = rows(U);
    vector[n] x = b; 
    real cs;
    array[n] int ids = sort_desc(linspaced_int_array(n, 1, n));
    for (i in ids){
      x[i] = x[i]/U[i,i];
      cs = 0;
      for (j in (i+1):n){
        cs = cs + x[j]*U[i,j];
        }
      x[i] = x[i] - cs/U[i,i];
    }
    return x;
  }


  real sparse_dot_prod(int l1, int u1, int l2, int u2, array[] int row_inds, vector cells){
  
  int l1star = l1;
  int l2star = l2;
  real result = 0.0;
  while(l1star<=u1 && l2star<=u2){
    //print("l1 ", l1star, " l2 ", l2star); # uncheck for consistentcy checking with C++ function
    if(row_inds[l1star]==row_inds[l2star]) {
      result = result + cells[l1star]*cells[l2star];
      l1star = l1star + 1; 
      l2star = l2star + 1;
    }
    else if(row_inds[l1star]<row_inds[l2star])
      l1star = l1star + 1;
    else
      l2star = l2star + 1;
  }
  return result;
 }
  
  //
  vector ic0(array[] int ptrs, array[] int inds, vector x){
   
   vector[size(x)] vals = x;
   int N = size(ptrs)-1;
   real dp;
   
   for(i in 1:N){
     for(j in ptrs[i]:(ptrs[i+1]-1)){
       
       int u1 = ptrs[i];
       int u2 = ptrs[inds[j]];
       ////print("i ", i, " j ", j, " u1 ",u1," u2 ",u2);
       ////print("i ", i, " j ", j, " ptrs[i+1]-2 ", ptrs[i+1]-2," u2 ",ptrs[inds[j]+1]-2);
       dp = sparse_dot_prod(u1, ptrs[i+1]-2, u2, ptrs[inds[j]+1]-2, inds, vals);
       //print("i ", i, " j ", j, " dp ",dp);
       ////print("i ", i, " j ", j, " vals[ptrs[inds[j]+1]-1] ", vals[ptrs[inds[j]+1]-1]);
       if(inds[j]<i){
         vals[j] = (vals[j]-dp)/vals[ptrs[inds[j]+1]-1];
       }
       if(inds[j]==i){
        vals[j] = sqrt(vals[j]-dp);
       }
      }
    }
      return vals;
    }


  matrix sparse_icholesky(array[] int ptrs, array[] int inds, vector x){
    int N = size(ptrs)-1;
    matrix[N,N] out;
    vector[size(x)] w;
    w = ic0(ptrs, inds, x);
    out = csr_to_dense_matrix(N, N, w, inds, ptrs);
    return out;
    }
  
   matrix icholesky_decompose(matrix A){
    
    int n = rows(A);
    matrix[rows(A),cols(A)] B = A;
    
    for(k in 1:n) {
      
      B[k,k] = sqrt(B[k,k]);
      
		for(i in (k+1):n) {
		  if(B[i,k] != 0){
		    B[i,k] = B[i,k]/B[k,k];            
		  }
		}
		
		for(j in (k+1):n){
		  for(i in j:n) {
		    if(B[i,j] != 0){
		      B[i,j] = B[i,j] - B[i,k]*B[j,k];
		    }
		  }
		}
		
		for(i in 1:n){
		  for(j in i+1:n){
		    B[i,j] = 0;
		    }
		  }
		}
		return B;
		}
		
		matrix hvcholesky_matern32_prec(array[] vector coords, real sigma, real lscale, real tau, array[,] int neiID, array[] int startingID, array[] int endingID, int N){
    matrix[N,N] L = rep_matrix(0, N, N);
    for(l in 1:N){
      if(l ==1){
        L[l,l] = inv_sqrt(square(sigma) + square(tau));
        } else {
          int dim = endingID[l] - startingID[l] + 1;
          array[dim] int index = linspaced_int_array(dim, startingID[l], endingID[l]);
          array[dim] int thisIDs = neiID[l,index];
          matrix[dim,dim] C = add_diag(gp_matern32_cov(coords[thisIDs,1:2], sigma, lscale), rep_vector(square(tau), dim));
          matrix[dim-1,dim-1] C1 = C[1:(dim-1),1:(dim-1)];
          real C2 = C[dim,dim];
          vector[(dim-1)] C12 = to_vector(C[1:(dim-1),dim]);
          L[l,l] = inv_sqrt(C2 - C12' * inverse(C1) * C12);
          L[l,thisIDs[1:(dim-1)]] = to_row_vector(-inverse(C1) * C12 * L[l,l]);
          }
        }
      return L;
      } 
      
      
  // Cholesky factor L where C^{-1} = L'L
  matrix nncholesky_matern32_prec(real sigma, real lscale, real tau, matrix site2neiDist, matrix neiDistMat, array[,] int neiID, int N, int K){
  
    int h;
    matrix[N,N] L = rep_matrix(0, N, N);
    
    real variance_ratio_plus_1 = square(tau*inv(sigma)) + 1;
    
    L[1,1] = inv(sigma)*inv_sqrt(variance_ratio_plus_1);
    
    for (i in 2:N) {
      int dim = (i < (K + 1))? (i - 1) : K;
      matrix[dim, dim] neiCorMat;
      matrix[dim, dim] neiCorChol;
      vector[dim] site2neiCor;
      vector[dim] v;
      
      if(dim == 1){
        neiCorMat[1, 1] = variance_ratio_plus_1;
      }else {
        h = 0;
        for (j in 1:(dim - 1)){
          for (k in (j + 1):dim){
            h = h + 1;
            neiCorMat[j, k] = (1 + sqrt(3) * neiDistMat[(i - 1), h] * inv(lscale)) * exp(-sqrt(3) * neiDistMat[(i - 1), h] * inv(lscale));
            neiCorMat[k, j] = neiCorMat[j, k];
          }
        }
        for(j in 1:dim){
          neiCorMat[j, j] = variance_ratio_plus_1;
        }
      }
    neiCorChol = cholesky_decompose(neiCorMat);
    site2neiCor = to_vector((1 + sqrt(3) * site2neiDist[(i - 1)][1:dim] * inv(lscale)) .* exp(-sqrt(3) * site2neiDist[(i - 1)][1:dim] * inv(lscale)));
    v = mdivide_left_tri_low(neiCorChol, site2neiCor);
    L[i,i] = inv(sigma)*inv_sqrt(variance_ratio_plus_1 - dot_self(v));
    L[i,neiID[(i-1), 1:dim]] = -mdivide_right_tri_low(v', neiCorChol)*L[i,i];
  }
  return L;
}  

// fitted values using nearest neighbor approximimation of data likelihood with matern 3/2 
    vector vecchia_matern32_fitted_rng(vector y, vector mu, real sigmasq, real tausq,
                             real lscale, matrix site2neiDist, matrix neiDistMat, 
                             array[,] int neiID, int N, int K) {
                               
          vector[N] yfitted;
          vector[N] cond_mu; // conditional mean
          vector[N] V; // conditional variances = sigmasq*V
          vector[N] resid = y - mu;
          //vector[N] U = resid;
          real variance_ratio_plus_1 = tausq*inv(sigmasq) + 1; // variance ratio plus 1
          V[1] = variance_ratio_plus_1;
          int dim;
          int h;
          yfitted[1] = normal_rng(mu[1], sqrt(sigmasq*V[1]));
          
          for (i in 2:N) {
            dim = (i < (K + 1))? (i - 1) : K;
            matrix[dim, dim] neiCorMat;
              matrix[dim, dim] neiCorChol;
              vector[dim] site2neiCor;
              vector[dim] v;
              row_vector[dim] v2;

              if(dim == 1){neiCorMat[1, 1] = variance_ratio_plus_1;}
              else{
                  h = 0;
                  for (j in 1:(dim - 1)){
                      for (k in (j + 1):dim){
                          h = h + 1;
                          neiCorMat[j, k] = (1 + sqrt(3) * neiDistMat[(i - 1), h] * inv(lscale)) * exp(-sqrt(3) * neiDistMat[(i - 1), h] * inv(lscale));
                          neiCorMat[k, j] = neiCorMat[j, k];
                      }
                  }
                  for(j in 1:dim){
                      neiCorMat[j, j] = variance_ratio_plus_1;
                  }
              }

              neiCorChol = cholesky_decompose(neiCorMat);
              site2neiCor = to_vector((1 + sqrt(3) * site2neiDist[(i - 1), 1: dim] * inv(lscale)) .* exp(-sqrt(3) * site2neiDist[(i - 1), 1: dim] * inv(lscale)));
             v = mdivide_left_tri_low(neiCorChol, site2neiCor);
             V[i] = variance_ratio_plus_1 - dot_self(v); // conditional variances
             v2 = mdivide_right_tri_low(v', neiCorChol);
             cond_mu[i] = mu[i] + v2 * resid[neiID[(i - 1), 1:dim]];
             yfitted[i] = normal_rng(cond_mu[i], sqrt(sigmasq*V[i]));
          }
          return yfitted;
      }

  array[] vector predict_vecchia_rng(vector y, matrix obsX, matrix predX, array[] vector obsCoords, array[] vector pred2obsDist, array[,] int pred2obsNeiID, array[] vector beta, vector sigma, vector lscale, vector tau, int nsize, int psize, int postsize){
    array[postsize] vector[psize] out;
    int nprint = postsize %/% 10;
    int m = dims(pred2obsDist)[2];
    int p = dims(obsX)[2];
    for(l in 1:postsize) {
      if(l%nprint == 0) print("Starts for prediction location : ", l);
      for(i in 1:psize) {
        vector[m] res = y[pred2obsNeiID[i,1:m]] - obsX[pred2obsNeiID[i,1:m],1:p]*beta[l];
        vector[m] c0 = matern32_vec(pred2obsDist[i], 1, lscale[l]);
        matrix[m,m] Ch = cholesky_decompose(add_diag(gp_matern32_cov(obsCoords[pred2obsNeiID[i,1:m]], 1, lscale[l]), rep_vector(square(tau[l])*inv_square(sigma[l]), m)));
        real mu =  predX[i,1:p]*beta[l] + c0'*mdivide_left_tri_upper(Ch',mdivide_left_tri_low(Ch, res));
        real v = square(sigma[l])*(1 + square(tau[l])*inv_square(sigma[l]) - dot_self(mdivide_left_tri_low(Ch, c0)));
        out[l][i] = normal_rng(mu, sqrt(v));
      }
    }
    return out;
  }
  
  //Vecchia (response) GP with Matern32
  real vgp_matern32_lpdf(vector y, vector Xbeta, real sigmasq, real tausq, real lscale, matrix site2neiDist, matrix neiDistMat, array[,] int neiID, int N, int K) {
    
    vector[N] V;
    vector[N] resid = y - Xbeta;
    vector[N] U = resid;
    real variance_ratio_plus_1 = tausq * inv(sigmasq) + 1; // variance ratio plus 1
    int h;
    for (i in 2:N) {
      int dim = (i < (K + 1))? (i - 1) : K;
      matrix[dim, dim] neiCorMat;
      matrix[dim, dim] neiCorChol;
      vector[dim] site2neiCor;
      vector[dim] v;
      row_vector[dim] v2;
      
      if(dim == 1){
        neiCorMat[1, 1] = variance_ratio_plus_1;
        } else {
          h = 0;
          for (j in 1:(dim - 1)){
            for (k in (j + 1):dim){
              h = h + 1;
              neiCorMat[j, k] = (1 + sqrt(3) * neiDistMat[(i - 1), h] * inv(lscale)) * exp(-sqrt(3) * neiDistMat[(i - 1), h] * inv(lscale));
              neiCorMat[k, j] = neiCorMat[j, k];
              }
            }
            for(j in 1:dim){
              neiCorMat[j, j] = variance_ratio_plus_1;
            }
        }

        neiCorChol = cholesky_decompose(neiCorMat);
        site2neiCor = to_vector((1 + sqrt(3) * site2neiDist[(i - 1), 1: dim] * inv(lscale)) .* exp(-sqrt(3) * site2neiDist[(i - 1), 1: dim] * inv(lscale)));
        v = mdivide_left_tri_low(neiCorChol, site2neiCor);
        V[i] = variance_ratio_plus_1 - dot_self(v);
        v2 = mdivide_right_tri_low(v', neiCorChol);
        U[i] = U[i] - v2 * resid[neiID[(i - 1), 1:dim]];
        }
        V[1] = variance_ratio_plus_1;
        return - 0.5 * ( 1 / sigmasq * dot_product(U, (U ./ V)) + sum(log(V)) + N * log(sigmasq));
      }
  
  real multi_normal_cholesky_prec_lpdf(vector y, vector Xbeta, matrix L) {
    real qterm = dot_self(L*(y - Xbeta));
    return sum(log(diagonal(L))) - 0.5 * qterm;
    }
    
  real nngp_lpdf(vector y, vector Xbeta, vector w, array[] int v, array[] int u, real logdet, int N) {
    real qterm = dot_self(csr_matrix_times_vector(N, N, w, v, u, (y - Xbeta)));
    return logdet - 0.5 * qterm;
  }
  
  real hvgp_lpdf(vector y, vector mu, matrix L, array[,] int neiID, array[] int startingID, array[] int endingID){
    real qterm = dot_self(nnforward_sub_solve(L, y-mu, neiID, startingID, endingID));
    real logdet =  sum(log(diagonal(L)));
    return -logdet - 0.5 * qterm; // constant term is ignored
  }
}

data {
  int<lower=0> N;
  int<lower=0> K;
  int<lower=0> nnZero;
  int<lower=0> P;
  vector[N] y;
  matrix[N, P] X;
  array[N] vector[2] coords;
  matrix[N,N] S;
  array[N,K+1] int neiID;
  array[N] int startingID;
  array[N] int endingID;
  array[N] int tauID;
  vector[nnZero] distVec;
  vector[P] mu_theta;
  matrix[P, P] V_theta;
  real<lower=0> sigma_rate;
  real<lower=0> tau_rate;
  real a;
  real b;
}

transformed data {
  cholesky_factor_cov[P] chol_V_theta;
  chol_V_theta = cholesky_decompose(V_theta);
  vector[rows(csr_extract_w(S))] Sw;
  array[rows(Sw)] int Sv;
  array[size(csr_extract_u(S))] int Su;
  // sparsity structure for K
  Sw = csr_extract_w(S);
  Sv = csr_extract_v(S);
  Su = csr_extract_u(S);
}

parameters{
  vector[P] theta_std;
  real<lower = 0> sigma;
  real<lower = 0> ell;
  real<lower = 0> tau;
}

transformed parameters{
  vector[P] theta = mu_theta + chol_V_theta * theta_std;
}

model {
  theta_std ~ std_normal();
  sigma ~ exponential(sigma_rate);
  ell ~ inv_gamma(a,b);
  tau ~ exponential(tau_rate);
  
  // Method 1: exact GP is slow
  //matrix[N,N] C = add_diag(gp_matern32_cov(coords, sigma, ell), rep_vector(square(tau), N));
  //y ~ multi_normal_cholesky(X*theta, cholesky_decompose(C));

  // Method 2: 7.001761 Minutes
  //matrix[N,N] L = hvcholesky_matern32_prec(coords, sigma, ell, tau, neiID, startingID, endingID, N); // chol of precision matrix
  //y ~ multi_normal_cholesky_prec(X*theta, L);

  //Method 3: 8.779619 Minutes
  //vector[size(Sw)] Cw = marginal_matern32_vec(distVec, sigma, ell, tau, tauID); 
  //matrix[N,N] L = sparse_icholesky(Su, Sv, Cw); // ichol of covariance matrix
  //y ~ multi_normal_cholesky(X*theta, L);

  //Method 4: 1.112345 Minutes
  vector[size(Sw)] Cw = marginal_matern32_vec(distVec, sigma, ell, tau, tauID);
  matrix[N,N] L = sparse_icholesky(Su, Sv, Cw); // ichol of covariance matrix
  y ~ hvgp(X*theta, L, neiID, startingID, endingID);
  
}

generated quantities {
 
}



