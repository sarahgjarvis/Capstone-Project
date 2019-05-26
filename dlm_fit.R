dlm.fit = function(y, d1, d2, m0, C0, l0, S0, F1, F2, G1, G2, just.lik)
{
  ### Gets the Time Series Length / Replicate number
  T = nrow(y)
  r = ncol(y)
  ### Gets the State Parameter dimension and submodels
  n = length(m0)
  n1 = length(F1)
  n2 = length(F2)
  ### Constructs F and G
  G2 = as.matrix(G2)
  F = matrix(rep(c(F1,F2),r),n,r)
  G = as.matrix(bdiag(G1,G2))
  ### Variable Saving
  ### Posterior Distribution
  m = matrix(0,T,n)
  C = array(0,c(T,n,n))
  ### Predictive State Distribution
  a = matrix(0,T,n)
  R = array(0,dim = c(T,n,n))
  P = array(0,dim = c(T,n,n))
  W = array(0,dim = c(T,n,n))
  ### One-Step Ahead Forecast
  f = matrix(0,T,r)
  Q = array(0,c(T,r,r))
  inv.Q = array(0,c(T,r,r))
  ### Regression Variables
  e = matrix(0,T,r)
  A = array(0,c(T,n,r))
  ### Sample Variance
  S = vector("numeric",T)
  l = vector("numeric",T)
  
  # prior dim check 
  m0 = matrix(m0,n,1)
  C0 = matrix(C0,n,n)
  
  ### First Update
  ### One-step state forecast
  a[1,]  = G %*% m0
  P[1,,] = G %*% C0 %*% t(G)
  W[1,1:n1,1:n1] = (1-d1) * G1 %*% C0[1:n1,1:n1] %*% t(G1) / d1
  for(s in 1:ceiling(n2/2)){
    index = (n1+1+(s-1)*2):min(n1+2+(s-1)*2,n)
    W[1,index,index] = (1-d2) * G2[index-n1,index-n1] %*% C0[index,index] %*% t(G2[index-n1,index-n1]) / d2
  }
  R[1,,] = P[1,,] + W[1,,]
  ### One-step ahead forecast
  f[1,] = t(F) %*% a[1,]
  Q[1,,] = as.matrix(diag(1,r) + t(F) %*% R[1,,] %*% F,r,r)
  inv.Q[1,,] = solve(Q[1,,])
  ### Auxilary Variables
  e[1,]  = as.matrix(y[1,] - f[1,],r,1)
  A[1,,] = R[1,,] %*% F %*% inv.Q[1,,]
  ### Variance update
  l[1] = l0 + 1
  S[1] = l0 * S0 / l[1] + (t(e[1,]) %*% inv.Q[1,,] %*% e[1,] / l[1])
  ### Posterior Distribution
  m[1,]  = a[1,] + as.matrix(A[1,,],n,r) %*% e[1,]
  C[1,,] = R[1,,] - as.matrix(A[1,,],n,r) %*% Q[1,,] %*% t(A[1,,]) 
  C[1,,] = (C[1,,] + t(C[1,,]))/2
  
  for(i in 2:T)
  {
    ### One-step state forecast
    a[i,]  = G %*% m[i-1,]
    P[i,,] = G %*% C[i-1,,] %*% t(G)
    W[i,1:n1,1:n1] = (1-d1) * G1 %*% C[i-1,1:n1,1:n1] %*% t(G1) / d1
    for(s in 1:ceiling(n2/2)){
      index = (n1+1+(s-1)*2):min(n1+2+(s-1)*2,n)
      W[i,index,index] = (1-d2) * G2[index-n1,index-n1] %*% C[i-1,index,index] %*% t(G2[index-n1,index-n1]) / d2
    }
    R[i,,] = P[i,,] + W[i,,]
    ### One-step ahead forecast
    f[i,] = t(F) %*% a[i,]
    Q[i,,] = matrix(diag(1,r) + t(F)%*% R[i,,]%*% F,r,r)
    inv.Q[i,,] = solve(Q[i,,])
    ### Auxilary Variables
    e[i,]  = as.matrix(y[i,] - f[i,],r,1)
    A[i,,] = as.matrix(R[i,,] %*% F %*% inv.Q[i,,],n,r)
    ### Variance update
    l[i] = l[i-1] + 1
    S[i] = l[i-1] * S[i-1] / l[i] + (t(e[i,]) %*% inv.Q[i,,] %*% e[i,] / l[i])
    ### Posterior Distribution
    m[i,]  = a[i,] + as.matrix(A[i,,],n,r) %*% e[i,]
    C[i,,] = R[i,,] - as.matrix(A[i,,],n,r) %*% Q[i,,] %*% t(as.matrix(A[i,,],n,r)) 
    C[i,,] = (C[i,,] + t(C[i,,]))/2
  }
  
  ### Adjust By Variance
  R[1,,] = S0 * R[1,,]
  Q[1,,]   = S0 * Q[1,,]
  C[1,,]   = S[1] * C[1,,]
  for(i in 2:T)
  {
    R[i,,] = S[i-1] * R[i,,]
    Q[i,,]   = S[i-1] * Q[i,,]
    C[i,,]   = S[i] * C[i,,]
  }
  
  # always calulate likelihood
  if(r == 1){ det.Q = log(abs(Q[1,,])) }else{ det.Q = log(det(as.matrix(Q[1,,],r,r))) }
  llik = lgamma((l0+r)/2)-lgamma(l0/2)-r*log(pi*l0)/2-det.Q/2-(l0+r)*log(1+t(e[1,])%*%inv.Q[1,,]%*%e[1,]/l0)/2
  for(t in 2:T){
    if(r == 1){ det.Q = log(abs(Q[t,,])) }else{ det.Q = log(det(as.matrix(Q[t,,],r,r))) }
    llik = llik + lgamma((l[t-1]+r)/2)-lgamma(l[t-1]/2)-r*log(pi*l[t-1])/2-det.Q/2-(l[t-1]+r)*log(1+t(e[t,])%*%inv.Q[t,,]%*%e[t,]/l[t-1])/2
    if(is.nan(llik)){stop("llik NaN at t=",t)}
  }
  if(just.lik){
    return(list(llik = llik))
  }else{
    ## SMOOTHING
    ### Initializes recursive relations
    sa = matrix(0,T,n)
    sR = array(0, dim = c(T,n,n))
    ### Runs the recursive equations
    sa[T,]  = m[T,]
    sR[T,,] = C[T,,]
    for(k in 1:(T-1))
    {
      ### Computes the Auxilary recursion Variable B
      B = C[T-k,,] %*% t(G) %*% solve(R[T-k+1,,])
      sa[T-k,] = m[T-k,] + B %*% (sa[T-k+1,] - a[T-k+1,])
      sR[T-k,,] = C[T-k,,] + B %*% (sR[T-k+1,,] - R[T-k+1,,]) %*% t(B)
    }
    ### Adjusts the variance update
    for(k in 1:T)
    {
      sR[T-k,,] = S[T] * sR[T-k,,] / S[T-k]
    }
    
    return(list(llik=llik,fm = m, fC = C, m = sa, C = sR, F = F, s = S, n = l, W = W))
  }
}

trend.dlm.mod = function(degree)
{
  ### Creates the Matrix G and vector F
  G = diag(rep(1,degree))
  F = vector("numeric", degree)
  ### Sets values
  if(degree > 1)
  {
    for(i in 2:degree)
    {
      G[i-1,i] = 1
    }
  }
  F[1] = 1
  ### Returns FF and GG
  return(list(F = F, G = G))
}

### Seasonal Form DLM
sea.dlm.mod = function(p, h)
{
  library(Matrix)
  ### Obtains the number of harmonic comonents
  nh = length(h)
  
  ### Obtains the frequencies for each harmonic
  w = h * 2 * pi / p
  ### Checks if the nyquist frequency was provided
  if( max(w) == pi)
  {
    ### Constructs the G matrices for each harmonic
    ### Creates the G sub-matrices
    G = array(0,c(nh-1,2,2))
    for(i in 1:(nh-1))
    {
      G[i,1,1] =  cos(w[i])
      G[i,1,2] =  sin(w[i])
      G[i,2,1] = -sin(w[i])
      G[i,2,2] =  cos(w[i])
    }
    ### Constructes the Seasonal Matrix
    GG = matrix(0,0,0)
    for(i in 1:(nh-1))
    {
      if(i ==1 ){ GG = G[1,,]}
      else{GG = bdiag(GG,G[i,,])}
    }
    GG = bdiag(GG,-1)
    ### Creates the Observation vector
    FF = vector("numeric", 2*nh -1)
    FF[1:(2*nh-1) %% 2 == 1] = 1
  } else {
    ### Constructs the G matrices for each harmonic
    ### Creates the G sub-matrices
    G = array(0,c(nh,2,2))
    for(i in 1:nh)
    {
      G[i,1,1] =  cos(w[i])
      G[i,1,2] =  sin(w[i])
      G[i,2,1] = -sin(w[i])
      G[i,2,2] =  cos(w[i])
    }
    ### Constructes the Seasonal Matrix
    GG = matrix(0,0,0)
    for(i in 1:nh)
    {
      if(i ==1 ){ GG = G[1,,]}
      else{GG = bdiag(GG,G[i,,])}
    }
    ### Creates the Observation vector
    FF = vector("numeric", 2*nh)
    FF[1:(2*nh) %% 2 == 1] = 1
  }
  ### Returns FF and GG
  return(list(F = FF, G = GG))
}
