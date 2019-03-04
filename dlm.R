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

#load data
y <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/tropical_MIROC5_control_30years.dat"))

#initial values
degree = 2
p = 12
#h = 1:4 
mod.trend = trend.dlm.mod(degree)
mod.sea = sea.dlm.mod(p,h)
# set priors
#m0 = c(285,0,rep(0,2*length(h)))
m0 = c(285, rep(0,2*length(h)))
#C0 = diag(c(5,2e-6,5,rep(1,2*length(h)-1)))
C0 = diag(c(5,2e-6,5,rep(1,2*length(h)-2)))
n0 = 1
S0 = 0.01
# set discount factors
#df.trend = 0.94
df.seas = 1

#------------------------------------------------------------------------------

#discount factor selection
h = 1:6
#maximize the loglikelihood over a grid of discount factors
df.trend <- seq(0.8, 1, 0.01)
ll <- matrix(NA, length(df.trend), ncol(y))
for(ts in 1:ncol(y)){
  for (i in 1:length(df.trend)){
    ll[i,ts] <- unlist(dlm.fit(as.matrix(y[,ts]), df.trend[i], df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, TRUE))
    
    print(i)
  }
}
row.names(ll) <- df.trend
which(ll == max(ll), arr.ind = TRUE)

#------------------------------------------------------------------------------

# fit DLM

#average the replicate rows/UNIVARIATE
y.mean <- as.matrix(rowMeans(y))

#global dlm
{global.hist.mean = dlm.fit(y.mean, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
  global.control.mean =  dlm.fit(y.mean, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
  global.era.mean = dlm.fit(y.mean, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
  global.ncep.mean = dlm.fit(y.mean, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
  global.deca.mean = dlm.fit(y.mean, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)}

#tropical dlm
{trop.hist.mean = dlm.fit(y.mean, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
  trop.control.mean = dlm.fit(y.mean, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
  trop.deca.mean = dlm.fit(y.mean, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
  trop.era.mean = dlm.fit(y.mean, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
  trop.ncep.mean = dlm.fit(y.mean, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)}

#northern dlm
{north.hist.mean = dlm.fit(y.mean, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
  north.control.mean = dlm.fit(y.mean, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
  north.era.mean = dlm.fit(y.mean, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
  north.ncep.mean = dlm.fit(y.mean, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
  north.deca.mean = dlm.fit(y.mean, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)}

#southern dlm
{south.hist.mean = dlm.fit(y.mean, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
  south.control.mean = dlm.fit(y.mean, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
  south.era.mean = dlm.fit(y.mean, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
  south.ncep.mean = dlm.fit(y.mean, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
  south.deca.mean = dlm.fit(y.mean, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)}



#------------------------------------------------------------------------------
#MULTIVARIATE
#global dlm
{global.hist = dlm.fit(y, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
global.control =  dlm.fit(y, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, TRUE)
global.era = dlm.fit(y, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
global.ncep = dlm.fit(y, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
global.deca = dlm.fit(y, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)}

#tropical dlm
{trop.hist = dlm.fit(y, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
trop.control = dlm.fit(y, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
trop.deca = dlm.fit(y, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
trop.era = dlm.fit(y, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
trop.ncep = dlm.fit(y, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)}

#northern dlm
{north.hist = dlm.fit(y, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
north.control = dlm.fit(y, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
north.era = dlm.fit(y, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
north.ncep = dlm.fit(y, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
north.deca = dlm.fit(y, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)}

#southern dlm
{south.hist = dlm.fit(y, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
south.control = dlm.fit(y, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
south.era = dlm.fit(y, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
south.ncep = dlm.fit(y, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
south.deca = dlm.fit(y, df.trend, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)}


#BASELINE COMPONENT
library(PearsonDS)

#which dlm?
dlm <- global.era.mean

trend.index = 1:degree
t.mean <- NULL
t.var <- NULL
dist <- NULL
quan <- matrix(NA, 360, 3)

#which dlm?
dlm <- global.deca.mean
for (t in 1:360){
  t.mean <- as.matrix(t(mod.trend$F)) %*% dlm$m[t,trend.index] 
  t.var <- t(as.matrix(mod.trend$F)) %*% dlm$C[t,trend.index,trend.index] %*% as.matrix(mod.trend$F)
  dist <- rpearsonVII(10000, dlm$n[dim(y)[1]], t.mean, sqrt(t.var))
  quan[t,] <- quantile(dist, c(0.025, 0.50, 0.975))
}

#save quantiles
#means/univariate
{# q.mean.global.hist <- quan
  # q.mean.global.control <- quan
  # q.mean.global.era <- quan
  # q.mean.global.ncep <- quan
  # q.mean.global.deca <- quan
  
  # q.mean.trop.hist <- quan
  # q.mean.trop.control <- quan
  # q.mean.trop.era <- quan
  # q.mean.trop.deca <- quan
  # q.mean.trop.ncep <- quan
  
  # q.mean.north.hist <- quan
  # q.mean.north.control <- quan
  # q.mean.north.era <- quan
  # q.mean.north.ncep <- quan
  # q.mean.north.deca <- quan
  
  # q.mean.south.hist <- quan
  # q.mean.south.control <- quan
  # q.mean.south.era <- quan
  # q.mean.south.ncep <- quan
  # q.mean.south.deca <- quan
}

x <- 1:360
par(mfrow=c(2,2))
#Global plot
{plot.ts(global.hist.mean$m[,1], col='blue', lwd=1.5, ylab = 'ave temp (K)',
         xlab = 'time (months)', ylim = c(287, 289), main='Global')
  polygon(c(rev(x), x), 
          c(rev(q.mean.global.hist[ ,3]), q.mean.global.hist[ ,1]), 
          col =rgb(0,0,225, max= 225, alpha=50), border = NA)
  lines(ts(global.control.mean$m[,1]), col='green', lwd=1.5)
  polygon(c(rev(x), x), 
          c(rev(q.mean.global.control[ ,3]), q.mean.global.control[ ,1]), 
          col =rgb(0,225,0, max= 225, alpha=50), border = NA)
  lines(ts(global.era.mean$m[,1]), col='purple', lwd=1.5)
  polygon(c(rev(x), x), 
          c(rev(q.mean.global.era[ ,3]), q.mean.global.era[ ,1]), 
          col =rgb(100,0,225, max= 225, alpha=50), border = NA)
  lines(ts(global.ncep.mean$m[,1]), col='black', lwd=1.5)
  polygon(c(rev(x), x), 
          c(rev(q.mean.global.ncep[ ,3]), q.mean.global.ncep[ ,1]), 
          col =rgb(0,0,0, max= 225, alpha=50), border = NA)
  lines(ts(global.deca.mean$m[,1]), col='red', lwd=1.5)
  polygon(c(rev(x), x), 
          c(rev(q.mean.global.deca[ ,3]), q.mean.global.deca[ ,1]), 
          col =rgb(225,0,0, max= 225, alpha=50), border = NA)}
#Tropical plot
{plot.ts(trop.hist.mean$m[,1], col='blue', lwd=1.5, ylab = 'ave temp (K)',
         xlab = 'time (months)', ylim = c(297, 300), main='Tropical')
  polygon(c(rev(x), x), 
          c(rev(q.mean.trop.hist[ ,3]), q.mean.trop.hist[ ,1]), 
          col =rgb(0,0,225, max= 225, alpha=50), border = NA)
  lines(ts(trop.control.mean$m[,1]), col='green', lwd=1.5)
  polygon(c(rev(x), x), 
          c(rev(q.mean.trop.control[ ,3]), q.mean.trop.control[ ,1]), 
          col =rgb(0,225,0, max= 225, alpha=50), border = NA)
  lines(ts(trop.era.mean$m[,1]), col='purple', lwd=1.5)
  polygon(c(rev(x), x), 
          c(rev(q.mean.trop.era[ ,3]), q.mean.trop.era[ ,1]), 
          col =rgb(100,0,225, max= 225, alpha=50), border = NA)
  lines(ts(trop.ncep.mean$m[,1]), col='black', lwd=1.5)
  polygon(c(rev(x), x), 
          c(rev(q.mean.trop.ncep[ ,3]), q.mean.trop.ncep[ ,1]), 
          col =rgb(0,0,0, max= 225, alpha=50), border = NA)
  lines(ts(trop.deca.mean$m[,1]), col='red', lwd=1.5)
  polygon(c(rev(x), x), 
          c(rev(q.mean.trop.deca[ ,3]), q.mean.trop.deca[ ,1]), 
          col =rgb(225,0,0, max= 225, alpha=50), border = NA)}
#Northern plot
{plot.ts(north.hist.mean$m[,1], col='blue', lwd=1.5, ylab = 'ave temp (K)',
         xlab = 'time (months)', ylim = c(287.5, 289), main='Northern')
  polygon(c(rev(x), x), 
          c(rev(q.mean.north.hist[ ,3]), q.mean.north.hist[ ,1]), 
          col =rgb(0,0,225, max= 225, alpha=50), border = NA)
  lines(ts(north.control.mean$m[,1]), col='green', lwd=1.5)
  polygon(c(rev(x), x), 
          c(rev(q.mean.north.control[ ,3]), q.mean.north.control[ ,1]), 
          col =rgb(0,225,0, max= 225, alpha=50), border = NA)
  lines(ts(north.era.mean$m[,1]), col='purple', lwd=1.5)
  polygon(c(rev(x), x), 
          c(rev(q.mean.north.era[ ,3]), q.mean.north.era[ ,1]), 
          col =rgb(100,0,225, max= 225, alpha=50), border = NA)
  lines(ts(north.ncep.mean$m[,1]), col='black', lwd=1.5)
  polygon(c(rev(x), x), 
          c(rev(q.mean.north.ncep[ ,3]), q.mean.north.ncep[ ,1]), 
          col =rgb(0,0,0, max= 225, alpha=50), border = NA)
  lines(ts(north.deca.mean$m[,1]), col='red', lwd=1.5)
  polygon(c(rev(x), x), 
          c(rev(q.mean.north.deca[ ,3]), q.mean.north.deca[ ,1]), 
          col =rgb(225,0,0, max= 225, alpha=50), border = NA)}
#Southern plot
{plot.ts(south.hist.mean$m[,1], col='blue', lwd=1.5, ylab = 'ave temp (K)',
         xlab = 'time (months)', ylim = c(286, 288.5), main='Southern')
  polygon(c(rev(x), x), 
          c(rev(q.mean.south.hist[ ,3]), q.mean.south.hist[ ,1]), 
          col =rgb(0,0,225, max= 225, alpha=50), border = NA)
  lines(ts(south.control.mean$m[,1]), col='green', lwd=1.5)
  polygon(c(rev(x), x), 
          c(rev(q.mean.south.control[ ,3]), q.mean.south.control[ ,1]), 
          col =rgb(0,225,0, max= 225, alpha=50), border = NA)
  lines(ts(south.era.mean$m[,1]), col='purple', lwd=1.5)
  polygon(c(rev(x), x), 
          c(rev(q.mean.south.era[ ,3]), q.mean.south.era[ ,1]), 
          col =rgb(100,0,225, max= 225, alpha=50), border = NA)
  lines(ts(south.ncep.mean$m[,1]), col='black', lwd=1.5)
  polygon(c(rev(x), x), 
          c(rev(q.mean.south.ncep[ ,3]), q.mean.south.ncep[ ,1]), 
          col =rgb(0,0,0, max= 225, alpha=50), border = NA)
  lines(ts(south.deca.mean$m[,1]), col='red', lwd=1.5)
  polygon(c(rev(x), x), 
          c(rev(q.mean.south.deca[ ,3]), q.mean.south.deca[ ,1]), 
          col =rgb(225,0,0, max= 225, alpha=50), border = NA)}



#all columns/multivariate
{# q.global.hist <- quan
# q.global.control <- quan
# q.global.era <- quan
# q.global.ncep <- quan
# q.global.deca <- quan
# q.trop.hist <- quan
# q.trop.control <- quan
# q.trop.era <- quan
# q.trop.deca <- quan
# q.trop.ncep <- quan
# q.north.hist <- quan
# q.north.control <- quan
# q.north.era <- quan
# q.north.ncep <- quan
# q.north.deca <- quan

# q.south.hist <- quan
# q.south.control <- quan
# q.south.era <- quan
# q.south.ncep <- quan
# q.south.deca <- quan
  }

x <- 1:360
par(mfrow=c(2,2))
#Global plot
{plot.ts(global.hist$m[,1], col='blue', lwd=1.5, ylab = 'ave temp (K)',
        xlab = 'time (months)', ylim = c(287, 289), main='Global')
polygon(c(rev(x), x), 
        c(rev(q.global.hist[ ,3]), q.global.hist[ ,1]), 
        col =rgb(0,0,225, max= 225, alpha=50), border = NA)
lines(ts(global.control$m[,1]), col='green', lwd=1.5)
polygon(c(rev(x), x), 
        c(rev(q.global.control[ ,3]), q.global.control[ ,1]), 
        col =rgb(0,225,0, max= 225, alpha=50), border = NA)
lines(ts(global.era$m[,1]), col='purple', lwd=1.5)
polygon(c(rev(x), x), 
       c(rev(q.global.era[ ,3]), q.global.era[ ,1]), 
       col =rgb(100,0,225, max= 225, alpha=50), border = NA)
lines(ts(global.ncep$m[,1]), col='black', lwd=1.5)
polygon(c(rev(x), x), 
        c(rev(q.global.ncep[ ,3]), q.global.ncep[ ,1]), 
        col =rgb(0,0,0, max= 225, alpha=50), border = NA)
lines(ts(global.deca$m[,1]), col='red', lwd=1.5)
polygon(c(rev(x), x), 
        c(rev(q.global.deca[ ,3]), q.global.deca[ ,1]), 
        col =rgb(225,0,0, max= 225, alpha=50), border = NA)}
#Tropical plot
{plot.ts(trop.hist$m[,1], col='blue', lwd=1.5, ylab = 'ave temp (K)',
         xlab = 'time (months)', ylim = c(297, 300), main='Tropical')
  polygon(c(rev(x), x), 
          c(rev(q.trop.hist[ ,3]), q.trop.hist[ ,1]), 
          col =rgb(0,0,225, max= 225, alpha=50), border = NA)
  lines(ts(trop.control$m[,1]), col='green', lwd=1.5)
  polygon(c(rev(x), x), 
          c(rev(q.trop.control[ ,3]), q.trop.control[ ,1]), 
          col =rgb(0,225,0, max= 225, alpha=50), border = NA)
  lines(ts(trop.era$m[,1]), col='purple', lwd=1.5)
  polygon(c(rev(x), x), 
          c(rev(q.trop.era[ ,3]), q.trop.era[ ,1]), 
          col =rgb(100,0,225, max= 225, alpha=50), border = NA)
  lines(ts(trop.ncep$m[,1]), col='black', lwd=1.5)
  polygon(c(rev(x), x), 
          c(rev(q.trop.ncep[ ,3]), q.trop.ncep[ ,1]), 
          col =rgb(0,0,0, max= 225, alpha=50), border = NA)
  lines(ts(trop.deca$m[,1]), col='red', lwd=1.5)
  polygon(c(rev(x), x), 
          c(rev(q.trop.deca[ ,3]), q.trop.deca[ ,1]), 
          col =rgb(225,0,0, max= 225, alpha=50), border = NA)}
#Northern plot
{plot.ts(north.hist$m[,1], col='blue', lwd=1.5, ylab = 'ave temp (K)',
         xlab = 'time (months)', ylim = c(287.5, 289), main='Northern')
  polygon(c(rev(x), x), 
          c(rev(q.north.hist[ ,3]), q.north.hist[ ,1]), 
          col =rgb(0,0,225, max= 225, alpha=50), border = NA)
  lines(ts(north.control$m[,1]), col='green', lwd=1.5)
  polygon(c(rev(x), x), 
          c(rev(q.north.control[ ,3]), q.north.control[ ,1]), 
          col =rgb(0,225,0, max= 225, alpha=50), border = NA)
  lines(ts(north.era$m[,1]), col='purple', lwd=1.5)
  polygon(c(rev(x), x), 
          c(rev(q.north.era[ ,3]), q.north.era[ ,1]), 
          col =rgb(100,0,225, max= 225, alpha=50), border = NA)
  lines(ts(north.ncep$m[,1]), col='black', lwd=1.5)
  polygon(c(rev(x), x), 
          c(rev(q.north.ncep[ ,3]), q.north.ncep[ ,1]), 
          col =rgb(0,0,0, max= 225, alpha=50), border = NA)
  lines(ts(north.deca$m[,1]), col='red', lwd=1.5)
  polygon(c(rev(x), x), 
          c(rev(q.north.deca[ ,3]), q.north.deca[ ,1]), 
          col =rgb(225,0,0, max= 225, alpha=50), border = NA)}
#Southern plot
{plot.ts(south.hist$m[,1], col='blue', lwd=1.5, ylab = 'ave temp (K)',
         xlab = 'time (months)', ylim = c(286, 288.5), main='Southern')
  polygon(c(rev(x), x), 
          c(rev(q.south.hist[ ,3]), q.south.hist[ ,1]), 
          col =rgb(0,0,225, max= 225, alpha=50), border = NA)
  lines(ts(south.control$m[,1]), col='green', lwd=1.5)
  polygon(c(rev(x), x), 
          c(rev(q.south.control[ ,3]), q.south.control[ ,1]), 
          col =rgb(0,225,0, max= 225, alpha=50), border = NA)
  lines(ts(south.era$m[,1]), col='purple', lwd=1.5)
  polygon(c(rev(x), x), 
          c(rev(q.south.era[ ,3]), q.south.era[ ,1]), 
          col =rgb(100,0,225, max= 225, alpha=50), border = NA)
  lines(ts(south.ncep$m[,1]), col='black', lwd=1.5)
  polygon(c(rev(x), x), 
          c(rev(q.south.ncep[ ,3]), q.south.ncep[ ,1]), 
          col =rgb(0,0,0, max= 225, alpha=50), border = NA)
  lines(ts(south.deca$m[,1]), col='red', lwd=1.5)
  polygon(c(rev(x), x), 
          c(rev(q.south.deca[ ,3]), q.south.deca[ ,1]), 
          col =rgb(225,0,0, max= 225, alpha=50), border = NA)}


#AMPLITUDE

#last time point not over time
seas.index = (degree+1):dim(dlm$m)[2]
seas.quan <- matrix(NA, 360, 3)

t <- 360
#change which dlm

#which dlm?
dlm <- trop.control
theta <- rmvt(10000, dlm$C[t,,], dlm$n[dim(y)[1]], dlm$m[t,])

#global
{
  #harmonic 1 (annual)
  global.ncep.amp.h1 <- sqrt(theta[,3]^2 + theta[,4]^2)
  global.era.amp.h1 <- sqrt(theta[,3]^2 + theta[,4]^2)
  global.deca.amp.h1 <- sqrt(theta[,3]^2 + theta[,4]^2)
  global.hist.amp.h1 <- sqrt(theta[,3]^2 + theta[,4]^2)
  global.control.amp.h1 <- sqrt(theta[,3]^2 + theta[,4]^2)
  global.har1 <- list(global.ncep.amp.h1, global.era.amp.h1, global.deca.amp.h1, global.hist.amp.h1, global.control.amp.h1)
  names(global.har1) <- c('ncep', 'era', 'dec', 'hist', 'cont')
  
  #harmonic 2 (semi-annual)
  global.ncep.amp.h2 <- sqrt(theta[,5]^2 + theta[,6]^2)
  global.era.amp.h2 <- sqrt(theta[,5]^2 + theta[,6]^2)
  global.deca.amp.h2 <- sqrt(theta[,5]^2 + theta[,6]^2)
  global.hist.amp.h2 <- sqrt(theta[,5]^2 + theta[,6]^2)
  global.control.amp.h2 <- sqrt(theta[,5]^2 + theta[,6]^2)
  global.har2 <- list(global.ncep.amp.h2, global.era.amp.h2, global.deca.amp.h2, global.hist.amp.h2, global.control.amp.h2)
  names(global.har2) <- c('ncep', 'era', 'dec', 'hist', 'cont')
}

#tropical
{
  #harmonic 1 (annual)
  trop.ncep.amp.h1 <- sqrt(theta[,3]^2 + theta[,4]^2)
  trop.era.amp.h1 <- sqrt(theta[,3]^2 + theta[,4]^2)
  trop.deca.amp.h1 <- sqrt(theta[,3]^2 + theta[,4]^2)
  trop.hist.amp.h1 <- sqrt(theta[,3]^2 + theta[,4]^2)
  trop.control.amp.h1 <- sqrt(theta[,3]^2 + theta[,4]^2)
  trop.har1 <- list(trop.ncep.amp.h1, trop.era.amp.h1, trop.deca.amp.h1, trop.hist.amp.h1, trop.control.amp.h1)
  names(trop.har1) <- c('ncep', 'era', 'dec', 'hist', 'cont')
  
  #harmonic 2 (semi-annual)
  trop.ncep.amp.h2 <- sqrt(theta[,5]^2 + theta[,6]^2)
  trop.era.amp.h2 <- sqrt(theta[,5]^2 + theta[,6]^2)
  trop.deca.amp.h2 <- sqrt(theta[,5]^2 + theta[,6]^2)
  trop.hist.amp.h2 <- sqrt(theta[,5]^2 + theta[,6]^2)
  trop.control.amp.h2 <- sqrt(theta[,5]^2 + theta[,6]^2)
  trop.har2 <- list(trop.ncep.amp.h2, trop.era.amp.h2, trop.deca.amp.h2, trop.hist.amp.h2, trop.control.amp.h2)
  names(trop.har2) <- c('ncep', 'era', 'dec', 'hist', 'cont')
}

#south
{
  #harmonic 1 (annual)
  south.ncep.amp.h1 <- sqrt(theta[,3]^2 + theta[,4]^2)
  south.era.amp.h1 <- sqrt(theta[,3]^2 + theta[,4]^2)
  south.deca.amp.h1 <- sqrt(theta[,3]^2 + theta[,4]^2)
  south.hist.amp.h1 <- sqrt(theta[,3]^2 + theta[,4]^2)
  south.control.amp.h1 <- sqrt(theta[,3]^2 + theta[,4]^2)
  south.har1 <- list(south.ncep.amp.h1, south.era.amp.h1, south.deca.amp.h1, south.hist.amp.h1, south.control.amp.h1)
  names(south.har1) <- c('ncep', 'era', 'dec', 'hist', 'cont')
  
  #harmonic 2 (semi-annual)
  south.ncep.amp.h2 <- sqrt(theta[,5]^2 + theta[,6]^2)
  south.era.amp.h2 <- sqrt(theta[,5]^2 + theta[,6]^2)
  south.deca.amp.h2 <- sqrt(theta[,5]^2 + theta[,6]^2)
  south.hist.amp.h2 <- sqrt(theta[,5]^2 + theta[,6]^2)
  south.control.amp.h2 <- sqrt(theta[,5]^2 + theta[,6]^2)
  south.har2 <- list(south.ncep.amp.h2, south.era.amp.h2, south.deca.amp.h2, south.hist.amp.h2, south.control.amp.h2)
  names(south.har2) <- c('ncep', 'era', 'dec', 'hist', 'cont')
}

#north
{
  #harmonic 1 (annual)
  north.ncep.amp.h1 <- sqrt(theta[,3]^2 + theta[,4]^2)
  north.era.amp.h1 <- sqrt(theta[,3]^2 + theta[,4]^2)
  north.deca.amp.h1 <- sqrt(theta[,3]^2 + theta[,4]^2)
  north.hist.amp.h1 <- sqrt(theta[,3]^2 + theta[,4]^2)
  north.control.amp.h1 <- sqrt(theta[,3]^2 + theta[,4]^2)
  north.har1 <- list(north.ncep.amp.h1, north.era.amp.h1, north.deca.amp.h1, north.hist.amp.h1, north.control.amp.h1)
  names(north.har1) <- c('ncep', 'era', 'dec', 'hist', 'cont')
  
  #harmonic 2 (semi-annual)
  north.ncep.amp.h2 <- sqrt(theta[,5]^2 + theta[,6]^2)
  north.era.amp.h2 <- sqrt(theta[,5]^2 + theta[,6]^2)
  north.deca.amp.h2 <- sqrt(theta[,5]^2 + theta[,6]^2)
  north.hist.amp.h2 <- sqrt(theta[,5]^2 + theta[,6]^2)
  north.control.amp.h2 <- sqrt(theta[,5]^2 + theta[,6]^2)
  north.har2 <- list(north.ncep.amp.h2, north.era.amp.h2, north.deca.amp.h2, north.hist.amp.h2, north.control.amp.h2)
  names(north.har2) <- c('ncep', 'era', 'dec', 'hist', 'cont')
}

#Amplitude plots
par(mfrow=c(2,4))
{boxplot(sapply(global.har1, function(x) quantile(x, c(0, 0.025, 0.5, 0.975, 1))), col=c('black', 'purple', 'red', 'blue', 'green'), main = 'Global, annual')

boxplot(sapply(global.har2, function(x) quantile(x, c(0, 0.025, 0.5, 0.975, 1))), col=c('black', 'purple', 'red', 'blue', 'green'), main = 'Global, semi-annual')

boxplot(sapply(trop.har1, function(x) quantile(x, c(0, 0.025, 0.5, 0.975, 1))), col=c('black', 'purple', 'red', 'blue', 'green'), main = 'Tropical, annual')

boxplot(sapply(trop.har2, function(x) quantile(x, c(0, 0.025, 0.5, 0.975, 1))), col=c('black', 'purple', 'red', 'blue', 'green'), main = 'Tropical, semi-annual')


boxplot(sapply(north.har1, function(x) quantile(x, c(0, 0.025, 0.5, 0.975, 1))), col=c('black', 'purple', 'red', 'blue', 'green'), main = 'Northern, annual')

boxplot(sapply(north.har2, function(x) quantile(x, c(0, 0.025, 0.5, 0.975, 1))), col=c('black', 'purple', 'red', 'blue', 'green'), main = 'Northern, semi-annual')

boxplot(sapply(south.har1, function(x) quantile(x, c(0, 0.025, 0.5, 0.975, 1))), col=c('black', 'purple', 'red', 'blue', 'green'), main = 'Southern, annual')

boxplot(sapply(south.har2, function(x) quantile(x, c(0, 0.025, 0.5, 0.975, 1))), col=c('black', 'purple', 'red', 'blue', 'green'), main = 'Southern, semi-annual')
}


