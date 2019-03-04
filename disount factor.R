
#initial values
degree = 2
p = 12
#h = 1:4 
mod.trend = trend.dlm.mod(degree)
mod.sea = sea.dlm.mod(p,h)
# set priors
m0 = c(285, rep(0,2*length(h)))
C0 = diag(c(5,2e-6,5,rep(1,2*length(h)-2)))
n0 = 1
S0 = 0.01
df.seas = 1

#discount factor selection
df <- function(data){
  h = 1:6
  #maximize the loglikelihood over a grid of discount factors
  df.trend <- seq(0.8, 1, 0.01)
  ll <- matrix(NA, length(df.trend), ncol(y))
  for(ts in 1:ncol(y)){
    for (i in 1:length(df.trend)){
      ll[i,ts] <- unlist(dlm.fit(as.matrix(data[,ts]), df.trend[i], df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, TRUE))

    }
    print(ts)
  }
  row.names(ll) <- df.trend
  return(which(ll == max(ll), arr.ind = TRUE))
} 

df(can.global.deca)
