# univariate DLMS
# run DLM for one replicate and save the mean

degree = 2
trend.index = 1:degree
p = 12
h = 1:4 
mod.trend = trend.dlm.mod(degree)
mod.sea = sea.dlm.mod(p,h)
# set priors
m0 = c(285,0,rep(0,2*length(h)))
C0 = diag(c(5,2e-6,5,rep(1,2*length(h)-1)))
n0 = 1
S0 = 0.01
# set discount factor
df.seas = 1

# observations
# Global
{univar.dlm.list.obs = list()
  for (i in 1:length(obs.ts$global)){
    dlm = dlm.fit(as.matrix(obs.ts$global[[i]][,1]), 0.96, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
    univar.dlm.list.obs$global[[i]] = dlm
  }
  names(univar.dlm.list.obs$global) = names(obs.ts$global)
}

# Tropical
{
  for (i in 1:length(obs.ts$trop)){
    dlm = dlm.fit(as.matrix(obs.ts$trop[[i]][,1]), 0.94, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
    univar.dlm.list.obs$trop[[i]] = dlm
  }
  names(univar.dlm.list.obs$trop) = names(obs.ts$trop)
}

# Northern
{
  for (i in 1:length(obs.ts$north)){
    dlm = dlm.fit(as.matrix(obs.ts$north[[i]][,1]), 0.99, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
    univar.dlm.list.obs$north[[i]] = dlm
  }
  names(univar.dlm.list.obs$north) = names(obs.ts$north)
}

# Southern
{
  for (i in 1:length(obs.ts$south)){
    dlm = dlm.fit(as.matrix(obs.ts$south[[i]][,1]), 0.99, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
    univar.dlm.list.obs$south[[i]] = dlm
  }
  names(univar.dlm.list.obs$south) = names(obs.ts$south)
}


# ----------------------------models-------------------------------


# Global
# fit dlm and save to list
{univar.dlm.list = list()
  for (i in 1:length(model.ts$global)){
    dlm = dlm.fit(as.matrix(model.ts$global[[i]][,1]), 0.94, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
    univar.dlm.list$global[[i]] = dlm
  }
  names(univar.dlm.list$global) = names(model.ts$global)
}

# Tropical
{
  for (i in 1:length(model.ts$trop)){
    dlm = dlm.fit(as.matrix(model.ts$trop[[i]][,1]), 0.91, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
    univar.dlm.list$trop[[i]] = dlm
  }
  names(univar.dlm.list$trop) = names(model.ts$trop)
}

# Northern
{
  for (i in 1:length(model.ts$north)){
    dlm = dlm.fit(as.matrix(model.ts$north[[i]][,1]), 0.99, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
    univar.dlm.list$north[[i]] = dlm
  }
  names(univar.dlm.list$north) = names(model.ts$north)
}

# Southern
{
  for (i in 1:length(model.ts$south)){
    dlm = dlm.fit(as.matrix(model.ts$south[[i]][,1]), 0.99, df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
    univar.dlm.list$south[[i]] = dlm
  }
  names(univar.dlm.list$south) = names(model.ts$south)
}




#check
plot.ts(model.ts$trop$miroc.control[,1])
lines(univar.dlm.list$trop$miroc.control$m[,1])


save(univar.dlm.list, file = 'univariate_dlm_list.RData')
