#calculate residual time series
z = apply(miroc.global.control, 2, function(x) x - t(t(F) %*% t(miroc.dlm.list.univar$miroc.global.control$m)))

F = matrix(c(1,0), nrow= 10)

#make list of means from dlms
{means.obs=list()
means.obs$global = list(obs.dlm.list$global.ncep$m, obs.dlm.list$global.era$m)
means.obs$trop = list(obs.dlm.list$trop.ncep$m, obs.dlm.list$trop.era$m)
means.obs$north = list(obs.dlm.list$north.ncep$m, obs.dlm.list$north.era$m)
means.obs$south = list(obs.dlm.list$south.ncep$m, obs.dlm.list$south.era$m)

{global = list(miroc.dlm.list$miroc.global.control$m, miroc.dlm.list$miroc.global.deca$m, miroc.dlm.list$miroc.global.hist$m, 
              can.dlm.list$can.global.control$m, can.dlm.list$can.global.deca$m, can.dlm.list$can.global.hist$m,
              had.dlm.list$had.global.control$m, had.dlm.list$had.global.deca$m, had.dlm.list$had.global.hist$m, 
              gfdl.dlm.list$gfdl.global.control$m, gfdl.dlm.list$gfdl.global.hist$m)
names(global) = c('miroc.control', 'miroc.deca', 'miroc.hist', 'can.control', 'can.deca', 'can.hist', 'had.control', 'had.deca', 'had.hist', 'gfdl.control', 'gfdl.hist')

trop = list(miroc.dlm.list$miroc.trop.control$m, miroc.dlm.list$miroc.trop.deca$m, miroc.dlm.list$miroc.trop.hist$m, 
              can.dlm.list$can.trop.control$m, can.dlm.list$can.trop.deca$m, can.dlm.list$can.trop.hist$m,
              had.dlm.list$had.trop.control$m, had.dlm.list$had.trop.deca$m, had.dlm.list$had.trop.hist$m, 
              gfdl.dlm.list$gfdl.trop.control$m, gfdl.dlm.list$gfdl.trop.hist$m)
names(trop) = c('miroc.control', 'miroc.deca', 'miroc.hist', 'can.control', 'can.deca', 'can.hist', 'had.control', 'had.deca', 'had.hist', 'gfdl.control', 'gfdl.hist')

north = list(miroc.dlm.list$miroc.north.control$m, miroc.dlm.list$miroc.north.deca$m, miroc.dlm.list$miroc.north.hist$m, 
              can.dlm.list$can.north.control$m, can.dlm.list$can.north.deca$m, can.dlm.list$can.north.hist$m,
              had.dlm.list$had.north.control$m, had.dlm.list$had.north.deca$m, had.dlm.list$had.north.hist$m, 
              gfdl.dlm.list$gfdl.north.control$m, gfdl.dlm.list$gfdl.north.hist$m)
names(north) = c('miroc.control', 'miroc.deca', 'miroc.hist', 'can.control', 'can.deca', 'can.hist', 'had.control', 'had.deca', 'had.hist', 'gfdl.control', 'gfdl.hist')

south = list(miroc.dlm.list$miroc.south.control$m, miroc.dlm.list$miroc.south.deca$m, miroc.dlm.list$miroc.south.hist$m, 
              can.dlm.list$can.south.control$m, can.dlm.list$can.south.deca$m, can.dlm.list$can.south.hist$m,
              had.dlm.list$had.south.control$m, had.dlm.list$had.south.deca$m, had.dlm.list$had.south.hist$m, 
              gfdl.dlm.list$gfdl.south.control$m, gfdl.dlm.list$gfdl.south.hist$m)
names(south) = c('miroc.control', 'miroc.deca', 'miroc.hist', 'can.control', 'can.deca', 'can.hist', 'had.control', 'had.deca', 'had.hist', 'gfdl.control', 'gfdl.hist')

means = list(global, trop, north, south)
names(means) = c('global', 'trop', 'north', 'south')
}
}



##############################save a list of residuals#######################
resid.obs = list()

#global
for (i in 1:length(obs.ts$global)){
    resid.obs$global[[i]] = obs.ts$global[[i]] - t(t(F) %*% t(means.obs$global[[i]]))
  }
names(resid.obs$global) = names(obs.ts$global)

#trop
for (i in 1:length(obs.ts$trop)){
  resid.obs$trop[[i]] = obs.ts$trop[[i]] - t(t(F) %*% t(means.obs$trop[[i]]))
}
names(resid.obs$trop) = names(obs.ts$trop)

#north
for (i in 1:length(obs.ts$north)){
  resid.obs$north[[i]] = obs.ts$north[[i]] - t(t(F) %*% t(means.obs$north[[i]]))
}
names(resid.obs$north) = names(obs.ts$north)

#south
for (i in 1:length(obs.ts$south)){
  resid.obs$south[[i]] = obs.ts$south[[i]] - t(t(F) %*% t(means.obs$south[[i]]))
}
names(resid.obs$south) = names(obs.ts$south)

##################-models-####################
resid = list()
#global
for (i in 1:length(model.ts$global)){
  resid$global[[i]] = apply(model.ts$global[[i]], 2, function(x) x - t(t(F) %*% t(means$global[[i]])))
}
names(resid$global) = names(model.ts$global)

#tropical
for (i in 1:length(model.ts$trop)){
  resid$trop[[i]] = apply(model.ts$trop[[i]], 2, function(x) x - t(t(F) %*% t(means$trop[[i]])))
}
names(resid$trop) = names(model.ts$trop)

#north
for (i in 1:length(model.ts$north)){
  resid$north[[i]] = apply(model.ts$north[[i]], 2, function(x) x - t(t(F) %*% t(means$north[[i]])))
}
names(resid$north) = names(model.ts$north)

#south
for (i in 1:length(model.ts$south)){
  resid$south[[i]] = apply(model.ts$south[[i]], 2, function(x) x - t(t(F) %*% t(means$south[[i]])))
}
names(resid$south) = names(model.ts$south)

  
  
