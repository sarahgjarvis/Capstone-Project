library(mvtnorm)

#Posterior Amplitudes

degree = 2
seas.index = (degree+1):dim(dlm$m)[2]
seas.quan <- matrix(NA, 360, 3)

t <- 360

#save theta
#obs theta
obs.theta.list <- list()
for (i in 1:length(obs.dlm.list)){
  theta <- rmvt(5000, obs.dlm.list[[i]]$C[t,,], obs.dlm.list[[i]]$n[360], obs.dlm.list[[i]]$m[t,])
  obs.theta.list[[i]] <- theta
}
names(obs.theta.list) <- names(obs.dlm.list)

#MIROC5 theta
miroc.theta.list <- list()
for (i in 1:length(miroc.dlm.list)){
  theta <- rmvt(5000, miroc.dlm.list[[i]]$C[t,,], miroc.dlm.list[[i]]$n[360], miroc.dlm.list[[i]]$m[t,])
  miroc.theta.list[[i]] <- theta
}
names(miroc.theta.list) <- names(miroc.dlm.list)

#CanCM4 theta
can.theta.list <- list()
for (i in 1:length(can.dlm.list)){
  theta <- rmvt(5000, can.dlm.list[[i]]$C[t,,], can.dlm.list[[i]]$n[360], can.dlm.list[[i]]$m[t,])
  can.theta.list[[i]] <- theta
}
names(can.theta.list) <- names(can.dlm.list)

#HadCM3 theta
had.theta.list <- list()
for (i in 1:length(had.dlm.list)){
  theta <- rmvt(5000, had.dlm.list[[i]]$C[t,,], had.dlm.list[[i]]$n[360], had.dlm.list[[i]]$m[t,])
  had.theta.list[[i]] <- theta
}
names(had.theta.list) <- names(had.dlm.list)

#GFDL theta
gfdl.theta.list <- list()
for (i in 1:length(gfdl.dlm.list)){
  theta <- rmvt(5000, gfdl.dlm.list[[i]]$C[t,,], gfdl.dlm.list[[i]]$n[360], gfdl.dlm.list[[i]]$m[t,])
  gfdl.theta.list[[i]] <- theta
}
names(gfdl.theta.list) <- names(gfdl.dlm.list)


#save amplitudes for both harmonics
#obs
{obs.amp.list <- list()
for (i in 1:length(obs.theta.list)){
  h1 <- sqrt(obs.theta.list[[i]][,3]^2 + obs.theta.list[[i]][,4]^2)
  h2 <- sqrt(obs.theta.list[[i]][,5]^2 + obs.theta.list[[i]][,6]^2)
  obs.amp.list[[i]] <- cbind(h1, h2)
}
names(obs.amp.list) <- names(obs.dlm.list)
}
#MIROC5
{miroc.amp.list <- list()
  for (i in 1:length(miroc.theta.list)){
    h1 <- sqrt(miroc.theta.list[[i]][,3]^2 + miroc.theta.list[[i]][,4]^2)
    h2 <- sqrt(miroc.theta.list[[i]][,5]^2 + miroc.theta.list[[i]][,6]^2)
    miroc.amp.list[[i]] <- cbind(h1, h2)
  }
  names(miroc.amp.list) <- names(miroc.dlm.list)
}
#CanCM4
{can.amp.list <- list()
  for (i in 1:length(can.theta.list)){
    h1 <- sqrt(can.theta.list[[i]][,3]^2 + can.theta.list[[i]][,4]^2)
    h2 <- sqrt(can.theta.list[[i]][,5]^2 + can.theta.list[[i]][,6]^2)
    can.amp.list[[i]] <- cbind(h1, h2)
  }
  names(can.amp.list) <- names(can.dlm.list)
}
#HadCM3
{had.amp.list <- list()
  for (i in 1:length(had.theta.list)){
    h1 <- sqrt(had.theta.list[[i]][,3]^2 + had.theta.list[[i]][,4]^2)
    h2 <- sqrt(had.theta.list[[i]][,5]^2 + had.theta.list[[i]][,6]^2)
    had.amp.list[[i]] <- cbind(h1, h2)
  }
  names(had.amp.list) <- names(had.dlm.list)
}
#CanCM3
{can.amp.list <- list()
  for (i in 1:length(can.theta.list)){
    h1 <- sqrt(can.theta.list[[i]][,3]^2 + can.theta.list[[i]][,4]^2)
    h2 <- sqrt(can.theta.list[[i]][,5]^2 + can.theta.list[[i]][,6]^2)
    can.amp.list[[i]] <- cbind(h1, h2)
  }
  names(can.amp.list) <- names(can.dlm.list)
}
#GFDL
{gfdl.amp.list <- list()
  for (i in 1:length(gfdl.theta.list)){
    h1 <- sqrt(gfdl.theta.list[[i]][,3]^2 + gfdl.theta.list[[i]][,4]^2)
    h2 <- sqrt(gfdl.theta.list[[i]][,5]^2 + gfdl.theta.list[[i]][,6]^2)
    gfdl.amp.list[[i]] <- cbind(h1, h2)
  }
  names(gfdl.amp.list) <- names(gfdl.dlm.list)
}

#save amplitude quantiles for plotting
#harmonic 1
{obs.amp.h1 <- lapply(obs.amp.list, function(x) quantile(x[,1], c(0, 0.025, 0.5, 0.975, 1)))
miroc.amp.h1 <- lapply(miroc.amp.list, function(x) quantile(x[,1], c(0, 0.025, 0.5, 0.975, 1)))
can.amp.h1 <- lapply(can.amp.list, function(x) quantile(x[,1], c(0, 0.025, 0.5, 0.975, 1)))
had.amp.h1 <- lapply(had.amp.list, function(x) quantile(x[,1], c(0, 0.025, 0.5, 0.975, 1)))
gfdl.amp.h1 <- lapply(gfdl.amp.list, function(x) quantile(x[,1], c(0, 0.025, 0.5, 0.975, 1)))}
#harmonic 2
{obs.amp.h2 <- lapply(obs.amp.list, function(x) quantile(x[,2], c(0, 0.025, 0.5, 0.975, 1)))
miroc.amp.h2 <- lapply(miroc.amp.list, function(x) quantile(x[,2], c(0, 0.025, 0.5, 0.975, 1)))
can.amp.h2 <- lapply(can.amp.list, function(x) quantile(x[,2], c(0, 0.025, 0.5, 0.975, 1)))
had.amp.h2 <- lapply(had.amp.list, function(x) quantile(x[,2], c(0, 0.025, 0.5, 0.975, 1)))
gfdl.amp.h2 <- lapply(gfdl.amp.list, function(x) quantile(x[,2], c(0, 0.025, 0.5, 0.975, 1)))}

#Amplitude boxplot for MIROC5
par(mfrow=c(2, 4),oma = c(0, 0, 2, 0))
{
{
#global annual
boxplot(obs.amp.h1$global.ncep, obs.amp.h1$global.era, miroc.amp.h1$miroc.global.deca, miroc.amp.h1$miroc.global.hist, miroc.amp.h1$miroc.global.control, col = c('black', 'purple', 'red', 'blue', 'green'), names = c('ncep', 'era', 'dec', 'hist', 'cont'), main='Global, annual (k=1)', ylab='ave temp (K)')
  #tropical annual
  boxplot(obs.amp.h1$trop.ncep, obs.amp.h1$trop.era, miroc.amp.h1$miroc.trop.deca, miroc.amp.h1$miroc.trop.hist, miroc.amp.h1$miroc.trop.control, col = c('black', 'purple', 'red', 'blue', 'green'), names = c('ncep', 'era', 'dec', 'hist', 'cont'), main='Tropical, annual (k=1)', ylab='ave temp (K)')
  #northern annual
  boxplot(obs.amp.h1$north.ncep, obs.amp.h1$north.era, miroc.amp.h1$miroc.north.deca, miroc.amp.h1$miroc.north.hist, miroc.amp.h1$miroc.north.control, col = c('black', 'purple', 'red', 'blue', 'green'), names = c('ncep', 'era', 'dec', 'hist', 'cont'), main='Northern, annual (k=1)', ylab='ave temp (K)')
  #southern annual
  boxplot(obs.amp.h1$south.ncep, obs.amp.h1$south.era, miroc.amp.h1$miroc.south.deca, miroc.amp.h1$miroc.south.hist, miroc.amp.h1$miroc.south.control, col = c('black', 'purple', 'red', 'blue', 'green'), names = c('ncep', 'era', 'dec', 'hist', 'cont'), main='Southern, annual (k=1)', ylab='ave temp (K)')
  
#semi-annual
boxplot(obs.amp.h2$global.ncep, obs.amp.h2$global.era, miroc.amp.h2$miroc.global.deca, miroc.amp.h2$miroc.global.hist, miroc.amp.h2$miroc.global.control, col = c('black', 'purple', 'red', 'blue', 'green'), names = c('ncep', 'era', 'dec', 'hist', 'cont'), main='Global, semi-annual (k=2)', ylab='ave temp (K)')
boxplot(obs.amp.h2$trop.ncep, obs.amp.h2$trop.era, miroc.amp.h2$miroc.trop.deca, miroc.amp.h2$miroc.trop.hist, miroc.amp.h2$miroc.trop.control, col = c('black', 'purple', 'red', 'blue', 'green'), names = c('ncep', 'era', 'dec', 'hist', 'cont'), main='Tropical, semi-annual (k=2)', ylab='ave temp (K)')
boxplot(obs.amp.h2$north.ncep, obs.amp.h2$north.era, miroc.amp.h2$miroc.north.deca, miroc.amp.h2$miroc.north.hist, miroc.amp.h2$miroc.north.control, col = c('black', 'purple', 'red', 'blue', 'green'), names = c('ncep', 'era', 'dec', 'hist', 'cont'), main='Northern, semi-annual (k=2)', ylab='ave temp (K)')
#semi-annual
boxplot(obs.amp.h2$south.ncep, obs.amp.h2$south.era, miroc.amp.h2$miroc.south.deca, miroc.amp.h2$miroc.south.hist, miroc.amp.h2$miroc.south.control, col = c('black', 'purple', 'red', 'blue', 'green'), names = c('ncep', 'era', 'dec', 'hist', 'cont'), main='Southern, semi-annual (k=2)', ylab='ave temp (K)')
}
}
mtext(c('NCEP', 'ERA', 'M-Decadal', 'M-Historical', 'M-Control'), adj = c(0.43, 0.46, 0.5, 0.56, 0.615), outer = TRUE, cex = 0.8, col = c('black', 'purple', 'red', 'blue', 'green'))


#Amplitude boxplot for CanCM4
#Global
{
  #annual
  boxplot(obs.amp.h1$global.ncep, obs.amp.h1$global.era, can.amp.h1$can.global.deca, can.amp.h1$can.global.hist, can.amp.h1$can.global.control, col = c('black', 'purple', 'red', 'blue', 'green'), names = c('ncep', 'era', 'dec', 'hist', 'cont'), main='Global, annual (k=1)', ylab='ave temp (K)')
  boxplot(obs.amp.h1$trop.ncep, obs.amp.h1$trop.era, can.amp.h1$can.trop.deca, can.amp.h1$can.trop.hist, can.amp.h1$can.trop.control, col = c('black', 'purple', 'red', 'blue', 'green'), names = c('ncep', 'era', 'dec', 'hist', 'cont'), main='Tropical, annual (k=1)', ylab='ave temp (K)')
  boxplot(obs.amp.h1$north.ncep, obs.amp.h1$north.era, can.amp.h1$can.north.deca, can.amp.h1$can.north.hist, can.amp.h1$can.north.control, col = c('black', 'purple', 'red', 'blue', 'green'), names = c('ncep', 'era', 'dec', 'hist', 'cont'), main='Northern, annual (k=1)', ylab='ave temp (K)')
  boxplot(obs.amp.h1$south.ncep, obs.amp.h1$south.era, can.amp.h1$can.south.deca, can.amp.h1$can.south.hist, can.amp.h1$can.south.control, col = c('black', 'purple', 'red', 'blue', 'green'), names = c('ncep', 'era', 'dec', 'hist', 'cont'), main='Southern, annual (k=1)', ylab='ave temp (K)')
  
  
  #semi-annual
  boxplot(obs.amp.h2$global.ncep, obs.amp.h2$global.era, can.amp.h2$can.global.deca, can.amp.h2$can.global.hist, can.amp.h2$can.global.control, col = c('black', 'purple', 'red', 'blue', 'green'), names = c('ncep', 'era', 'dec', 'hist', 'cont'), main='Global, semi-annual (k=2)', ylab='ave temp (K)')
  boxplot(obs.amp.h2$trop.ncep, obs.amp.h2$trop.era, can.amp.h2$can.trop.deca, can.amp.h2$can.trop.hist, can.amp.h2$can.trop.control, col = c('black', 'purple', 'red', 'blue', 'green'), names = c('ncep', 'era', 'dec', 'hist', 'cont'), main='Tropical, semi-annual (k=2)', ylab='ave temp (K)')
  boxplot(obs.amp.h2$north.ncep, obs.amp.h2$north.era, can.amp.h2$can.north.deca, can.amp.h2$can.north.hist, can.amp.h2$can.north.control, col = c('black', 'purple', 'red', 'blue', 'green'), names = c('ncep', 'era', 'dec', 'hist', 'cont'), main='Northern, semi-annual (k=2)', ylab='ave temp (K)')
  boxplot(obs.amp.h2$south.ncep, obs.amp.h2$south.era, can.amp.h2$can.south.deca, can.amp.h2$can.south.hist, can.amp.h2$can.south.control, col = c('black', 'purple', 'red', 'blue', 'green'), names = c('ncep', 'era', 'dec', 'hist', 'cont'), main='Southern, semi-annual (k=2)', ylab='ave temp (K)')
}
mtext(c('NCEP', 'ERA', 'M-Decadal', 'M-Historical', 'M-Control'), adj = c(0.43, 0.46, 0.5, 0.56, 0.615), outer = TRUE, cex = 0.8, col = c('black', 'purple', 'red', 'blue', 'green'))

#Amplitude boxplot for HadCM3
#Global
{
  #annual
  boxplot(obs.amp.h1$global.ncep, obs.amp.h1$global.era, had.amp.h1$had.global.deca, had.amp.h1$had.global.hist, had.amp.h1$had.global.control, col = c('black', 'purple', 'red', 'blue', 'green'), names = c('ncep', 'era', 'dec', 'hist', 'cont'), main='Global, annual (k=1)', ylab='ave temp (K)')
  boxplot(obs.amp.h1$trop.ncep, obs.amp.h1$trop.era, had.amp.h1$had.trop.deca, had.amp.h1$had.trop.hist, had.amp.h1$had.trop.control, col = c('black', 'purple', 'red', 'blue', 'green'), names = c('ncep', 'era', 'dec', 'hist', 'cont'), main='Tropical, annual (k=1)', ylab='ave temp (K)')
  boxplot(obs.amp.h1$north.ncep, obs.amp.h1$north.era, had.amp.h1$had.north.deca, had.amp.h1$had.north.hist, had.amp.h1$had.north.control, col = c('black', 'purple', 'red', 'blue', 'green'), names = c('ncep', 'era', 'dec', 'hist', 'cont'), main='Northern, annual (k=1)', ylab='ave temp (K)')
  #annual
  boxplot(obs.amp.h1$south.ncep, obs.amp.h1$south.era, had.amp.h1$had.south.deca, had.amp.h1$had.south.hist, had.amp.h1$had.south.control, col = c('black', 'purple', 'red', 'blue', 'green'), names = c('ncep', 'era', 'dec', 'hist', 'cont'), main='Southern, annual (k=1)', ylab='ave temp (K)')

  
  boxplot(obs.amp.h2$global.ncep, obs.amp.h2$global.era, had.amp.h2$had.global.deca, had.amp.h2$had.global.hist, had.amp.h2$had.global.control, col = c('black', 'purple', 'red', 'blue', 'green'), names = c('ncep', 'era', 'dec', 'hist', 'cont'), main='Global, semi-annual (k=2)', ylab='ave temp (K)')
  boxplot(obs.amp.h2$trop.ncep, obs.amp.h2$trop.era, had.amp.h2$had.trop.deca, had.amp.h2$had.trop.hist, had.amp.h2$had.trop.control, col = c('black', 'purple', 'red', 'blue', 'green'), names = c('ncep', 'era', 'dec', 'hist', 'cont'), main='Tropical, semi-annual (k=2)', ylab='ave temp (K)')
boxplot(obs.amp.h2$north.ncep, obs.amp.h2$north.era, had.amp.h2$had.north.deca, had.amp.h2$had.north.hist, had.amp.h2$had.north.control, col = c('black', 'purple', 'red', 'blue', 'green'), names = c('ncep', 'era', 'dec', 'hist', 'cont'), main='Northern, semi-annual (k=2)', ylab='ave temp (K)')
boxplot(obs.amp.h2$south.ncep, obs.amp.h2$south.era, had.amp.h2$had.south.deca, had.amp.h2$had.south.hist, had.amp.h2$had.south.control, col = c('black', 'purple', 'red', 'blue', 'green'), names = c('ncep', 'era', 'dec', 'hist', 'cont'), main='Southern, semi-annual (k=2)', ylab='ave temp (K)')
mtext(c('NCEP', 'ERA', 'M-Decadal', 'M-Historical', 'M-Control'), adj = c(0.43, 0.46, 0.5, 0.56, 0.615), outer = TRUE, cex = 0.8, col = c('black', 'purple', 'red', 'blue', 'green'))

}

#Amplitude boxplot for GFDL
{
  #annual
  boxplot(obs.amp.h1$global.ncep, obs.amp.h1$global.era, gfdl.amp.h1$gfdl.global.hist, gfdl.amp.h1$gfdl.global.control, col = c('black', 'purple', 'blue', 'green'), names = c('ncep', 'era', 'hist', 'cont'), main='Global, annual (k=1)', ylab='ave temp (K)')
  boxplot(obs.amp.h1$trop.ncep, obs.amp.h1$trop.era, gfdl.amp.h1$gfdl.trop.hist, gfdl.amp.h1$gfdl.trop.control, col = c('black', 'purple', 'blue', 'green'), names = c('ncep', 'era', 'hist', 'cont'), main='Tropical, annual (k=1)', ylab='ave temp (K)')
  #annual
  boxplot(obs.amp.h1$north.ncep, obs.amp.h1$north.era, gfdl.amp.h1$gfdl.north.hist, gfdl.amp.h1$gfdl.north.control, col = c('black', 'purple', 'blue', 'green'), names = c('ncep', 'era', 'hist', 'cont'), main='Northern, annual (k=1)', ylab='ave temp (K)')
  #annual
  boxplot(obs.amp.h1$south.ncep, obs.amp.h1$south.era, gfdl.amp.h1$gfdl.south.hist, gfdl.amp.h1$gfdl.south.control, col = c('black', 'purple', 'blue', 'green'), names = c('ncep', 'era', 'hist', 'cont'), main='Southern, annual (k=1)', ylab='ave temp (K)')
  
  #semi-annual
  boxplot(obs.amp.h2$global.ncep, obs.amp.h2$global.era, gfdl.amp.h2$gfdl.global.hist, gfdl.amp.h2$gfdl.global.control, col = c('black', 'purple', 'blue', 'green'), names = c('ncep', 'era', 'hist', 'cont'), main='Global, semi-annual (k=2)', ylab='ave temp (K)')
  #semi-annual
  boxplot(obs.amp.h2$trop.ncep, obs.amp.h2$trop.era, gfdl.amp.h2$gfdl.trop.hist, gfdl.amp.h2$gfdl.trop.control, col = c('black', 'purple', 'blue', 'green'), names = c('ncep', 'era', 'hist', 'cont'), main='Tropical, semi-annual (k=2)', ylab='ave temp (K)')
  #semi-annual
  boxplot(obs.amp.h2$north.ncep, obs.amp.h2$north.era, gfdl.amp.h2$gfdl.north.hist, gfdl.amp.h2$gfdl.north.control, col = c('black', 'purple', 'blue', 'green'), names = c('ncep', 'era', 'hist', 'cont'), main='Northern, semi-annual (k=2)', ylab='ave temp (K)')
  #semi-annual
  boxplot(obs.amp.h2$south.ncep, obs.amp.h2$south.era, gfdl.amp.h2$gfdl.south.hist, gfdl.amp.h2$gfdl.south.control, col = c('black', 'purple', 'blue', 'green'), names = c('ncep', 'era', 'hist', 'cont'), main='Southern, semi-annual (k=2)', ylab='ave temp (K)')
  mtext(c('NCEP', 'ERA', 'M-Decadal', 'M-Historical', 'M-Control'), adj = c(0.43, 0.46, 0.5, 0.56, 0.615), outer = TRUE, cex = 0.8, col = c('black', 'purple', 'red', 'blue', 'green'))
}


#annual
#Amplitude plot for comparing models
par(mar=c(4,4,1,1))
par(mfrow=c(4,3), oma = c(0, 0, 2, 0))
#Global
{
boxplot(miroc.amp.h1$miroc.global.deca, can.amp.h1$can.global.deca, had.amp.h1$had.global.deca, names = c('MIROC5', 'CanCM4', 'HadCM3'), col='red')

boxplot(miroc.amp.h1$miroc.global.hist, can.amp.h1$can.global.hist, had.amp.h1$had.global.hist, gfdl.amp.h1$gfdl.global.hist, names = c('MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col='blue', main='Global, annual')

boxplot(miroc.amp.h1$miroc.global.control, can.amp.h1$can.global.control, had.amp.h1$had.global.control, gfdl.amp.h1$gfdl.global.control, names = c('MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col='green')
}
#Tropical
{
  boxplot(miroc.amp.h1$miroc.trop.deca, can.amp.h1$can.trop.deca, had.amp.h1$had.trop.deca, names = c('MIROC5', 'CanCM4', 'HadCM3'), col='red')
  
  boxplot(miroc.amp.h1$miroc.trop.hist, can.amp.h1$can.trop.hist, had.amp.h1$had.trop.hist, gfdl.amp.h1$gfdl.trop.hist, names = c('MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col='blue', main='Tropical, annual')
  
  boxplot(miroc.amp.h1$miroc.trop.control, can.amp.h1$can.trop.control, had.amp.h1$had.trop.control, gfdl.amp.h1$gfdl.trop.control, names = c('MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col='green')
}
#Northern
{
  boxplot(miroc.amp.h1$miroc.north.deca, can.amp.h1$can.north.deca, had.amp.h1$had.north.deca, names = c('MIROC5', 'CanCM4', 'HadCM3'), col='red')
  
  boxplot(miroc.amp.h1$miroc.north.hist, can.amp.h1$can.north.hist, had.amp.h1$had.north.hist, gfdl.amp.h1$gfdl.north.hist, names = c('MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col='blue', main='Northern, annual')
  
  boxplot(miroc.amp.h1$miroc.north.control, can.amp.h1$can.north.control, had.amp.h1$had.north.control, gfdl.amp.h1$gfdl.north.control, names = c('MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col='green')
}
#Southern
{
  boxplot(miroc.amp.h1$miroc.south.deca, can.amp.h1$can.south.deca, had.amp.h1$had.south.deca, names = c('MIROC5', 'CanCM4', 'HadCM3'), col='red')
  
  boxplot(miroc.amp.h1$miroc.south.hist, can.amp.h1$can.south.hist, had.amp.h1$had.south.hist, gfdl.amp.h1$gfdl.south.hist, names = c('MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col='blue', main='Southern, annual')
  
  boxplot(miroc.amp.h1$miroc.south.control, can.amp.h1$can.south.control, had.amp.h1$had.south.control, gfdl.amp.h1$gfdl.south.control, names = c('MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col='green')
}




#Amplitude plot for comparing models and observations
par(mar=c(4,4,1,1))
par(mfrow=c(4,3), oma = c(0, 0, 2, 0))
#Global
{
  boxplot(obs.amp.h1$global.ncep, obs.amp.h1$global.era, miroc.amp.h1$miroc.global.deca, can.amp.h1$can.global.deca, had.amp.h1$had.global.deca, names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3'), col=c('black', 'purple', rep('red', 3)), cex.axis=0.9, ylim = c(1.8, 2.3), ylab = 'ave temp (K)')
  
  boxplot(obs.amp.h1$global.ncep, obs.amp.h1$global.era, miroc.amp.h1$miroc.global.hist, can.amp.h1$can.global.hist, had.amp.h1$had.global.hist, gfdl.amp.h1$gfdl.global.hist, names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col=c('black', 'purple', rep('blue', 4)), main='Global, annual', cex.axis=0.9, ylim = c(1.8, 2.3), ylab = 'ave temp (K)')
  
  boxplot(obs.amp.h1$global.ncep, obs.amp.h1$global.era, miroc.amp.h1$miroc.global.control, can.amp.h1$can.global.control, had.amp.h1$had.global.control, gfdl.amp.h1$gfdl.global.control, names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col=c('black', 'purple', rep('green', 4)), cex.axis=0.9, ylim = c(1.8, 2.3), ylab = 'ave temp (K)')
}
#Tropical
{
  boxplot(obs.amp.h1$trop.ncep, obs.amp.h1$trop.era, miroc.amp.h1$miroc.trop.deca, can.amp.h1$can.trop.deca, had.amp.h1$had.trop.deca, names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3'), col=c('black', 'purple', rep('red', 3)), cex.axis=0.9, ylim = c(0.1, 0.5), ylab = 'ave temp (K)')
  
  boxplot(obs.amp.h1$trop.ncep, obs.amp.h1$trop.era, miroc.amp.h1$miroc.trop.hist, can.amp.h1$can.trop.hist, had.amp.h1$had.trop.hist, gfdl.amp.h1$gfdl.trop.hist, names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col=c('black', 'purple', rep('blue', 4)), main='Tropical, annual', cex.axis=0.9, ylim = c(0.1, 0.5), ylab = 'ave temp (K)')
  
  boxplot(obs.amp.h1$trop.ncep, obs.amp.h1$trop.era, miroc.amp.h1$miroc.trop.control, can.amp.h1$can.trop.control, had.amp.h1$had.trop.control, gfdl.amp.h1$gfdl.trop.control, names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col=c('black', 'purple', rep('green', 4)), cex.axis=0.9, ylim = c(0.1, 0.5), ylab = 'ave temp (K)')
}
#Northern
{
  boxplot(obs.amp.h1$north.ncep, obs.amp.h1$north.era, miroc.amp.h1$miroc.north.deca, can.amp.h1$can.north.deca, had.amp.h1$had.north.deca, names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3'), col=c('black', 'purple', rep('red', 3)), cex.axis=0.9, ylim = c(6.3, 7.5), ylab = 'ave temp (K)')
  
  boxplot(obs.amp.h1$north.ncep, obs.amp.h1$north.era, miroc.amp.h1$miroc.north.hist, can.amp.h1$can.north.hist, had.amp.h1$had.north.hist, gfdl.amp.h1$gfdl.north.hist, names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col=c('black', 'purple', rep('blue', 4)), main='Northern, annual', cex.axis=0.9, ylim = c(6.3, 7.5), ylab = 'ave temp (K)')
  
  boxplot(obs.amp.h1$north.ncep, obs.amp.h1$north.era, miroc.amp.h1$miroc.north.control, can.amp.h1$can.north.control, had.amp.h1$had.north.control, gfdl.amp.h1$gfdl.north.control, names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col=c('black', 'purple', rep('green', 4)), cex.axis=0.9, ylim = c(6.3, 7.5), ylab = 'ave temp (K)')
}
#Southern
{
  boxplot(obs.amp.h1$south.ncep, obs.amp.h1$south.era, miroc.amp.h1$miroc.south.deca, can.amp.h1$can.south.deca, had.amp.h1$had.south.deca, names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3'), col=c('black', 'purple', rep('red', 3)), cex.axis=0.9, ylim = c(2.3, 3.3), ylab = 'ave temp (K)')
  
  boxplot(obs.amp.h1$south.ncep, obs.amp.h1$south.era, miroc.amp.h1$miroc.south.hist, can.amp.h1$can.south.hist, had.amp.h1$had.south.hist, gfdl.amp.h1$gfdl.south.hist, names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col=c('black', 'purple', rep('blue', 4)), main='Southern, annual', cex.axis=0.9, ylim = c(2.3, 3.3), ylab = 'ave temp (K)')
  
  boxplot(obs.amp.h1$south.ncep, obs.amp.h1$south.era, miroc.amp.h1$miroc.south.control, can.amp.h1$can.south.control, had.amp.h1$had.south.control, gfdl.amp.h1$gfdl.south.control, names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col=c('black', 'purple', rep('green', 4)), cex.axis=0.9, ylim = c(2.3, 3.3), ylab = 'ave temp (K)')
}
mtext(c('NCEP', 'ERA', 'M-Decadal', 'M-Historical', 'M-Control'), adj = c(0.43, 0.46, 0.5, 0.56, 0.615), outer = TRUE, cex = 0.8, col = c('black', 'purple', 'red', 'blue', 'green'))


#semi-annual
#Amplitude plot for comparing models
par(mfrow=c(4,3))
#Global
{
  boxplot(miroc.amp.h2$miroc.global.deca, can.amp.h2$can.global.deca, had.amp.h2$had.global.deca, gfdl.amp.h2$gfdl.global.deca, names = c('MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col='red')
  
  boxplot(miroc.amp.h2$miroc.global.hist, can.amp.h2$can.global.hist, had.amp.h2$had.global.hist, gfdl.amp.h2$gfdl.global.hist, names = c('MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col='blue', main='Global, semi-annual')
  
  boxplot(miroc.amp.h2$miroc.global.control, can.amp.h2$can.global.control, had.amp.h2$had.global.control, gfdl.amp.h2$gfdl.global.control, names = c('MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col='green')
}
#Tropical
{
  boxplot(miroc.amp.h2$miroc.trop.deca, can.amp.h2$can.trop.deca, had.amp.h2$had.trop.deca, gfdl.amp.h2$gfdl.trop.deca, names = c('MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col='red')
  
  boxplot(miroc.amp.h2$miroc.trop.hist, can.amp.h2$can.trop.hist, had.amp.h2$had.trop.hist, gfdl.amp.h2$gfdl.trop.hist, names = c('MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col='blue', main='Tropical, semi-annual')
  
  boxplot(miroc.amp.h2$miroc.trop.control, can.amp.h2$can.trop.control, had.amp.h2$had.trop.control, gfdl.amp.h2$gfdl.trop.control, names = c('MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col='green')
}
#Northern
{
  boxplot(miroc.amp.h2$miroc.north.deca, can.amp.h2$can.north.deca, had.amp.h2$had.north.deca, gfdl.amp.h2$gfdl.north.deca, names = c('MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col='red')
  
  boxplot(miroc.amp.h2$miroc.north.hist, can.amp.h2$can.north.hist, had.amp.h2$had.north.hist, gfdl.amp.h2$gfdl.north.hist, names = c('MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col='blue', main='Northern, semi-annual')
  
  boxplot(miroc.amp.h2$miroc.north.control, can.amp.h2$can.north.control, had.amp.h2$had.north.control, gfdl.amp.h2$gfdl.north.control, names = c('MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col='green')
}
#Southern
{
  boxplot(miroc.amp.h2$miroc.south.deca, can.amp.h2$can.south.deca, had.amp.h2$had.south.deca, gfdl.amp.h2$gfdl.south.deca, names = c('MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col='red')
  
  boxplot(miroc.amp.h2$miroc.south.hist, can.amp.h2$can.south.hist, had.amp.h2$had.south.hist, gfdl.amp.h2$gfdl.south.hist, names = c('MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col='blue', main='Southern, semi-annual')
  
  boxplot(miroc.amp.h2$miroc.south.control, can.amp.h2$can.south.control, had.amp.h2$had.south.control, gfdl.amp.h2$gfdl.south.control, names = c('MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col='green')
}

#Amplitude plot for comparing models and observations
par(mar=c(4,4,1,1))
par(mfrow=c(4,3), oma = c(0, 0, 2, 0))
#Global
{
  boxplot(obs.amp.h2$global.ncep, obs.amp.h2$global.era, miroc.amp.h2$miroc.global.deca, can.amp.h2$can.global.deca, had.amp.h2$had.global.deca, names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3'), col=c('black', 'purple', rep('red', 3)), cex.axis=0.9, ylim = c(0, 0.25))
  
  boxplot(obs.amp.h2$global.ncep, obs.amp.h2$global.era, miroc.amp.h2$miroc.global.hist, can.amp.h2$can.global.hist, had.amp.h2$had.global.hist, gfdl.amp.h2$gfdl.global.hist, names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col=c('black', 'purple', rep('blue', 4)), main='Global, semi-annual', cex.axis=0.9, ylim = c(0, 0.25))
  
  boxplot(obs.amp.h2$global.ncep, obs.amp.h2$global.era, miroc.amp.h2$miroc.global.control, can.amp.h2$can.global.control, had.amp.h2$had.global.control, gfdl.amp.h2$gfdl.global.control, names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col=c('black', 'purple', rep('green', 4)), cex.axis=0.9, ylim = c(0, 0.25))
}
#Tropical
{
  boxplot(obs.amp.h2$trop.ncep, obs.amp.h2$trop.era, miroc.amp.h2$miroc.trop.deca, can.amp.h2$can.trop.deca, had.amp.h2$had.trop.deca, names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3'), col=c('black', 'purple', rep('red', 3)), cex.axis=0.9, ylim = c(0.2, 0.6))
  
  boxplot(obs.amp.h2$trop.ncep, obs.amp.h2$trop.era, miroc.amp.h2$miroc.trop.hist, can.amp.h2$can.trop.hist, had.amp.h2$had.trop.hist, gfdl.amp.h2$gfdl.trop.hist, names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col=c('black', 'purple', rep('blue', 4)), main='Tropical, semi-annual', cex.axis=0.9, ylim = c(0.2, 0.6))
  
  boxplot(obs.amp.h2$trop.ncep, obs.amp.h2$trop.era, miroc.amp.h2$miroc.trop.control, can.amp.h2$can.trop.control, had.amp.h2$had.trop.control, gfdl.amp.h2$gfdl.trop.control, names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col=c('black', 'purple', rep('green', 4)), cex.axis=0.9, ylim = c(0.2, 0.6))
}
#Northern
{
  boxplot(obs.amp.h2$north.ncep, obs.amp.h2$north.era, miroc.amp.h2$miroc.north.deca, can.amp.h2$can.north.deca, had.amp.h2$had.north.deca, names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3'), col=c('black', 'purple', rep('red', 3)), cex.axis=0.9, ylim = c(0.05, 0.4))
  
  boxplot(obs.amp.h2$north.ncep, obs.amp.h2$north.era, miroc.amp.h2$miroc.north.hist, can.amp.h2$can.north.hist, had.amp.h2$had.north.hist, gfdl.amp.h2$gfdl.north.hist, names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col=c('black', 'purple', rep('blue', 4)), main='Northern, semi-annual', cex.axis=0.9, ylim = c(0.05, 0.4))
  
  boxplot(obs.amp.h2$north.ncep, obs.amp.h2$north.era, miroc.amp.h2$miroc.north.control, can.amp.h2$can.north.control, had.amp.h2$had.north.control, gfdl.amp.h2$gfdl.north.control, names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col=c('black', 'purple', rep('green', 4)), cex.axis=0.9, ylim = c(0.05, 0.4))
}
#Southern
{
  boxplot(obs.amp.h2$south.ncep, obs.amp.h2$south.era, miroc.amp.h2$miroc.south.deca, can.amp.h2$can.south.deca, had.amp.h2$had.south.deca, names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3'), col=c('black', 'purple', rep('red', 3)), cex.axis=0.9, ylim = c(0.0, 0.3))
  
  boxplot(obs.amp.h2$south.ncep, obs.amp.h2$south.era, miroc.amp.h2$miroc.south.hist, can.amp.h2$can.south.hist, had.amp.h2$had.south.hist, gfdl.amp.h2$gfdl.south.hist, names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col=c('black', 'purple', rep('blue', 4)), main='Southern, semi-annual', cex.axis=0.9, ylim = c(0.0, 0.3))
  
  boxplot(obs.amp.h2$south.ncep, obs.amp.h2$south.era, miroc.amp.h2$miroc.south.control, can.amp.h2$can.south.control, had.amp.h2$had.south.control, gfdl.amp.h2$gfdl.south.control, names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), col=c('black', 'purple', rep('green', 4)), cex.axis=0.9, ylim = c(0.0, 0.3))
}


plot.ts(miroc.global.deca[,1], ylim=c(284, 291))
lines(miroc.dlm.list$miroc.global.deca$m[,1])
lines(gfdl.dlm.list$gfdl.global.deca$m[,1], col='red')
lines(ts(gfdl.global.deca[,1]), col='red')
