#CanCM4 Analysis
source("~/Desktop/Climate Model Data/dlm.R")

#initial values
degree = 2
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



#------------------------observational------------------------
data = list(global.ncep=global.ncep, global.era=global.era, trop.ncep=trop.ncep, trop.era=trop.era, north.ncep=north.ncep, north.era=north.era, south.ncep=south.ncep, south.era=south.era)

df.trend <- c(rep(0.96, 2), rep(0.94, 2), rep(0.99, 4))

#fit dlm and save to list
{obs.dlm.list = list()
for (i in 1:length(data)){
  dlm = dlm.fit(data[[i]], df.trend[i], df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
  obs.dlm.list[[i]] = dlm
}
names(obs.dlm.list) = names(data)
}
#credible intervals
{t.mean <- NULL
t.var <- NULL
dist <- NULL
quan <- matrix(NA, 360, 3)
obs.quan.list <- list()
for(i in 1:length(data)){
  for (t in 1:360){
    t.mean <- t(mod.trend$F) %*% obs.dlm.list[[i]]$m[t,trend.index] 
    t.var <- t(mod.trend$F) %*% obs.dlm.list[[i]]$C[t,trend.index,trend.index] %*% as.matrix(mod.trend$F)
    dist <- rpearsonVII(10000, dlm$n[dim(y)[1]], t.mean, sqrt(t.var))
    quan[t,] <- quantile(dist, c(0.025, 0.50, 0.975))
  }
  obs.quan.list[[i]] <- quan
}
names(obs.quan.list) <- names(data)
}

#--------------------------MIROC5--------------------------

data <- list(miroc.global.control=miroc.global.control, miroc.global.deca=miroc.global.deca, miroc.global.hist=miroc.global.hist, miroc.trop.control=miroc.trop.control, miroc.trop.deca=miroc.trop.deca, miroc.trop.hist=miroc.trop.hist, miroc.north.control=miroc.north.control, miroc.north.deca=miroc.north.deca, miroc.north.hist=miroc.north.hist, miroc.south.control=miroc.south.control, miroc.south.deca=miroc.south.deca, miroc.south.hist=miroc.south.hist)

df.trend = list(0.94, 0.94, 0.94, 0.91, 0.91, 0.91, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99)

#fit dlm and save to list
{miroc.dlm.list = list()
  for (i in 1:length(data)){
    dlm = dlm.fit(data[[i]], df.trend[[i]], df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
    miroc.dlm.list[[i]] = dlm
  }
  names(miroc.dlm.list) = names(data)
  
  
  
  #create credible intervals
  library(PearsonDS)
  trend.index = 1:degree
  t.mean <- NULL
  t.var <- NULL
  dist <- NULL
  quan <- matrix(NA, 360, 3)
  miroc.quan.list <- list()
  for(i in 1:length(data)){
    for (t in 1:360){
      t.mean <- t(mod.trend$F) %*% miroc.dlm.list[[i]]$m[t,trend.index] 
      t.var <- t(mod.trend$F) %*% miroc.dlm.list[[i]]$C[t,trend.index,trend.index] %*% as.matrix(mod.trend$F)
      dist <- rpearsonVII(10000, dlm$n[dim(y)[1]], t.mean, sqrt(t.var))
      quan[t,] <- quantile(dist, c(0.025, 0.50, 0.975))
    }
    miroc.quan.list[[i]] <- quan
  }
  names(miroc.quan.list) <- names(data)
}

#plots
{x <- 1:360
  par(mfrow=c(2,2))
  #Global
  {#plot.ts(miroc.global.control[,1])
    plot.ts(miroc.dlm.list$miroc.global.control$m[,1], ylim = c(287,288.7), lwd=1.5, col='green', main='Global', ylab = 'ave temp (K)', xlab = 'time (years)', xaxt="n")
    axis(1, c(1, 120, 240, 360), labels = c("1981", "1991", '2001', '2011'))
    polygon(c(rev(x), x), c(rev(miroc.quan.list$miroc.global.control[ ,3]), miroc.quan.list$miroc.global.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)
    
    #plot.ts(miroc.global.deca[,1])
    lines(miroc.dlm.list$miroc.global.deca$m[,1], lwd=1.5, col='red')
    polygon(c(rev(x), x), c(rev(miroc.quan.list$miroc.global.deca[ ,3]), miroc.quan.list$miroc.global.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)
    
    #plot.ts(miroc.global.hist[,1])
    lines(miroc.dlm.list$miroc.global.hist$m[,1], lwd=1.5, col='blue')
    polygon(c(rev(x), x), c(rev(miroc.quan.list$miroc.global.hist[ ,3]), miroc.quan.list$miroc.global.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)
    
    #plot.ts(miroc.global.hist[,1])
    lines(obs.quan.list$global.ncep[,1], lwd=1.5, col='black')
    polygon(c(rev(x), x), c(rev(obs.quan.list$global.ncep[ ,3]), obs.quan.list$global.ncep[ ,1]), col =rgb(0,0,0, max= 225, alpha=50), border = NA)
    
    lines(ts(obs.dlm.list$global.era$m[,1]), col='purple', lwd=1.5)
    polygon(c(rev(x), x), c(rev(obs.quan.list$global.era[ ,3]), obs.quan.list$global.era[ ,1]), col =rgb(100,0,225, max= 225, alpha=50), border = NA)}
  
  #Tropical
  {#plot.ts(miroc.trop.control[,1])
    plot.ts(miroc.dlm.list$miroc.trop.control$m[,1], lwd=1.5, col='green', ylim = c(297.5,300), main='Tropical', ylab = 'ave temp (K)', xlab = 'time (years)', xaxt="n")
    axis(1, c(1, 120, 240, 360), labels = c("1981", "1991", '2001', '2011'))
    polygon(c(rev(x), x), c(rev(miroc.quan.list$miroc.trop.control[ ,3]), miroc.quan.list$miroc.trop.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)
    
    #plot.ts(miroc.trop.deca[,1])
    lines(miroc.dlm.list$miroc.trop.deca$m[,1], ylim = c(286.5,287.5), lwd=1.5, col='red')
    polygon(c(rev(x), x), c(rev(miroc.quan.list$miroc.trop.deca[ ,3]), miroc.quan.list$miroc.trop.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)
    
    #plot.ts(miroc.trop.hist[,1])
    lines(miroc.dlm.list$miroc.trop.hist$m[,1], ylim = c(286.5,287.5), lwd=1.5, col='blue')
    polygon(c(rev(x), x), c(rev(miroc.quan.list$miroc.trop.hist[ ,3]), miroc.quan.list$miroc.trop.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)
    
    #plot.ts(miroc.trop.ncep[,1])
    lines(obs.dlm.list$trop.ncep$m[,1], ylim = c(286.5,287.5), lwd=1.5, col='black')
    polygon(c(rev(x), x), c(rev(obs.quan.list$trop.ncep[ ,3]), obs.quan.list$trop.ncep[ ,1]), col =rgb(0,0,0, max= 225, alpha=50), border = NA)
    
    lines(ts(obs.dlm.list$trop.era$m[,1]), col='purple', lwd=1.5)
    polygon(c(rev(x), x), c(rev(obs.quan.list$trop.era[ ,3]), obs.quan.list$trop.era[ ,1]), col =rgb(100,0,225, max= 225, alpha=50), border = NA)}
  
  #Northern
  {#plot.ts(miroc.north.control[,1])
    plot.ts(miroc.dlm.list$miroc.north.control$m[,1], ylim = c(287.5,289), lwd=1.5, col='green', main='Northern', ylab = 'ave temp (K)', xlab = 'time (years)', xaxt="n")
    axis(1, c(1, 120, 240, 360), labels = c("1981", "1991", '2001', '2011'))
    polygon(c(rev(x), x), c(rev(miroc.quan.list$miroc.north.control[ ,3]), miroc.quan.list$miroc.north.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)
    
    #plot.ts(miroc.north.deca[,1])
    lines(miroc.dlm.list$miroc.north.deca$m[,1], ylim = c(286.5,287.5), lwd=1.5, col='red')
    polygon(c(rev(x), x), c(rev(miroc.quan.list$miroc.north.deca[ ,3]), miroc.quan.list$miroc.north.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)
    
    #plot.ts(miroc.north.hist[,1])
    lines(miroc.dlm.list$miroc.north.hist$m[,1], ylim = c(286.5,287.5), lwd=1.5, col='blue')
    polygon(c(rev(x), x), c(rev(miroc.quan.list$miroc.north.hist[ ,3]), miroc.quan.list$miroc.north.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)
    
    #plot.ts(miroc.north.hist[,1])
    lines(obs.dlm.list$north.ncep$m[,1], ylim = c(286.5,287.5), lwd=1.5, col='black')
    polygon(c(rev(x), x), c(rev(obs.quan.list$north.ncep[ ,3]), obs.quan.list$north.ncep[ ,1]), col =rgb(0,0,0, max= 225, alpha=50), border = NA)
    
    lines(ts(obs.dlm.list$north.era$m[,1]), col='purple', lwd=1.5)
    polygon(c(rev(x), x), c(rev(obs.quan.list$north.era[ ,3]), obs.quan.list$north.era[ ,1]), col =rgb(100,0,225, max= 225, alpha=50), border = NA)}
  
  #Southern
  {#plot.ts(miroc.south.control[,1])
    plot.ts(miroc.dlm.list$miroc.south.control$m[,1], ylim = c(286, 288.3), lwd=1.5, col='green', main='Southern', ylab = 'ave temp (K)', xlab = 'time (years)', xaxt="n")
    axis(1, c(1, 120, 240, 360), labels = c("1981", "1991", '2001', '2011'))
    polygon(c(rev(x), x), c(rev(miroc.quan.list$miroc.south.control[ ,3]), miroc.quan.list$miroc.south.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)
    
    #plot.ts(miroc.south.deca[,1])
    lines(miroc.dlm.list$miroc.south.deca$m[,1], ylim = c(286.5,287.5), lwd=1.5, col='red')
    polygon(c(rev(x), x), c(rev(miroc.quan.list$miroc.south.deca[ ,3]), miroc.quan.list$miroc.south.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)
    
    #plot.ts(miroc.south.hist[,1])
    lines(miroc.dlm.list$miroc.south.hist$m[,1], ylim = c(286.5,287.5), lwd=1.5, col='blue')
    polygon(c(rev(x), x), c(rev(miroc.quan.list$miroc.south.hist[ ,3]), miroc.quan.list$miroc.south.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)
    
    #plot.ts(miroc.south.hist[,1])
    lines(obs.dlm.list$south.ncep$m[,1], ylim = c(286.5,287.5), lwd=1.5, col='black')
    polygon(c(rev(x), x), c(rev(obs.quan.list$south.ncep[ ,3]), obs.quan.list$south.ncep[ ,1]), col =rgb(0,0,0, max= 225, alpha=50), border = NA)
    
    lines(ts(obs.dlm.list$south.era$m[,1]), col='purple', lwd=1.5)
    polygon(c(rev(x), x), c(rev(obs.quan.list$south.era[ ,3]), obs.quan.list$south.era[ ,1]), col =rgb(100,0,225, max= 225, alpha=50), border = NA)}
}
#--------------------------CanCM4--------------------------

data <- list(can.global.control=can.global.control, can.global.deca=can.global.deca, can.global.hist=can.global.hist, can.trop.control=can.trop.control, can.trop.deca=can.trop.deca, can.trop.hist=can.trop.hist, can.north.control=can.north.control, can.north.deca=can.north.deca, can.north.hist=can.north.hist, can.south.control=can.south.control, can.south.deca=can.south.deca, can.south.hist=can.south.hist)

df.trend = list(0.94, 0.94, 0.94, 0.87, 0.87, 0.87, 0.99, 0.99, 0.99, 0.97, 0.97, 0.97)

#fit dlm and save to list
{can.dlm.list = list()
for (i in 1:length(data)){
  dlm = dlm.fit(data[[i]], df.trend[[i]], df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
  can.dlm.list[[i]] = dlm
}
names(can.dlm.list) = names(data)



#create credible intervals
library(PearsonDS)
trend.index = 1:degree
t.mean <- NULL
t.var <- NULL
dist <- NULL
quan <- matrix(NA, 360, 3)
can.quan.list <- list()
for(i in 1:length(data)){
  for (t in 1:360){
    t.mean <- t(mod.trend$F) %*% can.dlm.list[[i]]$m[t,trend.index] 
    t.var <- t(mod.trend$F) %*% can.dlm.list[[i]]$C[t,trend.index,trend.index] %*% as.matrix(mod.trend$F)
    dist <- rpearsonVII(10000, dlm$n[dim(y)[1]], t.mean, sqrt(t.var))
    quan[t,] <- quantile(dist, c(0.025, 0.50, 0.975))
  }
  can.quan.list[[i]] <- quan
}
names(can.quan.list) <- names(data)
}

#plots
{x <- 1:360
par(mfrow=c(2,2))
#Global
{#plot.ts(can.global.control[,1])
plot.ts(can.dlm.list$can.global.control$m[,1], ylim = c(286.5,288), lwd=1.5, col='green', main='Global', ylab = 'ave temp (K)', xlab = 'time (years)', xaxt="n")
axis(1, c(1, 120, 240, 360), labels = c("1981", "1991", '2001', '2011'))
polygon(c(rev(x), x), c(rev(can.quan.list$can.global.control[ ,3]), can.quan.list$can.global.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)

#plot.ts(can.global.deca[,1])
lines(can.dlm.list$can.global.deca$m[,1], ylim = c(286.5,287.5), lwd=1.5, col='red')
polygon(c(rev(x), x), c(rev(can.quan.list$can.global.deca[ ,3]), can.quan.list$can.global.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)

#plot.ts(can.global.hist[,1])
lines(can.dlm.list$can.global.hist$m[,1], ylim = c(286.5,287.5), lwd=1.5, col='blue')
polygon(c(rev(x), x), c(rev(can.quan.list$can.global.hist[ ,3]), can.quan.list$can.global.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)

#plot.ts(can.global.hist[,1])
lines(obs.quan.list$global.ncep[,1], ylim = c(286.5,287.5), lwd=1.5, col='black')
polygon(c(rev(x), x), c(rev(obs.quan.list$global.ncep[ ,3]), obs.quan.list$global.ncep[ ,1]), col =rgb(0,0,0, max= 225, alpha=50), border = NA)

lines(ts(obs.dlm.list$global.era$m[,1]), col='purple', lwd=1.5)
polygon(c(rev(x), x), c(rev(obs.quan.list$global.era[ ,3]), obs.quan.list$global.era[ ,1]), col =rgb(100,0,225, max= 225, alpha=50), border = NA)}

#Tropical
{#plot.ts(can.trop.control[,1])
  plot.ts(can.dlm.list$can.trop.control$m[,1], lwd=1.5, col='green', ylim = c(298,300), main='Tropical', ylab = 'ave temp (K)', xlab = 'time (years)', xaxt="n")
  axis(1, c(1, 120, 240, 360), labels = c("1981", "1991", '2001', '2011'))
  polygon(c(rev(x), x), c(rev(can.quan.list$can.trop.control[ ,3]), can.quan.list$can.trop.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)
  
  #plot.ts(can.trop.deca[,1])
  lines(can.dlm.list$can.trop.deca$m[,1], ylim = c(286.5,287.5), lwd=1.5, col='red')
  polygon(c(rev(x), x), c(rev(can.quan.list$can.trop.deca[ ,3]), can.quan.list$can.trop.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)
  
  #plot.ts(can.trop.hist[,1])
  lines(can.dlm.list$can.trop.hist$m[,1], ylim = c(286.5,287.5), lwd=1.5, col='blue')
  polygon(c(rev(x), x), c(rev(can.quan.list$can.trop.hist[ ,3]), can.quan.list$can.trop.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)
  
  #plot.ts(can.trop.ncep[,1])
  lines(obs.dlm.list$trop.ncep$m[,1], ylim = c(286.5,287.5), lwd=1.5, col='black')
  polygon(c(rev(x), x), c(rev(obs.quan.list$trop.ncep[ ,3]), obs.quan.list$trop.ncep[ ,1]), col =rgb(0,0,0, max= 225, alpha=50), border = NA)
  
  lines(ts(obs.dlm.list$trop.era$m[,1]), col='purple', lwd=1.5)
  polygon(c(rev(x), x), c(rev(obs.quan.list$trop.era[ ,3]), obs.quan.list$trop.era[ ,1]), col =rgb(100,0,225, max= 225, alpha=50), border = NA)}

#Northern
{#plot.ts(can.north.control[,1])
  plot.ts(can.dlm.list$can.north.control$m[,1], ylim = c(287,289), lwd=1.5, col='green', main='Northern', ylab = 'ave temp (K)', xlab = 'time (years)', xaxt="n")
  axis(1, c(1, 120, 240, 360), labels = c("1981", "1991", '2001', '2011'))
  polygon(c(rev(x), x), c(rev(can.quan.list$can.north.control[ ,3]), can.quan.list$can.north.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)
  
  #plot.ts(can.north.deca[,1])
  lines(can.dlm.list$can.north.deca$m[,1], ylim = c(286.5,287.5), lwd=1.5, col='red')
  polygon(c(rev(x), x), c(rev(can.quan.list$can.north.deca[ ,3]), can.quan.list$can.north.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)
  
  #plot.ts(can.north.hist[,1])
  lines(can.dlm.list$can.north.hist$m[,1], ylim = c(286.5,287.5), lwd=1.5, col='blue')
  polygon(c(rev(x), x), c(rev(can.quan.list$can.north.hist[ ,3]), can.quan.list$can.north.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)
  
  #plot.ts(can.north.hist[,1])
  lines(obs.dlm.list$north.ncep$m[,1], ylim = c(286.5,287.5), lwd=1.5, col='black')
  polygon(c(rev(x), x), c(rev(obs.quan.list$north.ncep[ ,3]), obs.quan.list$north.ncep[ ,1]), col =rgb(0,0,0, max= 225, alpha=50), border = NA)
  
  lines(ts(obs.dlm.list$north.era$m[,1]), col='purple', lwd=1.5)
  polygon(c(rev(x), x), c(rev(obs.quan.list$north.era[ ,3]), obs.quan.list$north.era[ ,1]), col =rgb(100,0,225, max= 225, alpha=50), border = NA)}

#Southern
{#plot.ts(can.south.control[,1])
  plot.ts(can.dlm.list$can.south.control$m[,1], ylim = c(285.5, 287.5), lwd=1.5, col='green', main='Southern', ylab = 'ave temp (K)', xlab = 'time (years)', xaxt="n")
  axis(1, c(1, 120, 240, 360), labels = c("1981", "1991", '2001', '2011'))
  polygon(c(rev(x), x), c(rev(can.quan.list$can.south.control[ ,3]), can.quan.list$can.south.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)
  
  #plot.ts(can.south.deca[,1])
  lines(can.dlm.list$can.south.deca$m[,1], ylim = c(286.5,287.5), lwd=1.5, col='red')
  polygon(c(rev(x), x), c(rev(can.quan.list$can.south.deca[ ,3]), can.quan.list$can.south.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)
  
  #plot.ts(can.south.hist[,1])
  lines(can.dlm.list$can.south.hist$m[,1], ylim = c(286.5,287.5), lwd=1.5, col='blue')
  polygon(c(rev(x), x), c(rev(can.quan.list$can.south.hist[ ,3]), can.quan.list$can.south.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)
  
  #plot.ts(can.south.hist[,1])
  lines(obs.dlm.list$south.ncep$m[,1], ylim = c(286.5,287.5), lwd=1.5, col='black')
  polygon(c(rev(x), x), c(rev(obs.quan.list$south.ncep[ ,3]), obs.quan.list$south.ncep[ ,1]), col =rgb(0,0,0, max= 225, alpha=50), border = NA)
  
  lines(ts(obs.dlm.list$south.era$m[,1]), col='purple', lwd=1.5)
  polygon(c(rev(x), x), c(rev(obs.quan.list$south.era[ ,3]), obs.quan.list$south.era[ ,1]), col =rgb(100,0,225, max= 225, alpha=50), border = NA)}
}
#--------------------------HadCM3--------------------------

data <- list(had.global.control=had.global.control, had.global.deca=had.global.deca, had.global.hist=had.global.hist, had.trop.control=had.trop.control, had.trop.deca=had.trop.deca, had.trop.hist=had.trop.hist, had.north.control=had.north.control, had.north.deca=had.north.deca, had.north.hist=had.north.hist, had.south.control=had.south.control, had.south.deca=had.south.deca, had.south.hist=had.south.hist)

df.trend = list(0.94, 0.94, 0.94, 0.87, 0.87, 0.87, 0.99, 0.99, 0.99, 0.97, 0.97, 0.97)

#fit dlm and save to list
{had.dlm.list = list()
for (i in 1:length(data)){
  dlm = dlm.fit(data[[i]], df.trend[[i]], df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
  had.dlm.list[[i]] = dlm
}
names(had.dlm.list) = names(data)
}

#create credible intervals
{library(PearsonDS)
trend.index = 1:degree
t.mean <- NULL
t.var <- NULL
dist <- NULL
quan <- matrix(NA, 360, 3)
had.quan.list <- list()
for(i in 1:length(data)){
  for (t in 1:360){
    t.mean <- t(mod.trend$F) %*% had.dlm.list[[i]]$m[t,trend.index] 
    t.var <- t(mod.trend$F) %*% had.dlm.list[[i]]$C[t,trend.index,trend.index] %*% as.matrix(mod.trend$F)
    dist <- rpearsonVII(10000, dlm$n[dim(y)[1]], t.mean, sqrt(t.var))
    quan[t,] <- quantile(dist, c(0.025, 0.50, 0.975))
  }
  had.quan.list[[i]] <- quan
}
names(had.quan.list) <- names(data)
}

#plots
{x <- 1:360
par(mfrow=c(2,2))
#Global
{#plot.ts(had.global.control[,1])
  plot.ts(had.dlm.list$had.global.control$m[,1], ylim = c(286,288), lwd=1.5, col='green', main='Global', ylab = 'ave temp (K)', xlab = 'time (years)', xaxt="n")
  axis(1, c(1, 120, 240, 360), labels = c("1981", "1991", '2001', '2011'))
  polygon(c(rev(x), x), c(rev(had.quan.list$had.global.control[ ,3]), had.quan.list$had.global.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)
  
  #plot.ts(had.global.deca[,1])
  lines(had.dlm.list$had.global.deca$m[,1], ylim = c(286.5,287.5), lwd=1.5, col='red')
  polygon(c(rev(x), x), c(rev(had.quan.list$had.global.deca[ ,3]), had.quan.list$had.global.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)
  
  #plot.ts(had.global.hist[,1])
  lines(had.dlm.list$had.global.hist$m[,1], ylim = c(286.5,287.5), lwd=1.5, col='blue')
  polygon(c(rev(x), x), c(rev(had.quan.list$had.global.hist[ ,3]), had.quan.list$had.global.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)
  
  #plot.ts(had.global.hist[,1])
  lines(obs.quan.list$global.ncep[,1], ylim = c(286.5,287.5), lwd=1.5, col='black')
  polygon(c(rev(x), x), c(rev(obs.quan.list$global.ncep[ ,3]), obs.quan.list$global.ncep[ ,1]), col =rgb(0,0,0, max= 225, alpha=50), border = NA)
  
  lines(ts(obs.dlm.list$global.era$m[,1]), col='purple', lwd=1.5)
  polygon(c(rev(x), x), c(rev(obs.quan.list$global.era[ ,3]), obs.quan.list$global.era[ ,1]), col =rgb(100,0,225, max= 225, alpha=50), border = NA)}

#Tropical
{#plot.ts(had.trop.control[,1])
  plot.ts(had.dlm.list$had.trop.control$m[,1], lwd=1.5, col='green', ylim = c(298,300.5), main='Tropical', ylab = 'ave temp (K)', xlab = 'time (years)', xaxt="n")
  axis(1, c(1, 120, 240, 360), labels = c("1981", "1991", '2001', '2011'))
  polygon(c(rev(x), x), c(rev(had.quan.list$had.trop.control[ ,3]), had.quan.list$had.trop.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)
  
  #plot.ts(had.trop.deca[,1])
  lines(had.dlm.list$had.trop.deca$m[,1], ylim = c(286.5,287.5), lwd=1.5, col='red')
  polygon(c(rev(x), x), c(rev(had.quan.list$had.trop.deca[ ,3]), had.quan.list$had.trop.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)
  
  #plot.ts(had.trop.hist[,1])
  lines(had.dlm.list$had.trop.hist$m[,1], ylim = c(286.5,287.5), lwd=1.5, col='blue')
  polygon(c(rev(x), x), c(rev(had.quan.list$had.trop.hist[ ,3]), had.quan.list$had.trop.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)
  
  #plot.ts(had.trop.ncep[,1])
  lines(obs.dlm.list$trop.ncep$m[,1], ylim = c(286.5,287.5), lwd=1.5, col='black')
  polygon(c(rev(x), x), c(rev(obs.quan.list$trop.ncep[ ,3]), obs.quan.list$trop.ncep[ ,1]), col =rgb(0,0,0, max= 225, alpha=50), border = NA)
  
  lines(ts(obs.dlm.list$trop.era$m[,1]), col='purple', lwd=1.5)
  polygon(c(rev(x), x), c(rev(obs.quan.list$trop.era[ ,3]), obs.quan.list$trop.era[ ,1]), col =rgb(100,0,225, max= 225, alpha=50), border = NA)}

#Northern
{#plot.ts(had.north.control[,1])
  plot.ts(had.dlm.list$had.north.control$m[,1], ylim = c(286,289), lwd=1.5, col='green', main='Northern', ylab = 'ave temp (K)', xlab = 'time (years)', xaxt="n")
  axis(1, c(1, 120, 240, 360), labels = c("1981", "1991", '2001', '2011'))
  polygon(c(rev(x), x), c(rev(had.quan.list$had.north.control[ ,3]), had.quan.list$had.north.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)
  
  #plot.ts(had.north.deca[,1])
  lines(had.dlm.list$had.north.deca$m[,1], lwd=1.5, col='red')
  polygon(c(rev(x), x), c(rev(had.quan.list$had.north.deca[ ,3]), had.quan.list$had.north.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)
  
  #plot.ts(had.north.hist[,1])
  lines(had.dlm.list$had.north.hist$m[,1], lwd=1.5, col='blue')
  polygon(c(rev(x), x), c(rev(had.quan.list$had.north.hist[ ,3]), had.quan.list$had.north.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)
  
  #plot.ts(had.north.hist[,1])
  lines(obs.dlm.list$north.ncep$m[,1], lwd=1.5, col='black')
  polygon(c(rev(x), x), c(rev(obs.quan.list$north.ncep[ ,3]), obs.quan.list$north.ncep[ ,1]), col =rgb(0,0,0, max= 225, alpha=50), border = NA)
  
  lines(ts(obs.dlm.list$north.era$m[,1]), col='purple', lwd=1.5)
  polygon(c(rev(x), x), c(rev(obs.quan.list$north.era[ ,3]), obs.quan.list$north.era[ ,1]), col =rgb(100,0,225, max= 225, alpha=50), border = NA)}

#Southern
{#plot.ts(had.south.control[,1])
  plot.ts(had.dlm.list$had.south.control$m[,1], ylim = c(285.5, 287.5), lwd=1.5, col='green', main='Southern', ylab = 'ave temp (K)', xlab = 'time (years)', xaxt="n")
  axis(1, c(1, 120, 240, 360), labels = c("1981", "1991", '2001', '2011'))
  polygon(c(rev(x), x), c(rev(had.quan.list$had.south.control[ ,3]), had.quan.list$had.south.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)
  
  #plot.ts(had.south.deca[,1])
  lines(had.dlm.list$had.south.deca$m[,1], ylim = c(286.5,287.5), lwd=1.5, col='red')
  polygon(c(rev(x), x), c(rev(had.quan.list$had.south.deca[ ,3]), had.quan.list$had.south.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)
  
  #plot.ts(had.south.hist[,1])
  lines(had.dlm.list$had.south.hist$m[,1], ylim = c(286.5,287.5), lwd=1.5, col='blue')
  polygon(c(rev(x), x), c(rev(had.quan.list$had.south.hist[ ,3]), had.quan.list$had.south.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)
  
  #plot.ts(had.south.hist[,1])
  lines(obs.dlm.list$south.ncep$m[,1], ylim = c(286.5,287.5), lwd=1.5, col='black')
  polygon(c(rev(x), x), c(rev(obs.quan.list$south.ncep[ ,3]), obs.quan.list$south.ncep[ ,1]), col =rgb(0,0,0, max= 225, alpha=50), border = NA)
  
  lines(ts(obs.dlm.list$south.era$m[,1]), col='purple', lwd=1.5)
  polygon(c(rev(x), x), c(rev(obs.quan.list$south.era[ ,3]), obs.quan.list$south.era[ ,1]), col =rgb(100,0,225, max= 225, alpha=50), border = NA)}
}

#--------------------------GFDL--------------------------

data <- list(gfdl.global.control=gfdl.global.control, gfdl.global.deca=gfdl.global.deca, gfdl.global.hist=gfdl.global.hist, gfdl.trop.control=gfdl.trop.control, gfdl.trop.deca=gfdl.trop.deca, gfdl.trop.hist=gfdl.trop.hist, gfdl.north.control=gfdl.north.control, gfdl.north.deca=gfdl.north.deca, gfdl.north.hist=gfdl.north.hist, gfdl.south.control=gfdl.south.control, gfdl.south.deca=gfdl.south.deca, gfdl.south.hist=gfdl.south.hist)

df.trend = list(0.88, 0.88, 0.88, 0.91, 0.91, 0.91, 0.99, 0.99, 0.99, 0.94, 0.94, 0.94)

#fit dlm and save to list
{gfdl.dlm.list = list()
  for (i in 1:length(data)){
    dlm = dlm.fit(data[[i]], df.trend[[i]], df.seas, m0, C0, n0, S0, mod.trend$F, mod.sea$F, mod.trend$G, mod.sea$G, FALSE)
    gfdl.dlm.list[[i]] = dlm
  }
  names(gfdl.dlm.list) = names(data)
}

#create credible intervals
{library(PearsonDS)
  trend.index = 1:degree
  t.mean <- NULL
  t.var <- NULL
  dist <- NULL
  quan <- matrix(NA, 360, 3)
  gfdl.quan.list <- list()
  for(i in 1:length(data)){
    for (t in 1:360){
      t.mean <- t(mod.trend$F) %*% gfdl.dlm.list[[i]]$m[t,trend.index] 
      t.var <- t(mod.trend$F) %*% gfdl.dlm.list[[i]]$C[t,trend.index,trend.index] %*% as.matrix(mod.trend$F)
      dist <- rpearsonVII(10000, dlm$n[dim(y)[1]], t.mean, sqrt(t.var))
      quan[t,] <- quantile(dist, c(0.025, 0.50, 0.975))
    }
    gfdl.quan.list[[i]] <- quan
  }
  names(gfdl.quan.list) <- names(data)
}

#plots
{x <- 1:360
  par(mfrow=c(2,2))
  #Global
  {#plot.ts(gfdl.global.control[,1])
    plot.ts(gfdl.dlm.list$gfdl.global.control$m[,1], ylim = c(286,288), lwd=1.5, col='green', main='Global', ylab = 'ave temp (K)', xlab = 'time (years)', xaxt="n")
    axis(1, c(1, 120, 240, 360), labels = c("1981", "1991", '2001', '2011'))
    polygon(c(rev(x), x), c(rev(gfdl.quan.list$gfdl.global.control[ ,3]), gfdl.quan.list$gfdl.global.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)
    
    #plot.ts(gfdl.global.deca[,1])
    lines(gfdl.dlm.list$gfdl.global.deca$m[,1], lwd=1.5, col='red')
    polygon(c(rev(x), x), c(rev(gfdl.quan.list$gfdl.global.deca[ ,3]), gfdl.quan.list$gfdl.global.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)
    
    #plot.ts(gfdl.global.hist[,1])
    lines(gfdl.dlm.list$gfdl.global.hist$m[,1], lwd=1.5, col='blue')
    polygon(c(rev(x), x), c(rev(gfdl.quan.list$gfdl.global.hist[ ,3]), gfdl.quan.list$gfdl.global.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)
    
    #plot.ts(gfdl.global.ncep[,1])
    lines(obs.quan.list$global.ncep[,1], lwd=1.5, col='black')
    polygon(c(rev(x), x), c(rev(obs.quan.list$global.ncep[ ,3]), obs.quan.list$global.ncep[ ,1]), col =rgb(0,0,0, max= 225, alpha=50), border = NA)
    
    lines(ts(obs.dlm.list$global.era$m[,1]), col='purple', lwd=1.5)
    polygon(c(rev(x), x), c(rev(obs.quan.list$global.era[ ,3]), obs.quan.list$global.era[ ,1]), col =rgb(100,0,225, max= 225, alpha=50), border = NA)}
  
  #Tropical
  {#plot.ts(gfdl.trop.control[,1])
    plot.ts(gfdl.dlm.list$gfdl.trop.control$m[,1], lwd=1.5, col='green', ylim = c(297.5,299.7), main='Tropical', ylab = 'ave temp (K)', xlab = 'time (years)', xaxt="n")
    axis(1, c(1, 120, 240, 360), labels = c("1981", "1991", '2001', '2011'))
    polygon(c(rev(x), x), c(rev(gfdl.quan.list$gfdl.trop.control[ ,3]), gfdl.quan.list$gfdl.trop.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)
    
    #plot.ts(gfdl.trop.deca[,1])
    lines(gfdl.dlm.list$gfdl.trop.deca$m[,1], lwd=1.5, col='red')
    polygon(c(rev(x), x), c(rev(gfdl.quan.list$gfdl.trop.deca[ ,3]), gfdl.quan.list$gfdl.trop.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)
    
    #plot.ts(gfdl.trop.hist[,1])
    lines(gfdl.dlm.list$gfdl.trop.hist$m[,1], lwd=1.5, col='blue')
    polygon(c(rev(x), x), c(rev(gfdl.quan.list$gfdl.trop.hist[ ,3]), gfdl.quan.list$gfdl.trop.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)
    
    #plot.ts(gfdl.trop.ncep[,1])
    lines(obs.dlm.list$trop.ncep$m[,1], lwd=1.5, col='black')
    polygon(c(rev(x), x), c(rev(obs.quan.list$trop.ncep[ ,3]), obs.quan.list$trop.ncep[ ,1]), col =rgb(0,0,0, max= 225, alpha=50), border = NA)
    
    lines(ts(obs.dlm.list$trop.era$m[,1]), col='purple', lwd=1.5)
    polygon(c(rev(x), x), c(rev(obs.quan.list$trop.era[ ,3]), obs.quan.list$trop.era[ ,1]), col =rgb(100,0,225, max= 225, alpha=50), border = NA)}
  
  #Northern
  {#plot.ts(gfdl.north.control[,1])
    plot.ts(gfdl.dlm.list$gfdl.north.control$m[,1], ylim = c(286,289), lwd=1.5, col='green', main='Northern', ylab = 'ave temp (K)', xlab = 'time (years)', xaxt="n")
    axis(1, c(1, 120, 240, 360), labels = c("1981", "1991", '2001', '2011'))
    polygon(c(rev(x), x), c(rev(gfdl.quan.list$gfdl.north.control[ ,3]), gfdl.quan.list$gfdl.north.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)
    
    #plot.ts(gfdl.north.deca[,1])
    lines(gfdl.dlm.list$gfdl.north.deca$m[,1], lwd=1.5, col='red')
    polygon(c(rev(x), x), c(rev(gfdl.quan.list$gfdl.north.deca[ ,3]), gfdl.quan.list$gfdl.north.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)
    
    #plot.ts(gfdl.north.hist[,1])
    lines(gfdl.dlm.list$gfdl.north.hist$m[,1], lwd=1.5, col='blue')
    polygon(c(rev(x), x), c(rev(gfdl.quan.list$gfdl.north.hist[ ,3]), gfdl.quan.list$gfdl.north.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)
    
    #plot.ts(gfdl.north.hist[,1])
    lines(obs.dlm.list$north.ncep$m[,1], lwd=1.5, col='black')
    polygon(c(rev(x), x), c(rev(obs.quan.list$north.ncep[ ,3]), obs.quan.list$north.ncep[ ,1]), col =rgb(0,0,0, max= 225, alpha=50), border = NA)
    
    lines(ts(obs.dlm.list$north.era$m[,1]), col='purple', lwd=1.5)
    polygon(c(rev(x), x), c(rev(obs.quan.list$north.era[ ,3]), obs.quan.list$north.era[ ,1]), col =rgb(100,0,225, max= 225, alpha=50), border = NA)}
  
  #Southern
  {#plot.ts(gfdl.south.control[,1])
    plot.ts(gfdl.dlm.list$gfdl.south.control$m[,1], ylim = c(286.2, 288), lwd=1.5, col='green', main='Southern', ylab = 'ave temp (K)', xlab = 'time (years)', xaxt="n")
    axis(1, c(1, 120, 240, 360), labels = c("1981", "1991", '2001', '2011'))
    polygon(c(rev(x), x), c(rev(gfdl.quan.list$gfdl.south.control[ ,3]), gfdl.quan.list$gfdl.south.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)
    
    #plot.ts(gfdl.south.deca[,1])
    lines(gfdl.dlm.list$gfdl.south.deca$m[,1], lwd=1.5, col='red')
    polygon(c(rev(x), x), c(rev(gfdl.quan.list$gfdl.south.deca[ ,3]), gfdl.quan.list$gfdl.south.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)
    
    #plot.ts(gfdl.south.hist[,1])
    lines(gfdl.dlm.list$gfdl.south.hist$m[,1], lwd=1.5, col='blue')
    polygon(c(rev(x), x), c(rev(gfdl.quan.list$gfdl.south.hist[ ,3]), gfdl.quan.list$gfdl.south.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)
    

    lines(obs.dlm.list$south.ncep$m[,1], lwd=1.5, col='black')
    polygon(c(rev(x), x), c(rev(obs.quan.list$south.ncep[ ,3]), obs.quan.list$south.ncep[ ,1]), col =rgb(0,0,0, max= 225, alpha=50), border = NA)
    
    lines(ts(obs.dlm.list$south.era$m[,1]), col='purple', lwd=1.5)
    polygon(c(rev(x), x), c(rev(obs.quan.list$south.era[ ,3]), obs.quan.list$south.era[ ,1]), col =rgb(100,0,225, max= 225, alpha=50), border = NA)}
}

#---------------Summary Plot-----------------
par(mar=c(3,3,1,1))
par(mfrow=c(4,3))
#Global
{
  #Decadal
  {
    plot.ts(miroc.dlm.list$miroc.global.deca$m[,1], lwd=1.5, col='red', ylim=c(286, 289), main='', xaxt="n", ylab="ave temp (K)")
    axis(1, c(1, 120, 240, 360), labels = c("1981", "1991", '2001', '2011'))
    polygon(c(rev(x), x), c(rev(miroc.quan.list$miroc.global.deca[ ,3]), miroc.quan.list$miroc.global.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)
    
    lines(can.dlm.list$can.global.deca$m[,1], lwd=1.5, col='red', lty=2)
    polygon(c(rev(x), x), c(rev(can.quan.list$can.global.deca[ ,3]), can.quan.list$can.global.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)
    
    lines(had.dlm.list$had.global.deca$m[,1], lwd=1.5, col='red', lty=3)
    polygon(c(rev(x), x), c(rev(had.quan.list$had.global.deca[ ,3]), had.quan.list$had.global.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)
    
    lines(gfdl.dlm.list$gfdl.global.deca$m[,1], lwd=1.5, col='red', lty=4)
    polygon(c(rev(x), x), c(rev(gfdl.quan.list$gfdl.global.deca[ ,3]), gfdl.quan.list$gfdl.global.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)
    
    # lines(obs.quan.list$global.ncep[,1], lwd=1.5, col='black')
    # polygon(c(rev(x), x), c(rev(obs.quan.list$global.ncep[ ,3]), obs.quan.list$global.ncep[ ,1]), col =rgb(0,0,0, max= 225, alpha=50), border = NA)
    # 
    # lines(ts(obs.dlm.list$global.era$m[,1]), col='purple', lwd=1.5)
    # polygon(c(rev(x), x), c(rev(obs.quan.list$global.era[ ,3]), obs.quan.list$global.era[ ,1]), col =rgb(100,0,225, max= 225, alpha=50), border = NA)
  }
  #historical
  {
    plot.ts(miroc.dlm.list$miroc.global.hist$m[,1], lwd=1.5, col='blue', ylim=c(286, 289), main='Global', xaxt="n", ylab="ave temp (K)")
    axis(1, c(1, 120, 240, 360), labels = c("1981", "1991", '2001', '2011'))
    
    polygon(c(rev(x), x), c(rev(miroc.quan.list$miroc.global.hist[ ,3]), miroc.quan.list$miroc.global.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)
    
    lines(can.dlm.list$can.global.hist$m[,1], lwd=1.5, col='blue', lty=2)
    polygon(c(rev(x), x), c(rev(can.quan.list$can.global.hist[ ,3]), can.quan.list$can.global.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)
    
    lines(had.dlm.list$had.global.hist$m[,1], lwd=1.5, col='blue', lty=3)
    polygon(c(rev(x), x), c(rev(had.quan.list$had.global.hist[ ,3]), had.quan.list$had.global.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)
    
    lines(gfdl.dlm.list$gfdl.global.hist$m[,1], lwd=1.5, col='blue', lty=4)
    polygon(c(rev(x), x), c(rev(gfdl.quan.list$gfdl.global.hist[ ,3]), gfdl.quan.list$gfdl.global.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)
    
    # lines(obs.quan.list$global.ncep[,1], lwd=1.5, col='black')
    # polygon(c(rev(x), x), c(rev(obs.quan.list$global.ncep[ ,3]), obs.quan.list$global.ncep[ ,1]), col =rgb(0,0,0, max= 225, alpha=50), border = NA)
    # 
    # lines(ts(obs.dlm.list$global.era$m[,1]), col='purple', lwd=1.5)
    # polygon(c(rev(x), x), c(rev(obs.quan.list$global.era[ ,3]), obs.quan.list$global.era[ ,1]), col =rgb(100,0,225, max= 225, alpha=50), border = NA)
  }
  #control
  {
    plot.ts(miroc.dlm.list$miroc.global.control$m[,1], lwd=1.5, col='green', ylim=c(286, 289), main='', xaxt="n", ylab="ave temp (K)")
    axis(1, c(1, 120, 240, 360), labels = c("1981", "1991", '2001', '2011'))
    
    polygon(c(rev(x), x), c(rev(miroc.quan.list$miroc.global.control[ ,3]), miroc.quan.list$miroc.global.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)
    
    lines(can.dlm.list$can.global.control$m[,1], lwd=1.5, col='green', lty=2)
    polygon(c(rev(x), x), c(rev(can.quan.list$can.global.control[ ,3]), can.quan.list$can.global.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)
    
    lines(had.dlm.list$had.global.control$m[,1], lwd=1.5, col='green', lty=3)
    polygon(c(rev(x), x), c(rev(had.quan.list$had.global.control[ ,3]), had.quan.list$had.global.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)
    
    lines(gfdl.dlm.list$gfdl.global.control$m[,1], lwd=1.5, col='green', lty=4)
    polygon(c(rev(x), x), c(rev(gfdl.quan.list$gfdl.global.control[ ,3]), gfdl.quan.list$gfdl.global.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)
    
  #   lines(obs.quan.list$global.ncep[,1], lwd=1.5, col='black')
  #   polygon(c(rev(x), x), c(rev(obs.quan.list$global.ncep[ ,3]), obs.quan.list$global.ncep[ ,1]), col =rgb(0,0,0, max= 225, alpha=50), border = NA)
  #   
  #   lines(ts(obs.dlm.list$global.era$m[,1]), col='purple', lwd=1.5)
  #   polygon(c(rev(x), x), c(rev(obs.quan.list$global.era[ ,3]), obs.quan.list$global.era[ ,1]), col =rgb(100,0,225, max= 225, alpha=50), border = NA)
  # }
  }
}
#Tropical
{
  #Decadal
  {
    plot.ts(miroc.dlm.list$miroc.trop.deca$m[,1], lwd=1.5, col='red', main='', ylim=c(297.5, 300.5), xaxt="n", ylab="ave temp (K)")
    axis(1, c(1, 120, 240, 360), labels = c("1981", "1991", '2001', '2011'))
    
    polygon(c(rev(x), x), c(rev(miroc.quan.list$miroc.trop.deca[ ,3]), miroc.quan.list$miroc.trop.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)
    
    lines(can.dlm.list$can.trop.deca$m[,1], lwd=1.5, col='red', lty=2)
    polygon(c(rev(x), x), c(rev(can.quan.list$can.trop.deca[ ,3]), can.quan.list$can.trop.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)
    
    lines(had.dlm.list$had.trop.deca$m[,1], lwd=1.5, col='red', lty=3)
    polygon(c(rev(x), x), c(rev(had.quan.list$had.trop.deca[ ,3]), had.quan.list$had.trop.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)
    
    lines(gfdl.dlm.list$gfdl.trop.deca$m[,1], lwd=1.5, col='red', lty=4)
    polygon(c(rev(x), x), c(rev(gfdl.quan.list$gfdl.trop.deca[ ,3]), gfdl.quan.list$gfdl.trop.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)
    
    # lines(obs.quan.list$trop.ncep[,1], lwd=1.5, col='black')
    # polygon(c(rev(x), x), c(rev(obs.quan.list$trop.ncep[ ,3]), obs.quan.list$trop.ncep[ ,1]), col =rgb(0,0,0, max= 225, alpha=50), border = NA)
    # 
    # lines(ts(obs.dlm.list$trop.era$m[,1]), col='purple', lwd=1.5)
    # polygon(c(rev(x), x), c(rev(obs.quan.list$trop.era[ ,3]), obs.quan.list$trop.era[ ,1]), col =rgb(100,0,225, max= 225, alpha=50), border = NA)
  }
  #historical
  {
    plot.ts(miroc.dlm.list$miroc.trop.hist$m[,1], lwd=1.5, col='blue', ylim=c(297.5, 300.5), main='Tropical', xaxt="n", ylab="ave temp (K)")
    axis(1, c(1, 120, 240, 360), labels = c("1981", "1991", '2001', '2011'))
    
    polygon(c(rev(x), x), c(rev(miroc.quan.list$miroc.trop.hist[ ,3]), miroc.quan.list$miroc.trop.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)
    
    lines(can.dlm.list$can.trop.hist$m[,1], lwd=1.5, col='blue', lty=2)
    polygon(c(rev(x), x), c(rev(can.quan.list$can.trop.hist[ ,3]), can.quan.list$can.trop.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)
    
    lines(had.dlm.list$had.trop.hist$m[,1], lwd=1.5, col='blue', lty=3)
    polygon(c(rev(x), x), c(rev(had.quan.list$had.trop.hist[ ,3]), had.quan.list$had.trop.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)
    
    lines(gfdl.dlm.list$gfdl.trop.hist$m[,1], lwd=1.5, col='blue', lty=4)
    polygon(c(rev(x), x), c(rev(gfdl.quan.list$gfdl.trop.hist[ ,3]), gfdl.quan.list$gfdl.trop.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)
    
    # lines(obs.quan.list$trop.ncep[,1], lwd=1.5, col='black')
    # polygon(c(rev(x), x), c(rev(obs.quan.list$trop.ncep[ ,3]), obs.quan.list$trop.ncep[ ,1]), col =rgb(0,0,0, max= 225, alpha=50), border = NA)
    # 
    # lines(ts(obs.dlm.list$trop.era$m[,1]), col='purple', lwd=1.5)
    # polygon(c(rev(x), x), c(rev(obs.quan.list$trop.era[ ,3]), obs.quan.list$trop.era[ ,1]), col =rgb(100,0,225, max= 225, alpha=50), border = NA)
  }
  #control
  {
    plot.ts(miroc.dlm.list$miroc.trop.control$m[,1], lwd=1.5, col='green',ylim=c(297.5, 300.5), main='', xaxt="n", ylab="ave temp (K)")
    axis(1, c(1, 120, 240, 360), labels = c("1981", "1991", '2001', '2011'))
    
    polygon(c(rev(x), x), c(rev(miroc.quan.list$miroc.trop.control[ ,3]), miroc.quan.list$miroc.trop.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)
    
    lines(can.dlm.list$can.trop.control$m[,1], lwd=1.5, col='green', lty=2)
    polygon(c(rev(x), x), c(rev(can.quan.list$can.trop.control[ ,3]), can.quan.list$can.trop.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)
    
    lines(had.dlm.list$had.trop.control$m[,1], lwd=1.5, col='green', lty=3)
    polygon(c(rev(x), x), c(rev(had.quan.list$had.trop.control[ ,3]), had.quan.list$had.trop.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)
    
    lines(gfdl.dlm.list$gfdl.trop.control$m[,1], lwd=1.5, col='green', lty=4)
    polygon(c(rev(x), x), c(rev(gfdl.quan.list$gfdl.trop.control[ ,3]), gfdl.quan.list$gfdl.trop.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)
    
    # lines(obs.quan.list$trop.ncep[,1], lwd=1.5, col='black')
    # polygon(c(rev(x), x), c(rev(obs.quan.list$trop.ncep[ ,3]), obs.quan.list$trop.ncep[ ,1]), col =rgb(0,0,0, max= 225, alpha=50), border = NA)
    # 
    # lines(ts(obs.dlm.list$trop.era$m[,1]), col='purple', lwd=1.5)
    # polygon(c(rev(x), x), c(rev(obs.quan.list$trop.era[ ,3]), obs.quan.list$trop.era[ ,1]), col =rgb(100,0,225, max= 225, alpha=50), border = NA)
  }
}
#Northern
{
  #Decadal
  {
    plot.ts(miroc.dlm.list$miroc.north.deca$m[,1], lwd=1.5, col='red', main='', ylim=c(286, 289), xaxt="n", ylab="ave temp (K)")
    axis(1, c(1, 120, 240, 360), labels = c("1981", "1991", '2001', '2011'))
    
    polygon(c(rev(x), x), c(rev(miroc.quan.list$miroc.north.deca[ ,3]), miroc.quan.list$miroc.north.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)
    
    lines(can.dlm.list$can.north.deca$m[,1], lwd=1.5, col='red', lty=2)
    polygon(c(rev(x), x), c(rev(can.quan.list$can.north.deca[ ,3]), can.quan.list$can.north.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)
    
    lines(had.dlm.list$had.north.deca$m[,1], lwd=1.5, col='red', lty=3)
    polygon(c(rev(x), x), c(rev(had.quan.list$had.north.deca[ ,3]), had.quan.list$had.north.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)
    
    lines(gfdl.dlm.list$gfdl.north.deca$m[,1], lwd=1.5, col='red', lty=4)
    polygon(c(rev(x), x), c(rev(gfdl.quan.list$gfdl.north.deca[ ,3]), gfdl.quan.list$gfdl.north.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)
    
    # lines(obs.quan.list$north.ncep[,1], lwd=1.5, col='black')
    # polygon(c(rev(x), x), c(rev(obs.quan.list$north.ncep[ ,3]), obs.quan.list$north.ncep[ ,1]), col =rgb(0,0,0, max= 225, alpha=50), border = NA)
    # 
    # lines(ts(obs.dlm.list$north.era$m[,1]), col='purple', lwd=1.5)
    # polygon(c(rev(x), x), c(rev(obs.quan.list$north.era[ ,3]), obs.quan.list$north.era[ ,1]), col =rgb(100,0,225, max= 225, alpha=50), border = NA)
  }
  #historical
  {
    plot.ts(miroc.dlm.list$miroc.north.hist$m[,1], lwd=1.5, col='blue', ylim=c(286, 289), main='Northern', xaxt="n", ylab="ave temp (K)")
    axis(1, c(1, 120, 240, 360), labels = c("1981", "1991", '2001', '2011'))
    
    polygon(c(rev(x), x), c(rev(miroc.quan.list$miroc.north.hist[ ,3]), miroc.quan.list$miroc.north.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)
    
    lines(can.dlm.list$can.north.hist$m[,1], lwd=1.5, col='blue', lty=2)
    polygon(c(rev(x), x), c(rev(can.quan.list$can.north.hist[ ,3]), can.quan.list$can.north.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)
    
    lines(had.dlm.list$had.north.hist$m[,1], lwd=1.5, col='blue', lty=3)
    polygon(c(rev(x), x), c(rev(had.quan.list$had.north.hist[ ,3]), had.quan.list$had.north.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)
    
    lines(gfdl.dlm.list$gfdl.north.hist$m[,1], lwd=1.5, col='blue', lty=4)
    polygon(c(rev(x), x), c(rev(gfdl.quan.list$gfdl.north.hist[ ,3]), gfdl.quan.list$gfdl.north.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)
    
    # lines(obs.quan.list$north.ncep[,1], lwd=1.5, col='black')
    # polygon(c(rev(x), x), c(rev(obs.quan.list$north.ncep[ ,3]), obs.quan.list$north.ncep[ ,1]), col =rgb(0,0,0, max= 225, alpha=50), border = NA)
    # 
    # lines(ts(obs.dlm.list$north.era$m[,1]), col='purple', lwd=1.5)
    # polygon(c(rev(x), x), c(rev(obs.quan.list$north.era[ ,3]), obs.quan.list$north.era[ ,1]), col =rgb(100,0,225, max= 225, alpha=50), border = NA)
  }
  #control
  {
    plot.ts(miroc.dlm.list$miroc.north.control$m[,1], lwd=1.5, col='green',ylim=c(286.5, 289), main='', xaxt="n", ylab="ave temp (K)")
    axis(1, c(1, 120, 240, 360), labels = c("1981", "1991", '2001', '2011'))
    
    polygon(c(rev(x), x), c(rev(miroc.quan.list$miroc.north.control[ ,3]), miroc.quan.list$miroc.north.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)
    
    lines(can.dlm.list$can.north.control$m[,1], lwd=1.5, col='green', lty=2)
    polygon(c(rev(x), x), c(rev(can.quan.list$can.north.control[ ,3]), can.quan.list$can.north.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)
    
    lines(had.dlm.list$had.north.control$m[,1], lwd=1.5, col='green', lty=3)
    polygon(c(rev(x), x), c(rev(had.quan.list$had.north.control[ ,3]), had.quan.list$had.north.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)
    
    lines(gfdl.dlm.list$gfdl.north.control$m[,1], lwd=1.5, col='green', lty=4)
    polygon(c(rev(x), x), c(rev(gfdl.quan.list$gfdl.north.control[ ,3]), gfdl.quan.list$gfdl.north.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)
    
    # lines(obs.quan.list$north.ncep[,1], lwd=1.5, col='black')
    # polygon(c(rev(x), x), c(rev(obs.quan.list$north.ncep[ ,3]), obs.quan.list$north.ncep[ ,1]), col =rgb(0,0,0, max= 225, alpha=50), border = NA)
    # 
    # lines(ts(obs.dlm.list$north.era$m[,1]), col='purple', lwd=1.5)
    # polygon(c(rev(x), x), c(rev(obs.quan.list$north.era[ ,3]), obs.quan.list$north.era[ ,1]), col =rgb(100,0,225, max= 225, alpha=50), border = NA)
  }
}
#Southern
{
  #Decadal
  {
    plot.ts(miroc.dlm.list$miroc.south.deca$m[,1], lwd=1.5, col='red', main='', ylim=c(286, 288.5), xaxt="n", ylab="ave temp (K)")
    axis(1, c(1, 120, 240, 360), labels = c("1981", "1991", '2001', '2011'))
    polygon(c(rev(x), x), c(rev(miroc.quan.list$miroc.south.deca[ ,3]), miroc.quan.list$miroc.south.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)
    
    lines(can.dlm.list$can.south.deca$m[,1], lwd=1.5, col='red', lty=2)
    polygon(c(rev(x), x), c(rev(can.quan.list$can.south.deca[ ,3]), can.quan.list$can.south.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)
    
    lines(had.dlm.list$had.south.deca$m[,1], lwd=1.5, col='red', lty=3)
    polygon(c(rev(x), x), c(rev(had.quan.list$had.south.deca[ ,3]), had.quan.list$had.south.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)
    
    lines(gfdl.dlm.list$gfdl.south.deca$m[,1], lwd=1.5, col='red', lty=4)
    polygon(c(rev(x), x), c(rev(gfdl.quan.list$gfdl.south.deca[ ,3]), gfdl.quan.list$gfdl.south.deca[ ,1]), col =rgb(225,0,0, max= 225, alpha=50), border = NA)
    
    # lines(obs.quan.list$south.ncep[,1], lwd=1.5, col='black')
    # polygon(c(rev(x), x), c(rev(obs.quan.list$south.ncep[ ,3]), obs.quan.list$south.ncep[ ,1]), col =rgb(0,0,0, max= 225, alpha=50), border = NA)
    # 
    # lines(ts(obs.dlm.list$south.era$m[,1]), col='purple', lwd=1.5)
    # polygon(c(rev(x), x), c(rev(obs.quan.list$south.era[ ,3]), obs.quan.list$south.era[ ,1]), col =rgb(100,0,225, max= 225, alpha=50), border = NA)
  }
  #historical
  {
    plot.ts(miroc.dlm.list$miroc.south.hist$m[,1], lwd=1.5, col='blue', main='Southern', xaxt="n", ylab="ave temp (K)", ylim=c(286, 288.2))
    axis(1, c(1, 120, 240, 360), labels = c("1981", "1991", '2001', '2011'))
    
    polygon(c(rev(x), x), c(rev(miroc.quan.list$miroc.south.hist[ ,3]), miroc.quan.list$miroc.south.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)
    
    lines(can.dlm.list$can.south.hist$m[,1], lwd=1.5, col='blue', lty=2)
    polygon(c(rev(x), x), c(rev(can.quan.list$can.south.hist[ ,3]), can.quan.list$can.south.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)
    
    lines(had.dlm.list$had.south.hist$m[,1], lwd=1.5, col='blue', lty=3)
    polygon(c(rev(x), x), c(rev(had.quan.list$had.south.hist[ ,3]), had.quan.list$had.south.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)
    
    lines(gfdl.dlm.list$gfdl.south.hist$m[,1], lwd=1.5, col='blue', lty=4)
    polygon(c(rev(x), x), c(rev(gfdl.quan.list$gfdl.south.hist[ ,3]), gfdl.quan.list$gfdl.south.hist[ ,1]), col =rgb(0,0,225, max= 225, alpha=50), border = NA)
    
    # lines(obs.quan.list$south.ncep[,1], lwd=1.5, col='black')
    # polygon(c(rev(x), x), c(rev(obs.quan.list$south.ncep[ ,3]), obs.quan.list$south.ncep[ ,1]), col =rgb(0,0,0, max= 225, alpha=50), border = NA)
    # 
    # lines(ts(obs.dlm.list$south.era$m[,1]), col='purple', lwd=1.5)
    # polygon(c(rev(x), x), c(rev(obs.quan.list$south.era[ ,3]), obs.quan.list$south.era[ ,1]), col =rgb(100,0,225, max= 225, alpha=50), border = NA)
  }
  #control
  {
    plot.ts(miroc.dlm.list$miroc.south.control$m[,1], lwd=1.5, col='green', main='', xaxt="n", ylab="ave temp (K)", ylim=c(285, 288.5))
    axis(1, c(1, 120, 240, 360), labels = c("1981", "1991", '2001', '2011'))
    
    polygon(c(rev(x), x), c(rev(miroc.quan.list$miroc.south.control[ ,3]), miroc.quan.list$miroc.south.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)
    
    lines(can.dlm.list$can.south.control$m[,1], lwd=1.5, col='green', lty=2)
    polygon(c(rev(x), x), c(rev(can.quan.list$can.south.control[ ,3]), can.quan.list$can.south.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)
    
    lines(had.dlm.list$had.south.control$m[,1], lwd=1.5, col='green', lty=3)
    polygon(c(rev(x), x), c(rev(had.quan.list$had.south.control[ ,3]), had.quan.list$had.south.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)
    
    lines(gfdl.dlm.list$gfdl.south.control$m[,1], lwd=1.5, col='green', lty=4)
    polygon(c(rev(x), x), c(rev(gfdl.quan.list$gfdl.south.control[ ,3]), gfdl.quan.list$gfdl.south.control[ ,1]), col =rgb(0,225,0, max= 225, alpha=50), border = NA)
    
    # lines(obs.quan.list$south.ncep[,1], lwd=1.5, col='black')
    # polygon(c(rev(x), x), c(rev(obs.quan.list$south.ncep[ ,3]), obs.quan.list$south.ncep[ ,1]), col =rgb(0,0,0, max= 225, alpha=50), border = NA)
    # 
    # lines(ts(obs.dlm.list$south.era$m[,1]), col='purple', lwd=1.5)
    # polygon(c(rev(x), x), c(rev(obs.quan.list$south.era[ ,3]), obs.quan.list$south.era[ ,1]), col =rgb(100,0,225, max= 225, alpha=50), border = NA)
  }
}

#legend
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", legend =c('MIROC5', 'CanCM4', 'HadCM3', 'GFDL'), lty = c(1,2,3,4), pt.cex=3, lwd=3, cex=1.5, bty='n')
mtext("Species", at=0.2, cex=2)