source('/Users/SarJar/Desktop/Climate\ Model\ Data/Spectral\ Analysis/mv.tvar.lik.R')
source('/Users/SarJar/Desktop/Climate\ Model\ Data/Spectral\ Analysis/mv.tvar.fit.R')
source('/Users/SarJar/Desktop/Climate\ Model\ Data/Spectral\ Analysis/uni.spec.dens.intervals.R')
source('/Users/SarJar/Desktop/Climate\ Model\ Data/Spectral\ Analysis/uni.spec.dens.R')
source('/Users/SarJar/Desktop/Climate\ Model\ Data/Spectral\ Analysis/mv.tvar.spec.samp.R')
source('/Users/SarJar/Desktop/Climate\ Model\ Data/Spectral\ Analysis/TVD.dist.R')



# likelihood to select order 

# select the order q using the univariate time series of residuals for the each simulation type, each domain, and each individual replicate

o = lapply(resid.obs$south, function(x) {
  apply(x, 2, function(y) which.max(mv.tvar.lik(y, c(1,max.order), c(1,1), 0.001, 0.01)$llik))})

lapply(o, max)
  
order = c(4, 7, 5, 5)
                  


# fit AR dfs = c(1,1) make it non-time-varying
# multivariate data


s0 = 0.01
n0 = 1



col = c(rep(c('green', 'red', 'blue'), 3), 'green', 'blue')

# plot spectra

#MIROC5
# AR order (global, tropical, northern, southern)
{
par(mfrow=c(2,2), oma = c(0, 0, 2, 2))
par(mar=c(4,4,1,1))

#global
{#obs
#ncep
m0 = rep(0, order[1])
C0 = diag(1, order[1])
ncep.fit = mv.tvar.fit(resid.obs$global$ncep, order[1], c(1,1), m0, C0, s0, n0)
uni.spec.dens.intervals(ncep.fit, 'black' , FALSE, "Global", TRUE, c(-0.5, 3))
axis(1, at = c(2*pi/(12*5), 2*pi/(12*2), 2*pi/(12*1), 2*pi/(12*0.5)), labels = c(5, 2, 1, 0.5))

#era
m0 = rep(0, order[1])
C0 = diag(1, order[1])
era.fit = mv.tvar.fit(resid.obs$global$era, order[1], c(1,1), m0, C0, s0, n0)
uni.spec.dens.intervals(era.fit, 'purple' , TRUE, "", TRUE, c(-0.5, 3))

#miroc5
m0 = rep(0, order[1])
C0 = diag(1, order[1])
for(i in 1:3){
  fit = mv.tvar.fit(resid$global[[i]], order[1], c(1,1), m0, C0, s0, n0)
  uni.spec.dens.intervals(fit, col[i], TRUE, "", TRUE, c(-0.5, 3))
}

}

#trop
{#obs
  #ncep
  m0 = rep(0, order[2])
  C0 = diag(1, order[2])
  fit = mv.tvar.fit(resid.obs$trop$ncep, order[2], c(1,1), m0, C0, s0, n0)
  uni.spec.dens.intervals(fit, 'black' , FALSE, "Tropical", TRUE, c(-0.5, 3))
  axis(1, at = c(2*pi/(12*5), 2*pi/(12*2), 2*pi/(12*1), 2*pi/(12*0.5)), labels = c(5, 2, 1, 0.5))
  
  #era
  m0 = rep(0, order[2])
  C0 = diag(1, order[2])
  fit = mv.tvar.fit(resid.obs$trop$era, order[2], c(1,1), m0, C0, s0, n0)
  uni.spec.dens.intervals(fit, 'purple' , TRUE, "", TRUE, c(-0.5, 3))
  
  #miroc5
  m0 = rep(0, order[2])
  C0 = diag(1, order[2])
  for(i in 1:3){
    fit = mv.tvar.fit(resid$trop[[i]], order[2], c(1,1), m0, C0, s0, n0)
    uni.spec.dens.intervals(fit, col[i], TRUE, "", TRUE, c(-0.5, 3))
  }
  
}

#north
{#obs
  #ncep
  m0 = rep(0, order[3])
  C0 = diag(1, order[3])
  fit = mv.tvar.fit(resid.obs$north$ncep, order[3], c(1,1), m0, C0, s0, n0)
  uni.spec.dens.intervals(fit, 'black' , FALSE, "Northern", TRUE, c(-0.5, 3))
  axis(1, at = c(2*pi/(12*5), 2*pi/(12*2), 2*pi/(12*1), 2*pi/(12*0.5)), labels = c(5, 2, 1, 0.5))
  
  #era
  m0 = rep(0, order[3])
  C0 = diag(1, order[3])
  fit = mv.tvar.fit(resid.obs$north$era, order[3], c(1,1), m0, C0, s0, n0)
  uni.spec.dens.intervals(fit, 'purple' , TRUE, "", TRUE, c(-0.5, 3))
  
  #miroc5
  m0 = rep(0, order[3])
  C0 = diag(1, order[3])
  for(i in 1:3){
    fit = mv.tvar.fit(resid$north[[i]], order[3], c(1,1), m0, C0, s0, n0)
    uni.spec.dens.intervals(fit, col[i], TRUE, "", TRUE, c(-0.5, 3))
  }
  
}

#south
{#obs
  #ncep
  m0 = rep(0, order[4])
  C0 = diag(1, order[4])
  fit = mv.tvar.fit(resid.obs$south$ncep, order[4], c(1,1), m0, C0, s0, n0)
  uni.spec.dens.intervals(fit, 'black' , FALSE, "Southern", TRUE, c(-0.5, 3))
  axis(1, at = c(2*pi/(12*5), 2*pi/(12*2), 2*pi/(12*1), 2*pi/(12*0.5)), labels = c(5, 2, 1, 0.5))
  #era
  m0 = rep(0, order[4])
  C0 = diag(1, order[4])
  fit = mv.tvar.fit(resid.obs$south$era, order[4], c(1,1), m0, C0, s0, n0)
  uni.spec.dens.intervals(fit, 'purple' , TRUE, "", TRUE, c(-0.5, 3))
  
  #miroc5
  m0 = rep(0, order[4])
  C0 = diag(1, order[4])
  for(i in 1:3){
    fit = mv.tvar.fit(resid$south[[i]], order[4], c(1,1), m0, C0, s0, n0)
    uni.spec.dens.intervals(fit, col[i], TRUE, "", TRUE, c(-0.5, 3))
  }
  mtext(c('NCEP', 'ERA', 'M-Decadal', 'M-Historical', 'M-Control'), adj = c(0.43, 0.46, 0.5, 0.56, 0.615), outer = TRUE, cex = 0.8, col = c('black', 'purple', 'red', 'blue', 'green'))
}
}

 #CanCM4
{
  par(mfrow=c(2,2), oma = c(0, 0, 2, 2))
  par(mar=c(4,4,1,1))
  
  #global
  {#obs
    #ncep
    m0 = rep(0, order[1])
    C0 = diag(1, order[1])
    fit = mv.tvar.fit(resid.obs$global$ncep, order[1], c(1,1), m0, C0, s0, n0)
    uni.spec.dens.intervals(fit, 'black' , FALSE, "Global", TRUE, c(-0.5, 3))
    axis(1, at = c(2*pi/(12*5), 2*pi/(12*2), 2*pi/(12*1), 2*pi/(12*0.5)), labels = c(5, 2, 1, 0.5))
    
    #era
    fit = mv.tvar.fit(resid.obs$global$era, order[1], c(1,1), m0, C0, s0, n0)
    uni.spec.dens.intervals(fit, 'purple' , TRUE, "", TRUE, c(-0.5, 3))
    
    #can
    for(i in 4:6){
      fit = mv.tvar.fit(resid$global[[i]], order[1], c(1,1), m0, C0, s0, n0)
      uni.spec.dens.intervals(fit, col[i], TRUE, "", TRUE, c(-0.5, 3))
    }
    
  }
  
  #trop
  m0 = rep(0, order[2])
  C0 = diag(1, order[2])
  {#obs
    #ncep
    fit = mv.tvar.fit(resid.obs$trop$ncep, order[2], c(1,1), m0, C0, s0, n0)
    uni.spec.dens.intervals(fit, 'black' , FALSE, "Tropical", TRUE, c(-0.5, 3))
    axis(1, at = c(2*pi/(12*5), 2*pi/(12*2), 2*pi/(12*1), 2*pi/(12*0.5)), labels = c(5, 2, 1, 0.5))
    
    #era
    fit = mv.tvar.fit(resid.obs$trop$era, order[2], c(1,1), m0, C0, s0, n0)
    uni.spec.dens.intervals(fit, 'purple' , TRUE, "", TRUE, c(-0.5, 3))
    
    #can5
    for(i in 4:6){
      fit = mv.tvar.fit(resid$trop[[i]], order[2], c(1,1), m0, C0, s0, n0)
      uni.spec.dens.intervals(fit, col[i], TRUE, "", TRUE, c(-0.5, 3))
    }
    
  }
  
  #north
  m0 = rep(0, order[3])
  C0 = diag(1, order[3])
  {#obs
    #ncep
    fit = mv.tvar.fit(resid.obs$north$ncep, order[3], c(1,1), m0, C0, s0, n0)
    uni.spec.dens.intervals(fit, 'black' , FALSE, "Northern", TRUE, c(-0.5, 3))
    axis(1, at = c(2*pi/(12*5), 2*pi/(12*2), 2*pi/(12*1), 2*pi/(12*0.5)), labels = c(5, 2, 1, 0.5))
    
    #era
    fit = mv.tvar.fit(resid.obs$north$era, order[3], c(1,1), m0, C0, s0, n0)
    uni.spec.dens.intervals(fit, 'purple' , TRUE, "", TRUE, c(-0.5, 3))
    
    #can5
    for(i in 4:6){
      fit = mv.tvar.fit(resid$north[[i]], order[3], c(1,1), m0, C0, s0, n0)
      uni.spec.dens.intervals(fit, col[i], TRUE, "", TRUE, c(-0.5, 3))
    }
    
  }
  
  #south
  m0 = rep(0, order[4])
  C0 = diag(1, order[4])
  {#obs
    #ncep
    fit = mv.tvar.fit(resid.obs$south$ncep, order[4], c(1,1), m0, C0, s0, n0)
    uni.spec.dens.intervals(fit, 'black' , FALSE, "Southern", TRUE, c(-0.5, 3))
    axis(1, at = c(2*pi/(12*5), 2*pi/(12*2), 2*pi/(12*1), 2*pi/(12*0.5)), labels = c(5, 2, 1, 0.5))
    #era
    fit = mv.tvar.fit(resid.obs$south$era, order[4], c(1,1), m0, C0, s0, n0)
    uni.spec.dens.intervals(fit, 'purple' , TRUE, "", TRUE, c(-0.5, 3))
    
    #can5
    for(i in 4:6){
      fit = mv.tvar.fit(resid$south[[i]], order[4], c(1,1), m0, C0, s0, n0)
      uni.spec.dens.intervals(fit, col[i], TRUE, "", TRUE, c(-0.5, 3))
    }
    
  }
  mtext(c('NCEP', 'ERA', 'M-Decadal', 'M-Historical', 'M-Control'), adj = c(0.43, 0.46, 0.5, 0.56, 0.615), outer = TRUE, cex = 0.8, col = c('black', 'purple', 'red', 'blue', 'green'))
}

#hadCM4
{
  par(mfrow=c(2,2), oma = c(0, 0, 2, 2))
  par(mar=c(4,4,1,1))
  
  #global
  m0 = rep(0, order[1])
  C0 = diag(1, order[1])
  {#obs
    #ncep
    fit = mv.tvar.fit(resid.obs$global$ncep, order[1], c(1,1), m0, C0, s0, n0)
    uni.spec.dens.intervals(fit, 'black' , FALSE, "Global", TRUE, c(-0.5, 3))
    axis(1, at = c(2*pi/(12*5), 2*pi/(12*2), 2*pi/(12*1), 2*pi/(12*0.5)), labels = c(5, 2, 1, 0.5))
    
    #era
    fit = mv.tvar.fit(resid.obs$global$era, order[1], c(1,1), m0, C0, s0, n0)
    uni.spec.dens.intervals(fit, 'purple' , TRUE, "", TRUE, c(-0.5, 3))
    
    #had
    for(i in 7:9){
      fit = mv.tvar.fit(resid$global[[i]], order[1], c(1,1), m0, C0, s0, n0)
      uni.spec.dens.intervals(fit, col[i], TRUE, "", TRUE, c(-0.5, 3))
    }
    
  }
  
  #trop
  m0 = rep(0, order[2])
  C0 = diag(1, order[2])
  {#obs
    #ncep
    fit = mv.tvar.fit(resid.obs$trop$ncep, order[2], c(1,1), m0, C0, s0, n0)
    uni.spec.dens.intervals(fit, 'black' , FALSE, "Tropical", TRUE, c(-0.5, 3))
    axis(1, at = c(2*pi/(12*5), 2*pi/(12*2), 2*pi/(12*1), 2*pi/(12*0.5)), labels = c(5, 2, 1, 0.5))
    
    #era
    fit = mv.tvar.fit(resid.obs$trop$era, order[2], c(1,1), m0, C0, s0, n0)
    uni.spec.dens.intervals(fit, 'purple' , TRUE, "", TRUE, c(-0.5, 3))
    
    #had5
    for(i in 7:9){
      fit = mv.tvar.fit(resid$trop[[i]], order[2], c(1,1), m0, C0, s0, n0)
      uni.spec.dens.intervals(fit, col[i], TRUE, "", TRUE, c(-0.5, 3))
    }
    
  }
  
  #north
  m0 = rep(0, order[3])
  C0 = diag(1, order[3])
  {#obs
    #ncep
    fit = mv.tvar.fit(resid.obs$north$ncep, order[3], c(1,1), m0, C0, s0, n0)
    uni.spec.dens.intervals(fit, 'black' , FALSE, "Northern", TRUE, c(-0.5, 3))
    axis(1, at = c(2*pi/(12*5), 2*pi/(12*2), 2*pi/(12*1), 2*pi/(12*0.5)), labels = c(5, 2, 1, 0.5))
    
    #era
    fit = mv.tvar.fit(resid.obs$north$era, order[3], c(1,1), m0, C0, s0, n0)
    uni.spec.dens.intervals(fit, 'purple' , TRUE, "", TRUE, c(-0.5, 3))
    
    #had5
    for(i in 7:9){
      fit = mv.tvar.fit(resid$north[[i]], order[3], c(1,1), m0, C0, s0, n0)
      uni.spec.dens.intervals(fit, col[i], TRUE, "", TRUE, c(-0.5, 3))
    }
    
  }
  
  #south
  m0 = rep(0, order[4])
  C0 = diag(1, order[4])
  {#obs
    #ncep
    fit = mv.tvar.fit(resid.obs$south$ncep, order[4], c(1,1), m0, C0, s0, n0)
    uni.spec.dens.intervals(fit, 'black' , FALSE, "Southern", TRUE, c(-0.5, 3))
    axis(1, at = c(2*pi/(12*5), 2*pi/(12*2), 2*pi/(12*1), 2*pi/(12*0.5)), labels = c(5, 2, 1, 0.5))
    #era
    fit = mv.tvar.fit(resid.obs$south$era, order[4], c(1,1), m0, C0, s0, n0)
    uni.spec.dens.intervals(fit, 'purple' , TRUE, "", TRUE, c(-0.5, 3))
    
    #had5
    for(i in 7:9){
      fit = mv.tvar.fit(resid$south[[i]], order[4], c(1,1), m0, C0, s0, n0)
      uni.spec.dens.intervals(fit, col[i], TRUE, "", TRUE, c(-0.5, 3))
    }
    
  }
  mtext(c('NCEP', 'ERA', 'M-Decadal', 'M-Historical', 'M-Control'), adj = c(0.43, 0.46, 0.5, 0.56, 0.615), outer = TRUE, cex = 0.8, col = c('black', 'purple', 'red', 'blue', 'green'))
}

#GFDL
{
  par(mfrow=c(2,2), oma = c(0, 0, 2, 2))
  par(mar=c(4,4,1,1))
  
  #global
  m0 = rep(0, order[1])
  C0 = diag(1, order[1])
  {#obs
    #ncep
    fit = mv.tvar.fit(resid.obs$global$ncep, order[1], c(1,1), m0, C0, s0, n0)
    uni.spec.dens.intervals(fit, 'black' , FALSE, "Global", TRUE, c(-0.5, 3))
    axis(1, at = c(2*pi/(12*5), 2*pi/(12*2), 2*pi/(12*1), 2*pi/(12*0.5)), labels = c(5, 2, 1, 0.5))
    
    #era
    fit = mv.tvar.fit(resid.obs$global$era, order[1], c(1,1), m0, C0, s0, n0)
    uni.spec.dens.intervals(fit, 'purple' , TRUE, "", TRUE, c(-0.5, 3))
    
    #had
    for(i in 10:11){
      fit = mv.tvar.fit(resid$global[[i]], order[1], c(1,1), m0, C0, s0, n0)
      uni.spec.dens.intervals(fit, col[i], TRUE, "", TRUE, c(-0.5, 3))
    }
    
  }
  
  #trop
  m0 = rep(0, order[2])
  C0 = diag(1, order[2])
  {#obs
    #ncep
    fit = mv.tvar.fit(resid.obs$trop$ncep, order[2], c(1,1), m0, C0, s0, n0)
    uni.spec.dens.intervals(fit, 'black' , FALSE, "Tropical", TRUE, c(-0.5, 3))
    axis(1, at = c(2*pi/(12*5), 2*pi/(12*2), 2*pi/(12*1), 2*pi/(12*0.5)), labels = c(5, 2, 1, 0.5))
    
    #era
    fit = mv.tvar.fit(resid.obs$trop$era, order[2], c(1,1), m0, C0, s0, n0)
    uni.spec.dens.intervals(fit, 'purple' , TRUE, "", TRUE, c(-0.5, 3))
    
    #had5
    for(i in 10:11){
      fit = mv.tvar.fit(resid$trop[[i]], order[2], c(1,1), m0, C0, s0, n0)
      uni.spec.dens.intervals(fit, col[i], TRUE, "", TRUE, c(-0.5, 3))
    }
    
  }
  
  #north
  m0 = rep(0, order[3])
  C0 = diag(1, order[3])
  {#obs
    #ncep
    fit = mv.tvar.fit(resid.obs$north$ncep, order[3], c(1,1), m0, C0, s0, n0)
    uni.spec.dens.intervals(fit, 'black' , FALSE, "Northern", TRUE, c(-0.5, 3))
    axis(1, at = c(2*pi/(12*5), 2*pi/(12*2), 2*pi/(12*1), 2*pi/(12*0.5)), labels = c(5, 2, 1, 0.5))
    
    #era
    fit = mv.tvar.fit(resid.obs$north$era, order[3], c(1,1), m0, C0, s0, n0)
    uni.spec.dens.intervals(fit, 'purple' , TRUE, "", TRUE, c(-0.5, 3))
    
    #had5
    for(i in 10:11){
      fit = mv.tvar.fit(resid$north[[i]], order[3], c(1,1), m0, C0, s0, n0)
      uni.spec.dens.intervals(fit, col[i], TRUE, "", TRUE, c(-0.5, 3))
    }
    
  }
  
  #south
  m0 = rep(0, order[4])
  C0 = diag(1, order[4])
  {#obs
    #ncep
    fit = mv.tvar.fit(resid.obs$south$ncep, order[4], c(1,1), m0, C0, s0, n0)
    uni.spec.dens.intervals(fit, 'black' , FALSE, "Southern", TRUE, c(-0.5, 3))
    axis(1, at = c(2*pi/(12*5), 2*pi/(12*2), 2*pi/(12*1), 2*pi/(12*0.5)), labels = c(5, 2, 1, 0.5))
    #era
    fit = mv.tvar.fit(resid.obs$south$era, order[4], c(1,1), m0, C0, s0, n0)
    uni.spec.dens.intervals(fit, 'purple' , TRUE, "", TRUE, c(-0.5, 3))
    
    #had5
    for(i in 10:11){
      fit = mv.tvar.fit(resid$south[[i]], order[4], c(1,1), m0, C0, s0, n0)
      uni.spec.dens.intervals(fit, col[i], TRUE, "", TRUE, c(-0.5, 3))
    }
    
  }
  mtext(c('NCEP', 'ERA', 'M-Decadal', 'M-Historical', 'M-Control'), adj = c(0.43, 0.46, 0.5, 0.56, 0.615), outer = TRUE, cex = 0.8, col = c('black', 'purple', 'red', 'blue', 'green'))
}



#compare models
par(mar=c(4,4,1,1))
par(mfrow=c(4,3), oma = c(0, 0, 2, 0))
{
#Global
#decadal
{  m0 = rep(0, order[1])
  C0 = diag(1, order[1])
#ncep
fit = mv.tvar.fit(resid.obs$global$ncep, order[1], c(1,1), m0, C0, s0, n0)
uni.spec.dens.intervals(fit, 'black' , FALSE, "", TRUE, c(-0.5, 3), linty = 1)
axis(1, at = c(2*pi/(12*5), 2*pi/(12*2), 2*pi/(12*1), 2*pi/(12*0.5)), labels = c(5, 2, 1, 0.5))
#era
fit = mv.tvar.fit(resid.obs$global$era, order[1], c(1,1), m0, C0, s0, n0)
uni.spec.dens.intervals(fit, 'purple' , TRUE, "", TRUE, c(-0.5, 3), linty = 1)
#models
l=0
for (i in 1:length(resid$global)){
  if(grepl('deca',names(resid$global[i]))){
    l=l+1
    fit = mv.tvar.fit(resid$global[[i]], order[1], c(1,1), m0, C0, s0, n0)
    uni.spec.dens.intervals(fit, 'red', TRUE, "", TRUE, c(-0.5, 3), linty = l)}
}
}
#historical
{#ncep
fit = mv.tvar.fit(resid.obs$global$ncep, order[1], c(1,1), m0, C0, s0, n0)
uni.spec.dens.intervals(fit, 'black' , FALSE, "Global", TRUE, c(-0.5, 3), linty = 1)
axis(1, at = c(2*pi/(12*5), 2*pi/(12*2), 2*pi/(12*1), 2*pi/(12*0.5)), labels = c(5, 2, 1, 0.5))
#era
fit = mv.tvar.fit(resid.obs$global$era, order[1], c(1,1), m0, C0, s0, n0)
uni.spec.dens.intervals(fit, 'purple' , TRUE, "", TRUE, c(-0.5, 3), linty = 1)
#models
l=0
for (i in 1:length(resid$global)){
  if(grepl('hist',names(resid$global[i]))){
    l=l+1
    fit = mv.tvar.fit(resid$global[[i]], order[1], c(1,1), m0, C0, s0, n0)
    uni.spec.dens.intervals(fit, 'blue', TRUE, "", TRUE, c(-0.5, 3), linty = l)}
}
}
#control
{#ncep
fit = mv.tvar.fit(resid.obs$global$ncep, order[1], c(1,1), m0, C0, s0, n0)
uni.spec.dens.intervals(fit, 'black' , FALSE, "", TRUE, c(-0.5, 3), linty = 1)
axis(1, at = c(2*pi/(12*5), 2*pi/(12*2), 2*pi/(12*1), 2*pi/(12*0.5)), labels = c(5, 2, 1, 0.5))
#era
fit = mv.tvar.fit(resid.obs$global$era, order[1], c(1,1), m0, C0, s0, n0)
uni.spec.dens.intervals(fit, 'purple' , TRUE, "", TRUE, c(-0.5, 3), linty = 1)
#models
l=0
for (i in 1:length(resid$global)){
  if(grepl('control',names(resid$global[i]))){
    l=l+1
  fit = mv.tvar.fit(resid$global[[i]], order[1], c(1,1), m0, C0, s0, n0)
  uni.spec.dens.intervals(fit, 'green', TRUE, "", TRUE, c(-0.5, 3), linty = l)}
}
}

#trop
#decadal
{
  m0 = rep(0, order[2])
  C0 = diag(1, order[2])
  #ncep
  fit = mv.tvar.fit(resid.obs$trop$ncep, order[2], c(1,1), m0, C0, s0, n0)
  uni.spec.dens.intervals(fit, 'black' , FALSE, "", TRUE, c(-0.5, 3), linty = 1)
  axis(1, at = c(2*pi/(12*5), 2*pi/(12*2), 2*pi/(12*1), 2*pi/(12*0.5)), labels = c(5, 2, 1, 0.5))
  #era
  fit = mv.tvar.fit(resid.obs$trop$era, order[2], c(1,1), m0, C0, s0, n0)
  uni.spec.dens.intervals(fit, 'purple' , TRUE, "", TRUE, c(-0.5, 3), linty = 1)
  #models
  l=0
  for (i in 1:length(resid$trop)){
    if(grepl('deca',names(resid$trop[i]))){
      l=l+1
      fit = mv.tvar.fit(resid$trop[[i]], order[2], c(1,1), m0, C0, s0, n0)
      uni.spec.dens.intervals(fit, 'red', TRUE, "", TRUE, c(-0.5, 3), linty = l)}
  }
}
#historical
{#ncep
  fit = mv.tvar.fit(resid.obs$trop$ncep, order[2], c(1,1), m0, C0, s0, n0)
  uni.spec.dens.intervals(fit, 'black' , FALSE, "Tropical", TRUE, c(-0.5, 3), linty = 1)
  axis(1, at = c(2*pi/(12*5), 2*pi/(12*2), 2*pi/(12*1), 2*pi/(12*0.5)), labels = c(5, 2, 1, 0.5))
  #era
  fit = mv.tvar.fit(resid.obs$trop$era, order[2], c(1,1), m0, C0, s0, n0)
  uni.spec.dens.intervals(fit, 'purple' , TRUE, "", TRUE, c(-0.5, 3), linty = 1)
  #models
  l=0
  for (i in 1:length(resid$trop)){
    if(grepl('hist',names(resid$trop[i]))){
      l=l+1
      fit = mv.tvar.fit(resid$trop[[i]], order[2], c(1,1), m0, C0, s0, n0)
      uni.spec.dens.intervals(fit, 'blue', TRUE, "", TRUE, c(-0.5, 3), linty = l)}
  }
}
#control
{#ncep
  fit = mv.tvar.fit(resid.obs$trop$ncep, order[2], c(1,1), m0, C0, s0, n0)
  uni.spec.dens.intervals(fit, 'black' , FALSE, "", TRUE, c(-0.5, 3), linty = 1)
  axis(1, at = c(2*pi/(12*5), 2*pi/(12*2), 2*pi/(12*1), 2*pi/(12*0.5)), labels = c(5, 2, 1, 0.5))
  #era
  fit = mv.tvar.fit(resid.obs$trop$era, order[2], c(1,1), m0, C0, s0, n0)
  uni.spec.dens.intervals(fit, 'purple' , TRUE, "", TRUE, c(-0.5, 3), linty = 1)
  #models
  for (i in 1:length(resid$trop)){
    if(grepl('control',names(resid$trop[i]))){
      l=l+1
      fit = mv.tvar.fit(resid$trop[[i]], order[2], c(1,1), m0, C0, s0, n0)
      uni.spec.dens.intervals(fit, 'green', TRUE, "", TRUE, c(-0.5, 3), linty = l)}
  }
}

#north
{#decadal
{ m0 = rep(0, order[3])
  C0 = diag(1, order[3])
  #ncep
  fit = mv.tvar.fit(resid.obs$north$ncep, order[3], c(1,1), m0, C0, s0, n0)
  uni.spec.dens.intervals(fit, 'black' , FALSE, "", TRUE, c(-0.5, 3), linty = 1)
  axis(1, at = c(2*pi/(12*5), 2*pi/(12*2), 2*pi/(12*1), 2*pi/(12*0.5)), labels = c(5, 2, 1, 0.5))
  #era
  fit = mv.tvar.fit(resid.obs$north$era, order[3], c(1,1), m0, C0, s0, n0)
  uni.spec.dens.intervals(fit, 'purple' , TRUE, "", TRUE, c(-0.5, 3), linty = 1)
  #models
  l=0
  for (i in 1:length(resid$north)){
    if(grepl('deca',names(resid$north[i]))){
      l=l+1
      fit = mv.tvar.fit(resid$north[[i]], order[3], c(1,1), m0, C0, s0, n0)
      uni.spec.dens.intervals(fit, 'red', TRUE, "", TRUE, c(-0.5, 3), linty = l)}
  }
}
#historical
{#ncep
  fit = mv.tvar.fit(resid.obs$north$ncep, order[3], c(1,1), m0, C0, s0, n0)
  uni.spec.dens.intervals(fit, 'black' , FALSE, "Northern", TRUE, c(-0.5, 3), linty = 1)
  axis(1, at = c(2*pi/(12*5), 2*pi/(12*2), 2*pi/(12*1), 2*pi/(12*0.5)), labels = c(5, 2, 1, 0.5))
  #era
  fit = mv.tvar.fit(resid.obs$north$era, order[3], c(1,1), m0, C0, s0, n0)
  uni.spec.dens.intervals(fit, 'purple' , TRUE, "", TRUE, c(-0.5, 3), linty = 1)
  #models
  l=0
  for (i in 1:length(resid$north)){
    if(grepl('hist',names(resid$north[i]))){
      l=l+1
      fit = mv.tvar.fit(resid$north[[i]], order[3], c(1,1), m0, C0, s0, n0)
      uni.spec.dens.intervals(fit, 'blue', TRUE, "", TRUE, c(-0.5, 3), linty = l)}
  }
}
#control
{#ncep
  fit = mv.tvar.fit(resid.obs$north$ncep, order[3], c(1,1), m0, C0, s0, n0)
  uni.spec.dens.intervals(fit, 'black' , FALSE, "", TRUE, c(-0.5, 3), linty = 1)
  axis(1, at = c(2*pi/(12*5), 2*pi/(12*2), 2*pi/(12*1), 2*pi/(12*0.5)), labels = c(5, 2, 1, 0.5))
  #era
  fit = mv.tvar.fit(resid.obs$north$era, order[3], c(1,1), m0, C0, s0, n0)
  uni.spec.dens.intervals(fit, 'purple' , TRUE, "", TRUE, c(-0.5, 3), linty = 1)
  #models
  l=0
  for (i in 1:length(resid$north)){
    if(grepl('control',names(resid$north[i]))){
      l=l+1
      fit = mv.tvar.fit(resid$north[[i]], order[3], c(1,1), m0, C0, s0, n0)
      uni.spec.dens.intervals(fit, 'green', TRUE, "", TRUE, c(-0.5, 3), linty = l)}
  }
}}

#south
{#decadal
{  m0 = rep(0, order[4])
  C0 = diag(1, order[4])
  #ncep
  fit = mv.tvar.fit(resid.obs$south$ncep, order[4], c(1,1), m0, C0, s0, n0)
  uni.spec.dens.intervals(fit, 'black' , FALSE, "", TRUE, c(-0.5, 3), linty = 1)
  axis(1, at = c(2*pi/(12*5), 2*pi/(12*2), 2*pi/(12*1), 2*pi/(12*0.5)), labels = c(5, 2, 1, 0.5))
  #era
  fit = mv.tvar.fit(resid.obs$south$era, order[4], c(1,1), m0, C0, s0, n0)
  uni.spec.dens.intervals(fit, 'purple' , TRUE, "", TRUE, c(-0.5, 3), linty = 1)
  #models
  l=0
  for (i in 1:length(resid$south)){
    if(grepl('deca',names(resid$south[i]))){
      l=l+1
      fit = mv.tvar.fit(resid$south[[i]], order[4], c(1,1), m0, C0, s0, n0)
      uni.spec.dens.intervals(fit, 'red', TRUE, "", TRUE, c(-0.5, 3), linty = l)}
  }
}
#historical
{#ncep
  fit = mv.tvar.fit(resid.obs$south$ncep, order[4], c(1,1), m0, C0, s0, n0)
  uni.spec.dens.intervals(fit, 'black' , FALSE, "Southern", TRUE, c(-0.5, 3), linty = 1)
  axis(1, at = c(2*pi/(12*5), 2*pi/(12*2), 2*pi/(12*1), 2*pi/(12*0.5)), labels = c(5, 2, 1, 0.5))
  #era
  fit = mv.tvar.fit(resid.obs$south$era, order[4], c(1,1), m0, C0, s0, n0)
  uni.spec.dens.intervals(fit, 'purple' , TRUE, "", TRUE, c(-0.5, 3), linty = 1)
  #models
  l=0
  for (i in 1:length(resid$south)){
    if(grepl('hist',names(resid$south[i]))){
      l=l+1
      fit = mv.tvar.fit(resid$south[[i]], order[4], c(1,1), m0, C0, s0, n0)
      uni.spec.dens.intervals(fit, 'blue', TRUE, "", TRUE, c(-0.5, 3), linty = l)}
  }
}
#control
{#ncep
  fit = mv.tvar.fit(resid.obs$south$ncep, order[4], c(1,1), m0, C0, s0, n0)
  uni.spec.dens.intervals(fit, 'black' , FALSE, "", TRUE, c(-0.5, 3), linty = 1)
  axis(1, at = c(2*pi/(12*5), 2*pi/(12*2), 2*pi/(12*1), 2*pi/(12*0.5)), labels = c(5, 2, 1, 0.5))
  #era
  fit = mv.tvar.fit(resid.obs$south$era, order[4], c(1,1), m0, C0, s0, n0)
  uni.spec.dens.intervals(fit, 'purple' , TRUE, "", TRUE, c(-0.5, 3), linty = 1)
  #models
  l=0
  for (i in 1:length(resid$south)){
    if(grepl('control',names(resid$south[i]))){
      l=l+1
      fit = mv.tvar.fit(resid$south[[i]], order[4], c(1,1), m0, C0, s0, n0)
      uni.spec.dens.intervals(fit, 'green', TRUE, "", TRUE, c(-0.5, 3), linty = l)}
  }
}}
mtext(c('NCEP', 'ERA', 'M-Decadal', 'M-Historical', 'M-Control'), adj = c(0.43, 0.46, 0.5, 0.56, 0.615), outer = TRUE, cex = 0.8, col = c('black', 'purple', 'red', 'blue', 'green'))
}



#make list of autoregressive fit
#observatinal

fit.obs = list()
# Global
m0 = rep(0, order[1])
C0 = diag(1, order[1])
fit.obs$global = lapply(resid.obs$global, function(x) {mv.tvar.fit(x, order[1], c(1,1), m0, C0, s0, n0)})
# Tropical
m0 = rep(0, order[2])
C0 = diag(1, order[2])
fit.obs$trop = lapply(resid.obs$trop, function(x) {mv.tvar.fit(x, order[2], c(1,1), m0, C0, s0, n0)})

# North/South
m0 = rep(0, order[3])
C0 = diag(1, order[3])
fit.obs$north = lapply(resid.obs$north, function(x) {mv.tvar.fit(x, order[3], c(1,1), m0, C0, s0, n0)})
fit.obs$south = lapply(resid.obs$south, function(x) {mv.tvar.fit(x, order[4], c(1,1), m0, C0, s0, n0)})





fit = list()
for(i in 1:length(resid$global)){
  m0 = rep(0, order[1])
  C0 = diag(1, order[1])
  fit$global[[i]] = mv.tvar.fit(resid$global[[i]], order[1], c(1,1), m0, C0, s0, n0)
}

for(i in 1:length(resid$trop)){
  m0 = rep(0, order[2])
  C0 = diag(1, order[2])
  fit$trop[[i]] = mv.tvar.fit(resid$trop[[i]], order[2], c(1,1), m0, C0, s0, n0)
}

for(i in 1:length(resid$north)){
  m0 = rep(0, order[3])
  C0 = diag(1, order[3])
  fit$north[[i]] = mv.tvar.fit(resid$north[[i]], order[3], c(1,1), m0, C0, s0, n0)
}

for(i in 1:length(resid$south)){
  m0 = rep(0, order[4])
  C0 = diag(1, order[4])
  fit$south[[i]] = mv.tvar.fit(resid$south[[i]], order[4], c(1,1), m0, C0, s0, n0)
}

names(fit$global) = names(resid$global)
names(fit$trop) = names(resid$trop)
names(fit$north) = names(resid$north)
names(fit$south) = names(resid$south)


# white noise tvd comparison
lim = floor(T/2)
# example
# ncep.spec = uni.spec.dens(ncep.fit,lim,wo.wn=TRUE,TRUE)$spec
# wn.spec = rep(1,length(ncep.spec))
# ncep.spec.samp = mv.tvar.spec.samp(ncep.fit,lim,wo.wn=TRUE)$spec.samps
# ncep.wn.tvd = TVD.dist(wn.spec,ncep.spec.samp)$tvd.samp
# 


spec.obs = list()
spec.obs$global = lapply(fit.obs$global, function(x){uni.spec.dens(x,lim,wo.wn=TRUE,TRUE)$spec})
spec.obs$trop = lapply(fit.obs$trop, function(x){uni.spec.dens(x,lim,wo.wn=TRUE,TRUE)$spec})
spec.obs$north = lapply(fit.obs$north, function(x){uni.spec.dens(x,lim,wo.wn=TRUE,TRUE)$spec})
spec.obs$south = lapply(fit.obs$south, function(x){uni.spec.dens(x,lim,wo.wn=TRUE,TRUE)$spec})

spec = list()
spec$global = lapply(fit$global, function(x){uni.spec.dens(x,lim,wo.wn=TRUE,TRUE)$spec})
spec$trop = lapply(fit$trop, function(x){uni.spec.dens(x,lim,wo.wn=TRUE,TRUE)$spec})
spec$north = lapply(fit$north, function(x){uni.spec.dens(x,lim,wo.wn=TRUE,TRUE)$spec})
spec$south = lapply(fit$south, function(x){uni.spec.dens(x,lim,wo.wn=TRUE,TRUE)$spec})

wn.spec = rep(1,length(ncep.spec))

spec.samp.obs = list()
spec.samp.obs$global = lapply(fit.obs$global, function(x){mv.tvar.spec.samp(x,lim,wo.wn=TRUE)$spec.samps})
spec.samp.obs$trop = lapply(fit.obs$trop, function(x){mv.tvar.spec.samp(x,lim,wo.wn=TRUE)$spec.samps})
spec.samp.obs$north = lapply(fit.obs$north, function(x){mv.tvar.spec.samp(x,lim,wo.wn=TRUE)$spec.samps})
spec.samp.obs$south = lapply(fit.obs$south, function(x){mv.tvar.spec.samp(x,lim,wo.wn=TRUE)$spec.samps})

spec.samp = list()
spec.samp$global = lapply(fit$global, function(x){mv.tvar.spec.samp(x,lim,wo.wn=TRUE)$spec.samps})
spec.samp$trop = lapply(fit$trop, function(x){mv.tvar.spec.samp(x,lim,wo.wn=TRUE)$spec.samps})
spec.samp$north = lapply(fit$north, function(x){mv.tvar.spec.samp(x,lim,wo.wn=TRUE)$spec.samps})
spec.samp$south = lapply(fit$south, function(x){mv.tvar.spec.samp(x,lim,wo.wn=TRUE)$spec.samps})

tvd.wn.obs = list()
tvd.wn.obs$global = lapply(spec.samp.obs$global, function(x) {TVD.dist(wn.spec, x)$tvd.samp})
tvd.wn.obs$trop = lapply(spec.samp.obs$trop, function(x) {TVD.dist(wn.spec, x)$tvd.samp})
tvd.wn.obs$north = lapply(spec.samp.obs$north, function(x) {TVD.dist(wn.spec, x)$tvd.samp})
tvd.wn.obs$south = lapply(spec.samp.obs$south, function(x) {TVD.dist(wn.spec, x)$tvd.samp})


tvd.wn = list()
tvd.wn$global = lapply(spec.samp$global, function(x) {TVD.dist(wn.spec, x)$tvd.samp})
tvd.wn$trop = lapply(spec.samp$trop, function(x) {TVD.dist(wn.spec, x)$tvd.samp})
tvd.wn$north = lapply(spec.samp$north, function(x) {TVD.dist(wn.spec, x)$tvd.samp})
tvd.wn$south = lapply(spec.samp$south, function(x) {TVD.dist(wn.spec, x)$tvd.samp})

#miroc plot
{par(mfrow=c(1,4), oma = c(0, 0, 2, 0))
par(mar=c(3,3,1,1))
#global
boxplot(c(q.fn(tvd.wn.obs$global), q.fn(tvd.wn$global[2]), q.fn(tvd.wn$global[3]), q.fn(tvd.wn$global[1])), col = c('black', 'purple', 'red','blue', 'green'), names = c('NCEP', 'ERA', 'dec', 'hist', 'cont'))
#trop
boxplot(c(q.fn(tvd.wn.obs$trop), q.fn(tvd.wn$trop[2]), q.fn(tvd.wn$trop[3]), q.fn(tvd.wn$trop[1])), col = c('black', 'purple', 'red','blue', 'green'), names = c('NCEP', 'ERA', 'dec', 'hist', 'cont'))
#north
boxplot(c(q.fn(tvd.wn.obs$north), q.fn(tvd.wn$north[2]), q.fn(tvd.wn$north[3]), q.fn(tvd.wn$north[1])), col = c('black', 'purple', 'red','blue', 'green'), names = c('NCEP', 'ERA', 'dec', 'hist', 'cont'))
#south
boxplot(c(q.fn(tvd.wn.obs$south), q.fn(tvd.wn$south[2]), q.fn(tvd.wn$south[3]), q.fn(tvd.wn$south[1])), col = c('black', 'purple', 'red','blue', 'green'), names = c('NCEP', 'ERA', 'dec', 'hist', 'cont'))
mtext(c('NCEP', 'ERA', 'Decadal', 'Historical', 'Control'), adj = c(0.443, 0.470, 0.50, 0.543, 0.585), outer = TRUE, cex = 0.8, col = c('black', 'purple', 'red', 'blue', 'green'))}
#can plot
{par(mfrow=c(1,4), oma = c(0, 0, 2, 0))
par(mar=c(3,3,1,1))
#global
boxplot(c(q.fn(tvd.wn.obs$global), q.fn(tvd.wn$global[5]), q.fn(tvd.wn$global[6]), q.fn(tvd.wn$global[4])), col = c('black', 'purple', 'red','blue', 'green'), names = c('NCEP', 'ERA', 'dec', 'hist', 'cont'))
#trop
boxplot(c(q.fn(tvd.wn.obs$trop), q.fn(tvd.wn$trop[5]), q.fn(tvd.wn$trop[6]), q.fn(tvd.wn$trop[4])), col = c('black', 'purple', 'red','blue', 'green'), names = c('NCEP', 'ERA', 'dec', 'hist', 'cont'))
#north
boxplot(c(q.fn(tvd.wn.obs$north), q.fn(tvd.wn$north[5]), q.fn(tvd.wn$north[6]), q.fn(tvd.wn$north[4])), col = c('black', 'purple', 'red','blue', 'green'), names = c('NCEP', 'ERA', 'dec', 'hist', 'cont'))
#south
boxplot(c(q.fn(tvd.wn.obs$south), q.fn(tvd.wn$south[5]), q.fn(tvd.wn$south[6]), q.fn(tvd.wn$south[4])), col = c('black', 'purple', 'red','blue', 'green'), names = c('NCEP', 'ERA', 'dec', 'hist', 'cont'))
mtext(c('NCEP', 'ERA', 'Decadal', 'Historical', 'Control'), adj = c(0.443, 0.470, 0.50, 0.543, 0.585), outer = TRUE, cex = 0.8, col = c('black', 'purple', 'red', 'blue', 'green'))}
#had plot
{par(mfrow=c(1,4), oma = c(0, 0, 2, 0))
par(mar=c(3,3,1,1))
#global
boxplot(c(q.fn(tvd.wn.obs$global), q.fn(tvd.wn$global[8]), q.fn(tvd.wn$global[9]), q.fn(tvd.wn$global[7])), col = c('black', 'purple', 'red','blue', 'green'), names = c('NCEP', 'ERA', 'dec', 'hist', 'cont'))
#trop
boxplot(c(q.fn(tvd.wn.obs$trop), q.fn(tvd.wn$trop[8]), q.fn(tvd.wn$trop[9]), q.fn(tvd.wn$trop[7])), col = c('black', 'purple', 'red','blue', 'green'), names = c('NCEP', 'ERA', 'dec', 'hist', 'cont'))
#north
boxplot(c(q.fn(tvd.wn.obs$north), q.fn(tvd.wn$north[8]), q.fn(tvd.wn$north[9]), q.fn(tvd.wn$north[7])), col = c('black', 'purple', 'red','blue', 'green'), names = c('NCEP', 'ERA', 'dec', 'hist', 'cont'))
#south
boxplot(c(q.fn(tvd.wn.obs$south), q.fn(tvd.wn$south[8]), q.fn(tvd.wn$south[9]), q.fn(tvd.wn$south[7])), col = c('black', 'purple', 'red','blue', 'green'), names = c('NCEP', 'ERA', 'dec', 'hist', 'cont'))
mtext(c('NCEP', 'ERA', 'Decadal', 'Historical', 'Control'), adj = c(0.443, 0.470, 0.50, 0.543, 0.585), outer = TRUE, cex = 0.8, col = c('black', 'purple', 'red', 'blue', 'green'))}
#gfdl plot
{par(mfrow=c(1,4), oma = c(0, 0, 2, 0))
par(mar=c(3,3,1,1))
#global
boxplot(c(q.fn(tvd.wn.obs$global), q.fn(tvd.wn$global[11]), q.fn(tvd.wn$global[10])), col = c('black', 'purple','blue', 'green'), names = c('NCEP', 'ERA', 'hist', 'cont'))
#trop
boxplot(c(q.fn(tvd.wn.obs$trop), q.fn(tvd.wn$trop[11]), q.fn(tvd.wn$trop[10])), col = c('black', 'purple','blue', 'green'), names = c('NCEP', 'ERA', 'hist', 'cont'))
#north
boxplot(c(q.fn(tvd.wn.obs$north), q.fn(tvd.wn$north[11]), q.fn(tvd.wn$north[10])), col = c('black', 'purple','blue', 'green'), names = c('NCEP', 'ERA', 'hist', 'cont'))
#south
boxplot(c(q.fn(tvd.wn.obs$south), q.fn(tvd.wn$south[11]), q.fn(tvd.wn$south[10])), col = c('black', 'purple','blue', 'green'), names = c('NCEP', 'ERA', 'hist', 'cont'))
mtext(c('NCEP', 'ERA', 'Decadal', 'Historical', 'Control'), adj = c(0.443, 0.470, 0.50, 0.543, 0.585), outer = TRUE, cex = 0.8, col = c('black', 'purple', 'red', 'blue', 'green'))}

#plot to compare models
which(grepl('deca', names(tvd.wn$global)))
which(grepl('hist', names(tvd.wn$global)))
which(grepl('cont', names(tvd.wn$global)))


par(mfrow=c(4,3), oma = c(0, 0, 2, 0))
par(mar=c(3,3,1,1))

#global
{boxplot(c(q.fn(tvd.wn.obs$global), q.fn(tvd.wn$global[2]), q.fn(tvd.wn$global[5]), q.fn(tvd.wn$global[8])), col = c('black', 'purple', 'red', 'red', 'red'), ylim = c(0.2, 0.9), names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3'))
boxplot(c(q.fn(tvd.wn.obs$global), q.fn(tvd.wn$global[3]), q.fn(tvd.wn$global[6]), q.fn(tvd.wn$global[9]), q.fn(tvd.wn$global[11])), col = c('black', 'purple', rep('blue', 4)), main = 'Global', ylim = c(0.2, 0.9), names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3', 'GFDL'))
boxplot(c(q.fn(tvd.wn.obs$global), q.fn(tvd.wn$global[1]), q.fn(tvd.wn$global[4]), q.fn(tvd.wn$global[7]), q.fn(tvd.wn$global[10])), col = c('black', 'purple', rep('green', 4)), ylim = c(0.2, 0.9), names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3', 'GFDL'))}
#tropical
{boxplot(c(q.fn(tvd.wn.obs$trop), q.fn(tvd.wn$trop[2]), q.fn(tvd.wn$trop[5]), q.fn(tvd.wn$trop[8])), col = c('black', 'purple', 'red', 'red', 'red'), ylim = c(0.5, 0.9), names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3'))
boxplot(c(q.fn(tvd.wn.obs$trop), q.fn(tvd.wn$trop[3]), q.fn(tvd.wn$trop[6]), q.fn(tvd.wn$trop[9]), q.fn(tvd.wn$trop[11])), col = c('black', 'purple', rep('blue', 4)), main = 'Tropical', ylim = c(0.5, 0.9), names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3', 'GFDL'))
boxplot(c(q.fn(tvd.wn.obs$trop), q.fn(tvd.wn$trop[1]), q.fn(tvd.wn$trop[4]), q.fn(tvd.wn$trop[7]), q.fn(tvd.wn$trop[10])), col = c('black', 'purple', rep('green', 4)), ylim = c(0.5, 0.9), names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3', 'GFDL'))}
#north
{boxplot(c(q.fn(tvd.wn.obs$north), q.fn(tvd.wn$north[2]), q.fn(tvd.wn$north[5]), q.fn(tvd.wn$north[8])), col = c('black', 'purple', 'red', 'red', 'red'), ylim = c(0.2, 0.9), names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3'))
boxplot(c(q.fn(tvd.wn.obs$north), q.fn(tvd.wn$north[3]), q.fn(tvd.wn$north[6]), q.fn(tvd.wn$north[9]), q.fn(tvd.wn$north[11])), col = c('black', 'purple', rep('blue', 4)), main = 'Northern', ylim = c(0.2, 0.9), names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3', 'GFDL'))
boxplot(c(q.fn(tvd.wn.obs$north), q.fn(tvd.wn$north[1]), q.fn(tvd.wn$north[4]), q.fn(tvd.wn$north[7]), q.fn(tvd.wn$north[10])), col = c('black', 'purple', rep('green', 4)), ylim = c(0.2, 0.9), names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3', 'GFDL'))}
#south
{boxplot(c(q.fn(tvd.wn.obs$south), q.fn(tvd.wn$south[2]), q.fn(tvd.wn$south[5]), q.fn(tvd.wn$south[8])), col = c('black', 'purple', 'red', 'red', 'red'), ylim = c(0.2, 0.9), names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3'))
boxplot(c(q.fn(tvd.wn.obs$south), q.fn(tvd.wn$south[3]), q.fn(tvd.wn$south[6]), q.fn(tvd.wn$south[9]), q.fn(tvd.wn$south[11])), col = c('black', 'purple', rep('blue', 4)), main = 'Southern', ylim = c(0.2, 0.9), names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3', 'GFDL'))
boxplot(c(q.fn(tvd.wn.obs$south), q.fn(tvd.wn$south[1]), q.fn(tvd.wn$south[4]), q.fn(tvd.wn$south[7]), q.fn(tvd.wn$south[10])), col = c('black', 'purple', rep('green', 4)), ylim = c(0.2, 0.9), names = c('NCEP', 'ERA', 'MIROC5', 'CanCM4', 'HadCM3', 'GFDL'))}
mtext(c('NCEP', 'ERA', 'Decadal', 'Historical', 'Control'), adj = c(0.443, 0.470, 0.50, 0.543, 0.585), outer = TRUE, cex = 0.8, col = c('black', 'purple', 'red', 'blue', 'green'))

lapply(fit$south, function(x){mean(x$s)})






