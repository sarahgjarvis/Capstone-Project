#saves averaged data sets and makes time series plots
#observational-Global
{global.ncep <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/Observational\ Data/global_ncep2_obs_jan81dec10.dat"))
global.era <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/Observational\ Data/global_erainterm_t2m_jan81dec10.dat"))}
#observational-Tropical
{trop.ncep <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/Observational\ Data/tropical_ncep2_obs_jan81dec10.dat"))
  trop.era <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/Observational\ Data/tropical_erainterm_t2m_jan81dec10.dat"))}
#observational-Northern
{north.ncep <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/Observational\ Data/upperhemisphere_ncep2_obs_jan81dec10.dat"))
  north.era <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/Observational\ Data/upperhemisphere_erainterm_t2m_jan81dec10.dat"))}
#observational-southern
{south.ncep <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/Observational\ Data/lowerhemisphere_ncep2_obs_jan81dec10.dat"))
  south.era <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/Observational\ Data/lowerhemisphere_erainterm_t2m_jan81dec10.dat"))}
obs.ts = list()
obs.ts$global = list(global.ncep, global.era)
obs.ts$trop = list(trop.ncep, trop.era)
obs.ts$north = list(north.ncep, north.era)
obs.ts$south = list(south.ncep, south.era)
names(obs.ts$global) = names(obs.ts$trop) = names(obs.ts$north) = names(obs.ts$south) = c('ncep', 'era')

#model-Global

#MIROC5
{
  miroc.global.deca <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/MIROC5/global_MIROC5_decadal_jan81dec10_i1.dat"))
  miroc.global.hist <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/MIROC5/global_MIROC5_historical_jan81dec10_i1.dat"))
  miroc.global.control <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/MIROC5/global_MIROC5_control_30years.dat"))
}
#CanCM4
{can.global.deca <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/CanCM4/global_CanCM4_decadal_jan1981dec2010_i1.dat"))
c.global.hist <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/CanCM4/global_CanCM4_historical_jan1981dec2005_i1.dat"))
can.global.control <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/CanCM4/global_CanESM2_control_30years.dat"))
can.global.rcp <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/CanCM4/global_CanCM4_rcp_jan2006dec2035_i1.dat"))
#add rcp onto CanCM4
dim(c.global.hist)
dim(can.global.rcp)
can.global.hist <- rbind(c.global.hist, can.global.rcp[1:60,])
#jump?
plot.ts(can.global.hist[,4])
}
#HadCM3
{had.global.deca <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/HadCM3/global_HadCM3_decadal_nov1980dec2010_i2.dat"))[3:362,]
had.global.hist <- rbind(as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/HadCM3/global_HadCM3_historical_jan1981nov1984_i1.dat")), as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/HadCM3/global_HadCM3_historical_dec1984dec2005_i1.dat")), as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/HadCM3/global_HadCM3_rcp_jan2006dec2011_i1.dat")))
had.global.control <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/HadCM3/global_CanESM2_control_30years.dat"))}
#GFDL
{gfdl.global.deca <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/GFDL-CM2p1/global_GFDL-CM2p1_decadal_jan1981dec2010_i1.dat"))
gfdl.global.hist <- rbind(as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/GFDL-CM2p1/global_GFDL-CM2p1_historical_jan1981dec2005_i1.dat")), as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/GFDL-CM2p1/global_HadCM3_rcp_jan2006dec2010_i1.dat")))
gfdl.global.control <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/GFDL-CM2p1/global_CanESM2_control_30years.dat"))}
global = list(miroc.global.control, miroc.global.deca, miroc.global.hist, can.global.control, can.global.deca, can.global.hist, had.global.control, had.global.deca, had.global.hist, gfdl.global.control, gfdl.global.hist)
names(global) = c('miroc.control', 'miroc.deca', 'miroc.hist', 'can.control', 'can.deca', 'can.hist','had.control', 'had.deca', 'had.hist', 'gfdl.control', 'gfld.hist')

#model-Tropical
#MIROC5
{
  miroc.trop.deca <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/MIROC5/tropical_MIROC5_decadal_jan81dec10_i1.dat"))
  miroc.trop.hist <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/MIROC5/tropical_MIROC5_historical_jan81dec10_i1.dat"))
  miroc.trop.control <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/MIROC5/tropical_MIROC5_control_30years.dat"))
}
#CanCM4
{can.trop.deca <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/CanCM4/tropical_CanCM4_decadal_jan1981dec2010_i1.dat"))
c.trop.hist <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/CanCM4/tropical_CanCM4_historical_jan1981dec2005_i1.dat"))
can.trop.control <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/CanCM4/tropical_CanESM2_control_30years.dat"))
can.trop.rcp <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/CanCM4/tropical_CanCM4_rcp_jan2006dec2035_i1.dat"))
dim(c.trop.hist)
dim(can.trop.rcp)
can.trop.hist <- rbind(c.trop.hist, can.trop.rcp[1:60,])
#jump?
plot.ts(can.trop.hist[,1])
}
#HadCM3
{had.trop.deca <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/HadCM3/tropical_HadCM3_decadal_nov1980dec2010_i2.dat"))[3:362,]
had.trop.hist <- rbind(as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/HadCM3/tropical_HadCM3_historical_jan1981nov1984_i1.dat")), as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/HadCM3/tropical_HadCM3_historical_dec1984dec2005_i1.dat")), as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/HadCM3/tropical_HadCM3_rcp_jan2006dec2011_i1.dat")))
had.trop.control <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/HadCM3/tropical_CanESM2_control_30years.dat"))
}
#GFDL
{
  gfdl.trop.deca <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/GFDL-CM2p1/tropical_GFDL-CM2p1_decadal_jan1981dec2010_i1.dat"))
  gfdl.trop.hist <- rbind(as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/GFDL-CM2p1/tropical_GFDL-CM2p1_historical_jan1981dec2005_i1.dat")), as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/GFDL-CM2p1/tropical_HadCM3_rcp_jan2006dec2010_i1.dat")))
  gfdl.trop.control <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/GFDL-CM2p1/tropical_CanESM2_control_30years.dat"))
}
trop = list(miroc.trop.control, miroc.trop.deca, miroc.trop.hist, can.trop.control, can.trop.deca, can.trop.hist, had.trop.control, had.trop.deca, had.trop.hist, gfdl.trop.control, gfdl.trop.hist)
names(trop) = c('miroc.control', 'miroc.deca', 'miroc.hist', 'can.control', 'can.deca', 'can.hist','had.control', 'had.deca', 'had.hist', 'gfdl.control', 'gfld.hist')

#model-Northern
#MIROC5
{
  miroc.north.deca <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/MIROC5/upperhemisphere_MIROC5_decadal_jan81dec10_i1.dat"))
  miroc.north.hist <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/MIROC5/upperhemisphere_MIROC5_historical_jan81dec10_i1.dat"))
  miroc.north.control <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/MIROC5/upperhemisphere_MIROC5_control_30years.dat"))
}
#CanCM4
{can.north.deca <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/CanCM4/upperhemisphere_CanCM4_decadal_jan61dec90_i1.dat"))
  c.north.hist <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/CanCM4/upperhemisphere_CanCM4_historical_jan1981dec2005_i1.dat"))
  can.north.control <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/CanCM4/upperhemisphere_CanESM2_control_30years.dat"))
  can.north.rcp <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/CanCM4/upperhemisphere_CanCM4_rcp_jan2006dec2035_i1.dat"))
  dim(c.north.hist)
  dim(can.north.rcp)
  can.north.hist <- rbind(c.north.hist, can.north.rcp[1:60,])
  plot.ts(can.north.hist[,1])  
}
#HadCM3
{had.north.deca <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/HadCM3/upperhemisphere_HadCM3_decadal_nov1980dec2010_i2.dat"))[3:362,]
had.north.hist <- rbind(as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/HadCM3/upperhemisphere_HadCM3_historical_jan1981nov1984_i1.dat")), as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/HadCM3/upperhemisphere_HadCM3_historical_dec1984dec2005_i1.dat")), as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/HadCM3/upperhemisphere_HadCM3_rcp_jan2006dec2011_i1.dat")))
had.north.control <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/HadCM3/upperhemisphere_CanESM2_control_30years.dat"))
}
#GFDL
{
  gfdl.north.deca <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/GFDL-CM2p1/upperhemisphere_GFDL-CM2p1_decadal_jan1981dec2010_i1.dat"))
  gfdl.north.hist <- rbind(as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/GFDL-CM2p1/upperhemisphere_GFDL-CM2p1_historical_jan1981dec2005_i1.dat")), as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/GFDL-CM2p1/upperhemisphere_HadCM3_rcp_jan2006dec2010_i1.dat")))
  gfdl.north.control <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/GFDL-CM2p1/upperhemisphere_CanESM2_control_30years.dat"))
  }
north = list(miroc.north.control, miroc.north.deca, miroc.north.hist, can.north.control, can.north.deca, can.north.hist, had.north.control, had.north.deca, had.north.hist, gfdl.north.control, gfdl.north.hist)
names(north) = c('miroc.control', 'miroc.deca', 'miroc.hist', 'can.control', 'can.deca', 'can.hist','had.control', 'had.deca', 'had.hist', 'gfdl.control', 'gfld.hist')

#model-southern
#MIROC5
{
  miroc.south.deca <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/MIROC5/lowerhemisphere_MIROC5_decadal_jan81dec10_i1.dat"))
  miroc.south.hist <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/MIROC5/lowerhemisphere_MIROC5_historical_jan81dec10_i1.dat"))
  miroc.south.control <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/MIROC5/lowerhemisphere_MIROC5_control_30years.dat"))
}
#CanCM4
{can.south.deca <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/CanCM4/lowerhemisphere_CanCM4_decadal_jan1981dec2010_i1.dat"))
  c.south.hist <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/CanCM4/lowerhemisphere_CanCM4_historical_jan1981dec2005_i1.dat"))
  can.south.control <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/CanCM4/lowerhemisphere_CanESM2_control_30years.dat"))
  can.south.rcp <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/CanCM4/lowerhemisphere_CanCM4_rcp_jan2006dec2035_i1.dat"))
  
  dim(c.south.hist)
  dim(can.south.rcp)
  can.south.hist <- rbind(c.south.hist, can.south.rcp[1:60,])
  plot.ts(can.south.hist[,1])
  }
#HadCM3
{had.south.deca <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/HadCM3/lowerhemisphere_HadCM3_decadal_nov1980dec2010_i2.dat"))[3:362,]
  had.south.hist <- rbind(as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/HadCM3/lowerhemisphere_HadCM3_historical_jan1981nov1984_i1.dat")), as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/HadCM3/lowerhemisphere_HadCM3_historical_dec1984dec2005_i1.dat")), as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/HadCM3/lowerhemisphere_HadCM3_rcp_jan2006dec2011_i1.dat")))
  had.south.control <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/HadCM3/lowerhemisphere_CanESM2_control_30years.dat"))
}
#GFDL
{
  gfdl.south.deca <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/GFDL-CM2p1/lowerhemisphere_GFDL-CM2p1_decadal_jan1981dec2010_i1.dat"))
  gfdl.south.hist <- rbind(as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/GFDL-CM2p1/lowerhemisphere_GFDL-CM2p1_historical_jan1981dec2005_i1.dat")), as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/GFDL-CM2p1/lowerhemisphere_HadCM3_rcp_jan2006dec2010_i1.dat")))
  gfdl.south.control <- as.matrix(read.table("/Users/SarJar/Desktop/Climate\ Model\ Data/GFDL-CM2p1/lowerhemisphere_CanESM2_control_30years.dat"))
}
south = list(miroc.south.control, miroc.south.deca, miroc.south.hist, can.south.control, can.south.deca, can.south.hist, had.south.control, had.south.deca, had.south.hist, gfdl.south.control, gfdl.south.hist)
names(south) = c('miroc.control', 'miroc.deca', 'miroc.hist', 'can.control', 'can.deca', 'can.hist','had.control', 'had.deca', 'had.hist', 'gfdl.control', 'gfld.hist')


model.ts = list(global, trop, north, south)
names(model.ts) = c('global', 'trop', 'north', 'south')

par(mar=c(2,4,1,1))
par(mfrow=c(4,1))
#MIROC5 Plot
{#CanCM4 Global
  {plot.ts(global.ncep[1:120,], main = "Global", ylab='ave temp (K)', ylim=c(min(global.era), max(global.ncep)), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    lines(global.era[1:120,], col='purple')
    
    plot.ts(global.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(miroc.global.deca[1:120,]),max(miroc.global.deca[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(miroc.global.deca[1:120,i]), lty=i, col='red')
    }
    
    plot.ts(global.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(miroc.global.deca[1:120,]),max(miroc.global.hist[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(miroc.global.hist[1:120,i]), lty=i, col='blue')
    }
    
    plot.ts(global.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(global.ncep[1:120,]),max(miroc.global.control[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(miroc.global.control[1:120,i]), lty=i, col='green')
    }
  }
  
  #mirocCM4 Tropical
  {plot.ts(trop.ncep[1:120,], main = "Tropical", ylab='ave temp (K)', ylim=c(min(trop.era), max(trop.ncep)), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    lines(trop.era[1:120,], col='purple')
    
    plot.ts(trop.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(miroc.trop.deca[1:120,]),max(miroc.trop.deca[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(miroc.trop.deca[1:120,i]), lty=i, col='red')
    }
    
    plot.ts(trop.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(miroc.trop.deca[1:120,]),max(miroc.trop.hist[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(miroc.trop.hist[1:120,i]), lty=i, col='blue')
    }
    
    plot.ts(trop.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(miroc.trop.control[1:120,]),max(miroc.trop.control[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(miroc.trop.control[1:120,i]), lty=i, col='green')
    }
  }
  
  #mirocCM4 Northern
  {plot.ts(north.ncep[1:120,], main = "Northern", ylab='ave temp (K)', ylim=c(min(north.era), max(north.ncep)), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    lines(north.era[1:120,], col='purple')
    
    plot.ts(north.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(miroc.north.deca[1:120,]),max(miroc.north.deca[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(miroc.north.deca[1:120,i]), lty=i, col='red')
    }
    
    plot.ts(north.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(miroc.north.deca[1:120,]),max(miroc.north.hist[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(miroc.north.hist[1:120,i]), lty=i, col='blue')
    }
    
    plot.ts(north.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(miroc.north.control[1:120,]),max(miroc.north.control[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(miroc.north.control[1:120,i]), lty=i, col='green')
    }
  }
  
  #mirocCM4 Southern
  {plot.ts(south.ncep[1:120,], main = "Southern", ylab='ave temp (K)', ylim=c(min(south.era), max(south.ncep)), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    lines(south.era[1:120,], col='purple')
    
    plot.ts(south.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(south.ncep[1:120,]),max(miroc.south.deca[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(miroc.south.deca[1:120,i]), lty=i, col='red')
    }
    
    plot.ts(south.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(south.ncep[1:120,]),max(miroc.south.hist[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(miroc.south.hist[1:120,i]), lty=i, col='blue')
    }
    
    plot.ts(south.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(south.ncep[1:120,]),max(miroc.south.control[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(miroc.south.control[1:120,i]), lty=i, col='green')
    }
  }
}
#CanCM4 Plot
{#CanCM4 Global
{plot.ts(global.ncep[1:120,], main = "Global", ylab='ave temp (K)', ylim=c(min(global.era), max(global.ncep)), xaxt="n")
  axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
  lines(global.era[1:120,], col='purple')
  
  plot.ts(global.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(can.global.deca[1:120,]),max(can.global.deca[1:120,])), xaxt="n")
  axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
  for (i in 1:3){
    lines(ts(can.global.deca[1:120,i]), lty=i, col='red')
  }
  
  plot.ts(global.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(can.global.deca[1:120,]),max(can.global.hist[1:120,])), xaxt="n")
  axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
  for (i in 1:3){
    lines(ts(can.global.hist[1:120,i]), lty=i, col='blue')
  }
  
  plot.ts(global.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(can.global.control[1:120,]),max(can.global.control[1:120,])), xaxt="n")
  axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
  for (i in 1:3){
    lines(ts(can.global.control[1:120,i]), lty=i, col='green')
  }
}

#CanCM4 Tropical
{plot.ts(trop.ncep[1:120,], main = "Tropical", ylab='ave temp (K)', ylim=c(min(trop.era), max(trop.ncep)), xaxt="n")
  axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
lines(trop.era[1:120,], col='purple')

plot.ts(trop.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(can.trop.deca[1:120,]),max(can.trop.deca[1:120,])), xaxt="n")
axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
for (i in 1:3){
  lines(ts(can.trop.deca[1:120,i]), lty=i, col='red')
}

plot.ts(trop.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(can.trop.deca[1:120,]),max(can.trop.hist[1:120,])), xaxt="n")
axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
for (i in 1:3){
  lines(ts(can.trop.hist[1:120,i]), lty=i, col='blue')
}

plot.ts(trop.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(can.trop.control[1:120,]),max(can.trop.control[1:120,])), xaxt="n")
axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
for (i in 1:3){
  lines(ts(can.trop.control[1:120,i]), lty=i, col='green')
}
}

#CanCM4 Northern
{plot.ts(north.ncep[1:120,], main = "Northern", ylab='ave temp (K)', ylim=c(min(north.era), max(north.ncep)), xaxt="n")
  axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
  lines(north.era[1:120,], col='purple')
  
  plot.ts(north.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(can.north.deca[1:120,]),max(can.north.deca[1:120,])), xaxt="n")
  axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
  for (i in 1:3){
    lines(ts(can.north.deca[1:120,i]), lty=i, col='red')
  }
  
  plot.ts(north.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(can.north.deca[1:120,]),max(can.north.hist[1:120,])), xaxt="n")
  axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
  for (i in 1:3){
    lines(ts(can.north.hist[1:120,i]), lty=i, col='blue')
  }
  
  plot.ts(north.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(can.north.control[1:120,]),max(can.north.control[1:120,])), xaxt="n")
  axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
  for (i in 1:3){
    lines(ts(can.north.control[1:120,i]), lty=i, col='green')
  }
}

#CanCM4 Southern
{plot.ts(south.ncep[1:120,], main = "Southern", ylab='ave temp (K)', ylim=c(min(south.era), max(south.ncep)), xaxt="n")
  axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
  lines(south.era[1:120,], col='purple')
  
  plot.ts(south.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(can.south.deca[1:120,]),max(can.south.deca[1:120,])), xaxt="n")
  axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
  for (i in 1:3){
    lines(ts(can.south.deca[1:120,i]), lty=i, col='red')
  }
  
  plot.ts(south.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(can.south.deca[1:120,]),max(can.south.hist[1:120,])), xaxt="n")
  axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
  for (i in 1:3){
    lines(ts(can.south.hist[1:120,i]), lty=i, col='blue')
  }
  
  plot.ts(south.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(can.south.control[1:120,]),max(can.south.control[1:120,])), xaxt="n")
  axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
  for (i in 1:3){
    lines(ts(can.south.control[1:120,i]), lty=i, col='green')
  }
}}
#HadCM3 Plot
{#HadCM3 Global
  {plot.ts(global.ncep[1:120,], main = "Global", ylab='ave temp (K)', ylim=c(min(global.era), max(global.ncep)), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    lines(global.era[1:120,], col='purple')
    
    plot.ts(global.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(had.global.deca[1:120,]),max(had.global.control[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(had.global.deca[1:120,i]), lty=i, col='red')
    }
    
    plot.ts(global.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(had.global.deca[1:120,]),max(had.global.hist[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(had.global.hist[1:120,i]), lty=i, col='blue')
    }
    
    plot.ts(global.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(had.global.control[1:120,]),max(had.global.control[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(had.global.control[1:120,i]), lty=i, col='green')
    }
  }
  
  #HadCM3 Tropical
  {plot.ts(trop.ncep[1:120,], main = "Tropical", ylab='ave temp (K)', ylim=c(min(trop.era), max(trop.ncep)), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    lines(trop.era[1:120,], col='purple')
    
    plot.ts(trop.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(had.trop.deca[1:120,]),max(had.trop.deca[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(had.trop.deca[1:120,i]), lty=i, col='red')
    }
    
    plot.ts(trop.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(had.trop.deca[1:120,]),max(had.trop.hist[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(had.trop.hist[1:120,i]), lty=i, col='blue')
    }
    
    plot.ts(trop.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(had.trop.control[1:120,]),max(had.trop.control[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(had.trop.control[1:120,i]), lty=i, col='green')
    }
  }
  
  #HadCM3 Northern
  {plot.ts(north.ncep[1:120,], main = "Northern", ylab='ave temp (K)', ylim=c(min(north.era), max(north.ncep)), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    lines(north.era[1:120,], col='purple')
    
    plot.ts(north.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(had.north.deca[1:120,]),max(had.north.deca[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(had.north.deca[1:120,i]), lty=i, col='red')
    }
    
    plot.ts(north.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(had.north.deca[1:120,]),max(had.north.hist[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(had.north.hist[1:120,i]), lty=i, col='blue')
    }
    
    plot.ts(north.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(had.north.control[1:120,]),max(had.north.control[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(had.north.control[1:120,i]), lty=i, col='green')
    }
  }
  
  #HadCM3 Southern
  {plot.ts(south.ncep[1:120,], main = "Southern", ylab='ave temp (K)', ylim=c(min(south.era), max(south.ncep)), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    lines(south.era[1:120,], col='purple')
    
    plot.ts(south.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(had.south.deca[1:120,]),max(had.south.deca[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(had.south.deca[1:120,i]), lty=i, col='red')
    }
    
    plot.ts(south.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(had.south.deca[1:120,]),max(had.south.hist[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(had.south.hist[1:120,i]), lty=i, col='blue')
    }
    
    plot.ts(south.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(had.south.control[1:120,]),max(had.south.control[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(had.south.control[1:120,i]), lty=i, col='green')
    }
  }}
#GFDL Plot
{#gfdl Global
  {plot.ts(global.ncep[1:120,], main = "Global", ylab='ave temp (K)', ylim=c(min(global.era), max(global.ncep)), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    lines(global.era[1:120,], col='purple')
    
    plot.ts(global.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(gfdl.global.deca[1:120,]),max(gfdl.global.control[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(gfdl.global.deca[1:120,i]), lty=i, col='red')
    }
    
    plot.ts(global.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(gfdl.global.deca[1:120,]),max(gfdl.global.hist[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(gfdl.global.hist[1:120,i]), lty=i, col='blue')
    }
    
    plot.ts(global.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(gfdl.global.control[1:120,]),max(gfdl.global.control[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(gfdl.global.control[1:120,i]), lty=i, col='green')
    }
  }
  
  #gfdl Tropical
  {plot.ts(trop.ncep[1:120,], main = "Tropical", ylab='ave temp (K)', ylim=c(min(trop.era), max(trop.ncep)), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    lines(trop.era[1:120,], col='purple')
    
    plot.ts(trop.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(gfdl.trop.deca[1:120,]),max(gfdl.trop.deca[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(gfdl.trop.deca[1:120,i]), lty=i, col='red')
    }
    
    plot.ts(trop.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(gfdl.trop.deca[1:120,]),max(gfdl.trop.hist[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(gfdl.trop.hist[1:120,i]), lty=i, col='blue')
    }
    
    plot.ts(trop.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(gfdl.trop.control[1:120,]),max(gfdl.trop.control[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(gfdl.trop.control[1:120,i]), lty=i, col='green')
    }
  }
  
  #gfdl Northern
  {plot.ts(north.ncep[1:120,], main = "Northern", ylab='ave temp (K)', ylim=c(min(north.era), max(north.ncep)), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    lines(north.era[1:120,], col='purple')
    
    plot.ts(north.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(gfdl.north.deca[1:120,]),max(gfdl.north.deca[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(gfdl.north.deca[1:120,i]), lty=i, col='red')
    }
    
    plot.ts(north.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(gfdl.north.deca[1:120,]),max(gfdl.north.hist[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(gfdl.north.hist[1:120,i]), lty=i, col='blue')
    }
    
    plot.ts(north.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(gfdl.north.control[1:120,]),max(gfdl.north.control[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(gfdl.north.control[1:120,i]), lty=i, col='green')
    }
  }
  
  #gfdl Southern
  {plot.ts(south.ncep[1:120,], main = "Southern", ylab='ave temp (K)', ylim=c(min(south.era), max(south.ncep)), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    lines(south.era[1:120,], col='purple')
    
    plot.ts(south.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(gfdl.south.deca[1:120,]),max(gfdl.south.deca[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(gfdl.south.deca[1:120,i]), lty=i, col='red')
    }
    
    plot.ts(south.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(gfdl.south.deca[1:120,]),max(gfdl.south.hist[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(gfdl.south.hist[1:120,i]), lty=i, col='blue')
    }
    
    plot.ts(south.ncep[1:120,], ylab='ave temp (K)', ylim=c(min(gfdl.south.control[1:120,]),max(gfdl.south.control[1:120,])), xaxt="n")
    axis(1, c(1, 60, 120), labels = c("1981", "1986", "1991"))
    for (i in 1:3){
      lines(ts(gfdl.south.control[1:120,i]), lty=i, col='green')
    }
  }}

