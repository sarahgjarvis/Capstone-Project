TVD.dist = function(one.spec,multi.spec){
  source("/Users/SarJar/Desktop/Climate\ Model\ Data/Spectral\ Analysis/TVD.R")
  
  tvds = rep(NA,dim(multi.spec)[1])
  for(i in 1:length(tvds)){
    tvds[i] = TVD(one.spec,multi.spec[i,])
  }
  tvd.int = c(quantile(tvds,0.025),mean(tvds),quantile(tvds,0.975))
  
  return(list(tvd.samp=tvds,tvd.int=tvd.int))
}
