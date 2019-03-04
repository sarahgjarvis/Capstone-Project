averager = function(var,lon,lat,time,region){
  # var[lon,lat,time]
  # returns global ave over time and areas used for weighting
  
  if( range(lon)[1]<0 || range(lat)[2]>90 ){ stop("lat/lon not within usual bounds") }
  library(geosphere)
  
  if(region=="global"){
    lothresh = -90
    hithresh = 90 }
  if(region=="nearglobal"){
    lothresh = -60
    hithresh = 60 }
  if(region=="tropical"){
    lothresh = -20
    hithresh = 20 }
  if(region=="upperhemisphere"){
    lothresh = 0
    hithresh = 90 }
  if(region=="lowerhemisphere"){
    lothresh = -90
    hithresh = 0 }
  
  lat.index = which(lat < hithresh & lat > lothresh)
  var = var[,lat.index,]
  lat = lat[lat.index]
  
  lon = lon + (180-max(lon))
  
  ex.lon = c(-180,lon)
  ex.lat = lat
  
  mid.lon = ex.lon[-length(ex.lon)] + diff(ex.lon)/2
  ind.dec = all(diff(ex.lat) <= 0)
  mid.lat = sort(ex.lat[-length(ex.lat)] + diff(ex.lat)/2,decreasing=TRUE)
  
  mid.lon = c(mid.lon,mid.lon[1])
  hiedge = min(hithresh - mid.lat[1],mid.lat[1]-mid.lat[2])
  loedge = min(rev(mid.lat)[1] - lothresh,rev(mid.lat)[2] - rev(mid.lat)[1])
  mid.lat = sort(c(mid.lat[1]+hiedge,mid.lat,rev(mid.lat)[1]-loedge),decreasing = ind.dec)
  
  areas = matrix(NA,length(lon),length(lat))
  
  for(i in 1:length(lon)){
    for(j in 1:length(lat)){
      A = c(mid.lon[i],mid.lat[j])
      B = c(mid.lon[(i+1)],mid.lat[j])
      C = c(mid.lon[(i+1)],mid.lat[(j+1)])
      D = c(mid.lon[i],mid.lat[(j+1)])
      areas[i,j] = areaPolygon( matrix(c(A,B,C,D,A),ncol=2,byrow=TRUE) )
    }
  }
  # library(fields)
  # image.plot(areas)
  
  total.area = sum(areas)
  
 aves = rep(NA,length(time))
  for(t in 1:length(time)){
    aves[t] = sum(areas*var[,,t])/total.area
  }
 #plot.ts(aves)
  
  return(list(aves=aves,areas=areas)) 
}
