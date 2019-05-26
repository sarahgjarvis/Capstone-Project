TVD = function(spec.a,spec.b){
  
  spec.a = spec.a/sum(spec.a)
  spec.b = spec.b/sum(spec.b)
  
  tvd = sum( abs(spec.a - spec.b) )/ 2 

  return(tvd)  
}
