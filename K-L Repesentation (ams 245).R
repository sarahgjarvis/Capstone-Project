# Use the K-L representation to approximate the exponential correlation for range parameter equal to 1. Plot the approximation for several orders and compare to the actual correlation.

library(rootSolve)

# function to combine odd and even elements of a vector
combine = function(x, y){
  stopifnot(length(x) == length(y))
  res = rep(NA, length(x))
  for (i in 1:length(x)) {
    res[2*i - 1] = x[i] # odd
    res[2*i] = y[i] #even
  }
  return(res)
}

KL.exponential = function(s1, s2, J=5, L=1, phi=1){
  #uniroot.all searches for several roots
  w1 = uniroot.all(function(w) tan(w*L) - phi/w, interval = c(0,L), n=10000)
  w2 = uniroot.all(function(w) tan(w*L) + w/phi, interval = c(0,L), n=10000)
 
  
  # to get finite solutions
  w1 = head(sort(w1[which(abs(tan(w1/L)) < Inf)]), J)
  w2 = head(sort(w2[which(abs(tan(w2/L)) < Inf)]), J)
  
  # stop of w1 and w2 are different lengths
  stopifnot(length(w1) == length(w2))
  
  #create odd/even lambda
  # lambda1 = (2*phi) / (w1^2 + phi^2) #odd
  # lambda2 = (2*phi) / (w2^2 + phi^2) #even
  lambda1 = (2*phi) / (w1 + phi)^2 #odd
  lambda2 = (2*phi) / (w2 + phi)^2 #even
  #combine lambda
  lambda = combine(lambda1, lambda2)

  #create odd/even psi
  psi1 = function(s) {cos(w1*s) / sqrt(L + sin(2*w1*L)/(2*w1+.Machine$double.eps))}
  psi2 = function(s) {sin(w2*s) / sqrt(L - sin(2*w2*L)/(2*w2+.Machine$double.eps))}
  #combine psi
  psi = function(s) {combine(psi1(s), psi2(s))}
  
  #return the covariance function C(s,s')
  cov = sum(lambda * psi(s1) * Conj(psi(s2))) 
  return(list(cov=cov, 'lambda'=lambda, 'psi1'=psi1, 'psi2'=psi2, 'psi'=psi))
}

L=3

#distances
d = seq(0, 4, length = 100)

actual <- cov_fn$powered.exp(tau = d, phi=1, sig2=1, nu=1)
approx <- sapply(d, function(s) KL.exponential(0, s, J=2, L=L)$cov)
approx2 <- sapply(d, function(s) KL.exponential(0, s, J=4, L=L)$cov)
approx3 <- sapply(d, function(s) KL.exponential(0, s, J=6, L=L)$cov)

plot(d, actual, type = 'l', xlab = 'distance' ~ (tau), ylab = 'correlation', main = 'K-L Approximate of Exponential Correlation', ylim = c(-0.2, 1), lwd = 2)
lines(d, approx, lwd = 3, lty = 2, col = 5)
lines(d, approx2, lwd = 3, lty = 3, col = 6)
lines(d, approx3, lwd = 3, lty = 4, col = 4)
legend('topright', legend = c('actual', 'order 2 approx', 'order 4 approx', 'order 6 approx'), col = c(1,5,6,4), lwd = c(2,3,3,3), bty = 'n')


#Generate 100 realizations of a univariate Gaussian process with exponential correlation with range parameter 1. Compare the empirically estimated eigenvalues and eigenfunctions to the ones given by the K-L and the approximation on Page 12.
library(geoR)
range = seq(-4, 4, length = 100)

#function grf generates simulations of Gaussian random fields
gp = grf(n=length(range), grid = cbind(range, 0), cov.model = 'matern', kappa = 0.5 , cov.pars = c(1,1), nsim = 100)$data

plot(range, gp[,1], type = 'l', ylim = c(-4,4))
for(i in 2:ncol(gp)){
  lines(range, gp[,i], col = i, lty = i)
}

#get lambda
l = sapply(1, function(s) KL.exponential(0, s, J=50, L=12)$lambda)
p = sapply(1, function(s) KL.exponential(0, s, J=50, L=12)$psi)
p1 = sapply(1, function(s) KL.exponential(0, s, J=50, L=12)$psi1)
p2 = sapply(1, function(s) KL.exponential(0, s, J=50, L=12)$psi2)



#find eigen values with SVD
E = svd(gp)
J = 50
E$d = E$d[1:J]^2
psi1 = psi(range, 1:J)

#spectral
f = function(k) 1/(1 + k^2)
lamda2 = f(1:100 * pi / (2 * L))


plot(E$d/max(E$d), pch=20, main = 'Eigenvalues', ylab = 'scaled' ~ lambda)
points(l/max(l), pch = 20, col=20)
points(lamda2/max(lamda2), pch=20, col=30)
legend('topright', legend = c('empirical', 'K-L', 'spectral'), col = c(1,20, 30), pch = 20, bty='n')

#eigen functions
range = seq(-50, 50, length = 100)

psi1 = psi(range, 1:J)
get.psi2 = function(s, j){
  out = matrix(0, length(s), j)
  for (i in 1:length(s))
    out[i,] = cos(1:j * pi * s[i] / (2*L))
  return (out)
}

psi2 = get.psi2(range, J)

try = E$u[,1:J]
for (i in 1:6){
  plot(range, psi1[,i], ylim = c(-1, 1), type='l', bty = 'n', main = paste0("Eigenvector"),
       xlab = "", ylab = expression(psi))
  lines(range, psi2[,i], col = '13')
  lines(range, try[1:100,i], col = '30')
}
