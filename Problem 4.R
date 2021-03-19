library(MASS)
library(spatstat)
library(spatial)

set.seed(12345)

# Load data
ppregion(xl = 0, xu = 1, yl = 0, yu = 1)
window = owin(xrange = c(0, 1), yrange = c(0, 1))

cells <- read.table("cells.dat")
cells <- as.ppp(cells, window)

# We gestimate???
k = length(cells$x)
phi.1 = 
phi.2 = 
tau.0 = 

  
  
Stauss_sim = function(k, phi.1, phi.2, tau.0){
  x_vec = list(x = c(), y = c())
  
}

