library(MASS)
library(spatstat)
library(spatial)
library(parallel)

set.seed(12345)
num.cores = 4

# Load data
ppregion(xl = 0, xu = 1, yl = 0, yu = 1)
window = owin(xrange = c(0, 1), yrange = c(0, 1))

cells <- read.table("cells.dat")
cells <- as.ppp(cells, window)

# We gestimate???
k = length(cells$x)
tau.0 = 0.1
phi.0 = 0.1
phi.1 = 0.5

  
phi <- function(tau, tau.0, phi.0, phi.1) {
  if (tau <= tau.0) {
    return(phi.0)
  } else {
    return(phi.0*exp(-phi.1 * (tau - tau.0)))
  }
}

pdf <- function(x, y, tau.0, phi.0, phi.1) {
  res <- 1
  k <- length(x)
  for (i in 1:k) {
    for (j in 1:k) {
      tau <- norm(c(x[i], y[i]) - c(x[j], y[j]), type = "2")
      res = res * exp(-phi(tau, tau.0, phi.0, phi.1))
    }
  }
  return(res)
}
  
Strauss_sim = function(k, tau.0, phi.0, phi.1, n.iter = 1000){
  num_accept = 0
  x_vec = list(x = runif(k), y = runif(k))
  while (pdf(x_vec$x, x_vec$y, tau.0, phi.0, phi.1) == 0) {
    x_vec <- list(x = runif(k), y = runif(k))
  }
  
  for (i in 1:n.iter) {
    #if(i %% 100 == 0) print(i)
    u <- sample(1:k, 1)
    x.prop <- runif(2)
    
    alpha.tmp <- sum(
      sapply(1:k, function(j) {
        x1 <- x.prop - c(x_vec$x[j], x_vec$y[j])
        x2 <- c(x_vec$x[u], x_vec$y[u]) - c(x_vec$x[j], x_vec$y[j])
        return(phi(norm(x1, type = "2"), tau.0, phi.0, phi.1) - phi(norm(x2, type = "2"), tau.0, phi.0, phi.1))
        })
      )
    
    alpha <- min(1, exp(-alpha.tmp))
    
    if(runif(1) <= alpha) {
      num_accept = num_accept + 1
      x_vec$x[u] = x.prop[1]
      x_vec$y[u] = x.prop[2]
    }
  }
  print(paste0("num_accept: ", num_accept))
  return(x_vec)
}

compute_quantiles = function(s, k, tau.0, phi.0, phi.1){
  window = owin(xrange = c(0, 1), yrange = c(0, 1))
  sim <- as.ppp(Strauss_sim(k, tau.0, phi.0, phi.1), window)
  L_0 <- Kfn(sim, fs = 1.4, k = 100)
  Ls <- mclapply(2:s, function(i)  {
    print(i)
    sim = as.ppp(Strauss_sim(k, tau.0, phi.0, phi.1), window)
    return(Kfn(sim, fs = 1.4, k = 100)$y)
  },
    mc.cores = num.cores, mc.silent = FALSE)
  L <- matrix(unlist(Ls), byrow = FALSE, ncol = s - 1)
  L <- cbind(L_0$y, L)
  
  quantiles <- apply(L, 1, function(x) quantile(x, probs = c(0.95, 0.05)))
  Q_mat = data.frame(x = L_0$x, upper = quantiles["95%",], lower = quantiles["5%",])
  Q_mat
}

#compute_quantiles = function(s, k, tau.0, phi.0, phi.1){
#  window = owin(xrange = c(0, 1), yrange = c(0, 1))
#  sim = as.ppp(Strauss_sim(k, tau.0, phi.0, phi.1), window)
#  L_0 = Kfn(sim, fs = 1.4, k = 100)
#  L = matrix(ncol = s, nrow = length(L_0$x))
#  L[,1] = L_0$y
#  for (i in 2:s){
#    print(i)
#    sim = as.ppp(Strauss_sim(k, tau.0, phi.0, phi.1), window)
#    L[,i] = Kfn(sim, fs = 1.4, k = 100)$y
#  }
#  upper_q = c()
#  lower_q = c()
#  for (i in 1:length(L_0$x)){
#    upper_q = c(upper_q, quantile(L[i,], probs = 0.95))
#    lower_q = c(lower_q, quantile(L[i,], probs = 0.05))
#  }
#  Q_mat = data.frame(x = L_0$x, upper = upper_q, lower = lower_q)
#  Q_mat
#}

MC_test = function(s, dataset, k, tau.0, phi.0, phi.1, ...){  # input dataset as ppp
  quantiles = compute_quantiles(s, k, tau.0, phi.0, phi.1)
  plot(Kfn(dataset, fs = 1.4, k = 100), type = "l",  ...)
  lines(x = quantiles$x, y = quantiles$upper, col = "red", ...)
  lines(x = quantiles$x, y = quantiles$lower, col = "red", ...)
}


plot(cells$x, cells$y)
  
tau.0 = 0.1
phi.0 = 0.1
phi.1 = 0.5
res <- Strauss_sim(k, tau.0, phi.0, phi.1)
plot(res$x, res$y, type = "p")

res <- Strauss_sim(k, tau.0 = 0.1, phi.0 = 0.1, phi.1 = 5, n.iter = 1000)
plot(res$x, res$y, type = "p")
MC_test(100, cells, k, tau.0 = 0.1, phi.0 = 0.1, phi.1 = 5)

res <- Strauss_sim(k, tau.0 = 0.1, phi.0 = 0.1, phi.1 = 20, n.iter = 1000)
plot(res$x, res$y, type = "p")

res <- Strauss_sim(k, tau.0 = 0.1, phi.0 = 5, phi.1 = 20, n.iter = 1000)
plot(res$x, res$y, type = "p")

res <- Strauss_sim(k, tau.0 = 0, phi.0 = 5, phi.1 = 20, n.iter = 1000)
plot(res$x, res$y, type = "p")

res <- Strauss_sim(k, tau.0 = 0.1, phi.0 = 10, phi.1 = 2, n.iter = 5000)
plot(res$x, res$y, type = "p")

res <- Strauss_sim(k, tau.0 = 0.05, phi.0 = 5, phi.1 = 20, n.iter = 1000)
plot(res$x, res$y, type = "p")
MC_test(100, cells, k, tau.0 = 0.05, phi.0 = 5, phi.1 = 20)

MC_test(100, cells, k, tau.0 = 0.05, phi.0 = 9, phi.1 = 20)
