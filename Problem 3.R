library(MASS)

# We gestimate lambda_M to be 5

lambda_M = 8
sigma = 0.07 
lambda_k = 8 # antall rundt hver mor

kM = rpois(1, lambda_M)  # antall mødre

Neuman_sim = function(lambda_M, sigma, lambda_k){
  x_vec = list(x = c(), y = c())
  for (i in 1:kM){
    x = runif(1)
    y = runif(1)
    kcs = rpois(1, lambda_k)
    for (j in 1:kcs){
      xps = mvrnorm(1, c(x, y), diag(2)*sigma^2)
      x_vec$x = c(x_vec$x, xps[1] %% 1)
      x_vec$y = c(x_vec$y, xps[2] %% 1)
    }
  }
  x_vec
}



compute_quantiles = function(s, lambda_M, sigma, lambda_k){
  window = owin(xrange = c(0, 1), yrange = c(0, 1))
  sim = as.ppp(Neuman_sim(lambda_M, sigma, lambda_k), window)
  L_0 = Kfn(sim, fs = 1.4, k = 100)
  L = matrix(ncol = s, nrow = length(L_0$x))
  L[,1] = L_0$y
  for (i in 2:s){
    sim = as.ppp(Neuman_sim(lambda_M, sigma, lambda_k), window)
    L[,i] = Kfn(sim, fs = 1.4, k = 100)$y
  }
  upper_q = c()
  lower_q = c()
  for (i in 1:length(L_0$x)){
    upper_q = c(upper_q, quantile(L[i,], probs = 0.95))
    lower_q = c(lower_q, quantile(L[i,], probs = 0.05))
  }
  Q_mat = data.frame(x = L_0$x, upper = upper_q, lower = lower_q)
  Q_mat
}

MC_test = function(s, dataset, lambda_M, sigma, lambda_k){  # input dataset as ppp
  quantiles = compute_quantiles(s, lambda_M, sigma, lambda_k)
  plot(Kfn(dataset, fs = 1.4, k = 100), type = "l")
  lines(x = quantiles$x, y = quantiles$upper, col = "red")
  lines(x = quantiles$x, y = quantiles$lower, col = "red")
}

length(redwood$x) # = 62

MC_test(100, redwood, lambda_M = 8, sigma = 0.07, lambda_k = 8)
MC_test(100, redwood, lambda_M = 8, sigma = 0.08, lambda_k = 8)
MC_test(100, redwood, lambda_M = 8, sigma = 0.1, lambda_k = 8)
MC_test(100, redwood, lambda_M = 5, sigma = 0.1, lambda_k = 12)
MC_test(100, redwood, lambda_M = 5, sigma = 0.15, lambda_k = 12)
MC_test(100, redwood, lambda_M = 6, sigma = 0.1, lambda_k = 10)
MC_test(100, redwood, lambda_M = 10, sigma = 0.1, lambda_k = 6) # best
plot(Neuman_sim(lambda_M = 10, sigma = 0.1, lambda_k = 6)) # maybe not clustered

MC_test(100, redwood, lambda_M = 10, sigma = 0.07, lambda_k = 6)
plot(Neuman_sim(lambda_M = 10, sigma = 0.07, lambda_k = 6))


