library(MASS)
library(spatial)
library(spatstat)

set.seed(123)

cells <- read.table("cells.dat")
pines <- read.table("pines.dat")
redwood <- read.table("redwood.dat")

#pines <- ppinit("pines.dat")

names(cells) <- c("x", "y")
names(pines) <- c("x", "y")
names(redwood) <- c("x", "y")

data <- list(cells = cells, pines = pines, redwood = redwood)
names <- c("Cells", "Pines", "Redwood")

# Problem 1a
for(i in 1:length(data)) {
  plot(x = data[[i]]$x, y = data[[i]]$y, xlim = c(0, 1), ylim = c(0, 1), xlab = "x", ylab = "y", main = names[i])
}


# Problem 1b

ppregion(xl = 0, xu = 1, yl = 0, yu = 1)
window = owin(xrange = c(0, 1), yrange = c(0, 1))

cells = as.ppp(cells, window)
pines = as.ppp(pines, window)
redwood = as.ppp(redwood, window)

sim.unif <- function(n) {
  x <- runif(n, min = 0, max = 1)
  y <- runif(n, min = 0, max = 1)
  
  return(cbind(x, y))
}

sim = as.ppp(sim.unif(), window)
plot(Kfn(sim, fs = 1.4, k = 100))


plot(Kfn(cells, fs = 1.4, k = 100), type = "l")
plot(Kfn(pines, fs = 1.4, k = 100), type = "l")
plot(Kfn(redwood, fs = 1.4, k = 100), type = "l")


plot(Kfn(cells, fs = 1.4, k = 1000), type = "l")
abline(a = 0, b = 1, col = "red")
plot(Kfn(pines, fs = 1.4, k = 1000), type = "l")
abline(a = 0, b = 1, col = "red")
plot(Kfn(redwood, fs = 1.4, k = 1000), type = "l")
abline(a = 0, b = 1, col = "red")

# c)

compute_quantiles = function(s, n){
  window = owin(xrange = c(0, 1), yrange = c(0, 1))
  sim = as.ppp(sim.unif(n), window)
  L_0 = Kfn(sim, fs = 1.4, k = 100)
  L = matrix(ncol = s, nrow = length(L_0$x))
  L[,1] = L_0$y
  for (i in 2:s){
    sim = as.ppp(sim.unif(n), window)
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

MC_test = function(s, dataset){  # input dataset as ppp
  quantiles = compute_quantiles(s, length(dataset$x))
  plot(Kfn(dataset, fs = 1.4, k = 100), type = "l")
  lines(x = quantiles$x, y = quantiles$upper, col = "red")
  lines(x = quantiles$x, y = quantiles$lower, col = "red")
}
MC_test(100, cells)
MC_test(100, pines)
MC_test(100, redwood)

