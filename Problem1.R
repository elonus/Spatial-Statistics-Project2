library(MASS)
library(spatial)
library(spatstat)

set.seed(123)

# Load data
cells <- read.table("cells.dat")
pines <- read.table("pines.dat")
redwood <- read.table("redwood.dat")

names(cells) <- c("x", "y")
names(pines) <- c("x", "y")
names(redwood) <- c("x", "y")

data <- list(cells = cells, pines = pines, redwood = redwood)
names <- c("Cells", "Pines", "Redwood")

# Problem 1a
for(i in 1:length(data)) {
  pdf(paste0("images/", names[i], "_points.pdf"))
  op <- par(cex = 2, cex.lab = 2, cex.main = 2.3, mgp = c(2, 1, 0), mar = c(4.1, 3.1, 3.1, 1.1))
  #op <- par(cex.lab = 2, cex.main = 2, mgp = c(2, 1, 0))
  plot(x = data[[i]]$x, y = data[[i]]$y, xlim = c(0, 1), ylim = c(0, 1), xlab = "x", ylab = "y", main = names[i], pch = 16)
  par(op)
  dev.off()
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

data <- list(cells, pines, redwood)

for(i in 1:length(data)) {
  pdf(paste0("images/", names[i], "_L2hat.pdf"))
  op <- par(cex = 2, cex.lab = 2, cex.main = 2.3, mgp = c(2, 1, 0), mar = c(4.1, 4.6, 3.1, 1.1))
  plot(Kfn(data[[i]], fs = 1.4, k = 100), type = "p", main = names[i], xlab = "t", ylab = expression(hat(L[2])(t)))
  par(op)
  dev.off()
}

for(i in 1:length(data)) {
  pdf(paste0("images/", names[i], "_L2hat_with_line.pdf"))
  op <- par(cex = 2, cex.lab = 2, cex.main = 2.3, mgp = c(2, 1, 0), mar = c(4.1, 4.6, 3.1, 1.1))
  plot(Kfn(data[[i]], fs = 1.4, k = 100), type = "p", main = names[i], xlab = "t", ylab = expression(hat(L[2])(t)))
  abline(a = 0, b = 1, col = "red", lwd = 6)
  par(op)
  dev.off()
}

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

MC_test = function(s, dataset, ...){  # input dataset as ppp
  quantiles = compute_quantiles(s, length(dataset$x))
  plot(Kfn(dataset, fs = 1.4, k = 100), type = "l", ...)
  lines(x = quantiles$x, y = quantiles$upper, col = "red", lty = 2, ...)
  lines(x = quantiles$x, y = quantiles$lower, col = "red", lty = 2, ...)
  par(op)
}

MC_test(100, cells)
MC_test(100, pines)
MC_test(100, redwood)

for(i in 1:length(data)) {
  pdf(paste0("images/", names[i], "_L2hat_with_quantiles.pdf"))
  op <- par(cex = 2, cex.lab = 2, cex.main = 2.3, mgp = c(2, 1, 0), mar = c(4.1, 4.6, 3.1, 1.1))
  MC_test(1000, data[[i]], main = names[i], xlab = "t", ylab = expression(hat(L[2])(t)), lwd = 6)
  par(op)
  dev.off()
}
