library(MASS)
library(spatial)

set.seed(123)

cells <- read.table("cells.dat")
pines <- read.table("pines.dat")
redwood <- read.table("redwood.dat")

pines <- ppinit("pines.dat")

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
maxdist <- sqrt(2)
num.points <- 100

j.est.one.point <- function(t, data) {
  res <- 0
  for(i in 1:nrow(data)) {
    tmp <- 0
    for(j in 1:nrow(data)) {
      if(norm(data[i,] - data[j,], type = "2") <= t) {
        tmp <- tmp + 1
      }
    }
    res <- res + tmp - 1
  }
  res <- res / nrow(data)
  res <- res / (pi * t^2)
  return(res)
}

j.est <- function(data, num.points, xmin, xmax, ymin, ymax) {
  maxr <- norm(c(xmin, ymin) - c(xmax, ymax), type = "2")
  maxr <- 0.1
  res <- list(x = vector(mode = "numeric", length = num.points), 
              y = vector(mode = "numeric", length = num.points))
  for(i in 1:num.points) {
    res$x[i] <- i*maxr/num.points
    res$y[i] <- j.est.one.point(i*maxr/num.points, data)
  }
  return(res)
}

J <- j.est(pines, num.points = 20, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
plot(x = J$x, y = J$y, ylim = c(0, max(J$y)))

for(i in 1:length(data)){
  J <- j.est(data[[i]], num.points = 40, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
  plot(x = J$x, y = J$y, ylim = c(0, max(J$y)), main = names[i])
}
  
sim.unif <- function() {
  n <- rpois(1, lambda = 100)
  x <- runif(n, min = 0, max = 10)
  y <- runif(n, min = 0, max = 10)
  
  return(cbind(x, y))
}

tmp <- sim.unif()
tmp <- rpoispp(lambda = 100, win = square(1))
plot(x = tmp$x, y = tmp$y, xlim = c(0, 1), ylim = c(0, 1))
J <- j.est(tmp, num.points = 20, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
plot(x = J$x, y = J$y, ylim = c(0, max(J$y)))

tmp <- as.data.frame(tmp)
names(tmp) = c("x", "y")

plot(Kfn(tmp, fs = 1.4, k = 100))

tmp2 <- as.ppp(tmp, c(0, 10, 0, 10))

jest <- Jest(tmp2, r = seq(from = 0, to = 1.14, by = 0.001))

plot(Jest(tmp2, r = seq(from = 0, to = 1.14, by = 0.001)))

plot(Jest(as.ppp(pines, c(0, 1, 0, 1))))
jest <- Jest(as.ppp(pines, c(0, 1, 0, 1)))

#for(i in 1:length(data)) {
#  plot(Kfn(data[[i]], fs = maxdist, k = num.points), type = "l", xlab = "t", ylab = expression(L[2](t)), 
#       main = paste0("L(t) for ", names[i], " dataset"), ylim = c(0, 7))
#}


