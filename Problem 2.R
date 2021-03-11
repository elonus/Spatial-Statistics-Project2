library(MASS)
library(spatial)
library(fields)
library(numbers)

# a)
obs = read.table("obspines.txt", header = TRUE)
prob = read.table("obsprob.txt", header = TRUE)

obs.matrix = matrix(obs$N_obs, byrow = FALSE, nrow = sqrt(nrow(obs)))

image.plot(x = seq(5, 295, by = 10), y = seq(5, 295, by = 10), z = obs.matrix, 
           main = "The number of pine trees observed in each grid unit",
           xlab = "X", ylab = "Y")

prob.matrix = matrix(prob$alpha, byrow = FALSE, nrow = sqrt(nrow(prob)))

image.plot(x = seq(5, 295, by = 10), y = seq(5, 295, by = 10), z = prob.matrix,
           main = "Probability of observing a pine tree occuring in the grid",
           xlab = "X", ylab = "Y")


# kommenter: sammenheng mellom plots

# c)

est.lambda = sum(obs.matrix)/( 100*sum(prob.matrix) )

# simulate from prior

# create grid

simulations = matrix(rpois(900*6, est.lambda*100), ncol = 6)

# event locations
event.locs = list()
for (j in 1:6){
  event.loc = matrix(ncol = 2)
  for (i in 1:900){
    if (simulations[i, j] > 0){
      row = (i-1)%/%30 + 1
      col = mod(i-1, 30) + 1
      simx = runif(simulations[i, j], min = (row-1)*10, max =(row)*10)
      simy = runif(simulations[i, j], min = (col-1)*10, max =(col)*10)
      event.loc = rbind(event.loc, cbind(simx, simy))
    }
  }
  event.loc = event.loc[2:nrow(event.loc), ]
  event.locs[[j]] = event.loc
}
  
op = par(mfrow = c(2, 3), oma = c(1, 0, 2, 0))

plot(x = event.locs[[1]][, 1], y = event.locs[[1]][, 2], type = "p", 
     pch = 19, cex = 0.7, col = "aquamarine3", xlab = "X", ylab = "Y", 
     main = "Realization 1", xlim = c(0,300), ylim = c(0,300))  

plot(x = event.locs[[2]][, 1], y = event.locs[[2]][, 2], type = "p", 
     pch = 19, cex = 0.7, col = "aquamarine3", xlab = "X", ylab = "Y", 
     main = "Realization 2", xlim = c(0,300), ylim = c(0,300))  

plot(x = event.locs[[3]][, 1], y = event.locs[[3]][, 2], type = "p", 
     pch = 19, cex = 0.7, col = "aquamarine3", xlab = "X", ylab = "Y", 
     main = "Realization 3" , xlim = c(0,300), ylim = c(0,300))  

plot(x = event.locs[[4]][, 1], y = event.locs[[4]][, 2], type = "p", 
     pch = 19, cex = 0.7, col = "aquamarine3", xlab = "X", ylab = "Y", 
     main = "Realization 4", xlim = c(0,300), ylim = c(0,300))  

plot(x = event.locs[[5]][, 1], y = event.locs[[5]][, 2], type = "p", 
     pch = 19, cex = 0.7, col = "aquamarine3", xlab = "X", ylab = "Y", 
     main = "Realization 5", xlim = c(0,300), ylim = c(0,300))  

plot(x = event.locs[[6]][, 1], y = event.locs[[6]][, 2], type = "p", 
     pch = 19, cex = 0.7, col = "aquamarine3", xlab = "X", ylab = "Y", 
     main = "Realization 6", xlim = c(0,300), ylim = c(0,300))  
mtext("Prior Poisson eventlocation realizations", outer = TRUE, cex = 2)
par(op)

#d)
prob$alpha

post.event.locs = list()
for (j in (1:6)){
  event.loc = matrix(ncol = 2)
  for (i in (1:900)){
    lambda = (1-prob$alpha[i])*est.lambda*100
    count = rpois(1, lambda = lambda) + obs$N_obs[i]
    if (count>0){
      print(count)
      row = (i-1)%/%30 + 1
      col = mod(i-1, 30) + 1
      simx = runif(count, min = (row-1)*10, max =(row)*10)
      simy = runif(count, min = (col-1)*10, max =(col)*10)
      event.loc = rbind(event.loc, cbind(simx, simy))
    }
  }
  event.loc = event.loc[2:nrow(event.loc), ]
  post.event.locs[[j]] = event.loc
}

par(mfrow = c(2,3), oma = c(1, 0, 2, 0))

plot(x = post.event.locs[[1]][, 1], y = post.event.locs[[1]][, 2], type = "p", 
     pch = 19, cex = 0.7, col = "darkorchid2", xlab = "X", ylab = "Y", 
     main = "Realization 1", xlim = c(0,300), ylim = c(0,300)) 

plot(x = post.event.locs[[2]][, 1], y = post.event.locs[[2]][, 2], type = "p", 
     pch = 19, cex = 0.7, col = "darkorchid2", xlab = "X", ylab = "Y", 
     main = "Realization 2", xlim = c(0,300), ylim = c(0,300))  

plot(x = post.event.locs[[3]][, 1], y = post.event.locs[[3]][, 2], type = "p", 
     pch = 19, cex = 0.7, col = "darkorchid2", xlab = "X", ylab = "Y", 
     main = "Realization 3", xlim = c(0,300), ylim = c(0,300))  

plot(x = post.event.locs[[4]][, 1], y = post.event.locs[[4]][, 2], type = "p", 
     pch = 19, cex = 0.7, col = "darkorchid2", xlab = "X", ylab = "Y", 
     main = "Realization 4", xlim = c(0,300), ylim = c(0,300))  

plot(x = post.event.locs[[5]][, 1], y = post.event.locs[[5]][, 2], type = "p", 
     pch = 19, cex = 0.7, col = "darkorchid2", xlab = "X", ylab = "Y", 
     main = "Realization 5", xlim = c(0,300), ylim = c(0,300))  

plot(x = post.event.locs[[6]][, 1], y = post.event.locs[[6]][, 2], type = "p", 
     pch = 19, cex = 0.7, col = "darkorchid2", xlab = "X", ylab = "Y", 
     main = "Realization 6", xlim = c(0,300), ylim = c(0,300)) 
mtext("Posterior Poisson eventlocation realizations", outer = TRUE, cex = 2)
par(op)

#e)


prior = matrix(rpois(900*100, est.lambda*100), ncol = 100, nrow = 900)
prior = apply(prior, 1, mean)
posterior = rep(0, 900)
for (i in 1:900){
  lambda = (1-prob$alpha[i])*est.lambda*100
  count = rpois(100, lambda = lambda) + obs$N_obs[i]
  posterior[i] = mean(count)
}

prior.matrix = matrix(prior, byrow = FALSE, nrow = sqrt(length(prior)))
post.matrix = matrix(posterior, byrow = FALSE, nrow = sqrt(length(posterior)))



image.plot(x = seq(5, 295, by = 10), y = seq(5, 295, by = 10),
           z = prior.matrix, xlab = "X", ylab = "Y",
           main = "Prior", xlim = c(0,300),
           ylim = c(0,300))
image.plot(x = seq(5, 295, by = 10), y = seq(5, 295, by = 10),
           z = post.matrix, xlab = "X", ylab = "Y",
           main = "Posterior", xlim = c(0,300),
           ylim = c(0,300))
#mtext("Average of 100 realizations of thediscretized event-count model",
     # outer = TRUE, cex = 2)

           
