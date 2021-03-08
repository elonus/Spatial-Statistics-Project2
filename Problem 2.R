library(MASS)
library(spatial)
library(fields)

obs = read.table("obspines.txt", header = TRUE)
prob = read.table("obsprob.txt", header = TRUE)

obs.matrix = matrix(obs$N_obs, byrow = FALSE, nrow = sqrt(nrow(obs)))

image.plot(x = seq(5, 295, by = 10), y = seq(5, 295, by = 10), z = obs.matrix)

prob.matrix = matrix(prob$alpha, byrow = FALSE, nrow = sqrt(nrow(prob)))

image.plot(x = seq(5, 295, by = 10), y = seq(5, 295, by = 10), z = prob.matrix)

# kommenter: sammenheng mellom plots

