library(Matrix)
library(splines)

rm(list=ls())
rr2beta = function(x) (log((100+x)/100)/10);# x: Risk (% change) associated with 10 ug/m3 increase in PM
beta2rr = function(x) (exp(x*10)-1)*100;


source('functions.R')

ST = c(1, 0)
factor=.001

p = 117
t = 600


W = as.matrix(read.csv('Wcont.csv', F))
attributes(W)$dimnames = NULL
D=diag(colSums(W))
coords = read.csv('centroids.csv')

# load simulation data
load('example_dataset.RData')


# simulate beta
tempo <- (1:t)
normalize_range <- function(x, new_min, new_max) {
  min_x <- min(x)
  max_x <- max(x)
  scaled_x <- (x - min_x) / (max_x - min_x)
  new_min + scaled_x * (new_max - new_min)
}

coordsMX <- matrix(0, nrow(coords), 2)
coordsMX[, 1] <- normalize_range(coords[, 1], 0, 1)
coordsMX[, 2] <- normalize_range(coords[, 2], 0, 1)

rrval0 <- dnorm(seq(0, 1, length.out = t), mean = 0.07/2, sd = 0.05) +
  dnorm(seq(0, 1, length.out = t), mean = 0.6, sd = 0.08) * 1.5 +
  dnorm(seq(0, 1, length.out = t), mean = 1, sd = 0.15)^2

rrval0 <- normalize_range(rrval0, 1.5, 5)

dd1 <- as.matrix(dist(rbind(coordsMX, c(.4, .4))))
dd1 <- dd1[1:nrow(coordsMX), ncol(dd1)]

dd2 <- as.matrix(dist(rbind(coordsMX, c(.6, .6))))
dd2 <- dd2[1:nrow(coordsMX), ncol(dd2)]

rrval1 <- (- 6*(coordsMX[,1]-.1)^2 - 6*(coordsMX[,2]-.1)^2) * exp(-dd1/.5) +
  (- 6*(coordsMX[,1]-.9)^2 - 6*(coordsMX[,2]-.9)^2) * exp(-dd2/.5)

rrval1 <- rrval1 - mean(rrval1)
rm(dd1, dd2)


rrval3 <- kronecker(matrix(-4*(tempo-t/2)^2/t^2, nrow = 1), rrval1)
rrval3 <- rrval3 - rowMeans(rrval3) # Subtract row means
rrval3 <- rrval3 - colMeans(rrval3) # Subtract column means

beta <- rr2beta(2+mean(rrval0) + 1*(1*t(replicate(p, rrval0 - mean(rrval0))) + 4*replicate(t, rrval1)) + 4*rrval3*ST[2])
AB_it <- beta - replicate(t, rowMeans(beta)) - t(replicate(p, colMeans(beta))) + mean(beta)
beta <- beta + AB_it*(ST[2] - 1)
beta <- beta - (mean(beta) - round(mean(beta), 3))

matplot(t(beta), type = 'l', main = 'Betax')


beta0 <- 1.5


# offset
offset = as.matrix(read.csv('offset.csv', F))



# Preparing... ######
interaction = 4
nrep=10
nburn=10
ctuning = 1.5
thin=1


priors = list(
  a0=0.01, b0=0.01, 
  a1=10, b1=1, 
  s2_a = 1, s2_b = NA
)
# https://academic.oup.com/biostatistics/article/8/2/158/230741#2641848
priors$s2_b = log(1.02)^2*priors$s2_a/qt(.025,2*priors$s2_a)^2



Zmeasconf = vector('list', p)
for (i in 1:p) {
  Zmeasconf[[i]] = matrix(1, 1, t)
}


# FIT THE MODEL ####
res = fit_SDGLMC(Y, X*factor, Xmean*factor, Zmeasconf, 
                 W, offset, interaction, nrep, nburn, 
                 thin, ctuning, priors)
















