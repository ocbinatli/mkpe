library(pdist)
pdist <- pdist::pdist

source("mkpe_embedding_train.R")

# set the seed for random number generator used to initalize random variables
seed <- ???
# N_x x N_x drug compound structure similarity matrix
wdd <- sprintf("./%s_simmat_dc.txt", network)
# N_z x N_z protein sequence similarity matrix
wdp <- sprintf("./%s_simmat_dg.txt", network)
# N_x x N_z adjacency matrix containing drug-target interactions
cd <- sprintf("./%s_admat_dgc.txt", network)

K_c <- t(as.matrix(read.table(cd, header = TRUE)))
K_c[K_c == 0] <- NA
K_c <- K_c * 0.9

K_x <- t(as.matrix(read.table(wdd, header = TRUE)))
K_x <- (K_x + t(K_x)) / 2

K_z <- as.matrix(read.table(wdp, header = TRUE))
K_z <- (K_z + t(K_z)) / 2

#set the convergence threshold parameter
epsilon <- -Inf

#set the maximum number of iterations
iteration <- 100
#set the regularization parameter for cross-domain interactions
lambda_c <- 1.0

#set the regularization parameters for within-domain similarities
lambda_x <- 0.1
lambda_z <- 0.1

#determine whether you want to learn the kernel width used in the subspace
learn_sigma_e <- 1

#set the subspace dimensionality
R <- 2

#set the kernel width used in the subspace
sigma_e <- sqrt(R)

parameters <- data.frame(epsilon, iteration, lambda_c, lambda_x, lambda_z, learn_sigma_e, R, seed, sigma_e)

#perform training
state <- mkpe_embedding_train(K_c, K_x, K_z, parameters)

#display the embeddings for each domain
print(state$E_x)
print(state$E_z)

save("state", file = ???)