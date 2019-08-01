library(pdist)
library(XLConnect)
pdist <- pdist::pdist

code_path <- "./code"
data_path <- "./data"
result_path <- "./output"

source("mkpe_embedding_train.R")
# source(sprintf("%s/mkpe_embedding_train.R", code_path))
# 
# args <- commandArgs(trailingOnly = TRUE)
# R <- as.numeric(args[[1]])
# network <- args[[2]]
# seed <- as.numeric(args[[3]])

network <- "nr"
R <- 95

# if (network == "nr") {
#   cd <- sprintf("%s/drugtarget/nr_admat_dgc.txt", data_path)
#   wdd <- sprintf("%s/drugtarget/nr_simmat_dc.txt", data_path)
#   wdp <- sprintf("%s/drugtarget/nr_simmat_dg.txt", data_path)
# } else if (network == "gpcr") {
#   cd <- sprintf("%s/drugtarget/gpcr_admat_dgc.txt", data_path)
#   wdd <- sprintf("%s/drugtarget/gpcr_simmat_dc.txt", data_path)
#   wdp <- sprintf("%s/drugtarget/gpcr_simmat_dg.txt", data_path)
# } else if (network == "e") {
#   cd <- sprintf("%s/drugtarget/e_admat_dgc.txt", data_path)
#   wdd <- sprintf("%s/drugtarget/e_simmat_dc.txt", data_path)
#   wdp <- sprintf("%s/drugtarget/e_simmat_dg.txt", data_path)
# } else if (network == "ic") {
#   cd <- sprintf("%s/drugtarget/ic_admat_dgc.txt", data_path)
#   wdd <- sprintf("%s/drugtarget/ic_simmat_dc.txt", data_path)
#   wdp <- sprintf("%s/drugtarget/ic_simmat_dg.txt", data_path)
# }

if (network == "nr") {
  cd <- "../drugtarget/nr_admat_dgc.txt"
  wdd <- "../drugtarget/nr_simmat_dc.txt"
  wdp <- "../drugtarget/nr_simmat_dg.txt"
} else if (network == "gpcr") {
  cd <- "../drugtarget/gpcr_admat_dgc.txt"
  wdd <- "../drugtarget/gpcr_simmat_dc.txt"
  wdp <- "../drugtarget/gpcr_simmat_dg.txt"
} else if (network == "e") {
  cd <- "../drugtarget/e_admat_dgc.txt"
  wdd <- "../drugtarget/e_simmat_dc.txt"
  wdp <- "../drugtarget/e_simmat_dg.txt"
} else if (network == "ic") {
  cd <- "../drugtarget/ic_admat_dgc.txt"
  wdd <- "../drugtarget/ic_simmat_dc.txt"
  wdp <- "../drugtarget/ic_simmat_dg.txt"
}

rownum <- 10
#colnum <- 2
#objectives_matrix <- data.frame(matrix(NA, 90, colnum, dimnames = list(1:90, c("Objective", "ElapsedTime"))))

#set the seed for random number generator used to initalize random variables
#seed <- 1
for(seed in 1:rownum) {
  
  #cross-domain interaction score
  K_c <- as.matrix(read.table(cd, header = TRUE))
  K_c[K_c == 0] <- NA
  K_c <- K_c * 0.9

  #within-domain similarity score between proteins
  K_x <- as.matrix(read.table(wdp, header = TRUE))

  #within-domain similarity score between drugs
  K_z <- as.matrix(read.table(wdd, header = TRUE))

   
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
  #R <- 2
  
  #set the kernel width used in the subspace
  sigma_e <- sqrt(R)

#initialize the kernels
# xdim <- 50
# zdim <- 40
# K_c <- runif(xdim * zdim, min = 0, max = 1) #should be an N_x x N_z matrix containing cross-domain interactions between samples of domains X and Z
# K_x <- runif(xdim * xdim, min = 0, max = 1) #should be an N_x x N_x matrix containing within-domain similarities between samples of domain X
# K_z <- runif(zdim * zdim, min = 0, max = 1) #should be an N_z x N_z matrix containing within-domain similarities between samples of domain Z
# dim(K_c) <- c(xdim, zdim)
# dim(K_x) <- c(xdim, xdim)
# dim(K_z) <- c(zdim, zdim)

# K_c <- matrix(c(0.3398182, 0.88458175, 0.33904897, 0.01767709, 0.5320574,
#        0.1715670, 0.48588341, 0.94626936, 0.17343856, 0.4547357,
#        0.7536758, 0.71759808, 0.94159669, 0.78774892, 0.6617816,
#        0.9275573, 0.76630969, 0.02883306, 0.90231960, 0.5685501,
#        0.9395922, 0.17853935, 0.05513724, 0.95037413, 0.3058883,
#        0.7801624, 0.87337868, 0.51581468, 0.70330984, 0.3035004,
#        0.2270704, 0.04476505, 0.85510884, 0.64568789, 0.5203025,
#        0.2655142, 0.45434280, 0.63095834, 0.25620187, 0.6663581), nrow = 8, ncol = 5, byrow = TRUE)
# 
# K_x <- matrix(c(0.076784122, 0.1378996, 0.172022896, 0.03781349, 0.98904233, 0.7846628, 0.56901659, 0.7749361,
#        0.063177609, 0.8516822, 0.696318775, 0.22089039, 0.61352354, 0.2090009, 0.41595065, 0.8382859,
#        0.931916154, 0.1160602, 0.398394689, 0.83677730, 0.49690376, 0.3002289, 0.16752731, 0.7115377,
#        0.819686037, 0.3275561, 0.328347114, 0.51604206, 0.64003314, 0.5971587, 0.03479517, 0.3334558,
#        0.933109828, 0.3645221, 0.005380092, 0.18815976, 0.06142783, 0.1607764, 0.65422526, 0.5994234,
#        0.919302087, 0.6802923, 0.277949194, 0.60515441, 0.14376999, 0.6477999, 0.46533041, 0.4127936,
#        0.031248105, 0.9282517, 0.201141551, 0.38797903, 0.25693645, 0.2688659, 0.76590966, 0.4395602,
#        0.005249931, 0.9239280, 0.989641697, 0.31731845, 0.74399936, 0.4909254, 0.94899871, 0.6330180), nrow = 8, ncol = 8, byrow = TRUE)
# 
# K_z <- matrix(c(0.1622318, 0.38095114, 0.6432006, 0.3043789, 0.1094987,
#        0.1952338, 0.13307243, 0.9355000, 0.8101883, 0.1243233,
#        0.1042765, 0.27764779, 0.8323893, 0.7823773, 0.7589163,
#        0.7378679, 0.45444451, 0.4523470, 0.1308205, 0.6714693,
#        0.5529544, 0.08088534, 0.9914341, 0.1592636, 0.4249286), nrow = 5, ncol = 5, byrow = TRUE)


 #set the parameters data frame
  parameters <- data.frame(epsilon, iteration, lambda_c, lambda_x, lambda_z, learn_sigma_e, R, seed, sigma_e)
  
  #perform training
  state <- mkpe_embedding_train(K_c, K_x, K_z, parameters)
  
  #display the embeddings for each domain
  #print(state$E_x)
  #print(state$E_z)
  #if (file.exists(sprintf("%s/mkpe_%s_perplexity_%s_dataset_%s_seed.RData", result_path, R, network, seed)) == FALSE) {
  save("state", file = sprintf("%s/mkpe_%s_perplexity_%s_dataset_%s_seed.RData", result_path, R, network, seed))
  #objectives_matrix[seed, "Objective"] <- tail(state[["objective"]],1)
  #objectives_matrix[seed, "ElapsedTime"] <- state$time[[3]]
}

# book <- loadWorkbook("objectivesMKPE.xlsx")
# createSheet(book, name = "mkpe")
# writeWorksheet(book, objectives_matrix, sheet = "mkpe")
# saveWorkbook(book, file = "objectivesMKPE.xlsx")