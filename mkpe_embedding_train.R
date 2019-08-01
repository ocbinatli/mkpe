mkpe_embedding_train <- function(K_c, K_x, K_z, parameters) {

  N_x <- dim(K_x)[1]
  N_z <- dim(K_z)[1]
  R <- parameters$R
  sigma_e <- parameters$sigma_e
  epsilon <- parameters$epsilon
  lambda_c <- parameters$lambda_c
  lambda_x <- parameters$lambda_x
  lambda_z <- parameters$lambda_z
  iteration <- parameters$iteration
  learn_sigma_e <- parameters$learn_sigma_e
  
  set.seed(parameters$seed)

  indices_c <- which(is.na(K_c) == FALSE, arr.ind = T)
  count_c <- length(indices_c) / 2
  indices_x <- which(is.na(K_x) == FALSE, arr.ind = T)
  count_x <- length(indices_x) / 2
  indices_z <- which(is.na(K_z) == FALSE, arr.ind = T)
  count_z <- length(indices_z) / 2
  
  gamma_x <- 1
  gamma_z <- 1
  gamma_eta <- 1
  # P_x <- matrix(data = c(-1.53890033193870,	-0.773602845680918,
  #                             0.855129459504855,	-0.708036748219266,
  #                             -1.36657535786149,	0.864412969239083,
  #                             -1.76848564014379,	-0.0147542510685107,
  #                             1.32966647924592,	-0.193043562883580,
  #                             -0.449991936643678,	-1.03990751855869,
  #                             0.111448676021657,	-0.143056487916747,
  #                             -0.313294452779210,	0.543501639178857,
  #                             0.226919835123064,	0.681875759919915,
  #                             0.0650407997634434,	-0.299458280773562,
  #                             -0.642153510700408,	0.512463042891312,
  #                             0.976894558627171,	-1.00759603736923,
  #                             -1.54525975621760,	-0.0996245949315373,
  #                             -0.122850842348634,	-0.776030472241980,
  #                             0.842179336633990,	-0.552788028173069,
  #                             -2.38074571381066,	0.0790146698431145,
  #                             0.0511147726016121,	-0.660657127802827,
  #                             -0.487392978357807,	1.23130555940409,
  #                             -1.08237056501740,	-0.811602483767501,
  #                             0.950708238940756,	-0.593397071286631,
  #                             0.355116197539052,	0.525930896643522,
  #                             -0.427670870563882,	0.558499203559599,
  #                             0.763068554770551,	-0.647110587659478,
  #                             0.222418835620462,	-0.930569427147368,
  #                             0.256383006977626,	-0.794578978950119,
  #                             0.541583576031859,	-0.351755901806197), nrow = 26, ncol = 2, byrow = TRUE)
  # 
  # P_z <- matrix(data = c(-0.503226819050481,	0.177372323481206,
  #                      -1.21151877072153,	1.04652729129208,
  #                      1.62466055658946,	-0.166830905435649,
  #                      0.127024793054780,	-0.0404769739633148,
  #                      -1.15290088706912,	-0.111445595064485,
  #                      0.0682101654746600,	-1.05389239077415,
  #                      -0.266071678245521,	1.16576278466057,
  #                      -0.143819310110554,	-1.19288780123853,
  #                      0.486037926261697,	0.681460227698107,
  #                      -0.349997322499615,	-1.00456952416151,
  #                      1.13980169450248,	0.177720143089890,
  #                      -1.11131414630241,	-0.671171261523487,
  #                      1.37578965560780,	0.598993736963102,
  #                      -0.158970404728823,	-0.0554886026223039,
  #                      -0.731674839060567,	-0.147427711434326,
  #                      -1.28542120433262,	-0.670806768921277,
  #                      -0.132852124083494,	-0.352731962170597,
  #                      0.470682288358934,	-1.51335739744481,
  #                      1.32562367153053,	0.0484245435645341,
  #                      -2.38784201404015,	-0.552518306701820,
  #                      0.272517029325722,	2.77269410908998,
  #                      -0.0858764790485010,	-0.918493817866487,
  #                      0.0679431195276638,	0.627693868009014,
  #                      -0.696890907765267,	1.55446030564446,
  #                      0.485550926690108,	0.407948239275143,
  #                      -0.0708983963137185,	-0.0345005411497051,
  #                      0.551435811346911,	-0.363809536796243,
  #                      0.825085512685776,	1.39415026317769,
  #                      -0.782448082378367,	-0.426401607393573,
  #                      -0.607948535720416,	1.17179571701245,
  #                      -0.0598769673452117,	-0.0475655576218983,
  #                      -0.441943961123237,	-0.329504009863933,
  #                      -0.864170557961357,	0.197150527947774,
  #                      -0.0386719982948745,	-1.16767364226629,
  #                      -1.06123297942630,	-0.157609981836122,
  #                      0.391479624899552,	1.60245481610666,
  #                      0.657401171982278,	0.327120579502262,
  #                      1.43978940003184,	0.277453091269237,
  #                      -1.12042498662545,	1.70636518647107,
  #                      0.615304208468136,	1.61842973815650,
  #                      0.620361198172113,	2.02480416242137,
  #                      -0.441980384297731,	-0.344065574628753,
  #                      -0.346126814089098,	0.490137958550425,
  #                      -0.216582823920814,	-2.12180466899737,
  #                      0.0907110628586373,	-0.677772730987262,
  #                      1.06552144598826,	0.120612804840037,
  #                      0.377574051107238,	-1.48832913178564,
  #                      -0.695214013508033,	-1.05918103922983,
  #                      -1.33358526076681,	-1.67019646151308,
  #                      -0.371611133139236,	0.991840443709460,
  #                      -1.11171426888148,	-1.27915519766338,
  #                      -0.142469433590110,	0.729732751090863,
  #                      -0.235292317153440,	-0.349069024707421,
  #                      -0.374880147961807,	0.495074916294177), nrow = 54, ncol = 2, byrow = TRUE)
   # P_x <- matrix(c(0.02780379,  0.56666883, -1.01419299,  0.33316433, -0.94052644,  1.90717476, -1.17342867, -0.79201270, 2.03210402, 1.12236415, 0.89393422, -0.75422609,  1.13852185, -1.21158750, -0.29333122,  1.05149459), nrow = 8, ncol =2, byrow = TRUE)
   # P_z <- matrix(c(0.30771813, -0.04885156,  1.72654839, -0.49162719, -0.26494832, -0.41863172,  0.58786279,  1.28437716,  0.98469239, -0.02472406), nrow = 5, ncol = 2, byrow = TRUE)
  P_x <- rnorm(N_x * R)
  P_z <- rnorm(N_z * R)
  dim(P_x) <- c(N_x, R)
  dim(P_z) <- c(N_z, R)
  E_x <- project_to_stiefel_manifold(P_x)
  E_z <- project_to_stiefel_manifold(P_z)
  DE_c <- (as.matrix(pdist(E_x, E_z)))^2
  KE_c <- exp(-DE_c / sigma_e^2)
  #print(KE_c)
  DE_x <- (as.matrix(pdist(E_x, E_x)))^2
  KE_x <- exp(-DE_x / sigma_e^2)
  DE_z <- (as.matrix(pdist(E_z, E_z)))^2
  KE_z <- exp(-DE_z / sigma_e^2)
  objective_c <- (sum((KE_c[indices_c] - K_c[indices_c])^2)) / count_c
  objective_x <- (sum((KE_x[indices_x] - K_x[indices_x])^2)) / count_x
  objective_z <- (sum((KE_z[indices_z] - K_z[indices_z])^2)) / count_z
  objective <- lambda_c * objective_c + lambda_x * objective_x + lambda_z * objective_z
  
  iter <- 0
  print(paste("iteration=", iter, "objective=", objective))
  elapsed_time <- system.time({
  while (TRUE) {
    iter <- iter + 1
    return_x <- 0
    return_z <- 0
    return_eta <- 0
    
    E_x_gradient <- matrix(0L, nrow = N_x, ncol = R) 
    for (i in 1:N_x) {
      j <- !is.na(K_c[i,])
      if (sum(j) == 1) {
        E_x_gradient[i,] <- c(E_x_gradient[i,] - 4 * lambda_c * (as.vector(t(KE_c[i, j]) * t(t(KE_c[i, j] - t(K_c[i, j]))) ) * t(E_x[i,] - t(E_z[j,])) ) / sigma_e^2 / count_c)
      } else {
        E_x_gradient[i,] <- E_x_gradient[i,] - 4 * lambda_c * colSums(as.vector(t(KE_c[i, j]) * t(t(KE_c[i, j] - t(K_c[i, j]))) ) * t(E_x[i,] - t(E_z[j,])) ) / sigma_e^2 / count_c
      } 
    }

    if (lambda_x != 0) {
      for (i in 1:N_x) {
        j <- !is.na(K_x[i,])
        if (sum(j) == 1) {
          E_x_gradient[i,] <- c(E_x_gradient[i,] - 8 * lambda_x * (as.vector(t(KE_x[i, j]) * t(t(KE_x[i, j] - t(K_x[i, j]))) ) * t(E_x[i,] - t(E_x[j,])) ) / sigma_e^2 / count_x)
        } else {
          E_x_gradient[i,] <- E_x_gradient[i,] - 8 * lambda_x * colSums(as.vector(t(KE_x[i, j]) * t(t(KE_x[i, j] - t(K_x[i, j]))) ) * t(E_x[i,] - t(E_x[j,])) )  / sigma_e^2 / count_x
        }
      }
    }

    E_x_gradient <- (E_x %*% t(E_x_gradient) %*% E_x - E_x_gradient)
    E_x_gradient_norm <- sum(diag(t(E_x_gradient) %*% (diag(N_x) - 0.5 * (E_x %*% t(E_x))) %*% E_x_gradient))
    if (sqrt(E_x_gradient_norm) < epsilon) {
      return_x <- 1
      objective_c <- c(objective_c, objective_c[length(objective_c)])
      objective_x <- c(objective_x, objective_x[length(objective_x)])
      objective_z <- c(objective_z, objective_z[length(objective_z)])
      objective <- c(objective, objective[length(objective)])
    } else {
        while (TRUE) {
          E_x_new <- E_x + 2 * gamma_x * E_x_gradient
          E_x_new <- project_to_stiefel_manifold(E_x_new)
          DE_c <- (as.matrix(pdist(E_x_new, E_z)))^2
          KE_c <- exp(-DE_c / sigma_e^2)
          DE_x <- (as.matrix(pdist(E_x_new, E_x_new)))^2
          KE_x <- exp(-DE_x / sigma_e^2)
          objective_new <- lambda_c * sum((KE_c[indices_c] - K_c[indices_c])^2) / count_c + lambda_x * sum((KE_x[indices_x] - K_x[indices_x])^2) / count_x + lambda_z * sum((KE_z[indices_z] - K_z[indices_z])^2) / count_z
          if (gamma_x == 0) {
            break
          } else if (objective[length(objective)] - objective_new >= gamma_x * E_x_gradient_norm) {
            gamma_x <- 2 * gamma_x
          } else {
            break
          }
        }
    
        while (TRUE){
            E_x_new <- E_x + gamma_x * E_x_gradient
            E_x_new <- project_to_stiefel_manifold(E_x_new)
            DE_c <- (as.matrix(pdist(E_x_new, E_z)))^2
            KE_c <- exp(-DE_c / sigma_e^2)
            DE_x <- (as.matrix(pdist(E_x_new, E_x_new)))^2
            KE_x <- exp(-DE_x / sigma_e^2)
            objective_new <- lambda_c * sum((KE_c[indices_c] - K_c[indices_c])^2) / count_c + lambda_x * sum((KE_x[indices_x] - K_x[indices_x])^2) / count_x + lambda_z * sum((KE_z[indices_z] - K_z[indices_z])^2) / count_z
            #if ((objective[length(objective)] - objective_new < 0.5 * gamma_x * E_x_gradient_norm) || (gamma_x > 0)){
            if (gamma_x == 0) {
              break
            } else if (objective[length(objective)] - objective_new < 0.5 * gamma_x * E_x_gradient_norm) {
              gamma_x <- 0.5 * gamma_x
            } else {
              break
            }
        }
    
        E_x <- E_x + gamma_x * E_x_gradient
        E_x <- project_to_stiefel_manifold(E_x)
        DE_c <- (as.matrix(pdist(E_x, E_z)))^2
        KE_c <- exp(-DE_c / sigma_e^2)
        DE_x <- (as.matrix(pdist(E_x, E_x)))^2
        KE_x <- exp(-DE_x / sigma_e^2)
        objective_c_last <- sum((KE_c[indices_c] - K_c[indices_c])^2) / count_c
        objective_x_last <- sum((KE_x[indices_x] - K_x[indices_x])^2) / count_x
        objective_z_last <- sum((KE_z[indices_z] - K_z[indices_z])^2) / count_z
        objective_last <- lambda_c * objective_c_last + lambda_x * objective_x_last + lambda_z * objective_z_last
        objective_c <- c(objective_c, objective_c_last)
        objective_x <- c(objective_x, objective_x_last)
        objective_z <- c(objective_z, objective_z_last)
        objective <- c(objective, objective_last)
        print(paste("iteration=", iter, "objective=", objective_last, "norm_x=", sqrt(E_x_gradient_norm), "gamma_x=", gamma_x))
    }
    
    E_z_gradient <- matrix(0L, nrow = N_z, ncol = R) 
    for (j in 1:N_z) {
      i <- !is.na(K_c[, j])
      if (sum(i) == 1) {
        E_z_gradient[j,] <- c(E_z_gradient[j,] - 4 * lambda_c * (as.vector(t(KE_c[i, j]) * t(t(KE_c[i, j] - t(K_c[i, j]))) ) * t(E_z[j,] - t(E_x[i,])) ) / sigma_e^2 / count_c)
      } else {
        E_z_gradient[j,] <- E_z_gradient[j,] - 4 * lambda_c * colSums(as.vector(t(KE_c[i, j]) * t(t(KE_c[i, j] - t(K_c[i, j]))) ) * t(E_z[j,] - t(E_x[i,])) ) / sigma_e^2 / count_c
      }
    }
    if (lambda_z != 0){
      for (i in 1:N_z){
        j <- !is.na(K_z[i,])
        if (sum(j) == 1) {
          E_z_gradient[i,] <- c(E_z_gradient[i,] - 8 * lambda_z * (as.vector(t(KE_z[i, j]) * t(t(KE_z[i, j] - t(K_z[i, j]))) ) * t(E_z[i,] - t(E_z[j,])) ) / sigma_e^2 / count_z)
        } else {
          E_z_gradient[i,] <- E_z_gradient[i,] - 8 * lambda_z * colSums(as.vector(t(KE_z[i, j]) * t(t(KE_z[i, j] - t(K_z[i, j]))) ) * t(E_z[i,] - t(E_z[j,])) ) / sigma_e^2 / count_z
        }
      }
    }
    E_z_gradient <- (E_z %*% t(E_z_gradient) %*% E_z - E_z_gradient)
    E_z_gradient_norm <- sum(diag(t(E_z_gradient) %*% (diag(N_z) - 0.5 * (E_z %*% t(E_z))) %*% E_z_gradient))
    if (sqrt(E_z_gradient_norm) < epsilon) {
      return_z <- 1
      objective_c <- c(objective_c, objective_c[length(objective_c)])
      objective_x <- c(objective_x, objective_x[length(objective_x)])
      objective_z <- c(objective_z, objective_z[length(objective_z)])
      objective <- c(objective, objective[length(objective)])
    } else {
        while (TRUE) {
          E_z_new <- E_z + 2 * gamma_z * E_z_gradient
          E_z_new <- project_to_stiefel_manifold(E_z_new)
          DE_c <- (as.matrix(pdist(E_x, E_z_new)))^2
          KE_c <- exp(-DE_c / sigma_e^2)
          DE_z <- (as.matrix(pdist(E_z_new, E_z_new)))^2
          KE_z <- exp(-DE_z / sigma_e^2)
          objective_new <- lambda_c * sum((KE_c[indices_c] - K_c[indices_c])^2) / count_c + lambda_x * sum((KE_x[indices_x] - K_x[indices_x])^2) / count_x + lambda_z * sum((KE_z[indices_z] - K_z[indices_z])^2) / count_z
          if (gamma_z == 0) {
            break
          } else if (objective[length(objective)] - objective_new >= gamma_z * E_z_gradient_norm) {
            gamma_z <- 2 * gamma_z
          } else {
            break
          }
        }
      
        while (TRUE) {
          E_z_new <- E_z + gamma_z * E_z_gradient
          E_z_new <- project_to_stiefel_manifold(E_z_new)
          DE_c <- (as.matrix(pdist(E_x, E_z_new)))^2
          KE_c <- exp(-DE_c / sigma_e^2)
          DE_z <- (as.matrix(pdist(E_z_new, E_z_new)))^2
          KE_z <- exp(-DE_z / sigma_e^2)
          objective_new <- lambda_c * sum((KE_c[indices_c] - K_c[indices_c])^2) / count_c + lambda_x * sum((KE_x[indices_x] - K_x[indices_x])^2) / count_x + lambda_z * sum((KE_z[indices_z] - K_z[indices_z])^2) / count_z
          if (gamma_z == 0) {
            break
          } else if (objective[length(objective)] - objective_new < 0.5 * gamma_z * E_z_gradient_norm) {
            gamma_z <- 0.5 * gamma_z
          } else {
            break
          }
        }
      
      E_z <- E_z + gamma_z * E_z_gradient
      E_z <- project_to_stiefel_manifold(E_z)
      DE_c <- (as.matrix(pdist(E_x, E_z)))^2
      KE_c <- exp(-DE_c / sigma_e^2)
      DE_z <- (as.matrix(pdist(E_z, E_z)))^2
      KE_z <- exp(-DE_z / sigma_e^2)
      objective_c_last <- sum((KE_c[indices_c] - K_c[indices_c])^2) / count_c
      objective_x_last <- sum((KE_x[indices_x] - K_x[indices_x])^2) / count_x
      objective_z_last <- sum((KE_z[indices_z] - K_z[indices_z])^2) / count_z
      objective_last <- lambda_c * objective_c_last + lambda_x * objective_x_last + lambda_z * objective_z_last
      objective_c <- c(objective_c, objective_c_last)
      objective_x <- c(objective_x, objective_x_last)
      objective_z <- c(objective_z, objective_z_last)
      objective <- c(objective, objective_last)
      
      print(paste("iteration=", iter, "objective=", objective_last, "norm_z=", sqrt(E_z_gradient_norm), "gamma_z=", gamma_z))
    }
   
    if (learn_sigma_e == 1) {
      eta_gradient <- 0
      eta_gradient <- eta_gradient - 4 * lambda_c * sum(KE_c[indices_c] * (KE_c[indices_c] - K_c[indices_c]) * DE_c[indices_c]) / sigma_e^2 / count_c
      eta_gradient <- eta_gradient - 4 * lambda_x * sum(KE_x[indices_x] * (KE_x[indices_x] - K_x[indices_x]) * DE_x[indices_x]) / sigma_e^2 / count_x
      eta_gradient <- eta_gradient - 4 * lambda_z * sum(KE_z[indices_z] * (KE_z[indices_z] - K_z[indices_z]) * DE_z[indices_z]) / sigma_e^2 / count_z
      eta_gradient_norm <- eta_gradient^2
      if (sqrt(eta_gradient_norm) < epsilon){
        return_eta <- 1
        objective_c <- c(objective_c, objective_c[length(objective_c)])
        objective_x <- c(objective_x, objective_x[length(objective_x)])
        objective_z <- c(objective_z, objective_z[length(objective_z)])
        objective <- c(objective, objective[length(objective)])
      } else {
        while (TRUE){
          sigma_e_new <- exp(log(sigma_e) + 2 * gamma_eta * eta_gradient)
          KE_c <- exp(-DE_c / sigma_e_new^2)
          KE_x <- exp(-DE_x / sigma_e_new^2)
          KE_z <- exp(-DE_z / sigma_e_new^2)
          objective_new <- lambda_c * sum((KE_c[indices_c] - K_c[indices_c])^2) / count_c + lambda_x * sum((KE_x[indices_x] - K_x[indices_x])^2) / count_x + lambda_z * sum((KE_z[indices_z] - K_z[indices_z])^2) / count_z
          if (gamma_eta == 0) {
            break
          } else if (objective[length(objective)] - objective_new >= gamma_eta * eta_gradient_norm) {
            gamma_eta <- 2 * gamma_eta
          } else {
            break
          }
        }
        
        while (TRUE) {
          sigma_e_new <- exp(log(sigma_e) + gamma_eta * eta_gradient)
          KE_c <- exp(-DE_c / sigma_e_new^2)
          KE_x <- exp(-DE_x / sigma_e_new^2)
          KE_z <- exp(-DE_z / sigma_e_new^2)
          objective_new <- lambda_c * sum((KE_c[indices_c] - K_c[indices_c])^2) / count_c + lambda_x * sum((KE_x[indices_x] - K_x[indices_x])^2) / count_x + lambda_z * sum((KE_z[indices_z] - K_z[indices_z])^2) / count_z
          if (gamma_eta == 0) {
            break
          } else if (objective[length(objective)] - objective_new < 0.5 * gamma_eta * eta_gradient_norm) {
            gamma_eta <- 0.5 * gamma_eta
          } else {
            break
          }
        }
        
        sigma_e <- exp(log(sigma_e) + gamma_eta * eta_gradient)
        KE_c <- exp(-DE_c / sigma_e_new^2)
        KE_x <- exp(-DE_x / sigma_e_new^2)
        KE_z <- exp(-DE_z / sigma_e_new^2)
        objective_c_last <- sum((KE_c[indices_c] - K_c[indices_c])^2) / count_c
        objective_x_last <- sum((KE_x[indices_x] - K_x[indices_x])^2) / count_x
        objective_z_last <- sum((KE_z[indices_z] - K_z[indices_z])^2) / count_z
        objective_last <- lambda_c * objective_c_last + lambda_x * objective_x_last + lambda_z * objective_z_last
        objective_c <- c(objective_c, objective_c_last)
        objective_x <- c(objective_x, objective_x_last)
        objective_z <- c(objective_z, objective_z_last)
        objective <- c(objective, objective_last)
        print(paste("iteration=", iter, "objective=", objective_last, "norm_e=", sqrt(eta_gradient_norm), "gamma_e=", gamma_eta))
      }
    } else {
      return_eta <- 1
      objective_c <- c(objective_c, objective_c[length(objective_c)])
      objective_x <- c(objective_x, objective_x[length(objective_x)])
      objective_z <- c(objective_z, objective_z[length(objective_z)])
      objective <- c(objective, objective[length(objective)])
    }
    
    if (return_x == 1 && return_z == 1 && return_eta == 1){
      break
    }
    if (iter == iteration){
      break
    }
  }
  })
  #print(elapsed_time)
  #set the state list
  state <- list(E_x = E_x, E_z = E_z, sigma_e = sigma_e, objective_c = objective_c, objective_x = objective_x, objective_z = objective_z, objective = objective, time = elapsed_time)
}

project_to_stiefel_manifold <- function(Q) {
  s <- svd(Q, nu = nrow(Q), nv = ncol(Q))
  Q <- s$u %*% diag(nrow = dim(Q)[1], ncol = dim(Q)[2]) %*% t(s$v)
}