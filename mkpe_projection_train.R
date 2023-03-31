mkpe_projection_train <- function(K_c, K_x, K_z, X, Z, parameters) {
  start_time <- Sys.time()
  N_x <- dim(K_x)[1]
  D_x <- dim(K_x)[2]
  N_z <- dim(K_z)[1]
  D_z <- dim(K_z)[2]
  R <- parameters$R
  sigma_e <- parameters$sigma_e
  epsilon <- parameters$epsilon
  lambda_c <- parameters$lambda_c
  lambda_x <- parameters$lambda_x
  lambda_z <- parameters$lambda_z
  iteration <- parameters$iteration
  learn_sigma_e <- parameters$learn_sigma_e
  seed <- parameters$seed
  set.seed(seed)
  
  indices_c <- which(is.na(K_c) == FALSE, arr.ind = TRUE)
  count_c <- length(indices_c) / 2
  indices_x <- which(is.na(K_x) == FALSE, arr.ind = TRUE)
  count_x <- length(indices_x) / 2
  indices_z <- which(is.na(K_z) == FALSE, arr.ind = TRUE)
  count_z <- length(indices_z) / 2
  
  gamma_x <- 1
  gamma_z <- 1
  gamma_eta <- 1
  
  P_x <- rnorm(D_x * R)
  P_z <- rnorm(D_z * R)
  dim(P_x) <- c(D_x, R)
  dim(P_z) <- c(D_z, R)
  
  Q_x <- project_to_stiefel_manifold(P_x)
  E_x <- X %*% Q_x 
  Q_z <- project_to_stiefel_manifold(P_z)
  E_z <- Z %*% Q_z 
  
  DE_c <- (as.matrix(pdist(E_x, E_z)))^2
  KE_c <- exp(-DE_c / sigma_e^2)
  DE_x <- (as.matrix(pdist(E_x, E_x)))^2
  KE_x <- exp(-DE_x / sigma_e^2)
  DE_z <- (as.matrix(pdist(E_z, E_z)))^2
  KE_z <- exp(-DE_z / sigma_e^2)
  objective_c <- (sum((KE_c[indices_c] - K_c[indices_c])^2)) / count_c
  objective_x <- (sum((KE_x[indices_x] - K_x[indices_x])^2)) / count_x
  objective_z <- (sum((KE_z[indices_z] - K_z[indices_z])^2)) / count_z
  objective <- lambda_c * objective_c + lambda_x * objective_x + lambda_z * objective_z
  
  iter <- 0
  current_time <- Sys.time()
  # elapsing time in seconds
  timediff_sec <- as.double(difftime(current_time, start_time, units = "secs"))
  print(paste("iteration=", iter, "objective=", objective, "elapsed time=", timediff_sec))
  
  while (TRUE) {
    iter <- iter + 1
    return_x <- 0
    return_z <- 0
    return_eta <- 0
    
    Q_x_gradient <- matrix(0L, nrow = D_x, ncol = R) 
    for (i in 1:N_x) {
      j <- !is.na(K_c[i,])
      for (s in 1:R) {
        if (sum(j) < 2) {
          Q_x_gradient[, s] <- t(t(Q_x_gradient[, s])) - 4 * lambda_c * as.vector((KE_c[i, j] * (KE_c[i, j] - K_c[i, j]) ) * t(t(E_x[i, s] - t(E_z[j, s]))) ) * t(t(X[i, ])) / sigma_e^2 / count_c
        } else {
          Q_x_gradient[, s] <- t(t(Q_x_gradient[, s])) - 4 * lambda_c * rowSums(as.vector(KE_c[i, j] * (KE_c[i, j] - K_c[i, j]) ) * t(t(E_x[i, s] - t(E_z[j, s]))) ) * t(t(X[i, ])) / sigma_e^2 / count_c
        }
      }
    }
    
    if (lambda_x != 0) {
      for (i in 1:N_x) {
        j <- !is.na(K_x[i,])
        for (s in 1:R) {
          if (sum(j) == 1) {
            Q_x_gradient[, s] <- c(Q_x_gradient[, s] - 4 * lambda_x * (as.vector(KE_x[i, j] * (KE_x[i, j] - K_x[i, j]) * t(t(E_x[i, s] - t(E_x[j, s]))) ) * t(t(X[i, ] - t(X[j, ]))) )  / sigma_e^2 / count_x)
          } else {
            fp <- as.vector(KE_x[i, j] * (KE_x[i, j] - K_x[i, j]) * t(t(E_x[i, s] - t(E_x[j, s]))) )
            sp <- t(t(X[i, ] - t(X[j, ])))
            tot <- matrix(fp, nrow = nrow(sp), ncol = ncol(sp), byrow = TRUE)
            Q_x_gradient[, s] <- Q_x_gradient[, s] - 4 * lambda_x * rowSums(tot * sp) / sigma_e^2 / count_x
          }
        }
      }
    }
    
    Q_x_gradient <- (Q_x %*% t(Q_x_gradient) %*% Q_x - Q_x_gradient)
    Q_x_gradient_norm <- sum(diag(t(Q_x_gradient) %*% (diag(D_x) - 0.5 * (Q_x %*% t(Q_x))) %*% Q_x_gradient))
    if (sqrt(Q_x_gradient_norm) < epsilon) {
      return_x <- 1
      objective_c <- c(objective_c, objective_c[length(objective_c)])
      objective_x <- c(objective_x, objective_x[length(objective_x)])
      objective_z <- c(objective_z, objective_z[length(objective_z)])
      objective <- c(objective, objective[length(objective)])
    } else {
      while (TRUE) {
        Q_x_new <- Q_x + 2 * gamma_x * Q_x_gradient
        Q_x_new <- project_to_stiefel_manifold(Q_x_new)
        E_x_new <- X %*% Q_x_new
        DE_c <- (as.matrix(pdist(E_x_new, E_z)))^2
        KE_c <- exp(-DE_c / sigma_e^2)
        DE_x <- (as.matrix(pdist(E_x_new, E_x_new)))^2
        KE_x <- exp(-DE_x / sigma_e^2)
        objective_new <- lambda_c * sum((KE_c[indices_c] - K_c[indices_c])^2) / count_c + lambda_x * sum((KE_x[indices_x] - K_x[indices_x])^2) / count_x + lambda_z * sum((KE_z[indices_z] - K_z[indices_z])^2) / count_z
        if (gamma_x == 0) {
          break
        } else if (objective[length(objective)] - objective_new >= gamma_x * Q_x_gradient_norm) {
          gamma_x <- 2 * gamma_x
        } else {
          break
        }
      }
      
      while (TRUE){
        Q_x_new <- Q_x + gamma_x * Q_x_gradient
        Q_x_new <- project_to_stiefel_manifold(Q_x_new)
        E_x_new <- X %*% Q_x_new
        DE_c <- (as.matrix(pdist(E_x_new, E_z)))^2
        KE_c <- exp(-DE_c / sigma_e^2)
        DE_x <- (as.matrix(pdist(E_x_new, E_x_new)))^2
        KE_x <- exp(-DE_x / sigma_e^2)
        objective_new <- lambda_c * sum((KE_c[indices_c] - K_c[indices_c])^2) / count_c + lambda_x * sum((KE_x[indices_x] - K_x[indices_x])^2) / count_x + lambda_z * sum((KE_z[indices_z] - K_z[indices_z])^2) / count_z
        if (gamma_x == 0) {
          break
        } else if (objective[length(objective)] - objective_new < 0.5 * gamma_x * Q_x_gradient_norm) {
          gamma_x <- 0.5 * gamma_x
        } else {
          break
        }
      }
      
      Q_x <- Q_x + gamma_x * Q_x_gradient
      Q_x <- project_to_stiefel_manifold(Q_x)
      E_x <- X %*% Q_x
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
      
      current_time <- Sys.time()
      # elapsing time in seconds
      timediff_sec_new <- as.double(difftime(current_time, start_time, units = "secs"))
      timediff_sec <- c(timediff_sec, timediff_sec_new)
      print(paste("iteration=", iter, "objective=", objective_last, "norm_x=", sqrt(Q_x_gradient_norm), "gamma_x=", gamma_x, "elapsed time=", timediff_sec_new))
    }
    
    Q_z_gradient <- matrix(0L, nrow = D_z, ncol = R) 
    for (j in 1:N_z) {
      i <- !is.na(K_c[, j])
      for (s in 1:R) {
        if (sum(i) == 1) {
          Q_z_gradient[, s] <- c(Q_z_gradient[, s] - 4 * lambda_c * sum(as.vector(KE_c[i, j] * (KE_c[i, j] - K_c[i, j]) * (E_z[j, s] - E_x[i, s])) ) %*% t(Z[j, ]) / sigma_e^2 / count_c)
        }
        else {
          Q_z_gradient[, s] <- Q_z_gradient[, s] - 4 * lambda_c * sum(as.vector(KE_c[i, j] * (KE_c[i, j] - K_c[i, j]) * (E_z[j, s] - E_x[i, s])) ) %*% t(Z[j, ]) / sigma_e^2 / count_c
        }
      }
    }
    
    if (lambda_z != 0){
      for (i in 1:N_z){
        j <- !is.na(K_z[i,])
        for (s in 1:R) {
          if (sum(j) == 1) {
            Q_z_gradient[, s] <- c(Q_z_gradient[, s] - 4 * lambda_z * (as.vector(KE_z[i, j] * (KE_z[i, j] - K_z[i, j]) * t(t(E_z[i, s] - t(E_z[j, s]))) ) * t(t(Z[i, ] - t(Z[j, ]))) ) / sigma_e^2 / count_z)
          } else {
            fp <- as.vector(KE_z[i, j] * (KE_z[i, j] - K_z[i, j]) * t(t(E_z[i, s] - t(E_z[j, s]))) )
            sp <- t(t(Z[i, ] - t(Z[j, ])))
            tot <- matrix(fp, nrow = nrow(sp), ncol = ncol(sp), byrow = TRUE)
            Q_z_gradient[, s] <- Q_z_gradient[, s] - 4 * lambda_z * rowSums(tot * sp) / sigma_e^2 / count_z
          }
        }
      }
    }
    
    Q_z_gradient <- (Q_z %*% t(Q_z_gradient) %*% Q_z - Q_z_gradient)
    Q_z_gradient_norm <- sum(diag(t(Q_z_gradient) %*% (diag(D_z) - 0.5 * (Q_z %*% t(Q_z))) %*% Q_z_gradient))
    if (sqrt(Q_z_gradient_norm) < epsilon) {
      return_z <- 1
      objective_c <- c(objective_c, objective_c[length(objective_c)])
      objective_x <- c(objective_x, objective_x[length(objective_x)])
      objective_z <- c(objective_z, objective_z[length(objective_z)])
      objective <- c(objective, objective[length(objective)])
    } else {
      
      while (TRUE) {
        Q_z_new <- Q_z + 2 * gamma_z * Q_z_gradient
        Q_z_new <- project_to_stiefel_manifold(Q_z_new)
        E_z_new <- Z %*% Q_z_new
        DE_c <- (as.matrix(pdist(E_x, E_z_new)))^2
        KE_c <- exp(-DE_c / sigma_e^2)
        DE_z <- (as.matrix(pdist(E_z_new, E_z_new)))^2
        KE_z <- exp(-DE_z / sigma_e^2)
        objective_new <- lambda_c * sum((KE_c[indices_c] - K_c[indices_c])^2) / count_c + lambda_x * sum((KE_x[indices_x] - K_x[indices_x])^2) / count_x + lambda_z * sum((KE_z[indices_z] - K_z[indices_z])^2) / count_z
        if (gamma_z == 0) {
          break
        } else if (objective[length(objective)] - objective_new >= gamma_z * Q_z_gradient_norm) {
          gamma_z <- 2 * gamma_z
        } else {
          break
        }
      }
      
      while (TRUE) {
        Q_z_new <- Q_z + gamma_z * Q_z_gradient
        Q_z_new <- project_to_stiefel_manifold(Q_z_new)
        E_z_new <- Z %*% Q_z_new
        DE_c <- (as.matrix(pdist(E_x, E_z_new)))^2
        KE_c <- exp(-DE_c / sigma_e^2)
        DE_z <- (as.matrix(pdist(E_z_new, E_z_new)))^2
        KE_z <- exp(-DE_z / sigma_e^2)
        objective_new <- lambda_c * sum((KE_c[indices_c] - K_c[indices_c])^2) / count_c + lambda_x * sum((KE_x[indices_x] - K_x[indices_x])^2) / count_x + lambda_z * sum((KE_z[indices_z] - K_z[indices_z])^2) / count_z
        if (gamma_z == 0) {
          break
        } else if (objective[length(objective)] - objective_new < 0.5 * gamma_z * Q_z_gradient_norm) {
          gamma_z <- 0.5 * gamma_z
        } else {
          break
        }
      }
      
      Q_z <- Q_z + gamma_z * Q_z_gradient
      Q_z <- project_to_stiefel_manifold(Q_z)
      E_z <- Z %*% Q_z
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
      
      current_time <- Sys.time()
      # elapsing time in seconds
      timediff_sec_new <- as.double(difftime(current_time, start_time, units = "secs"))
      timediff_sec <- c(timediff_sec, timediff_sec_new)
      print(paste("iteration=", iter, "objective=", objective_last, "norm_z=", sqrt(Q_z_gradient_norm), "gamma_z=", gamma_z, "elapsed time=", timediff_sec_new))
      
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
  # set the state list
  state <- list(Q_x = Q_x, Q_z = Q_z, sigma_e = sigma_e, objective_c = objective_c, objective_x = objective_x, objective_z = objective_z, objective = objective, time = timediff_sec)
}

project_to_stiefel_manifold <- function(Q) {
  s <- svd(Q, nu = nrow(Q), nv = ncol(Q))
  Q <- s$u %*% diag(nrow = dim(Q)[1], ncol = dim(Q)[2]) %*% t(s$v)
}
