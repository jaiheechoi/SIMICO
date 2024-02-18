# Left Survival Function
surv_left <- function(l, d, temp_beta, dmats, lt, rt, r1, k){
  
  #left design matrix
  left_dmat <- dmats$left_dmat
  
  # get beta values for chosen outcome
  covcol <- ncol(left_dmat)
  sub_beta <- temp_beta[((l-1)*covcol + 1): (l *covcol)]
  
  # get sigma squared value from parameter list
  sigmasq <- temp_beta[k*covcol + 1]
  
  # get censoring values
  tpos_ind <- as.numeric(lt > 0)
  obs_ind <- as.numeric(rt != Inf)
  
  #Calculate left survival times
  hl1 <- as.numeric(exp(left_dmat %*% t(matrix(sub_beta, nrow = 1)) + sqrt(2 * sigmasq)* r1[d]))
  sl1 <- ifelse(tpos_ind == 0, 1, exp(-hl1))
  
  #return the survival terms
  return(sl1)
  
}

surv_right <- function(l, d, temp_beta, dmats, lt, rt, r1, k){
  
  # right design matrix
  right_dmat <- dmats$right_dmat
  
  # get beta values for chosen outcome
  covcol <- ncol(right_dmat)
  sub_beta <- temp_beta[((l-1)*covcol + 1): (l *covcol)]
  
  # get sigma squared value from parameter list
  sigmasq <- temp_beta[k*covcol + 1]
  
  # get censoring values
  tpos_ind <- as.numeric(lt > 0)
  obs_ind <- as.numeric(rt != Inf)
  
  #Calculate right survival times
  hr1 <- as.numeric(exp(right_dmat %*% t(matrix(sub_beta, nrow = 1)) + sqrt(2 * sigmasq)* r1[d]))
  sr1 <- ifelse(obs_ind == 0, 0, exp(-hr1))
  sr1[!is.finite(sr1)] <- 0
  
  return(sr1)
  
}

# Left Hazard Function
haz_left <- function(l, d, temp_beta, dmats, lt, rt, r1, k){
  
  # left design matrix
  left_dmat <- dmats$left_dmat
  
  # get beta values for chosen outcome
  covcol <- ncol(left_dmat)
  sub_beta <- temp_beta[((l-1)*covcol + 1): (l *covcol)]
  
  # get sigma squared value from parameter list
  sigmasq <- temp_beta[k*covcol + 1]
  
  
  # get censoring values
  tpos_ind <- as.numeric(lt > 0)
  obs_ind <- as.numeric(rt != Inf)
  
  # calculate left hazard terms
  hl1 <- as.numeric(exp(left_dmat %*% t(matrix(sub_beta, nrow = 1)) + sqrt(2 * sigmasq)* r1[d]))
  return(hl1)
}


# Right Hazard Function
haz_right <- function(l, d, temp_beta, dmats, lt, rt, r1, k){
  
  #right design matrix
  right_dmat <- dmats$right_dmat
  
  # get beta values for chosen outcome
  covcol <- ncol(right_dmat)
  sub_beta <- temp_beta[((l-1)*covcol + 1): (l *covcol)]
  
  # get sigma squared value from parameter list
  sigmasq <- temp_beta[k*covcol + 1]
  
  # get censoring values
  tpos_ind <- as.numeric(lt > 0)
  obs_ind <- as.numeric(rt != Inf)
  
  # calculate right hazard times
  hr1 <- as.numeric(exp(right_dmat %*% t(matrix(sub_beta, nrow = 1)) + sqrt(2 * sigmasq)* r1[d]))
  
  return(hr1)
  
}


# Surv diff
surv_diff <- function(l,d,temp_beta, dmats, lt, rt, r1, k){
  
  (surv_left(l, d, temp_beta, dmats, lt, rt, r1, k) - surv_right(l, d, temp_beta, dmats, lt, rt, r1, k) )
  
}


without_one_phen <- function(l, k, store){
  
  #if total two outcomes
  if(k == 2){
    
    out <- store[,-l,]
    
    return(out)
    
    # case for more than two outcomes
  } else {
    
    #make index of outcomes not including the one of interest
    idx <- (1:k)[-l]
    
    # subset the data
    sub <- store[,idx,]
    
    # multiply the survival differences across the other outcomes (idx)
    sub_prod <- apply(sub, c(1,3), prod)
    
    return(sub_prod)
  }
  
}


without_two_phen <- function(l,m, k, store, n, d){
  
  # get index of outcomes for all not equal to outcomes l and m
  idx <- (1:k)[-c(l,m)]
  
  # subset the array of differences
  sub <- array(store[,idx,], dim=c (n,length(idx),d))
  
  # multiply the survival differences across the other outcomes (idx)
  sub_prod <- apply(sub, c(1,3), prod)
  
}


get_A <- function(store, weights, d, n){
  
  #number of outcomes
  k <- ncol(store[,,d])
  
  # multiply the survival differences across number of outcomes
  mult_across_k <- array(apply(store, c(1,3), prod), dim = c(n,1,d))
  
  # add quadrature weights
  add_weights <- t(matrix(mult_across_k, nrow = n)) * weights
  
  # sum terms and divide by sqrt(pi)
  A_i <- colSums(add_weights)/sqrt(pi)
  
  return(A_i)
  
}

# Simulate gmat
sim_gMat <- function(n,q,rho, maxMP){
  
  ## Construct a binary correlation matrix
  cmat <- toeplitz(c(1, rep(rho, q - 1)))
  meanparam <- runif(q, .01, maxMP)
  x <- bindata::rmvbin(n, margprob = meanparam, bincorr = cmat)  + bindata::rmvbin(n, margprob = meanparam, bincorr = cmat)
  
  return(x)
}


################################################# F3 NR TERMS
fd_term <- function(l, temp_beta, phen,d, apply_diffs, A_i, no_l_all,HL_array, HR_array){
  
  # quadrature weights and roots
  ghDat <- fastGHQuad::gaussHermiteData(d)
  w1 <- ghDat$w
  r1 <- ghDat$x
  
  
  # Get the survival (exp(-exp(eta))) / hazard (exp(eta)) terms
  # for given l, for all 100 weights
  # Result is a list of 100, each n x 1
  hl_d <- HL_array[,,l]
  hr_d <- HR_array[,,l]
  
  
  # prod l not equal to k, SL - SR
  mult_k_without_l <- no_l_all[,,l]
  
  
  # just the first derivative term times the weight
  # weight_d*(S(L)(-H(L))U - S(R)(-H(R))V)
  first_deriv <- function(l, d, sl_d, sr_d, hl_d, hr_d, phen){
    
    
    # get design matrices
    left_dmat <- phen$dmats$left_dmat
    right_dmat <- phen$dmats$right_dmat
    
    #left and right times + censoring
    lt <- phen$lt
    rt <- phen$rt
    tpos_ind <- as.numeric(lt > 0)
    obs_ind <- as.numeric(rt != Inf)
    
    
    # first derivative terms
    U1 <- left_dmat * ifelse(tpos_ind == 0, 0, (exp(-hl_d[,d]) * -hl_d[,d]))
    U2 <- right_dmat * ifelse(obs_ind == 0, 0, (exp(-hr_d[,d]) * -hr_d[,d] ))
    
    # check to make sure there are no NAs
    U2[is.na(U2)] <- 0
    
    # the whole term is a difference between the left values and right values
    inside <- U1 - U2
    
    return(inside)
  }
  
  
  #Get apply for all 100 weights
  #Result is a list of 100, each 1000 x 5
  insides <- lapply(1:d, first_deriv, l = l, hl_d = hl_d, hr_d = hr_d, phen = phen)
  
  
  #Make function that multiplies the prod of surv-diffs and the inside
  mult_together <- function(d, arrayA, listB, weights){
    
    # multiply all the terms together
    arrayA[,d]* listB[[d]] * weights[d]
  }
  
  
  #make for all 100 weights
  deriv_prod <- lapply(1:d, mult_together, arrayA = mult_k_without_l, listB = insides, weights = w1)
  
  # Make the list to an array
  dp_array <- simplify2array(deriv_prod)
  
  
  # Sum over D
  # output is 1000 x 5
  sum_over_d <- apply(dp_array, c(1,2), sum)
  
  
  # Combine with other values and sum over n
  # result is 1 x 5
  fd <- apply(((1/sqrt(pi)) *(sum_over_d/A_i)), 2, sum)
  
  
  return(fd)
  
}



sd_on <- function(l, k, temp_beta, phen, d, apply_diffs, A_i, no_l_all, HL_array, HR_array)
{
  # quadrature weights and roots
  ghDat <- fastGHQuad::gaussHermiteData(d)
  w1 <- ghDat$w
  r1 <- ghDat$x
  
  
  # prod l not equal to k, SL - SR
  no_l <- no_l_all[,,l]
  
  
  # left and righ design matrices
  
  left_dmat <- phen$dmats$left_dmat
  right_dmat <- phen$dmats$right_dmat
  
  #left and right times + censoring terms
  lt <- phen$lt
  rt <- phen$rt
  tpos_ind <- as.numeric(lt > 0)
  obs_ind <- as.numeric(rt != Inf)
  
  
  # survival terms
  hl_d <- HL_array[,,l]
  hr_d <- HR_array[,,l]
  
  
  # second derivative term
  get_sd <- function(hl_d, hr_d, phen, d, no_l){
    
    # left and right design matrices
    left_dmat <- phen$dmats$left_dmat
    right_dmat <- phen$dmats$right_dmat
    
    
    #left and right times + censoring terms
    lt <- phen$lt
    rt <- phen$rt
    tpos_ind <- as.numeric(lt > 0)
    obs_ind <- as.numeric(rt != Inf)
    
    
    
    # first derivative terms
    ul_1 <- ifelse(tpos_ind == 0, 0, -hl_d[,d] * exp(-hl_d[,d]) + (hl_d[,d]^2 * exp(-hl_d[,d])))
    ur_1 <- ifelse(obs_ind == 0, 0, -hr_d[,d] * exp(-hr_d[,d]) + (hr_d[,d]^2 * exp(-hr_d[,d])))
    ur_1[which(is.na(ur_1))] <- 0
    
    # second derivative terms
    sd_term1 <- t(left_dmat) %*% ( (no_l[,d] * as.numeric(ul_1/A_i)) * left_dmat)
    sd_term2 <- t(right_dmat) %*% ( (no_l[,d]* as.numeric(ur_1)/A_i) * right_dmat )
    
    
    # difference between left and right term multiplied by GQ weights
    sd_5x5 <- (sd_term1 - sd_term2) * w1[d]
    
    return(sd_5x5)
  }
  
  # apply the derivatives to each node of the quadrature d
  derivs <- lapply(1:d, get_sd, hl_d = hl_d, hr_d = hr_d, phen = phen, no_l = no_l)
  
  # Make the list to an array
  derivs_array <- simplify2array(derivs)
  
  # Sum over n and divide by sqrt(pi)
  term1 <- apply(derivs_array, c(1,2), sum)/sqrt(pi)
  
  
  ################################################
  # term 2
  
  first_deriv <- function(l, d, hl_d, hr_d, phen){
    
    # left and right design matrices
    left_dmat <- phen$dmats$left_dmat
    right_dmat <- phen$dmats$right_dmat
    
    
    #left and right times + censoring terms
    lt <- phen$lt
    rt <- phen$rt
    tpos_ind <- as.numeric(lt > 0)
    obs_ind <- as.numeric(rt != Inf)
    
    #first derivative terms
    U1 <- left_dmat * ifelse(tpos_ind == 0, 0, (exp(-hl_d[,d]) * -hl_d[,d]/A_i))
    U2 <- right_dmat * ifelse(obs_ind == 0, 0, (exp(-hr_d[,d]) * -hr_d[,d]/A_i ))
    U2[is.na(U2)] <- 0
    
    # whole term is the difference between left and right terms
    inside <- U1 - U2
    
    return(inside)
    
  }
  
  #Get apply for all 100 weights
  #Result is a list of 100, each 1000 x 5
  insides <- lapply(1:d, first_deriv, l = l, hl_d = hl_d, hr_d = hr_d, phen = phen)
  
  
  #Make function that multiplies the prod of surv-diffs and the inside
  mult_together <- function(d, arrayA, listB, weights){
    arrayA[,d]* listB[[d]] * weights[d]
  }
  
  
  #make for all 100 weights
  deriv_prod <- lapply(1:d, mult_together, arrayA = no_l, listB = insides, weights = w1)
  
  
  # Make the list to an array
  dp_array <- simplify2array(deriv_prod)
  
  
  # Sum over D
  # output is 1000 x 5
  term2_a <- apply(dp_array, c(1,2), sum)/sqrt(pi)
  
  # full term is the term multiplied to itself
  term2 <- (t(term2_a) %*% term2_a)
  
  
  #subtract the two terms
  sd_on <- term1 - term2
  
  return(sd_on)
  
}



sd_off <- function(l, m, phen_l, phen_m, temp_beta, d, apply_diffs,A_i, HL_array, HR_array, no_l_all, no_two_all, tpos_all, obs_all,k){
  
  # quadrature weights and roots
  ghDat <- fastGHQuad::gaussHermiteData(d)
  w1 <- ghDat$w
  r1 <- ghDat$x
  
  # get the left and right hazard for observation l
  hl_d <- HL_array[,,l]
  hr_d <- HR_array[,,l]
  
  #left and right design matrices for observation l
  ld_l <- phen_l$dmats$left_dmat
  rd_l <- phen_l$dmats$right_dmat
  
  #left and right design matrices for observation m
  ld_m <- phen_m$dmats$left_dmat
  rd_m <- phen_m$dmats$right_dmat
  
  
  #First derivative terms
  U1 <- sweep((exp(-HL_array) * -HL_array), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  U2 <- sweep((exp(-HR_array) * -HR_array), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  U2[is.na(U2)] <- 0
  
  # product of diff of survival without l
  surv_no_l <- no_l_all
  
  # term looks different for k = 2 observations vs more than 2 observations
  if(k == 2){
    sg_sd <- function(d, l, m) {
      
      #second derivative term
      w1[d]* ( t( ld_l/A_i *  U1[,d,l] -  rd_l/A_i* U2[,d,l] ) %*%
                 ( ld_m* U1[,d,m] -  rd_m *U2[,d,m]  ) )
      
    }
    
    # apply to all quadrature nodes
    term1 <- apply(simplify2array(lapply(1:d, sg_sd, l = l, m = m)), c(1,2), sum)/sqrt(pi)
    
  } else {
    
    # Product of survival terms without two phenotypes
    # get combination of all the indices of outcomes
    combs <- combn(1:k, 2)
    
    # order the observation indices
    min_k <- min(l,m)
    max_k <- max(l,m)
    
    # get the product of the differences of survival without observation l and m
    no_l_m <- no_two_all[,,which(combs[1,] == min_k  & combs[2,] == max_k)]
    
    # Function for the second deriv
    sg_sd <- function(d, l, m) {
      
      # term for the second derivative
      w1[d]* ( t( ld_l * (no_l_m[,d]/A_i * U1[,d,l]) -  rd_l*(no_l_m[,d]/A_i* U2[,d,l]) ) %*%
                 ( ld_m* U1[,d,m] -  rd_m *U2[,d,m]  ) )
      
    }
    
    #apply first derivative term to all
    term1 <- apply(simplify2array(lapply(1:d, sg_sd, l = l, m = m)), c(1,2), sum)/sqrt(pi)
  }
  
  # the first derivative of observation l
  fd_term_l <- function(d, l){
    
    w1[d] * (ld_l*U1[,d,l] - rd_l*U2[,d,l]) * (no_l_all[,d,l]/A_i)
    
  }
  
  #the first derivative of observation m
  fd_term_m <- function(d, m){
    
    w1[d] * (ld_m*U1[,d,m] - rd_m*U2[,d,m]) * (no_l_all[,d,m]/A_i)
    
  }
  
  # apply the first derivative functions to all the quadrature nodes
  term2_a <- apply(simplify2array(lapply(1:d, fd_term_l, l = l)), c(1,2), sum)/sqrt(pi)
  term2_b <- apply(simplify2array(lapply(1:d, fd_term_m, m = m)), c(1,2), sum)/sqrt(pi)
  
  # Multiply the terms together, matching the index for weight d
  term2 <- (t(term2_a) %*% term2_b)
  
  # combine all the terms for the second derivative
  sd_off <- term1 - term2
  
  return(sd_off)
  
}


############################################### F4 Sigma function
# SIGMA



ss_fd <- function(l, phen, HL_array, HR_array, tpos_all, obs_all, apply_diffs, temp_beta, A_i, no_l_all, k, d){
  
  # quadrature weights and roots
  ghDat <- fastGHQuad::gaussHermiteData(d)
  w1 <- ghDat$w
  r1 <- ghDat$x
  
  # define sigma term
  nocol <- ncol(phen$dmats$right_dmat)
  sigmasq <- temp_beta[k*nocol + 1]
  
  #calculate the first derivative terms
  U1 <- sweep((exp(-HL_array) * -HL_array), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  U2 <- sweep((exp(-HR_array) * -HR_array), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  U2[is.na(U2)] <- 0
  
  # get the product of the differences of survival for all k not outcome of interest
  surv_no_l <- no_l_all
  
  # multiply the terms together and sum across observations
  fd_out <- apply((U1 - U2) * surv_no_l, c(1,2), sum)
  
  
  # get the total first derivative term
  deriv <- sum(rowSums(sweep(fd_out, 2, w1 * r1/sqrt(2*sigmasq), FUN = "*"))/A_i)/sqrt(pi)
  
  return(deriv)
  
}


ss_sd <- function(HL_array, HR_array, xAll, apply_diffs, temp_beta, A_i, no_l_all, no_two_all, k, d){
  
  # quadrature weights and roots
  ghDat <- fastGHQuad::gaussHermiteData(d)
  w1 <- ghDat$w
  r1 <- ghDat$x
  
  # define sigma term
  nocol <- ncol(xAll$xDats[[1]]$dmats$right_dmat)
  sigmasq <- temp_beta[k*nocol + 1]
  
  # get censoring terms
  tpos_all <- xAll$ts_all
  obs_all <- xAll$ob_all
  
  # get first derivative terms
  U1 <- sweep((exp(-HL_array) * -HL_array), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  U2 <- sweep((exp(-HR_array) * -HR_array), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  U2[is.na(U2)] <- 0
  
  
  # make function for the third term of the derivative
  # Looks different for k = 2 vs more than 2 outcomes
  term_3_rep <- function(l){
    
    # get index of all outcomes not including outcome of interest
    idx <- (1:k)[-l]
    
    # get combination of all possible pairs of outcomes
    combs <- combn(1:k, 2)
    
    # make function to apply to above indices
    to_idx <- function(m){
      
      # order outcomes
      min_k <- min(l,m)
      max_k <- max(l,m)
      
      if(k == 2){
        out <- (U1-U2)[,,l] * (U1-U2)[,,m]
      } else {
        out <- (U1-U2)[,,l] * (U1-U2)[,,m] * no_two_all[,,which(combs[1,] == min_k  & combs[2,] == max_k)]
      }
      return(out)
    }
    
    
    # apply function to all the index of other outcomes
    apply_idx <- simplify2array(lapply(idx, to_idx))
    
    #sum across all observations
    return(apply(apply_idx, c(1,2), sum))
    
  }
  
  # apply above function to all observations
  term_three <- simplify2array(lapply(1:k, term_3_rep))
  
  # get product of differences between survival terms for all but observation of interest
  surv_no_l <- no_l_all
  
  # calculate second derivative terms
  Y1 <- sweep(-HL_array * exp(-HL_array) + (HL_array^2 * exp(-HL_array)), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  Y2 <- sweep(-HR_array * exp(-HR_array) + (HR_array^2 * exp(-HR_array)), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  Y2[is.na(Y2)] <- 0
  
  # get the first derivative terms
  S1 <- sweep((exp(-HL_array) * -HL_array), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  S2 <- sweep((exp(-HR_array) * -HR_array), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  S2[is.na(S2)] <- 0
  
  # combine all the terms to get full second derivative
  score <- sweep((U1 - U2) * surv_no_l, 2, -r1/(2*sqrt(2*sigmasq)), FUN = "*") +
    sweep((Y1 - Y2) * surv_no_l, 2, (r1/sqrt(2 * sigmasq))^2, FUN = "*") +
    sweep(term_three,2, (r1/sqrt(2 * sigmasq))^2, FUN = "*")
  
  # sum across all observations
  out_sumk <- apply(score, c(1,2), sum)
  
  # sum across node number and multiply by weights
  sum_d <- apply((t(out_sumk) * w1), 2, sum)
  
  # sum across all subjects and divide by sqrt(pi)
  ss_sd_t1 <- sum(sum_d/A_i)/sqrt(pi)
  
  
  # get the second term in the derivative
  ss_sd_t2 <- sum((apply(t(apply((S1 - S2)* surv_no_l, c(1,2), sum))*w1 * r1/sqrt(2*sigmasq),2,sum)/A_i/sqrt(pi))^2)
  
  # combine the two terms
  return(ss_sd_t1 - ss_sd_t2)
}




st_off<- function(l, HL_array, HR_array, xAll, apply_diffs, temp_beta, A_i, no_l_all, no_two_all, k, d) {
  
  # quadrature weights and roots
  ghDat <- fastGHQuad::gaussHermiteData(d)
  w1 <- ghDat$w
  r1 <- ghDat$x
  
  # for phenotype l, get incides of other observations
  idx <- (1:k)[-l]
  
  # subset the data for just the observation of interest
  phen <- xAll$xDats[[l]]
  
  # get censoring terms
  tpos_all <- xAll$ts_all
  obs_all <- xAll$ob_all
  
  # left and right design matrices
  ldm <- phen$dmats$left_dmat
  rdm <- phen$dmats$right_dmat
  
  # define sigma term
  nocol <- ncol(ldm)
  sigmasq <- temp_beta[k*nocol + 1]
  
  #calculate the second derivative terms
  Y1 <- sweep(-HL_array * exp(-HL_array) + (HL_array^2 * exp(-HL_array)), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  Y2 <- sweep(-HR_array * exp(-HR_array) + (HR_array^2 * exp(-HR_array)), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  Y2[is.na(Y2)] <- 0
  
  # calculate the first derivative terms
  U1 <- sweep((exp(-HL_array) * -HL_array), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  U2 <- sweep((exp(-HR_array) * -HR_array), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  U2[is.na(U2)] <- 0
  
  # product of diff of survival without l
  surv_no_l <- no_l_all
  
  # second derivative function
  # term is different for k
  st_sd <- function(d) {
    if(k == 2){
      out <- (ldm *(surv_no_l[,d,l] * Y1[,d,l])) - (rdm * (surv_no_l[,d,l]* Y2[,d,l]) ) +
        ( (ldm *U1[,d,l]) - (rdm *U2[,d,l] ) ) * (U1 - U2)[,d,idx]
    } else{
      #Product of survival terms without two phenotypes
      combs <- combn(1:k, 2)
      with_l_idx <- which(apply(combs, 2, function(x) which(x == l)) == 1 | apply(combs, 2, function(x) which(x == l)) == 2)
      surv_diff_no_two <- no_two_all[,,with_l_idx]
      out <- (ldm *(surv_no_l[,d,l] * Y1[,d,l])) - (rdm * (surv_no_l[,d,l]* Y2[,d,l]) ) +
        ( (ldm *U1[,d,l]) - (rdm *U2[,d,l] ) ) * apply((U1 - U2)[,d,idx]* surv_diff_no_two[,d,], 1, sum)
      
    }
    return(out)
  }
  
  
  # Apply to all the node points
  to_d <- apply(sweep(simplify2array(lapply(1:d, st_sd)), 3, w1 * (r1/sqrt(2*sigmasq)), FUN = "*"), c(1,2), sum)
  
  
  # first term, sum over d and divide by sqrt(pi)
  st_t1 <- colSums(sweep(to_d, 1, A_i, FUN = "/"))/sqrt(pi)
  
  # First deriv times the data mat
  term3_func <- function(d){
    (ldm * U1[,d,l] * surv_no_l[,d,l] ) - (rdm * U2[,d,l] * surv_no_l[,d,l] )
  }
  
  
  # THe term that is like the first deriv of the theta term
  t2a <- sweep(apply(sweep(simplify2array(lapply(1:d, term3_func)), 3, w1, FUN = "*"), c(1,2), sum),1, A_i, FUN = "/")/sqrt(pi)
  
  # The term that is like the first deriv of the sigma term
  t2b <- rowSums(sweep(apply((U1 - U2) * surv_no_l/A_i, c(1,2), sum),2, w1*r1/sqrt(2*sigmasq), FUN = "*"))/sqrt(pi)
  
  
  # Sum over n after multiplying both terms
  st_t2 <- colSums(t2a *t2b)
  
  # combine both terms to get the full derivative
  deriv <- st_t1 - st_t2
  
  return(deriv)
}


################################################### F5 Gamma terms
gamma_fd <- function(l, HL_array, HR_array, tpos_all, obs_all, temp_beta, A_i, no_l_all, gMat, a1, a2, d){
  
  # quadrature weights and roots
  ghDat <- fastGHQuad::gaussHermiteData(d)
  w1 <- ghDat$w
  r1 <- ghDat$x
  
  #Get weights from MAF
  MAF <- apply(gMat, 2, function(x) mean(x, na.rm = T)/2)
  
  #Weight MAF using beta distribution
  beta <- dbeta(MAF, a1, a2, ncp = 0, log = FALSE)
  
  #Apply weights to matrix of genotypes
  Z_w = t(t(gMat) * (beta))
  
  # calculate first derivative term
  U1 <- sweep((exp(-HL_array) * -HL_array), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  U2 <- sweep((exp(-HR_array) * -HR_array), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  U2[is.na(U2)] <- 0
  
  # difference between survival terms multiplied across observations not including k
  no_l <- no_l_all[,,l]
  
  # mutliply elements together to get the first term of the derivative
  t1 <- t(Z_w) %*% ((U1 - U2)[,,l] * no_l/A_i)
  
  #multiply quadrature weights and sum across subjects
  deriv <- apply(sweep(t1, 2, w1, FUN = "*"), 1, sum)/sqrt(pi)
  
  return(deriv)
}


gamma_on <- function(l, HL_array, HR_array, tpos_all, obs_all, temp_beta, A_i, no_l_all, gMat, a1, a2, d){
  
  # quadrature weights and roots
  ghDat <- fastGHQuad::gaussHermiteData(d)
  w1 <- ghDat$w
  r1 <- ghDat$x
  
  #Get weights from MAF
  MAF <- apply(gMat, 2, function(x) mean(x, na.rm = T)/2)
  
  #Weight MAF using beta distribution
  beta <- dbeta(MAF, a1, a2, ncp = 0, log = FALSE)
  
  #Apply weights to matrix of genotypes
  Z_w = t(t(gMat) * (beta))
  
  # calculate second derivative terms
  Y1 <- sweep(-HL_array * exp(-HL_array) + (HL_array^2 * exp(-HL_array)), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  Y2 <- sweep(-HR_array * exp(-HR_array) + (HR_array^2 * exp(-HR_array)), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  Y2[is.na(Y2)] <- 0
  
  # function for the first term
  term1 <- function(d,l){
    (t(Z_w * (Y1 - Y2)[,d,l] *no_l_all[,d,l]/A_i)) %*% Z_w
  }
  
  # apply to the quadrature notes
  gt_t1 <- apply(sweep(simplify2array(lapply(1:d, term1, l = l)),3, w1, FUN = "*"), c(1,2), sum)/sqrt(pi)
  
  # calculate the first derivative term
  U1 <- sweep((exp(-HL_array) * -HL_array), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  U2 <- sweep((exp(-HR_array) * -HR_array), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  U2[is.na(U2)] <- 0
  
  # term multiplying the difference of the survival terms for all observations other than k
  no_l <- no_l_all[,,l]
  
  # function combining all the terms
  fd_t <- function(l, d){
    out <- sweep(Z_w, 1, (U1 - U2)[,d,l] * no_l[,d]/A_i, FUN = "*")
    return(out)
  }
  
  # apply to all the quadrature nodes
  to_d <- apply(sweep(simplify2array(lapply(1:d, fd_t, l = l)), 3, w1, FUN = "*"), c(1,2), sum)/sqrt(pi)
  
  #multiply term to itself
  gt_t2 <- t(t(to_d) %*% to_d)
  
  return(gt_t1 - gt_t2)
  
}



# dl/dgamma_k d_gamma_j
# 50 x 50
gamma_off <- function(l, m, HL_array, HR_array, tpos_all, obs_all, temp_beta, A_i, no_l_all, no_two_all, gMat, a1, a2, k, d){
  #Get weights from MAF
  MAF <- apply(gMat, 2, function(x) mean(x, na.rm = T)/2)
  
  # quadrature weights and roots
  ghDat <- fastGHQuad::gaussHermiteData(d)
  w1 <- ghDat$w
  r1 <- ghDat$x
  
  #Weight MAF using beta distribution
  beta <- dbeta(MAF, a1, a2, ncp = 0, log = FALSE)
  
  #Apply weights to matrix of genotypes
  Z_w = t(t(gMat) * (beta))
  
  # Product of survival terms without two phenotypes
  
  # first derivative term
  U1 <- sweep((exp(-HL_array) * -HL_array), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  U2 <- sweep((exp(-HR_array) * -HR_array), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  U2[is.na(U2)] <- 0
  
  # function for the first term combining all the elements
  # different for k = 2 vs more than 2 observations
  term1 <- function(d,l,m){
    if(k == 2){
      out <- (t(Z_w * (U1 - U2)[,d,l] *(U1 - U2)[,d,m]/A_i)) %*% Z_w
    } else{
      # combination of all the outcome indices
      combs <- combn(1:k, 2)
      
      # order the observation indices
      min_k <- min(l,m)
      max_k <- max(l,m)
      
      # product of survival differences without observation l and m
      no_l_m <- no_two_all[,,which(combs[1,] == min_k  & combs[2,] == max_k)]
      out <- (t(Z_w * (U1 - U2)[,d,l] *(U1 - U2)[,d,m] * no_l_m[,d]/A_i)) %*% Z_w
    }
    
    return(out)
  }
  
  # apply to all quadrature nodes
  gt_t1 <- apply(sweep(simplify2array(lapply(1:d, term1, l = l, m = m)),3, w1, FUN = "*"), c(1,2), sum)/sqrt(pi)
  
  
  # function same as first derivative term
  fd_t <- function(l, d){
    out <- sweep(Z_w, 1, (U1 - U2)[,d,l] * no_l_all[,d,l]/A_i, FUN = "*")
    return(out)
  }
  
  # apply to all quadrature nodes, multiply by weight, sum over all d
  # apply for both observation l and observation m
  to_d_l <- apply(sweep(simplify2array(lapply(1:d, fd_t, l = l)), 3, w1, FUN = "*"), c(1,2), sum)/sqrt(pi)
  to_d_m <- apply(sweep(simplify2array(lapply(1:d, fd_t, l = m)), 3, w1, FUN = "*"), c(1,2), sum)/sqrt(pi)
  
  # multiply together
  gt_t2 <- t(t(to_d_l) %*% to_d_m)
  
  return(gt_t1 - gt_t2)
}


gammatheta <- function(l, HL_array, HR_array, tpos_all, obs_all, apply_diffs, temp_beta, A_i, xDats, no_l_all, gMat, a1, a2, d){
  
  #Get weights from MAF
  MAF <- apply(gMat, 2, function(x) mean(x, na.rm = T)/2)
  
  #Weight MAF using beta distribution
  beta <- dbeta(MAF, a1, a2, ncp = 0, log = FALSE)
  
  #Apply weights to matrix of genotypes
  Z_w = t(t(gMat) * (beta))
  
  # product of differences of survival terms excluding observation l
  no_l <- no_l_all[,,l]
  
  # design matrices
  phen = xDats[[l]]
  left_dmat <- phen$dmats$left_dmat
  right_dmat <- phen$dmats$right_dmat
  
  # left and right survival terms
  lt <- phen$lt
  rt <- phen$rt
  
  # quadrature weights and roots
  ghDat <- fastGHQuad::gaussHermiteData(d)
  w1 <- ghDat$w
  r1 <- ghDat$x
  
  # calculate the second derivative terms
  Y1 <- sweep(-HL_array * exp(-HL_array) + (HL_array^2 * exp(-HL_array)), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  Y2 <- sweep(-HR_array * exp(-HR_array) + (HR_array^2 * exp(-HR_array)), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  Y2[is.na(Y2)] <- 0
  
  # function for total second derivative term
  get_sg <- function(l, HL_array, HR_array, phen, d, no_l_all, Z_w){
    left_dmat <- phen$dmats$left_dmat
    right_dmat <- phen$dmats$right_dmat
    lt <- phen$lt
    rt <- phen$rt
    
    #left term
    sg_t1 <-  left_dmat * (no_l_all[,d,l]/A_i * as.numeric(Y1[,d,l]))
    
    #right term
    sg_t2 <-  right_dmat * (no_l_all[,d,l]/A_i * as.numeric(Y2[,d,l]))
    
    # multiply difference of left and right term times weighted genetic matrix
    out <- t(Z_w) %*% (sg_t1 - sg_t2)
    
    return(out)
  }
  
  # apply to all quadrature nodes
  derivs <- sweep(simplify2array(lapply(1:d, get_sg, l = l, HL_array = HL_array, HR_array = HR_array, phen = phen, no_l_all = no_l_all, Z_w = Z_w)), 3, w1, FUN = "*")
  
  # Sum over D
  term1 <- apply(derivs, c(1,2), sum)/sqrt(pi)
  
  #### now calculate term 2
  
  # calculate first derivative terms
  U1 <- sweep((exp(-HL_array) * -HL_array), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  U2 <- sweep((exp(-HR_array) * -HR_array), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  U2[is.na(U2)] <- 0
  
  # calculate the first part of the second term
  term2_a <- function(d,l){
    sweep(Z_w, 1, (U1 - U2)[,d,l] * no_l_all[,d,l]/A_i, FUN = "*")
  }
  
  # apply to all quadrature nodes, multiply weights, and sum over d
  t2a_tod <- apply(sweep(simplify2array(lapply(1:d, term2_a, l = l)),3,w1,FUN = "*"), c(1,2), sum)/sqrt(pi)
  
  # function for the third term of the derivative
  term3_func <- function(d,l){
    phen = xDats[[l]]
    left_dmat <- phen$dmats$left_dmat
    right_dmat <- phen$dmats$right_dmat
    
    # combine left and right terms
    (left_dmat * U1[,d,l] * no_l_all[,d,l]) - (right_dmat * U2[,d,l] * no_l_all[,d,l] )
  }
  
  # The term that is like the first deriv of the theta term
  t2b_tod <- sweep(apply(sweep(simplify2array(lapply(1:d, term3_func, l = l)), 3, w1, FUN = "*"), c(1,2), sum), 1, A_i, FUN = "/")/sqrt(pi)
  
  # mulply the two parts together
  term2 <- (t(t2a_tod) %*% t2b_tod)
  
  return(t(term1 - term2))
  
}


gammatheta_off <- function(l,m, HL_array, HR_array, xAll, apply_diffs, temp_beta, A_i, no_l_all, no_two_all, gMat, a1, a2, k, d){
  
  # quadrature weights and roots
  ghDat <- fastGHQuad::gaussHermiteData(d)
  w1 <- ghDat$w
  r1 <- ghDat$x
  
  # design matrices
  phen = xAll$xDats[[l]]
  left_dmat <- phen$dmats$left_dmat
  right_dmat <- phen$dmats$right_dmat
  
  # censor terms
  tpos_all <- xAll$ts_all
  obs_all <- xAll$ob_all
  
  # first derivative term
  U1 <- sweep((exp(-HL_array) * -HL_array), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  U2 <- sweep((exp(-HR_array) * -HR_array), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  U2[is.na(U2)] <- 0
  
  #Get weights from MAF
  MAF <- apply(gMat, 2, function(x) mean(x, na.rm = T)/2)
  
  #Weight MAF using beta distribution
  beta <- dbeta(MAF, a1, a2, ncp = 0, log = FALSE)
  
  #Apply weights to matrix of genotypes
  Z_w = t(t(gMat) * (beta))
  
  # product of differences of left and right survival excluding observation l
  no_l <- no_l_all[,,l]
  
  # second derivative term
  # different for k = 2 vs more than k outcomes
  get_sg_off <- function(l, m, HL_array, HR_array, phen, d, Z_w){
    
    # design matrices and left and right times
    left_dmat <- phen$dmats$left_dmat
    right_dmat <- phen$dmats$right_dmat
    lt <- phen$lt
    rt <- phen$rt
    
    if(k == 2) {
      # left_term of the first term
      sg_t1a <-  left_dmat * (as.numeric(U1[,d,l])/A_i)
      
      # right term of the first term
      sg_t1b <-  right_dmat * (as.numeric(U2[,d,l])/A_i)
      
    } else{
      combs <- combn(1:k, 2)
      min_k <- min(l,m)
      max_k <- max(l,m)
      no_l_m <- no_two_all[,,which(combs[1,] == min_k  & combs[2,] == max_k)]
      
      #left_term of first term
      sg_t1a <-  left_dmat * (no_l_m[,d]/A_i * as.numeric(U1[,d,l]))
      
      #right term of first term
      sg_t1b <-  right_dmat * (no_l_m[,d]/A_i * as.numeric(U2[,d,l]))
    }
    
    # second term
    sg_t2 <- Z_w *(U1[,d,m] - U2[,d,m])
    
    # combine all the terms
    out <-t(sg_t2)%*% (sg_t1a - sg_t1b)
    
    return(out)
    
  }
  
  # apply first derivative function to all quadrature nodes
  derivs <- simplify2array(lapply(1:d, get_sg_off, l = l, m = m, HL_array = HL_array, HR_array = HR_array, phen = phen, Z_w = Z_w))
  
  # Sum over D
  term1 <- apply(sweep(derivs, 3, w1, FUN = '*'), c(1,2), sum)/sqrt(pi)
  
  #### now term 2
  
  # first derivative terms
  U1 <- sweep((exp(-HL_array) * -HL_array), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  U2 <- sweep((exp(-HR_array) * -HR_array), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  U2[is.na(U2)] <- 0
  
  
  # first part of second term
  term2_a <- function(d,l){
    
    # design matrices
    left_dmat <- phen$dmats$left_dmat
    right_dmat <- phen$dmats$right_dmat
    
    # combine left and right terms
    sweep((left_dmat * U1[,d,l] - right_dmat*U2[,d,l]), 1, no_l_all[,d,l]/A_i, FUN = "*")
  }
  
  # apply to all quadrature nodes, sum over d, multiply by weights
  t2a_tod <- apply(sweep(simplify2array(lapply(1:d, term2_a, l = l)),3,w1,FUN = "*"), c(1,2), sum)/sqrt(pi)
  
  # function for second part of second term
  term2_b <- function(d,m){
    sweep(Z_w, 1, (U1 - U2)[,d,m] * no_l_all[,d,m]/A_i, FUN = "*")
  }
  
  # apply to all quadrature nodes, sum over d, multiply by weights
  t2b_tod <- apply(sweep(simplify2array(lapply(1:d, term2_b, m = m)),3,w1,FUN = "*"), c(1,2), sum)/sqrt(pi)
  
  # combine both parts to get second term
  term2 <- t(t(t2a_tod) %*% t2b_tod)
  
  return(t(term1 - term2))
}


# dl/dgamma_k dsigma
# 50 x 1
gammasigma <- function(l, HL_array, HR_array, tpos_all, obs_all, apply_diffs, temp_beta, A_i, xDats, no_l_all, no_two_all, gMat, a1, a2, k, d) {
  
  # quadrature weights and roots
  ghDat <- fastGHQuad::gaussHermiteData(d)
  w1 <- ghDat$w
  r1 <- ghDat$x
  
  #Get weights from MAF
  MAF <- apply(gMat, 2, function(x) mean(x, na.rm = T)/2)
  
  #Weight MAF using beta distribution
  beta <- dbeta(MAF, a1, a2, ncp = 0, log = FALSE)
  
  #Apply weights to matrix of genotypes
  Z_w = t(t(gMat) * (beta))
  
  # for phenotype l
  idx <- (1:k)[-l]
  
  # design matrices
  phen <- xDats[[l]]
  left_dmat <- phen$dmats$left_dmat
  right_dmat <- phen$dmats$right_dmat
  
  # get sigma term
  nocol <- ncol(left_dmat)
  sigmasq <- temp_beta[k*nocol + 1]
  
  # calculate second derivative terms
  Y1 <- sweep(-HL_array * exp(-HL_array) + (HL_array^2 * exp(-HL_array)), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  Y2 <- sweep(-HR_array * exp(-HR_array) + (HR_array^2 * exp(-HR_array)), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  Y2[is.na(Y2)] <- 0
  
  # calculate first derivative terms
  U1 <- sweep((exp(-HL_array) * -HL_array), c(1,3), ifelse(tpos_all == 0, 0,1), FUN = "*")
  U2 <- sweep((exp(-HR_array) * -HR_array), c(1,3), ifelse(obs_all == 0, 0,1), FUN = "*")
  U2[is.na(U2)] <- 0
  
  # product of diff of survival without l
  surv_no_l <- no_l_all
  
  # Function for the second deriv
  sg_sd <- function(d, l) {
    
    if(k == 2){
      out <- Z_w*( (surv_no_l[,d,l] * Y1[,d,l]) - (surv_no_l[,d,l]* Y2[,d,l]) ) +
        ( Z_w* (U1[,d,l] - U2[,d,l] ) ) * (U1 - U2)[,d,idx]
    } else {
      # Product of survival terms without two phenotypes
      combs <- combn(1:k, 2)
      no_l_m <- no_two_all[,,-which(combs[1,] != l & combs[2,] != l)]
      out <- Z_w*( (surv_no_l[,d,l] * Y1[,d,l]) - (surv_no_l[,d,l]* Y2[,d,l]) ) +
        ( Z_w* (U1[,d,l] - U2[,d,l] ) ) * apply((U1 - U2)[,d,idx]* no_l_m[,d,], 1, sum)
    }
    return(out)
  }
  
  # Apply to all the node points
  to_d <- apply(sweep(simplify2array(lapply(1:d, sg_sd, l = l)), 3, w1 * (r1/sqrt(2*sigmasq)), FUN = "*"), c(1,2), sum)
  
  # sum over d to get first term
  st_t1 <- colSums(sweep(to_d, 1, A_i, FUN = "/"))/sqrt(pi)
  
  # function for first part of second term
  term2_a <- function(d,l){
    sweep(Z_w, 1, (U1 - U2)[,d,l] * no_l_all[,d,l]/A_i, FUN = "*")
  }
  
  # apply to all quadrature nodes and sum
  t2a <- apply(sweep(simplify2array(lapply(1:d, term2_a, l = l)),3,w1,FUN = "*"), c(1,2), sum)/sqrt(pi)
  
  # The term that is like the first deriv of the sigma term
  t2b <- rowSums(sweep(apply((U1 - U2) * surv_no_l/A_i, c(1,2), sum),2, w1*r1/sqrt(2*sigmasq), FUN = "*"))/sqrt(pi)
  
  # Sum over n after multiplying both terms
  st_t2 <- colSums(t2a *t2b)
  
  # combine the two terms
  deriv <- st_t1 - st_t2
  
  return(matrix(deriv, nrow = 1))
  
}


# newton raphson, get p-value, generate data

simico_fit_null <- function(init_beta, epsilon, xDats, lt_all, rt_all, k, d) {
  
  # number of observations
  n = nrow(xDats[[1]]$dmats$right_dmat)
  
  # number of covariates
  nocol <- ncol(xDats[[1]]$dmats$right_dmat)
  
  # create matrix of censoring
  tpos_all <- matrix(NA, nrow = n, ncol = k)
  obs_all <- matrix(NA, nrow = n, ncol = k)
  
  # get censoring terms
  for(j in 1:k){
    tpos_all[,j] <- as.numeric(lt_all[,j] > 0)
    obs_all[,j] <- as.numeric(rt_all[,j] != Inf)
  }
  
  # make data in correct form
  threedmat <- list()
  for(i in 1:k){
    threedmat[[i]] <- list(dmats = xDats[[i]]$dmats, lt = lt_all[,i], rt = rt_all[,i])
  }
  
  # compile all elements
  xAll <- list(xDats = threedmat, ts_all = tpos_all, ob_all = obs_all)
  
  # get values for format of function input
  xDats <- xAll$xDats
  t_all <- xAll$ts_all
  o_all <- xAll$ob_all
  
  # quadrature weights and roots
  ghDat <- fastGHQuad::gaussHermiteData(d)
  w1 <- ghDat$w
  r1 <- ghDat$x
  
  # conditions for convergence
  iter = 1
  diff = 1
  
  # inital beta values
  temp_beta <- matrix(init_beta, nrow = 1)
  
  # loop for newton raphson
  while(diff > epsilon & iter < 200){
    
    # make left and right survival and hazard values
    HL_array <- array (c (NA, NA), dim=c (n,d,k))
    HR_array <- array (c (NA, NA), dim=c (n,d,k))
    
    # calculate the values of the array
    for(i in 1:k){
      # get left times, right times, and dmats for each outcome
      lt = lt_all[,i]
      rt = rt_all[,i]
      dmats = xDats[[i]]$dmats
      
      hl_d <- lapply(1:d, haz_left, l = i, temp_beta = temp_beta, dmats = dmats, lt = lt, rt = rt, r1 = r1, k =  k)
      hr_d <- lapply(1:d, haz_right, l = i, temp_beta = temp_beta, dmats = dmats, lt = lt, rt = rt, r1 = r1, k = k)
      HL_array[,,i] <- simplify2array(hl_d)
      HR_array[,,i] <- simplify2array(hr_d)
    }
    
    # calculate the differences between survival values
    apply_diffs <- array (c (NA, NA), dim=c (n,k,d))
    
    for(i in 1:k){
      
      # subset data for observation number
      lt = lt_all[,i]
      rt = rt_all[,i]
      dmats = xDats[[i]]$dmats
      
      # cacluate the values and fill in the array
      apply_diffs[,i,] <- simplify2array(lapply(1:d,surv_diff, l = i, temp_beta = temp_beta, dmats = dmats, lt = lt, rt = rt, r1 = r1, k = k))
    }
    
    # get denominator of values
    A_i <- get_A(apply_diffs, w1, 100, n)
    
    # make sure there are no zeroes
    A_i[which(A_i == 0)]<- min(A_i[which(A_i!= 0)])
    
    # calculate the product of the difference between the survival terms excluding outcome l
    no_l_all <- simplify2array(lapply(1:k, without_one_phen, k = k, store=apply_diffs))
    
    # get the combination of all outcome indices
    combs <- combn(1:k, 2)
    
    # calculate the product of the difference between the survival terms excluding two outcomes
    if(k == 2){
      no_two_all <- 1
    } else {
      no_two_all <- array(data = NA, dim = c(n,d, choose(k,2)))
      for(i in 1:choose(k,2)){
        no_two_all[,,i] <- without_two_phen(combs[1,i], combs[2,i], k, apply_diffs, n, d)
        
      }
    }
    
    # generalizable way to gt the gradient
    grad <- c()
    
    # loop through all outcomes to get the full gradient
    for(i in 1:k){
      temp_grad <- fd_term(i, temp_beta, xDats[[i]], d, apply_diffs = apply_diffs, A_i =A_i, no_l_all = no_l_all, HL_array = HL_array, HR_array = HR_array)
      
      # combine all the gradients
      grad <- c(grad, temp_grad)
      
    }
    
    # calculate the gradient of the sigma squared term
    grad_ss <- ss_fd(1, xDats[[1]], HL_array, HR_array, t_all, o_all, apply_diffs = apply_diffs, temp_beta = temp_beta, A_i =A_i, no_l_all = no_l_all, k = k, d = d)
    
    #combine all the gradient terms
    grad <- matrix(c(grad, grad_ss), ncol = 1)
    
    # total number of rows and columns for the information matrix
    totd <- (nocol*k) + 1
    
    # start building the information matrix
    jmat <- matrix(NA, nrow = totd, ncol = totd)
    
    # start filling in the matrix
    for(i in 1:k){
      
      # get the correct indices
      d1 <- (nocol*(i-1))+ 1
      d2 <- nocol* i
      
      # calculate the second derivatives
      jmat[d1:d2, d1:d2] <- sd_on(i, k, temp_beta, xDats[[i]], d, apply_diffs = apply_diffs, A_i =A_i, no_l_all = no_l_all, HL_array = HL_array, HR_array = HR_array)
      jmat[d1:d2, totd] <- st_off(i, HL_array, HR_array, xAll, apply_diffs = apply_diffs, temp_beta = temp_beta, A_i =A_i, no_l_all = no_l_all, no_two_all = no_two_all, k = k, d = d)
      jmat[totd, d1:d2] <- t(st_off(i, HL_array, HR_array, xAll, apply_diffs = apply_diffs, temp_beta = temp_beta, A_i =A_i, no_l_all = no_l_all, no_two_all = no_two_all, k = k, d = d))
      
      # get the off diagonal terms
      idx <- (1:k)[-i]
      
      # fill in the matrix
      for(j in 1:length(idx))
      {
        od1 <- (nocol*(idx[j]-1))+ 1
        od2 <- nocol* idx[j]
        
        # second derivative terms
        jmat[d1:d2, od1:od2] <- sd_off(i,idx[j], phen_l = xDats[[i]], phen_m = xDats[[idx[j]]], temp_beta, d = d, apply_diffs = apply_diffs, A_i =A_i, HL_array = HL_array, HR_array = HR_array, no_l_all = no_l_all, no_two_all = no_two_all, tpos_all = t_all, obs_all = o_all, k = k)
        jmat[od1:od2, d1:d2] <- sd_off(idx[j],i, phen_l = xDats[[idx[j]]], phen_m = xDats[[i]], temp_beta, d = d, apply_diffs = apply_diffs, A_i =A_i, HL_array = HL_array, HR_array = HR_array, no_l_all = no_l_all, no_two_all = no_two_all, tpos_all = t_all, obs_all = o_all, k = k)
        
      }
    }
    
    # second derivative of the sigma term
    jmat[totd, totd] <- ss_sd(HL_array, HR_array, xAll, apply_diffs = apply_diffs, temp_beta = temp_beta, A_i =A_i, no_l_all = no_l_all, no_two_all = no_two_all, k = k, d = d)
    
    # newton raphson
    beta_new <- temp_beta - t(grad) %*% solve(jmat)
    
    # calculate the difference
    diff = (-t(grad) %*% solve(jmat)) %*% t(-t(grad) %*% solve(jmat))
    print(diff)
    # update the iteration
    iter = iter + 1
    
    # update the new beta
    temp_beta <- matrix(beta_new, nrow = 1)
    
  }
  
  return(list(beta_fit = temp_beta, iter = iter, diff = diff, jmat = jmat, grad = grad))
  
}


############### get p value
simico_out <- function(xMat, lt_all, rt_all, a1, a2, G, k, d, numKnots= 1, checkAIC = FALSE){
  if(!is.null(numKnots) & checkAIC == TRUE){
    message("Cannot force number of knots and run checkAIC. Please specify a correct option.")
    invokeRestart("abort")
  }

  p = ncol(xMat)
  n = nrow(xMat)
  
  # quadrature weights and roots
  ghDat <- fastGHQuad::gaussHermiteData(d)
  w1 <- ghDat$w
  r1 <- ghDat$x
  
  # create matrix of censoring
  ts_all <- matrix(NA, nrow = n, ncol = k)
  ob_all <- matrix(NA, nrow = n, ncol = k)
  
  # get censoring terms
  for(j in 1:k){
    ts_all[,j] <- as.numeric(lt_all[,j] > 0)
    ob_all[,j] <- as.numeric(rt_all[,j] != Inf)
  }
  
  if(checkAIC == TRUE){
    
    knotsOut <- simico_fit_spline(xMat, lt_all, rt_all, k, d)
    
    #find the minimum
    nullFit <- knotsOut$nullFits_list[[which.min(knotsOut$aicAll[,2])]]
    xDats <- knotsOut$xDats_list[[which.min(knotsOut$aicAll[,2])]]
    
    
  } else{
    
    # first, fit the spline to the data
    leftArray <- array(data=NA, dim=c(n, ncol(xMat)+2 +numKnots, k))
    rightArray <- array(data=NA, dim=c(n, ncol(xMat)+2 + numKnots, k))
    
    for (pheno_it in 1:k) {
      tempDmats <- ICSKAT::make_IC_dmat(xMat = xMat, lt = lt_all[, pheno_it],
                                        rt = rt_all[, pheno_it], obs_ind = ob_all[, pheno_it],
                                        tpos_ind = ts_all[, pheno_it], nKnots= numKnots)
      leftArray[, , pheno_it] <- tempDmats$left_dmat
      rightArray[, , pheno_it] <- tempDmats$right_dmat
    }
    
    # put all the variables in the correct form
    xDats <- list()
    for(j in 1:k){
      dmatph <- list(right_dmat = rightArray[,,j], left_dmat = leftArray[,,j])
      xDats[[j]] <- list(dmats = dmatph, lt = lt_all[,j], rt = rt_all[,j])
    }
    
    xAll <- list(xDats = xDats, ts_all = ts_all, ob_all = ob_all)
    
    # number of covariates
    nocol <- ncol(xDats[[1]]$dmats$right_dmat)
    
    # get values for format of function input
    
    init_beta <- c()
    # get a good initial value for the null model fit
    for(outc in 1:length(xDats)){
      outcomeDat <- xDats[[outc]]
      
      # get events as non-right censored events
      event <- ob_all[,outc]
      
      # get midpoint of times
      mdTime <- (rt_all[,outc] + lt_all[,outc])/2
      mdTime[which(event == 0)] <- lt_all[,outc][which(event == 0)]
      
      # fit the kaplan meier estimator
      kmFit <- survival::survfit(Surv(mdTime, event) ~ 1)
      if(kmFit$surv[1] == 1){
        kmFit$surv[which(kmFit$surv == 1)] <- .999
      }

      if(min(kmFit$time) == 0){
        kmFit$time[which(kmFit$time == 0)] <- .001
      }
      
      # fit the linear regression
      lmMod <- lm(log(-log(kmFit$surv)) ~ log(kmFit$time))
      
      # save the coefficients
      onebeta <- c(rep(0,p), as.numeric(lmMod$coefficients), rep(0, numKnots))
      init_beta <- c(init_beta, onebeta)
    }
    
    # add the initial value for the subject specific random variable
    init_beta <- c(init_beta, 1)
    
    #fit the null model
    skip_to_next <- FALSE
    nullFit <- NA
    tryCatch(nullFit <- simico_fit_null(init_beta = init_beta, epsilon = 10^-5,
                           xDats = exampleDat$fullDat$xDats,
                           lt_all = exampleDat$leftTimesMat, rt_all = exampleDat$rightTimesMat,
                           k = k, d = 100), error = function(e){skip_to_next <<- TRUE})

    if(skip_to_next){ next}
  }
  
    # get the beta fits, Itt, and design matrices from spline selection
    temp_beta <- nullFit$beta_fit
    Itt <- nullFit$jmat
    
    # make left and right survival and hazard values
    HL_array <- array (c (NA, NA), dim=c (n,d,k))
    HR_array <- array (c (NA, NA), dim=c (n,d,k))
    
    
    # calculate the values of the array
    for(i in 1:k){
      lt = lt_all[,i]
      rt = rt_all[,i]
      dmats = xDats[[i]]$dmats
      hl_d <- lapply(1:d, haz_left, l = i, temp_beta = temp_beta, dmats = dmats, lt = lt, rt= rt, r1 = r1, k =  k)
      hr_d <- lapply(1:d, haz_right, l = i, temp_beta = temp_beta, dmats = dmats, lt = lt, rt = rt, r1 = r1, k = k)
      HL_array[,,i] <- simplify2array(hl_d)
      HR_array[,,i] <- simplify2array(hr_d)
    }
    
    # calculate the differences between survival values
    apply_diffs <- array (c (NA, NA), dim=c (n,k,d))
    
    for(i in 1:k){
      
      # subset data for observation number
      lt = lt_all[,i]
      rt = rt_all[,i]
      dmats = xDats[[i]]$dmats
      
      # cacluate the values and fill in the array
      apply_diffs[,i,] <- simplify2array(lapply(1:d,surv_diff, l = i, temp_beta = temp_beta, dmats = dmats, lt = lt, rt = rt, r1 = r1, k = k))
  }
    
  # get denominator of values
  A_i <- get_A(apply_diffs, w1, 100, n)
    
  # make sure there are no zeroes
  A_i[which(A_i == 0)]<- min(A_i[which(A_i!= 0)])
    
  # calculate AIC
  # first get log-likelihood
  loglik <- log(sum(A_i))
    
  #get number of parameters
  nParam = (nocol * k) + 1
    
  #calculate AIC
  AIC = (2 * nParam) - (2 * loglik)
  
  # number of covariates
  nocol <- ncol(xDats[[1]]$dmats$right_dmat)
  
  # calculate the product of the difference between the survival terms excluding outcome l
  no_l_all <- simplify2array(lapply(1:k, without_one_phen, k = k, store=apply_diffs))
  
  # get the combination of all outcome indices
  combs <- combn(1:k, 2)
  
  # calculate the product of the difference between the survival terms excluding two outcomes
  if(k == 2){
    no_two_all <- 1
  } else {
    no_two_all <- array(data = NA, dim = c(n,d, choose(k,2)))
    for(i in 1:choose(k,2)){
      no_two_all[,,i] <- without_two_phen(combs[1,i], combs[2,i], k, apply_diffs, n, d)
      
    }
  }
  
  # build the test statistic
  
  # get correct dimension values
  dmatdim <- ncol(xDats[[1]]$dmats$right_dmat)
  gdim <- ncol(G)
  
  # information matrix of gamma and theta terms
  Igt <- matrix(NA, nrow = (dmatdim * k + 1),ncol = k*gdim)
  
  # loop through observations and run derivative functions
  for(i in 1:k){
    
    # gamma theta
    Igt[(((i - 1)*dmatdim) + 1): (i*dmatdim), (((i - 1)*gdim) + 1): (i*gdim)] <- gammatheta(i, HL_array, HR_array, ts_all, ob_all, apply_diffs, temp_beta, A_i, xDats, no_l_all, G, a1, a2, d)
    Igt[(dmatdim * k + 1), (((i - 1)*gdim) + 1): (i*gdim)] <- gammasigma(i, HL_array, HR_array, ts_all, ob_all, apply_diffs, temp_beta, A_i, xDats, no_l_all, no_two_all, G, a1, a2, k, d)
    
    # get indices for off diagonal terms
    idx <- (1:k)[-i]
    
    # loop through index
    for(j in 1:length(idx)){
      
      Igt[(((idx[j] - 1)*dmatdim) + 1): (idx[j]*dmatdim), (((i - 1)*gdim) + 1): (i*gdim)] <- gammatheta_off(idx[j],i, HL_array, HR_array, xAll, apply_diffs, temp_beta, A_i, no_l_all, no_two_all, G, a1, a2, k, d)
      Igt[(dmatdim * k + 1), (((idx[j] - 1)*gdim) + 1): (idx[j]*gdim)] <- gammasigma(idx[j], HL_array, HR_array, ts_all, ob_all, apply_diffs, temp_beta, A_i, xDats, no_l_all, no_two_all, G, a1, a2, k, d)
      
    }
  }
  
  # information matrix of gamma gamma
  Igg <- matrix(NA, nrow = k * gdim, ncol = k*gdim)
  
  for(i in 1:k){
    
    # get correct dimension values
    d1 <- (gdim*(i-1))+ 1
    d2 <- gdim* i
    
    # second derivative of gamma term on the diagonal
    Igg[d1:d2, d1:d2] <- gamma_on(i, HL_array, HR_array, ts_all, ob_all, temp_beta, A_i, no_l_all, G, a1, a2, d)
    
    # get indices for off diagonals
    idx <- (1:k)[-i]
    
    # loop through indices
    for(j in 1:length(idx))
    {
      od1 <- (ncol(G)*(idx[j]-1))+ 1
      od2 <- ncol(G)* idx[j]
      
      # second derivative for gamma terms on the off diagonals
      Igg[d1:d2, od1:od2] <- gamma_off(i, idx[j], HL_array, HR_array, ts_all, ob_all, temp_beta, A_i, no_l_all, no_two_all, G, a1, a2, k, d)
      Igg[od1:od2, d1:d2] <- gamma_off(idx[j], i, HL_array, HR_array, ts_all, ob_all, temp_beta, A_i, no_l_all, no_two_all, G, a1, a2, k, d)
    }
  }
  
  # get the full information matrix
  sigmat <- (-Igg) - t(-Igt) %*% solve(-Itt) %*% (-Igt)
  
  # get the gradient term
  U_g <- c()
  
  # loop through all the observation values
  for(i in 1:k){
    
    temp_grad <- gamma_fd(i, HL_array, HR_array, ts_all, ob_all, temp_beta, A_i, no_l_all, G, a1, a2, d)
    U_g <- c(U_g, temp_grad)
  }
  
  # compute the score statistic
  gamma_score <- t(U_g) %*% U_g
  
  burden_score <- (sum(U_g))^2
  
  # get the eigen values
  lams <- eigen(sigmat)$values
  
  # compute the p-value using davies method
  pval <- CompQuadForm::davies(q=gamma_score, lambda=lams)
  
  B_burden= burden_score / sum(sigmat);
  p_burden= 1 - stats::pchisq(B_burden, df = 1)
  
  # compute omnibus pvalue
  #evenly space out rho values
  prhos <- seq(from = .2, to = .8, by = .2)
  
  #create vector to store all p-values
  omni_pvals <- rep(NA, length(prhos) + 2)
  
  # add burden and multiple p-values to vector (rho = 0 and rho = 1)
  omni_pvals[1] <- pval$Qq
  omni_pvals[length(omni_pvals)] <- p_burden
  
  # vector of 1s
  J <- matrix(1, nrow = length(U_g))
  
  # Loop through all rho to calculate omnibus test statstic at different rho values
  for(rh in 1:length(prhos)){
    tempscore <- (1 - prhos[rh]) * gamma_score + prhos[rh] * burden_score
    half <- chol((1 - prhos[rh]) * diag(length(U_g)) + prhos[rh] * (J %*% t(J)))
    templam <- eigen(half %*% sigmat %*% t(half))$values
    omni_pvals[i + 1] <- CompQuadForm::davies(q=tempscore, lambda=templam)$Qq
  }
  
  # sometimes there are NAs
  omni_pvals <- omni_pvals[-which(is.na(omni_pvals))]
  
  # if at least one value, ACAT all the p-values together
  if(length(omni_pvals) ==0){omniP <- NA} else {omniP <- ACAT(omni_pvals)}
  
  return(list(multQ = gamma_score, multP = pval$Qq, burdQ = burden_score, burdP = p_burden, omniP = omniP, AIC = AIC))
}

#function for fitting the spline
simico_fit_spline <- function(xMat, lt_all, rt_all, k, d){
  
  # number of covariates
  p = ncol(xMat)
  
  # quadrature weights and roots
  ghDat <- fastGHQuad::gaussHermiteData(d)
  w1 <- ghDat$w
  r1 <- ghDat$x

    # loop through 1 through 4 number of internal knots
    knots <- 1:4
    
    # store the aic values
    aicAll <- matrix(NA, nrow = 4, ncol = 2)
    
    # save null fits and design matrices
    nullFits_list <-list()
    xDats_list <- list()
    
    # loop through number of knots
    for(i in 1:length(knots)){
      
      #set number of knots
      numKnots <- knots[i]
      
      # create matrix of censoring
      ts_all <- matrix(NA, nrow = n, ncol = k)
      ob_all <- matrix(NA, nrow = n, ncol = k)
      
      # get censoring terms
      for(j in 1:k){
        ts_all[,j] <- as.numeric(lt_all[,j] > 0)
        ob_all[,j] <- as.numeric(rt_all[,j] != Inf)
      }
      
      # first, fit the spline to the data
      leftArray <- array(data=NA, dim=c(n, ncol(xMat)+2 +numKnots, k))
      rightArray <- array(data=NA, dim=c(n, ncol(xMat)+2 + numKnots, k))
      
      for (pheno_it in 1:k) {
        tempDmats <- ICSKAT::make_IC_dmat(xMat = xMat, lt = lt_all[, pheno_it],
                                  rt = rt_all[, pheno_it], obs_ind = ob_all[, pheno_it],
                                  tpos_ind = ts_all[, pheno_it], nKnots= numKnots)
        leftArray[, , pheno_it] <- tempDmats$left_dmat
        rightArray[, , pheno_it] <- tempDmats$right_dmat
      }
      
      # put all the variables in the correct form
      xDats <- list()
      for(j in 1:k){
        dmatph <- list(right_dmat = rightArray[,,j], left_dmat = leftArray[,,j])
        xDats[[j]] <- list(dmats = dmatph, lt = lt_all[,j], rt = rt_all[,j])
      }
      
      xAll <- list(xDats = xDats, ts_all = ts_all, ob_all = ob_all)
      
      # save in list for this number of knots
      xDats_list[[i]] <- xDats
      
      # number of covariates
      nocol <- ncol(xDats[[1]]$dmats$right_dmat)
      
      # get values for format of function input
      
      init_beta <- c()
      # get a good initial value for the null model fit
      for(outc in 1:length(xDats)){
        outcomeDat <- xDats[[outc]]
        
        # get events as non-right censored events
        event <- ob_all[,outc]
        
        # get midpoint of times
        mdTime <- (rt_all[,outc] + lt_all[,outc])/2
        mdTime[which(event == 0)] <- lt_all[,outc][which(event == 0)]
        
        # fit the kaplan meier estimator
        kmFit <- survival::survfit(Surv(mdTime, event) ~ 1)
        if(max(kmFit$surv) == 1){
          kmFit$surv[which(kmFit$surv == 1)] <- .999
        }
        
         if(min(kmFit$surv) == 0){
           kmFit$surv[which(kmFit$surv == 0)] <- .001
         }
        # 
        # if(max(kmFit$time) == Inf){
        #   kmFit$time[which(kmFit$time == Inf)] <- 100
        # }
        
        # fit the linear regression
        lmMod <- lm(log(-log(kmFit$surv)) ~ log(kmFit$time))
        
        # save the coefficients
        onebeta <- c(rep(0,p), as.numeric(lmMod$coefficients), rep(0, numKnots))
        init_beta <- c(init_beta, onebeta)
      }
      
      # add the initial value for the subject specific random variable
      init_beta <- c(init_beta, 1)
      
      # fit null
      nullFit <- simico_fit_null(init_beta, epsilon = 10^-5, xDats, lt_all, rt_all, k, d)
      
      # save the nullFits
      nullFits_list[[i]] <- nullFit
      
      # get null fit from NR
      temp_beta <- nullFit$beta_fit
      
      # make left and right survival and hazard values
      HL_array <- array (c (NA, NA), dim=c (n,d,k))
      HR_array <- array (c (NA, NA), dim=c (n,d,k))
      
      # calculate the values of the array
      for(j in 1:k){
        lt = lt_all[,j]
        rt = rt_all[,j]
        dmats <- xDats[[j]]$dmats
        hl_d <- lapply(1:d, haz_left, l = j, temp_beta = temp_beta, dmats = dmats, lt = lt, rt = rt, r1 = r1, k =  k)
        hr_d <- lapply(1:d, haz_right, l = j, temp_beta = temp_beta, dmats = dmats, lt = lt, rt = rt, r1 = r1, k = k)
        HL_array[,,j] <- simplify2array(hl_d)
        HR_array[,,j] <- simplify2array(hr_d)
      }
      
      # calculate the differences between survival values
      apply_diffs <- array (c (NA, NA), dim=c (n,k,d))
      
      for(j in 1:k){
        
        dmats <- xDats[[j]]$dmats
        
        # subset data for observation number
        lt = lt_all[,j]
        rt = rt_all[,j]
        
        # cacluate the values and fill in the array
        apply_diffs[,j,] <- simplify2array(lapply(1:d,surv_diff, l = j, temp_beta = temp_beta,dmats = dmats, lt = lt, rt = rt, r1 = r1, k = k))
      }
      
      # get denominator of values
      A_i <- get_A(apply_diffs, w1, 100, n)
      
      # make sure there are no zeroes
      A_i[which(A_i == 0)]<- min(A_i[which(A_i!= 0)])
      
      # calculate AIC
      # first get log-likelihood
      loglik <- log(sum(A_i))
      
      #get number of parameters
      nParam = (nocol * k) + 1
      
      #calculate AIC
      AIC = (2 * nParam) - (2 * loglik)
      
      #record knot number
      aicAll[i,1] <- numKnots
      
      # record AIC value
      aicAll[i,2] <- AIC
    }
  return(list(aicAll = aicAll, nullFits_list = nullFits_list, 
              xDats_list = xDats_list))
}

################################################# F6
#Sample generation
simico_gen_dat <- function(bhFunInv, obsTimes = 1:3, windowHalf = 0.1, probMiss = .1, n, p, k, tauSq, gMatCausal, xMat, effectSizes, oppSign = FALSE, nKnots) {
  
  nocol = p + 2 + nKnots
  # number of subjects and outcomes
  # true model has nothing
  fixedMat <- matrix(data=0, nrow=n, ncol=k)
  
  # get genetic effects
    # get genetic effects
  geneticVec <- c()

  # Calculate the effect size
  for (es in 1:length(effectSizes))
  {
    # create vector for each SNP
    Bk <- rep(NA, ncol(gMatCausal))

    # loop through all SNPs
    for(j in 1:length(Bk)){

      # calculate minor allele frequency
      MAF <- apply(gMatCausal, 2, function(x) mean(x)/2)

      if(oppSign == FALSE){
      # multply effect size and genetic matrix
      Bk[j] = effectSizes[es]* abs(log10(MAF[j]))} else if(oppSign == TRUE) {
      Bk[j] = effectSizes[es]* abs(log10(MAF[j])) * ((-1)^(j + 1))}
    }
    geneticVec <- c(geneticVec, (gMatCausal %*% Bk))
  }
  
  # random unif vector for PIT
  unifVec <- runif(n=n * k)
  
  # vectorize fixed effects, all first phenotype, then all second phenotype, etc.
  fixedVec <- c(fixedMat)
  
  # random intercept
  randomInt <- rnorm(n=n, sd = sqrt(tauSq))
  
  # get full term
  randomIntRep <- rep(randomInt, k)
  
  # add random effect
  etaVec <- fixedVec + randomIntRep + geneticVec
  
  # probability integral transform - assumes PH model
  toBeInv <- -log(1 - unifVec) / exp(etaVec)
  
  # all n*K exact failure times
  exactTimesVec <- bhFunInv(toBeInv)
  
  # all exact failure times for all K phenotypes
  exactTimesMat <- matrix(data=exactTimesVec, nrow=n, ncol=k, byrow = FALSE)
  
  # hold left and right intervals data for all K phenotypes
  leftTimesMat <- matrix(data=NA, nrow=n, ncol=k)
  rightTimesMat <- matrix(data=NA, nrow=n, ncol=k)
  obsInd <- matrix(data=NA, nrow=n, ncol=k)
  tposInd <- matrix(data=NA, nrow=n, ncol=k)
  
    # do visits separately for each phenotype
  nVisits <- length(obsTimes)

  for (pheno_it in 1:k) {
            # 1 - probMiss is the chance of making it to the visit
  madeVisit <- matrix(data = stats::rbinom(n=n*nVisits, size=1, prob=(1 - probMiss)), nrow=n, ncol=nVisits)

  # make sure there is at least one visit for each subject
  nMadeVisits <- apply(madeVisit, 1, sum)
  zeroVisits <- which(nMadeVisits == 0)

  while (length(zeroVisits) > 0) {
    madeVisit[zeroVisits, ] <- matrix(data = stats::rbinom(n=length(zeroVisits) * nVisits, size=1,
                prob=(1 - probMiss)), nrow=length(zeroVisits), ncol=nVisits)

    nMadeVisits <- apply(madeVisit, 1, sum)
    zeroVisits <- which(nMadeVisits == 0)
  }

    # actual visit time is uniformly distributed around the intended obsTime, windowHalf on each side
    visitTime <- sweep(matrix(data = stats::runif(n=n*nVisits, min=-windowHalf, max=windowHalf), nrow=n, ncol=nVisits), MARGIN=2, STATS=obsTimes, FUN="+")

    # get all visits for each subject
    allVisits <- madeVisit * visitTime
    
    # make the interval for each subject
    allInts <- t(mapply(FUN=ICSKAT::createInt, obsTimes = data.frame(t(allVisits)), eventTime=exactTimesMat[, pheno_it]))
    
    leftTimesMat[, pheno_it] <- allInts[, 1]
    rightTimesMat[, pheno_it] <- allInts[, 2]
    
    # event time indicators
    obsInd[, pheno_it] <- ifelse(rightTimesMat[, pheno_it] == Inf, 0, 1)
    tposInd[, pheno_it] <- ifelse(leftTimesMat[, pheno_it] == 0, 0, 1)
  }
  
    leftArray <- array(data=NA, dim=c(n, p + nKnots + 2, k))

  rightArray <- array(data=NA, dim=c(n, p + nKnots + 2, k))


  for (pheno_it in 1:k) {

    tempDmats <- ICSKAT::make_IC_dmat(xMat = xMat, lt = leftTimesMat[, pheno_it],

                                      rt = rightTimesMat[, pheno_it], obs_ind = obsInd[, pheno_it],

                                      tpos_ind = tposInd[, pheno_it], nKnots=nKnots)

    leftArray[, , pheno_it] <- tempDmats$left_dmat

    rightArray[, , pheno_it] <- tempDmats$right_dmat

  }
  
  #make placeholder to change the data later
  # *** this is how the functions take the inputs
  dataph <- matrix(NA, nrow= n, ncol = nocol)
  vecN <- rep(NA, n)
  dmatph <- list(right_dmat = dataph, left_dmat = dataph)
  #xmatph <- list(dmats = dmatph, lt = vecN, rt = vecN)
  xmatph <- list(dmats = dmatph)
  allph <- matrix(NA, nrow = n, ncol = k)
  threedmat <- list()
  for(i in 1:k){
    threedmat[[i]] <- xmatph
  }
  samp <- list(xDats = threedmat, ts_all = allph, ob_all = allph)
  
  for(pheno in 1:k){
    samp$xDats[[pheno]]$dmats$right_dmat <- rightArray[,,pheno]
    samp$xDats[[pheno]]$dmats$left_dmat <- leftArray[,,pheno]
    #samp$xDats[[pheno]]$lt <- leftTimesMat[,pheno]
    #samp$xDats[[pheno]]$rt <- rightTimesMat[,pheno]
  }
  
  samp$ob_all <- obsInd
  samp$ts_all <- tposInd
  
  
  # return
  return(list(exactTimesMat = exactTimesMat, leftTimesMat = leftTimesMat,
              rightTimesMat = rightTimesMat, obsInd = obsInd, tposInd = tposInd, fullDat = samp))
  
}

# Function to subset the gmat by a number of causal SNPs
Get_CausalSNPs_bynum<-function(gMat, num, Causal.MAF.Cutoff){
  #Calculate MAF for the genotypes
  MAF <- apply(gMat, 2, function(x) mean(x)/2)
  
  IDX<-which(MAF < Causal.MAF.Cutoff)
  if(length(IDX) == 0){
    msg<-sprintf("No SNPs with MAF < %f",Causal.MAF.Cutoff)
    stop(msg)
  }
  
  N.causal<-num
  if(N.causal < 1){
    N.causal = 1
  }
  
  re<-sort(sample(IDX,N.causal))
  return(re)
  
}
