#' @keywords internal
tr <- function(A) sum(diag(A))


#' @import statmod
#' @keywords internal
exp_values <- function(Phi, Lambda_s, Psi_s, Psi_s1, cov_s, getdet = FALSE)
{
   k <- dim(Phi)[2]
   I_k <- diag(1, k)
   S <- length(Lambda_s)

   ###defining objects
   j_s <- numeric(S)
   I_j <- list()
   Sig_s <- list()
   ds_s <- list()
   I_tot <- list()
   LambTOT <- list()
   Sig_s1 <- list()
   delta_Lambda <- list()
   delta_Phi <- list()
   Delta_Lambda <- list()
   Delta_Phi <- list()
   Covfcfs <- list()
   Txsfs <- list()
   Txsfcs <- list()
   Tfsfs <- list()
   Tfcsfcs <- list()
   Tfcsfs <- list()

  for (s in 1:S){
   	ds_s[[s]] <- NULL
   	j_s[s] <- c(dim(Lambda_s[[s]])[[2]])
    I_j[[s]] <- diag(1, j_s[s])
    Sig_s[[s]] <- Phi %*% t(Phi) + Lambda_s[[s]] %*% t(Lambda_s[[s]]) + Psi_s[[s]]
    if (getdet)  {ds_s[[s]] <- det(Sig_s[[s]])}
    I_tot[[s]] <- diag(1, k + j_s[s])
    LambTOT[[s]] <- cbind(Phi, Lambda_s[[s]])
    Sig_s1[[s]] <- Psi_s1[[s]] - (statmod::vecmat(diag(Psi_s1[[s]]), LambTOT[[s]]) %*%
                   solve(I_tot[[s]] + (t(LambTOT[[s]]) %*% statmod::vecmat(diag(Psi_s1[[s]]),
                  LambTOT[[s]]))) %*% statmod::matvec(t(LambTOT[[s]]), diag(Psi_s1[[s]])))
    delta_Lambda[[s]] <- t(Lambda_s[[s]]) %*% Sig_s1[[s]]
    delta_Phi[[s]] <- t(Phi) %*% Sig_s1[[s]]
    Delta_Lambda[[s]] <- I_j[[s]] - (t(Lambda_s[[s]]) %*% Sig_s1[[s]] %*% Lambda_s[[s]])
    Delta_Phi[[s]] <- I_k - (t(Phi) %*% Sig_s1[[s]] %*% Phi)
    Covfcfs[[s]] <- -t(Phi) %*% Sig_s1[[s]] %*% Lambda_s[[s]]
    Txsfs[[s]] <- cov_s[[s]] %*% t(delta_Lambda[[s]])
    Txsfcs[[s]] <- cov_s[[s]] %*% t(delta_Phi[[s]])
    Tfsfs[[s]] <- delta_Lambda[[s]] %*% cov_s[[s]] %*% t(delta_Lambda[[s]]) + Delta_Lambda[[s]]
    Tfcsfcs[[s]] <- delta_Phi[[s]] %*% cov_s[[s]] %*% t(delta_Phi[[s]]) + Delta_Phi[[s]]
    Tfcsfs[[s]] <- delta_Phi[[s]] %*% cov_s[[s]] %*% t(delta_Lambda[[s]]) + Covfcfs[[s]]
   }
 return(list(Txsxs = cov_s, Txsfs = Txsfs, Txsfcs = Txsfcs, Tfsfs = Tfsfs,
             Tfcsfcs =  Tfcsfcs, Tfcsfs = Tfcsfs, ds_s=ds_s,  Sig_s1 = Sig_s1))
}




#' Provides some starting values for the parameters of a MSFA model
#'
#' This is a supporting function for \code{ecm_msfa}. The method employed is documented in the reference.
#'
#' The upper-triangular zero constraint is adopted to achieve identification,
#' as detailed in the reference, though the function can also be run without such constraint.
#' @param X_s List of lenght \eqn{S}{S}, corresponding to number of different studies considered.
#' Each element of the list contains a data matrix, with the same number of columns \eqn{P}{P} for all the studies.
#' No standardization is carried out by the function.
#' @param k Number of common factors.
#' @param j_s Number of study-specific factors. A vector of positive integers of length \eqn{S}{S}.
#' @param constraint  Constraint for ensuring identifiability. The default is "block_lower2", which
#' corresponds to the main proposal of De Vito et al. (2018). An alternative identification
#' strategy is triggered by  "block_lower1"; this is more restrictive but may work also with smaller
#' number of variables.
#' @param method Which method should be used to find the starting values? The two possibilities are \code{"adhoc"} for
#' the method described in De Vito et al. (2016), and \code{"fa"} for averaging over separate study-specific FA models.
#' Default is \code{"adhoc"}.
#' @param robust If \code{TRUE}, robust covariance matrix is used in place of the sample covariance. Default
#' is \code{FALSE}.
#' @param corr If \code{TRUE}, the analysis will employ the correlation matrix instead of the covariance matrix.
#' @param mcd If \code{TRUE}, the robust estimator used for the covariance is the same proposed in Pison et al. (2003),
#' otherwise the default value of the function \code{CovRob} of the \code{robust} library is employed. Default is
#' \code{FALSE}.
#' @return A list  containing  \code{Phi},\code{Lambda_s} and  \code{psi_s}, starting values for the model matrices.
#' @import psych
#' @export
#' @references De Vito, R., Bellio, R., Parmigiani, G. and Trippa, L. (2019). Multi-study Factor Analysis. Biometrics,  75, 337-346.
start_msfa <- function(X_s, k, j_s, constraint = "block_lower2", method = "adhoc", robust = FALSE, corr = FALSE, mcd = FALSE)
{
  X_used_s <- X_s
  S <- length(X_s)
  if(corr & !robust)
    for(s in 1:S)  X_used_s[[s]] <- scale(X_s[[s]])
  if(robust & corr & method=="adhoc"){
    for(s in 1:S){
      ogg_s <- if(mcd) covRob(X_s[[s]], estim = "mcd", quan = .75, ntrial = 1000) else covRob(X_s[[s]])
    }
  X_used_s[[s]] <- scale(X_s[[s]], center = ogg_s$center, scale = sqrt(diag(ogg_s$cov)))
  }
  p <- dim(X_s[[1]])[2]
  Phi <- matrix(0, nrow=p, ncol=k)
  Lambda_s <- psi_s <- list()
  if(method=="adhoc"){
     X <- Reduce(rbind, X_used_s)
     X.pcr <- prcomp(X)
     Phi <- matrix(X.pcr$rotation[,1:k], nrow=p, ncol=k, byrow=FALSE)

     if (constraint == "block_lower1") {
     Phi[upper.tri(Phi)] <- 0
     for(s in 1:S){
       iniLS <- array(prcomp(X_used_s[[s]])$rotation, dim=c(p, j_s[s]))
       iniTot <- cbind(Phi, iniLS)
       iniTot[upper.tri(iniTot)] <- 0
       Lambda_s[[s]] <-  matrix(iniTot[,(k+1):(k+j_s[s])], p , j_s[s])
       psi_s[[s]] <- fa(X_used_s[[s]], nfactors = k+j_s[s])$uniq
       		}
      	}

     if (constraint == "block_lower2") {
     Phi[upper.tri(Phi)] <- 0
     for(s in 1:S){
       Lambda_s[[s]] = array(prcomp(X_used_s[[s]])$rotation, dim=c(p, j_s[s]))
       Lambda_s[[s]][upper.tri(Lambda_s[[s]])] <- 0
       psi_s[[s]] <- fa(X_used_s[[s]], nfactors = k+j_s[s])$uniq
      		}
    		}

     if (constraint == "null") {
     Phi <- Phi
     for(s in 1:S){
       Lambda_s[[s]] = array(prcomp(X_used_s[[s]])$rotation, dim=c(p, j_s[s]))
       psi_s[[s]] <- fa(X_used_s[[s]], nfactors = k+j_s[s])$uniq
      		}
    	}
    }
 #### it is important to post-process the output for avoiding sign changes
 if(method=="fa"){
      est <- ecm_fa(X_s, tot_s = k + j_s, robust = robust, mcd = mcd, corr = corr, tol = 10^-5, nIt = 5000, trace = FALSE)
      Phi <- est$Omega_s[[1]][,1:k] / S
      Lambda_s[[1]] <-  est$Omega_s[[1]][,(k+1):(k+j_s[1])]
      psi_s[[1]] <- est$psi_s[[1]]
      for(s in 2:S){
        Phi <- Phi + est$Omega_s[[s]][,1:k] / S * sign(Phi) * sign(est$Omega_s[[s]][,1:k]) ###to avoid sign changes
        Lambda_s[[s]] <-  est$Omega_s[[s]][,(k+1):(k+j_s[s])]
        psi_s[[s]] <- est$psi_s[[s]]
       }
  }
  out <- list(Phi=Phi, Lambda_s=Lambda_s, psi_s=psi_s)
  return(out)
}



#' @keywords internal
loglik_ecm <- function(Sig_s1,  ds_s, n_s, cov_s)
{
   S <- length(n_s)
   #####log likelihood value for each study
   val_s <- c()
   for(s in 1:S){
	    val_s[s] <- - (n_s[s]/2) * log(ds_s[[s]]) - (n_s[s]/2) * tr(Sig_s1[[s]] %*% cov_s[[s]])
	    }
   #####sum of each study-likelihood
   val_tot <- sum(val_s)
   return(val_tot)
}

loglik_msfax = function(Phi,Lambda_s,Gamma,H_s,n_s,cov_s)
{
  S <- length(n_s)
  Sig_s = list()
  for(s in 1:S)
  {
    Sig_s[[s]] = Phi %*% t(Phi) + Lambda_s[[s]] %*% t(Lambda_s[[s]]) + Gamma + H_s[[s]]
  }

  # get all the inverses
  Sig_s1 = lapply(Sig_s,solve)

  # get all the determinants
  ds_s = lapply(Sig_s,det)

  #####log likelihood value for each study
  val_s <- c()
  for(s in 1:S){
    val_s[s] <- - (n_s[s]/2) * log(ds_s[[s]]) - (n_s[s]/2) * tr(Sig_s1[[s]] %*% cov_s[[s]])
  }
  #####sum of each study-likelihood
  val_tot <- sum(val_s)
  return(val_tot)
}

loglik_vanilla <- function(Sig_s1,  ds_s, n_s, cov_s)
{
  val_tot = - (n_s/2) * log(ds_s) - (n_s/2) * tr(Sig_s1 %*% cov_s)
  return(val_tot)
}

#' @keywords internal
param2vect <- function(param, constraint)
{
  p <- length(param$psi_s[[1]])
  S <- length(param$psi_s)
  phi_vals <- as.vector(param$Phi[lower.tri(param$Phi, diag = TRUE)])
  k <- ncol(param$Phi)
  lambda_vals <- psi_vals <- j_s <- c()

  if(constraint=="block_lower1"){
    for(s in 1:S){
      Lambda_s <- param$Lambda_s[[s]][-(1:k),]
      lambda_vals <- c(lambda_vals, as.vector(Lambda_s[lower.tri(Lambda_s, diag = TRUE)]))
      psi_vals <- c(psi_vals, param$psi_s[[s]])
      j_s[[s]] <- ncol(param$Lambda_s[[s]])}
  }
  if(constraint=="block_lower2"){
    for(s in 1:S){
      Lambda_s <- param$Lambda_s[[s]]
      lambda_vals <- c(lambda_vals, as.vector(Lambda_s[lower.tri(Lambda_s, diag = TRUE)]))
      psi_vals <- c(psi_vals, param$psi_s[[s]])
      j_s[[s]] <- ncol(param$Lambda_s[[s]])}
  }
  theta <- c(phi_vals, lambda_vals, psi_vals)
  return(theta)
}






#' @keywords internal
vect2param <- function(vect, param.struct, constraint, p, k, j_s)
{
  S <- length(j_s)
  nP <- k * p - k * ( k - 1) / 2
  if (constraint == "block_lower1")  { nL <- j_s * (p - k)  - j_s *  (j_s - 1) / 2}
  if (constraint == "block_lower2")  { nL <- j_s * p  - j_s *  (j_s - 1) / 2}
  phi_vals <- vect[1:nP]
  Phi <- matrix(0, nrow=p, ncol=k)
  Phi[lower.tri(Phi, diag = TRUE)] <- phi_vals
  Lambda_s <- param.struct$Lambda_s
  psi_s <- param.struct$psi_s
  for(s in 1:S){
    nL_s  <- if(s==1) 0 else sum(nL[1:(s-1)])
    ind <-  (nP + nL_s + 1):(nP + nL_s + nL[s])
    lambda_vals_s <-  vect[ind]
    Lambda_s[[s]][Lambda_s[[s]]!=0] <- lambda_vals_s
    ind_s <- (nP + sum(nL) + p * (s-1) + 1):(nP + sum(nL) + p * s)
    psi_vals_s <- vect[ind_s]
    psi_s[[s]] <-  psi_vals_s
  }
  return(list(Phi=Phi, Lambda_s=Lambda_s, psi_s=psi_s))
}


#' @keywords internal
#### loglikelihood function re-expressed as a function of the model parameters
#### theta: c(Phi, Lambda_1,..,Lambda_S,Psi_1,..,Psi_S)
loglik_int <- function(theta, n_s, cov_s, k, j_s, constraint)
{
  S <- length(n_s)
  p <- ncol(cov_s[[1]])
  nP <- k * p - k * ( k - 1) / 2
  if (constraint == "block_lower1")  { nL <- j_s * (p - k)  - j_s *  (j_s - 1) / 2}
  if (constraint == "block_lower2")  { nL <- j_s * p  - j_s *  (j_s - 1) / 2}
  phi_vals <- theta[1:nP]
  Phi <- matrix(0, p, k)
  Phi[lower.tri(Phi, diag = TRUE)] <- phi_vals
  out <- 0
  for(s in 1:S){
    nL_s  <- if(s==1) 0 else sum(nL[1:(s-1)])
    ind <-  (nP + nL_s + 1):(nP + nL_s + nL[s])
    omega_vals_s <- c(phi_vals, theta[ind])
    Omega_s <- matrix(0, p, k + j_s[s])
    if (constraint == "block_lower1")  Omega_s[lower.tri(Omega_s, diag = TRUE)] <- omega_vals_s
    if (constraint == "block_lower2")
    {
      Lambda_s <- matrix(0, p, j_s[s])
      Lambda_s[lower.tri(Lambda_s, diag = TRUE)] <- theta[ind]
      Omega_s <- cbind(Phi, Lambda_s)
      }
    ind_s <- (nP + sum(nL) + p * (s-1) + 1):(nP + sum(nL) + p * s)
    psi_vals_s <- theta[ind_s]
    Psi_s1 <- diag(1 / psi_vals_s)
    D1L_s <- statmod::vecmat(1 / sqrt(psi_vals_s), Omega_s)
    LDL_s <- crossprod(D1L_s)
    A <- diag(k + j_s[s]) + LDL_s
    A1 <- chol2inv(chol(A))
    D2L_s <- statmod::vecmat(1/psi_vals_s, Omega_s)
    Sig_s1 <- Psi_s1 - D2L_s  %*% A1 %*% t(D2L_s)
    log_ds_s <-  log(det(A)) + sum(log(psi_vals_s))
    out  <- out - (n_s[s]/2) * log_ds_s - (n_s[s]/2) * tr(Sig_s1 %*% cov_s[[s]])
    }
  return(out)
}




#' Variance matrix of MLE estimates for a MSFA model
#'
#' Computes the inverse observed information for a MSFA model
#'
#'
#' Numerical differentiation is employed to obtain the observed information matrix at a
#' given parameter values, so that when the parameter values equals the MLE the function
#' returns the estimated variance matrix of the fitted model. The method is rather inefficient, and
#' it may lead to long computations, though the function is designed to be called only once after the
#' estimation has been carried out. However, it would be relatively straightforward to employ analytical
#' differentiation at least for the log-likelihood gradient, and this may be implemented in future
#' releases of the code.
#'
#' @param X_s List of lenght \eqn{S}{S}, corresponding to number of different studies considered.
#' Each element of the list contains a data matrix, with the same number of columns \eqn{P}{P} for all
#' the studies.
#' @param mle The object returned by \code{ecm_msfa}.
#' @param getgrad Should the function return also the gradient at \code{mle}? Default is \code{FALSE}.
#' @return A list with exactly the same structure of the three slots \code{Phi}, \code{Lambda_s} and
#' \code{Psi_s} of \code{mle}, but containing the standard errors rather than the point estimates.
#' Furthemore, slots for the hessian matrix and the gradient at \code{mle} are included, the latter
#' is not NULL when \code{getgrad}  is \code{TRUE}.
#' @export
#' @import statmod
#' @importFrom pracma grad
#' @importFrom pracma hessian
vcov_msfa <- function(X_s, mle, getgrad = TRUE)
{
  constraint <- mle$constraint
  p <- ncol(X_s[[1]])
  k <- ncol(mle$Phi)
  S <- length(X_s)
  j_s <- c()
  for(s in 1:S) j_s[[s]] <- ncol(mle$Lambda_s[[s]])
  theta <- param2vect(mle, mle$constraint)
  gout <- NULL
  if(getgrad) gout <- pracma::grad(loglik_int, x = theta, n_s = mle$n_s,
                                  cov_s = mle$cov_s, k = k, j_s = j_s, constraint=constraint)
  hout <- pracma::hessian(loglik_int, x = theta, n_s = mle$n_s, cov_s = mle$cov_s, k = k, j_s = j_s,
                          constraint=constraint)
  seout <- sqrt(diag(chol2inv(chol(-hout))))
  ###### re-arrange the ses like in param
  param <- vect2param(seout, mle, constraint, p, k, j_s)
  out <- list(Phi = param$Phi, Lambda_s = param$Lambda_s, psi_s=param$psi_s,  grad = gout, hessian = hout)
  return(out)
}




#' Estimates the parameters of a MSFA model
#'
#' Maximum likelihood estimation of the MSFA model parameters via the ECM
#' algorithm.
#'
#' There are two different constraints for achieving model identification,
#' as detailed in the reference,
#' though the function can also be run without such constraints (not recommended).
#' No checking is done on the starting value for the various model matrices,
#' since a suitable value for them  is produced by the function \code{start_msfa}.
#' @param X_s List of lenght \eqn{S}{S}, corresponding to number of different studies considered.
#' Each element of the list contains a data matrix, with the same number of columns \eqn{P}{P} for all the studies.
#' @param start A list containing the slots \code{Phi}, \code{Lambda_s} and \code{Psi_s}, containing the starting
#' values for the matrix  \code{Phi} of common factor loadings, of size \eqn{P \times K}{P x K}, for
#' the matrices \code{Lambda_s} of study-specific factor loadings, a list of size \eqn{S}{S}  where each element
#' contains a matrix with \eqn{P \times J_s}{P x J_s}, and finally for the study-specific matrices of uniquenesses,
#' a list of size \eqn{S}{S}, where each element contains a vector of length \eqn{P}{P}.
#' Note that a suitable list of this kind is produced by \code{start_msfa}.
#' @param nIt Maximum number of iterations for the ECM algorithm. Default is 50000.
#' @param tol Tolerance for declaring convergence of the ECM algorithm. Default is 10^-7.
#' @param constraint  Constraint for ensuring identifiability. The default is "block_lower2", which
#' corresponds to the main proposal of De Vito et al. (2018). An alternative identification
#' strategy is triggered by  "block_lower1"; this is more restrictive but may work also with smaller
#' number of variables. Again, the latter strategy is mentioned in De Vito et al. (2018).
#' @param robust If \code{TRUE}, robust covariance matrix is used in place of the sample covariance. Default
#' is \code{FALSE}.
#' @param corr If \code{TRUE}, the analysis will employ the correlation matrix instead of the covariance matrix.
#' @param mcd If \code{TRUE}, the robust estimator used for the covariance is the same proposed in Pison et al. (2003),
#' otherwise the default value of the function \code{CovRob} of the \code{robust} library is employed. Default is
#' \code{FALSE}.
#' @param trace If \code{TRUE} then trace information is being printed every 1000 iterations of the ECM algorithm.
#' @return A list  containing the following components:
#' \item{\code{Phi},\code{Lambda_s}, \code{psi_s}}{the estimated model matrices.}
#' \item{loglik}{the value of the log likelihood function at the final estimates.}
#' \item{\code{AIC, BIC}}{model selection criteria at the estimate.}
#' \item{\code{npar}}{number of model parameters.}
#' \item{iter}{the number of ECM iterations performed.}
#' \item{constraint}{the identification constraint enforced.}
#' @export
#' @import robust
#' @importFrom stats cor cov factanal prcomp
#' @references De Vito, R., Bellio, R., Trippa, L. and Parmigiani, G. (2018). (2019). Multi-study Factor Analysis. Biometrics,  75, 337-346.
#' @references Pison, G., Rousseeuw, P.J., Filzmoser, P. and Croux, C. (2003). Robust factor analysis. Journal
#' of Multivariate Analysis, 84, 145-172.
ecm_msfa <- function(X_s, start, nIt = 50000, tol = 10^-7, constraint = "block_lower2", robust = FALSE,
                     corr = TRUE, mcd = FALSE, trace = TRUE, extend=FALSE)
{
  #######
  S <- length(X_s)
  j_s <- n_s <- numeric(S)
  ############
  Phi <- start$Phi
  p <- dim(Phi)[[1]]
  k <- dim(Phi)[[2]]
  Lambda_s <- start$Lambda_s
  psi_s <- start$psi_s
  theta <- param2vect(start, constraint)

  #######defining objects
  Psi_s1 <- Psi_s <- cov_s <- list()
  L_s <- list()

  ######1st round of cycle
  for(s in 1:S){
  	n_s[s] <-  dim(X_s[[s]])[[1]]
  	j_s[s] <-  dim(Lambda_s[[s]])[[2]]
  	Psi_s[[s]] <- diag(psi_s[[s]])
    Psi_s1[[s]] <-  diag(1/psi_s[[s]])
    if((!robust) & (!corr)) cov_s[[s]] <- cov(X_s[[s]])
    if((!robust) & corr) cov_s[[s]] <- cor(X_s[[s]])
    if(robust & mcd) cov_s[[s]] <- covRob(X_s[[s]], estim = "mcd", quan = .75, ntrial = 1000, corr = corr)$cov
    if(robust & (!mcd)) cov_s[[s]] <- covRob(X_s[[s]], corr = corr)$cov
    }
  ######E-step
  out <- exp_values(Phi, Lambda_s, Psi_s, Psi_s1, cov_s, getdet = TRUE)
  Sig_s1 <- out$Sig_s1
  ds_s <- out$ds_s
  l_stop0 <- 0
  lm1 <- 0
  l0 <- loglik_ecm(Sig_s1,  ds_s, n_s, cov_s)
  for (i in (1:nIt))
  {
   ###########CM1 ---------------------------------------------------------------------------------------

   ######expected values
   out <- exp_values(Phi, Lambda_s, Psi_s, Psi_s1, cov_s)
   Txsfs <- out$Txsfs; Txsfcs <- out$Txsfcs; Tfsfs <- out$Tfsfs; Tfcsfcs <- out$Tfcsfcs; Tfcsfs <- out$Tfcsfs
   ######update  of Phi_s
   Psi_new <- list()
   Psi_new1 <- list()
   psi_new <- list()

   for(s in 1:S){
   	psi_new[[s]]  <- diag(cov_s[[s]] + Phi %*% Tfcsfcs[[s]] %*% t(Phi) + Lambda_s[[s]] %*%
   	                 Tfsfs[[s]] %*% t(Lambda_s[[s]]) - 2*Txsfcs[[s]] %*% t(Phi) -  2*Txsfs[[s]] %*% t(Lambda_s[[s]]) +
   	                 2 * Phi %*% Tfcsfs[[s]] %*% t(Lambda_s[[s]]))
   	Psi_new[[s]] <- diag(psi_new[[s]])
   	##########inverse
   	Psi_new1[[s]] <- diag(1/diag(Psi_new[[s]]))
   	 }

   ###########CM2 ---------------------------------------------------------------------------------------

   ######expected values
   out<- exp_values(Phi, Lambda_s, Psi_new, Psi_new1, cov_s)
   Txsfs <- out$Txsfs; Txsfcs <- out$Txsfcs; Tfsfs <- out$Tfsfs;
   Tfcsfcs <- out$Tfcsfcs; Tfcsfs <- out$Tfcsfs

   ######update of Phi
   C_s <- list()
   kron_s <- list()
   for(s in 1:S){
      	C_s[[s]] <- n_s[s] * Psi_new1[[s]] %*% Txsfcs[[s]] - n_s[s] * Psi_new1[[s]] %*% Lambda_s[[s]] %*% t(Tfcsfs[[s]])
      	kron_s[[s]] <- kronecker(t(Tfcsfcs[[s]]), n_s[s] * Psi_new1[[s]])
      	}
    C <- Reduce('+', C_s)
    kron <- Reduce('+', kron_s)
    Phi_vec <- solve(kron) %*% matrix(as.vector(C))
    Phi_new <- matrix(Phi_vec, p, k)


   ########CM3 ---------------------------------------------------------------------------------------

   ######expected values
   out <- exp_values(Phi_new, Lambda_s, Psi_new, Psi_new1, cov_s)
   Txsfs <- out$Txsfs; Txsfcs <- out$Txsfcs; Tfsfs <- out$Tfsfs;
   Tfcsfcs <-  out$Tfcsfcs; Tfcsfs <- out$Tfcsfs

   ######update of Phi
   Lambda_new <- list()
   for(s in 1:S){
   	Lambda_new[[s]] <- matrix(((Txsfs[[s]] - Phi_new %*% Tfcsfs[[s]]) %*% solve(Tfsfs[[s]])), p, j_s[s])
   }

	######constraint

  if (constraint == "null")  {
    Phi_new <- Phi_new
   	for (s in 1:S) Lambda_new[[s]] <- Lambda_new[[s]]
    }
  if (constraint == "block_lower1")  {
    Phi_new[upper.tri(Phi_new)] <- 0
    for (s in 1:S){
   	  L_s[[s]] <- cbind(Phi_new, Lambda_new[[s]])
   	  L_s[[s]][upper.tri(L_s[[s]])] <- 0
   	  Phi_new <- matrix(L_s[[s]][,1:k], nrow=p, ncol = k)
      Lambda_new[[s]] <- L_s[[s]][,(k+1):(k+j_s[s])]
     }
	 }

   ###The following block ensures the full rank condition holds
   if (constraint == "block_lower2")  {
     lambda_vals <- c()
     psi_vals <- psi_new <- c()
     Phi_new[upper.tri(Phi_new)] <- 0
     phi_val <- as.vector(Phi_new[lower.tri(Phi_new, diag = TRUE)])
   	 for (s in 1:S){
        Lambda_new[[s]][upper.tri(Lambda_new[[s]])] <- 0
        lambda_vals <- c(lambda_vals, as.vector(Lambda_new[[s]][lower.tri(Lambda_new[[s]], diag = TRUE)]))
        psi_new[[s]] <- diag(Psi_new[[s]])
        psi_vals <- c(psi_vals, psi_new[[s]])
       }
    L_sTOT <- Reduce('cbind', Lambda_new)
    Omega <- cbind(Phi_new, L_sTOT)
    rank_tot <-  qr(Omega)$rank
    theta_new <- c(phi_val, lambda_vals, psi_vals)
    param.struct <- list(Phi = Phi_new, Lambda_s = Lambda_new, psi_s=psi_new)
    Delta <- theta_new - theta
    sh <- 0   ###no more than 20 step-halving rounds

    while( (rank_tot < k + sum(j_s)) & (sh<20))
    {
       Delta <- Delta / 2
       sh <- sh + 1
       theta_new <- theta + Delta
       param <- vect2param(theta_new, param.struct, constraint, p, k, j_s)
       Lambda_new <- c()
       psi_new <- param$psi_new
       for(s in 1:S)
         {
           Lambda_new[[s]] <- param$Lambda_s[[s]]
           Psi_new[[s]] <- diag(psi_new[[s]])
           Psi1_new[[s]] <- diag(1 / psi_new[[s]])
       }
       L_sTOT <- Reduce('cbind', Lambda_new)
       Phi_new <- param$Phi
       Omega <- cbind(Phi_new, L_sTOT)
       rank_tot <-  qr(Omega)$rank
      }

   if(sh==20) stop("The full rank condition does not hold\n")
   }

  ###########stopping rule
  out <- exp_values(Phi_new, Lambda_new, Psi_new, Psi_new1, cov_s, getdet = TRUE)
  Sig_s1 <- out$Sig_s1
  ds_s <- out$ds_s
  l1 <- loglik_ecm(Sig_s1,  ds_s, n_s, cov_s)
  a <- (l1 - l0)/ (l0-lm1)
  l_stop <- lm1 + (1/ (1-a)) * (l0-lm1)
  l0 <- loglik_ecm(Sig_s1,  ds_s, n_s, cov_s)
  if((trace) & (i %% 1000 == 0))  cat("i=", i, "Criterion for convergence ", abs(l_stop-l_stop0),  "\n")
  if((abs(l_stop-l_stop0)<tol) & i > 1 & l_stop != Inf) break

  Psi_s <- Psi_new
  psi_s <- psi_new
  Phi <- Phi_new
  Lambda_s <- Lambda_new
  if (constraint == "block_lower2") theta <- theta_new
  Psi_s1 <- Psi_new1
  lm1 <- l0
  l0 <- l1
  l_stop0 <- l_stop

  #####AIC and BIC computation

  if (constraint == "block_lower1") npar <- p * S + k * (p - ( k - 1) / 2) +  sum(j_s * (p - k - (j_s - 1) / 2))
  if (constraint == "block_lower2")  npar <- p * S + k * (p - ( k - 1) / 2) +  sum(j_s * (p  - (j_s - 1) / 2))
  n_tot <- sum(n_s)
  AIC <- -2 * l1 + npar * 2
  BIC <- -2 * l1 + npar * log(n_tot)
  }
  ############return output
  res <- list(Phi = Phi, Lambda_s = Lambda_s, psi_s = psi_s, loglik = l1,
              AIC = AIC, BIC = BIC, npar=npar,
              iter = i,  cov_s = cov_s,  n_s = n_s, constraint=constraint)
  if(extend)
  {
    Gamma_trivial = diag(p)
    for(pred in 1:p)
    {
      minVar = 10000
      for(s in 1:S)
      {
        if(psi_s[[s]][pred] < minVar)
          minVar = psi_s[[s]][pred]

        Gamma_trivial[pred,pred] = 0.5*minVar
        if(minVar == 10000) print("[MSFA-X] High variance issue; no study had variance less than 10000.")
      }
    }

    H_s_trivial = list()
    for(s in 1:S)
    {
      H_s_trivial[[s]] = diag(psi_s[[s]])-Gamma_trivial
    }

    res$Gamma = Gamma_trivial
    res$H_s = H_s_trivial
    res$SharedCov = Phi%*%t(Phi) + Gamma_trivial
    res$StudyCov = list()
    for(s in 1:S)
      res$StudyCov[[s]] = Lambda_s[[s]] %*% t(Lambda_s[[s]]) + H_s_trivial[[s]]
    res$SharedPrec = solve(res$SharedCov)
    res$StudyPrec = lapply(res$StudyCov,solve)

  }

  res$X_s = X_s

  return(res)
}

#' Estimates the parameters of a MSFA-X model
#'
#' Maximum likelihood estimation of the MSFA-X model parameters via the ECM
#' algorithm.
#'
#' This function parallels the approach in the \code{ecm_msfa} function; please
#' see \code{?ecm_msfa} for important details on identifiability. The difference
#' between this function and original MSFA is that here we break down the noise
#' term into a shared noise term parameterized by Gamma and a study-specific noise term
#' parameterized by H_s for study s.
#'
#' @param X_s List of length \eqn{S}{S}, corresponding to number of different studies considered.
#' Each element of the list contains a data matrix, with the same number of columns \eqn{P}{P} for all the studies.
#' @param start A list containing the slots \code{Phi}, \code{Lambda_s} and \code{Psi_s}, containing the starting
#' values for the matrix  \code{Phi} of common factor loadings, of size \eqn{P \times K}{P x K}, for
#' the matrices \code{Lambda_s} of study-specific factor loadings, a list of size \eqn{S}{S}  where each element
#' contains a matrix with \eqn{P \times J_s}{P x J_s}, and finally for the study-specific matrices of uniquenesses,
#' a list of size \eqn{S}{S}, where each element contains a vector of length \eqn{P}{P}.
#' Note that a suitable list of this kind is produced by \code{start_msfa}.
#' @param nIt Maximum number of iterations for the ECM algorithm. Default is 50000.
#' @param tol Tolerance for declaring convergence of the ECM algorithm. Default is 10^-7.
#' @param constraint  Constraint for ensuring identifiability. The default is "block_lower2", which
#' corresponds to the main proposal of De Vito et al. (2018). An alternative identification
#' strategy is triggered by  "block_lower1"; this is more restrictive but may work also with smaller
#' number of variables. Again, the latter strategy is mentioned in De Vito et al. (2018).
#' @param robust If \code{TRUE}, robust covariance matrix is used in place of the sample covariance. Default
#' is \code{FALSE}.
#' @param corr If \code{TRUE}, the analysis will employ the correlation matrix instead of the covariance matrix.
#' @param mcd If \code{TRUE}, the robust estimator used for the covariance is the same proposed in Pison et al. (2003),
#' otherwise the default value of the function \code{CovRob} of the \code{robust} library is employed. Default is
#' \code{FALSE}.
#' @param trace If \code{TRUE} then trace information is being printed every 1000 iterations of the ECM algorithm.
#' @return A list  containing the following components:
#' \item{\code{Phi},\code{Lambda_s}, \code{psi_s}}{the estimated model matrices.}
#' \item{loglik}{the value of the log likelihood function at the final estimates.}
#' \item{\code{AIC, BIC}}{model selection criteria at the estimate.}
#' \item{\code{npar}}{number of model parameters.}
#' \item{iter}{the number of ECM iterations performed.}
#' \item{constraint}{the identification constraint enforced.}
#' @export
#' @import robust
#' @importFrom stats cor cov factanal prcomp
#' @references De Vito, R., Bellio, R., Trippa, L. and Parmigiani, G. (2018). (2019). Multi-study Factor Analysis. Biometrics,  75, 337-346.
#' @references Pison, G., Rousseeuw, P.J., Filzmoser, P. and Croux, C. (2003). Robust factor analysis. Journal
#' of Multivariate Analysis, 84, 145-172.
#'
ecm_msfax = function (X_s, start, nIt = 50000, tol = 10^-7, constraint = "block_lower2",
                      robust = FALSE, corr = TRUE, mcd = FALSE, verbose = TRUE, noise_trace = FALSE)
{
  S <- length(X_s)
  j_s <- n_s <- numeric(S)
  Phi <- start$Phi
  p <- dim(Phi)[[1]]
  k <- dim(Phi)[[2]]
  Lambda_s <- start$Lambda_s
  psi_s <- start$psi_s # this is a vector of length P
  theta <- param2vect(start, constraint)
  Psi_s1 <- Psi_s <- cov_s <- list() # this will store P x P matrices later
  L_s <- list()
  for (s in 1:S) {
    n_s[s] <- dim(X_s[[s]])[[1]]
    j_s[s] <- dim(as.matrix(Lambda_s[[s]]))[[2]]
    Psi_s[[s]] <- diag(psi_s[[s]])
    Psi_s1[[s]] <- diag(1/psi_s[[s]])
    if ((!robust) & (!corr))
      cov_s[[s]] <- cov(X_s[[s]])
    if ((!robust) & corr)
      cov_s[[s]] <- cor(X_s[[s]])
    if (robust & mcd)
      cov_s[[s]] <- covRob(X_s[[s]], estim = "mcd", quan = 0.75,
                           ntrial = 1000, corr = corr)$cov
    if (robust & (!mcd))
      cov_s[[s]] <- covRob(X_s[[s]], corr = corr)$cov
  }
  out <- exp_values(Phi, Lambda_s, Psi_s, Psi_s1, cov_s, getdet = TRUE)
  Sig_s1 <- out$Sig_s1
  ds_s <- out$ds_s
  l_stop0 <- 0
  lm1 <- 0
  l0 <- loglik_ecm(Sig_s1, ds_s, n_s, cov_s)

  ## for starting point, we need a guess for Gamma and H_s.
  ## Initializing Gamma
  ## Let's just take half of the minimum observed variance for each predictor
  Gamma = diag(p)
  for(pred in 1:p)
  {
    minVar = 1000
    for(s in 1:S)
    {
      if(start$psi_s[[s]][pred] < minVar)
        minVar = start$psi_s[[s]][pred]
    }

    Gamma[pred,pred] = 0.5*minVar
    if(minVar == 1000) print("Oops")
  }
  ## Initializing H_s
  H_s = list()
  for(s in 1:S)
    H_s[[s]] = Psi_s[[s]]-Gamma

  gamma_trace = data.frame("x"=diag(Gamma))
  h_1_trace = data.frame("x"=diag(H_s[[1]]))
  h_2_trace = data.frame("x"=diag(H_s[[2]]))

  updatePhi = T
  updateLambda_s = T
  updateGamma = T
  updateH_s = T

  for (i in (1:nIt)) {
    Theta_curr = list("Phi"=Phi,"Lambda_s"=Lambda_s,
                      "Gamma"=Gamma,
                      "H_s"=H_s)
    out <- exp_values(Phi, Lambda_s, Psi_s, Psi_s1, cov_s)
    Txsxs <- out$Txsxs
    Txsfs <- out$Txsfs
    Txsfcs <- out$Txsfcs
    Tfsfs <- out$Tfsfs
    Tfcsfcs <- out$Tfcsfcs
    Tfcsfs <- out$Tfcsfs

    Psi_new <- list()
    Psi_new1 <- list()
    psi_new <- list()
    H_s_new <- list()

    ## CM1: Update H_s
    ## Doing this first seems to keep Gamma positive

    if(updateH_s)
    {
      for(s in 1:S){
        M = Txsxs[[s]] +
          Phi %*% Tfcsfcs[[s]] %*% t(Phi) +
          Lambda_s[[s]] %*% Tfsfs[[s]] %*% t(Lambda_s[[s]]) -
          2*Txsfcs[[s]] %*% t(Phi) -
          2*Txsfs[[s]] %*% t(Lambda_s[[s]])+
          2*Phi %*% Tfcsfs[[s]] %*% t(Lambda_s[[s]])
        H_s_new[[s]] = diag(diag(M - Gamma))
      }

      h_1_trace = cbind.data.frame(h_1_trace, diag(H_s_new[[1]]))
    }
    if(!updateH_s)
      H_s_new = H_s

    ## CM2: Update Gamma
    ## Question: Do we need to update Psi_s and the conditional expectations here?
    ## If I do this, it tends to lead to zero or negative estimates for Gamma
    ## If I don't do it, then things stay positive

    # for(s in 1:S){
    #   Psi_new[[s]] <- H_s_new[[s]] + Gamma
    #   psi_new[[s]]  <-  diag(H_s_new[[s]] + Gamma)
    #   ##########inverse
    #   Psi_new1[[s]] <- diag(1/diag(Psi_new[[s]]))
    # }
    #
    # out <- exp_values(Phi, Lambda_s, Psi_new, Psi_new1, cov_s)
    # Txsxs <- out$Txsxs
    # Txsfs <- out$Txsfs
    # Txsfcs <- out$Txsfcs
    # Tfsfs <- out$Tfsfs
    # Tfcsfcs <- out$Tfcsfcs
    # Tfcsfs <- out$Tfcsfs

    if(updateGamma)
    {
      gamma_new = list() # each element here is solved with the polynomial
      for(pred in 1:p)
      {
        # update gamma_new by solving polynomial
        myPoly = function(x)
        {
          q = 0
          for(s in 1:S)
          {
            M = Txsxs[[s]] +
              Phi %*% Tfcsfcs[[s]] %*% t(Phi) +
              Lambda_s[[s]] %*% Tfsfs[[s]] %*% t(Lambda_s[[s]]) -
              2*Txsfcs[[s]] %*% t(Phi) -
              2*Txsfs[[s]] %*% t(Lambda_s[[s]])+
              2*Phi %*% Tfcsfs[[s]] %*% t(Lambda_s[[s]])

            otherStudies = setdiff(1:S,s)

            prodTerm = 1
            for(sPrime in otherStudies)
            {
              prodTerm <- prodTerm*(x+diag(H_s_new[[sPrime]])[pred])^2
            }

            q <- q + n_s[[s]]/2*prodTerm*(M[pred,pred] - (x+diag(H_s_new[[s]])[pred]))
          }

          return(q)
        }

        mySol = uniroot(myPoly, lower = 0, upper = 10, extendInt = "downX",trace=1)
        gamma_new[[pred]] = ifelse(mySol$root>0, mySol$root,0)

      }

      Gamma_new = diag(gamma_new)

      gamma_trace = cbind.data.frame(gamma_trace, diag(Gamma_new))
    }
    if(!updateGamma)
      Gamma_new = Gamma

    ## if either gamma or H_s is still being updated
    ## Update Psi with new value of Gamma and H_s

    for(s in 1:S){
      if(updateGamma | updateH_s)
      {
        Psi_new[[s]] <- H_s_new[[s]] + Gamma_new
        psi_new[[s]]  <-  diag(H_s_new[[s]] + Gamma_new)
        ##########inverse
        Psi_new1[[s]] <- diag(1/diag(Psi_new[[s]]))
      }
      else
      {
        Psi_new[[s]] = Psi_s
        psi_new[[s]] = psi_s
        Psi_new1[[s]] = Psi_s1
      }
    }


    ## CM3: Update Phi
    if(updatePhi)
    {
      out <- exp_values(Phi, Lambda_s, Psi_new, Psi_new1, cov_s)
      Txsfs <- out$Txsfs
      Txsfcs <- out$Txsfcs
      Tfsfs <- out$Tfsfs
      Tfcsfcs <- out$Tfcsfcs
      Tfcsfs <- out$Tfcsfs
      C_s <- list()
      kron_s <- list()
      for (s in 1:S) {
        C_s[[s]] <- n_s[s] * Psi_new1[[s]] %*% Txsfcs[[s]] -
          n_s[s] * Psi_new1[[s]] %*% Lambda_s[[s]] %*%
          t(Tfcsfs[[s]])
        kron_s[[s]] <- kronecker(t(Tfcsfcs[[s]]), n_s[s] *
                                   Psi_new1[[s]])
      }
      C <- Reduce("+", C_s)
      kron <- Reduce("+", kron_s)
      Phi_vec <- solve(kron) %*% matrix(as.vector(C))
      Phi_new <- matrix(Phi_vec, p, k)
    }
    if(!updatePhi)
      Phi_new = Phi

    ## CM4: Update Lambda
    if(updateLambda_s)
    {
      out <- exp_values(Phi_new, Lambda_s, Psi_new, Psi_new1,
                        cov_s)
      Txsfs <- out$Txsfs
      Txsfcs <- out$Txsfcs
      Tfsfs <- out$Tfsfs
      Tfcsfcs <- out$Tfcsfcs
      Tfcsfs <- out$Tfcsfs
      Lambda_new <- list()
      for (s in 1:S) {
        Lambda_new[[s]] <- matrix(((Txsfs[[s]] - Phi_new %*%
                                      Tfcsfs[[s]]) %*% solve(Tfsfs[[s]])), p, j_s[s])
      }
    }
    if(!updateLambda_s)
      Lambda_new = Lambda_s

    if (constraint == "null") {
      Phi_new <- Phi_new
      for (s in 1:S) Lambda_new[[s]] <- Lambda_new[[s]]
    }
    if (constraint == "block_lower1") {
      Phi_new[upper.tri(Phi_new)] <- 0
      for (s in 1:S) {
        L_s[[s]] <- cbind(Phi_new, Lambda_new[[s]])
        L_s[[s]][upper.tri(L_s[[s]])] <- 0
        Phi_new <- matrix(L_s[[s]][, 1:k], nrow = p,
                          ncol = k)
        Lambda_new[[s]] <- L_s[[s]][, (k + 1):(k + j_s[s])]
      }
    }
    if (constraint == "block_lower2") {
      lambda_vals <- c()
      psi_vals <- psi_new <- c()
      Phi_new[upper.tri(Phi_new)] <- 0
      phi_val <- as.vector(Phi_new[lower.tri(Phi_new, diag = TRUE)])
      for (s in 1:S) {
        Lambda_new[[s]][upper.tri(Lambda_new[[s]])] <- 0
        lambda_vals <- c(lambda_vals, as.vector(Lambda_new[[s]][lower.tri(Lambda_new[[s]],
                                                                          diag = TRUE)]))
        psi_new[[s]] <- diag(Psi_new[[s]])
        psi_vals <- c(psi_vals, psi_new[[s]])
      }

      L_sTOT <- Reduce("cbind", Lambda_new)


      Omega <- cbind(Phi_new, L_sTOT)
      rank_tot <- qr(Omega)$rank
      theta_new <- c(phi_val, lambda_vals, psi_vals)
      param.struct <- list(Phi = Phi_new, Lambda_s = Lambda_new,
                           psi_s = psi_new)
      Delta <- theta_new - theta
      sh <- 0

      while ((rank_tot < k + sum(j_s)) & (sh < 20)) {
        print("hello")
        Delta <- Delta/2
        sh <- sh + 1
        theta_new <- theta + Delta
        param <- vect2param(theta_new, param.struct,
                            constraint, p, k, j_s)
        Lambda_new <- c()
        psi_new <- param$psi_new
        for (s in 1:S) {
          Lambda_new[[s]] <- param$Lambda_s[[s]]
          Psi_new[[s]] <- diag(psi_new[[s]])
          Psi1_new[[s]] <- diag(1/psi_new[[s]])
        }
        L_sTOT <- Reduce("cbind", Lambda_new)
        Phi_new <- param$Phi
        Omega <- cbind(Phi_new, L_sTOT)
        rank_tot <- qr(Omega)$rank
      }
      if (sh == 20)
        stop("The full rank condition does not hold\n")
    }

    out <- exp_values(Phi_new, Lambda_new, Psi_new, Psi_new1,
                      cov_s, getdet = TRUE)
    Sig_s1 <- out$Sig_s1
    ds_s <- out$ds_s

    ## old stopping criteria
    l1 <- loglik_ecm(Sig_s1, ds_s, n_s, cov_s)
    a <- (l1 - l0)/(l0 - lm1)
    l_stop <- lm1 + (1/(1 - a)) * (l0 - lm1)
    l0 <- loglik_ecm(Sig_s1, ds_s, n_s, cov_s)

    #if ((abs(l_stop - l_stop0) < tol) & i > 1 & l_stop !=
    #    Inf)
    #  break
    Gamma <- Gamma_new
    H_s <- H_s_new
    Psi_s <- Psi_new
    psi_s <- psi_new
    Phi <- Phi_new
    Lambda_s <- Lambda_new
    if (constraint == "block_lower2")
      theta <- theta_new
    Psi_s1 <- Psi_new1
    lm1 <- l0
    l0 <- l1

    ## new stopping criteria: need to update the parameters
    new_ll = loglik_msfax(Phi = Phi,
                          Lambda_s = Lambda_s,
                          Gamma = Gamma,
                          H_s = H_s,
                          n_s = n_s,
                          cov_s = cov_s)

    ## log likelihood of model with new Phi vs model with old Phi, all other params are new
    old_ll_phi = loglik_msfax(Phi = Theta_curr$Phi,
                              Lambda_s = Lambda_s,
                              Gamma = Gamma,
                              H_s = H_s,
                              n_s = n_s,
                              cov_s = cov_s)
    ## log likelihood of model with new Lambda_s vs model with old Lambda_s, all other params are new
    old_ll_lambda_s = loglik_msfax(Phi = Phi,
                                   Lambda_s = Theta_curr$Lambda_s,
                                   Gamma = Gamma,
                                   H_s = H_s,
                                   n_s = n_s,
                                   cov_s = cov_s)

    ## log likelihood of model with new Gamma vs model with old Gamma, all other params are new
    old_ll_gamma = loglik_msfax(Phi = Phi,
                                   Lambda_s = Lambda_s,
                                   Gamma = Theta_curr$Gamma,
                                   H_s = H_s,
                                   n_s = n_s,
                                   cov_s = cov_s)

    ## log likelihood of model with new H_s vs model with old H_s, all other params are new
    old_ll_h_s = loglik_msfax(Phi = Phi,
                                Lambda_s = Lambda_s,
                                Gamma = Gamma,
                                H_s = Theta_curr$H_s,
                                n_s = n_s,
                                cov_s = cov_s)

    ## overall old ll
    old_ll = loglik_msfax(Phi = Theta_curr$Phi,
                          Lambda_s = Theta_curr$Lambda_s,
                          Gamma = Theta_curr$Gamma,
                          H_s = Theta_curr$H_s,
                          n_s = n_s,
                          cov_s = cov_s)
    # # print out ll ratios
    # print("Phi ratio:")
    phi_ratio = (new_ll-old_ll_phi)/abs(old_ll_phi)
    # print("Lambda_s ratio:")
    lambda_s_ratio = (new_ll-old_ll_lambda_s)/abs(old_ll_lambda_s)
    # print("Gamma ratio:")
    gamma_ratio = (new_ll-old_ll_gamma)/abs(old_ll_gamma)
    # print("H_s ratio")
    h_s_ratio = (new_ll-old_ll_h_s)/abs(old_ll_h_s)

    print(paste0("Gamma ratio:",gamma_ratio))
    if(phi_ratio < 1e-5)
      updatePhi = F
    if(lambda_s_ratio < 1e-5)
      updateLambda = F
    #if(gamma_ratio < 1e-10)
    #  updateGamma = F
    #if(h_s_ratio < 1e-5)
    #  updateH_s = F

    print("update rules:")
    print(c(updatePhi,updateLambda_s,updateGamma,updateH_s))
    if ((verbose)) #& (i%% == 0))
      cat("i=", i, "Criterion for convergence ", abs((old_ll-new_ll)/old_ll), "\n")
    if (abs((old_ll-new_ll)/old_ll) < tol)
      break
    # takeaway from this: Gamma simply does not affect the likelihood very much
    # the only thing we can really do is just see how much gamma itself is still changing
    gamma_percent_change = (diag(Gamma_new)-diag(Theta_curr$Gamma))/diag(Theta_curr$Gamma)
    print("Gamma percent change:")
    print(round(max(abs(gamma_percent_change)),8))

    l_stop0 <- l_stop
    if (constraint == "block_lower1")
      npar <- p * S + k * (p - (k - 1)/2) + sum(j_s * (p -
                                                         k - (j_s - 1)/2))
    if (constraint == "block_lower2")
      npar <- p * S + k * (p - (k - 1)/2) + sum(j_s * (p -
                                                         (j_s - 1)/2))
    n_tot <- sum(n_s)
    AIC <- -2 * l1 + npar * 2
    BIC <- -2 * l1 + npar * log(n_tot)
  }

  res <- list(Phi = Phi, Lambda_s = Lambda_s, psi_s = psi_s, Gamma = Gamma, H_s = H_s,
              loglik = l1, AIC = AIC, BIC = BIC, npar = npar, iter = i,
              cov_s = cov_s, n_s = n_s, constraint = constraint,
              gammaTrace = gamma_trace, h_1_trace = h_1_trace)
  return(res)
}

#' Estimates the parameters of study-specific FA models
#'
#' Maximum likelihood estimation of study-specific FA models parameters via the ECM
#' algorithm, adopting the upper-triangular zero constraint to achieve identification
#' for each loading matrix. Note: the function can also estimate a FA model for a single
#' study, by specifiyng \code{X_s = list(data)}, where \code{data} is the data matrix.
#' @param X_s List of lenght \eqn{S}{S}, corresponding to number of different studies considered.
#' Each element of the list contains a data matrix, with the same number of columns \eqn{P}{P} for all the studies.
#' @param tot_s Number of latent factors for each study. A vector of positive integers of length \eqn{S}{S}.
#' @param nIt Maximum number of iterations for the ECM algorithm. Default is 50000.
#' @param tol Tolerance for declaring convergence of the ECM algorithm. Default is 10^-7.
#' @param block_lower Should the upper-triangular zero constraint be enforced? Default is \code{TRUE}
#' (strongly suggested).
#' @param robust If \code{TRUE}, robust covariance matrix is used in place of the sample covariance. Default
#' is \code{FALSE}.
#' @param corr If \code{TRUE}, the analysis will employ the correlation matrix instead of the covariance matrix.
#' @param mcd If \code{TRUE}, the robust estimator used for the covariance is the same proposed in Pison et al. (2003),
#' otherwise the default value of the function \code{CovRob} of the \code{robust} library is employed. Default is
#' \code{FALSE}.
#' @param trace If \code{TRUE} then trace information is being printed every \code{traceIT} iterations of the ECM algorithm.
#' @param traceIT Frequency of tracing information.
#' @return A list  containing the following components:
#' \item{\code{Omega_s}, \code{Psi_s}}{the estimated model matrices.}
#' \item{loglik}{the value of the log likelihood function at the final estimates.}
#' \item{\code{AIC, BIC}}{model selection criteria at the estimate.}
#' \item{\code{npar}}{number of model parameters.}
#' \item{iter}{the number of ECM iterations performed.}
#' @export
#' @import robust psych
#' @references De Vito, R., Bellio, R., Trippa, L. and Parmigiani, G. (2019). Multi-study Factor Analysis. Biometrics,  75, 337-346.
#' @references Pison, G., Rousseeuw, P.J., Filzmoser, P. and Croux, C. (2003). Robust factor analysis. Journal
#' Multivariate Analysis, 84, 145-172.
ecm_fa <- function(X_s, tot_s, nIt = 50000, tol = 10^-7, block_lower = TRUE, robust = FALSE, corr = TRUE, mcd = FALSE, trace = TRUE, traceIT = 1000)
{
  Omega_s <- list()
  Psi_s <- psi_s <- list()
  #######
  p <- ncol(X_s[[1]])
  S <- length(X_s)
  n_s <- numeric(S)
  #######defining objects
  Psi_s1 <- list()
  cov_s <- list()
  Phi <- matrix(0, nrow=p, ncol=1)
  ######1st round of cycle
  for(s in 1:S){
    n_s[s] <-  dim(X_s[[s]])[[1]]
    if((!robust) & (!corr)) cov_s[[s]] <- cov(X_s[[s]])
    if((!robust) & corr) cov_s[[s]] <- cor(X_s[[s]])
    if(robust & mcd) cov_s[[s]] <- covRob(X_s[[s]], estim = "mcd", quan = .75, ntrial = 1000, corr = corr)$cov
    if(robust & (!mcd)) cov_s[[s]] <- covRob(X_s[[s]], corr = corr)$cov
    FA.s <- factanal(X_s[[s]], factors = tot_s[[s]], covmat = cov_s[[s]],  n.obs=nrow(X_s[[s]]), rotation = "none")
    Omega_s[[s]] <- FA.s$loadings
    Psi_s[[s]] <- diag(FA.s$uniq)
    Psi_s1[[s]] <-  diag(1/diag(Psi_s[[s]]))
  }
  ######E-step
  out <- exp_values(Phi, Omega_s, Psi_s, Psi_s1, cov_s, getdet = TRUE)
  Sig_s1 <- out$Sig_s1; ds_s= out$ds_s;
  l_stop0 <- 0
  lm1 <- 0
  l0 <- loglik_ecm(Sig_s1,  ds_s, n_s, cov_s)
  for (i in (1:nIt))
  {
    ###########CM1 ---------------------------------------------------------------------------------------

    ######expected values
    out <- exp_values(Phi, Omega_s, Psi_s, Psi_s1, cov_s)
    Txsfs <- out$Txsfs; Txsfcs <- out$Txsfcs; Tfsfs <- out$Tfsfs; Tfcsfcs <- out$Tfcsfcs; Tfcsfs <- out$Tfcsfs
    ######update  of Phi_s
    Psi_new <- list()
    Psi_new1 <- list()

    for(s in 1:S){
      Psi_new[[s]]  <- diag(cov_s[[s]] + Omega_s[[s]] %*% Tfsfs[[s]] %*% t(Omega_s[[s]]) -  2*Txsfs[[s]] %*% t(Omega_s[[s]]) )
      Psi_new[[s]] <- diag(Psi_new[[s]])
      ##########inverse
      Psi_new1[[s]] <- diag(1/diag(Psi_new[[s]]))
    }

    ###########CM2 ---------------------------------------------------------------------------------------

    ######expected values
    out<- exp_values(Phi, Omega_s, Psi_new, Psi_new1, cov_s)
    Txsfs <- out$Txsfs; Txsfcs <- out$Txsfcs; Tfsfs <- out$Tfsfs;
    Tfcsfcs <- out$Tfcsfcs; Tfcsfs <- out$Tfcsfs

    ######update of Phi: not needed

    ########CM3 ---------------------------------------------------------------------------------------

    ######expected values
    out <- exp_values(Phi, Omega_s, Psi_new, Psi_new1, cov_s)
    Txsfs <- out$Txsfs; Txsfcs <- out$Txsfcs; Tfsfs <- out$Tfsfs;
    Tfcsfcs <-  out$Tfcsfcs; Tfcsfs <- out$Tfcsfs

    ######update of Phi
    Omega_new <- list()
    for(s in 1:S) {
      Omega_new[[s]] <- matrix((Txsfs[[s]] %*% solve(Tfsfs[[s]])), p, tot_s[s])
      Omega_new[[s]][upper.tri(Omega_new[[s]])] <- 0
    }

    ###########stopping rule
    out <- exp_values(Phi, Omega_new, Psi_new, Psi_new1, cov_s, getdet = TRUE)

    Sig_s1 <- out$Sig_s1
    ds_s <- out$ds_s
    l1 <- loglik_ecm(Sig_s1,  ds_s, n_s, cov_s)
    a <- (l1 - l0)/ (l0-lm1)
    l_stop <- lm1 + (1/ (1-a)) * (l0-lm1)
    l0 <- loglik_ecm(Sig_s1,  ds_s, n_s, cov_s)

    if((trace) & (i %% 100 == 0))  cat("i=", i, "Criterion for convergence ", abs(l_stop-l_stop0), "\n")
    if( (abs(l_stop-l_stop0)<tol) & i > 1 & l_stop != Inf) break
    Psi_s <- Psi_new
    Omega_s <- Omega_new
    Psi_s1 <- Psi_new1
    lm1 <- l0
    l0 <- l1
    l_stop0 <- l_stop
  }
  ############return output
  for(s in 1:S) psi_s[[s]] <- diag(Psi_s[[s]])
  npar <- p * S + sum(tot_s * (p - (tot_s - 1) / 2))
  n_tot <- sum(n_s)
  AIC <- -2 * l1 + npar * 2
  BIC <- -2 * l1 + npar * log(n_tot)
  res <- list(Omega_s = Omega_s, psi_s = psi_s, loglik = l1, AIC = AIC,
              BIC = BIC, npar=npar, iter = i)
  return(res)
}

#' @export
check_constraint = function(p,k,j_s,s)
{
  lhs = p*k-k*(k-1)/2

  for(i in 1:s)
  {
    lhs = lhs + p*j_s[[i]]-j_s[[i]]*(j_s[[i]]-1)/2
  }

  lhs = lhs + s*p
  rhs = s*p*(p+1)/2

  return(lhs <= rhs)
}

# This function is a wrapper for various ways of getting the number of factors from a single
# dataset (i.e., not in the multi-study setting)
get_n_factors_vanilla = function(X, method = "cng")
{
  if(method == "bartlett")
    nTot = nFactors::nBartlett(data.frame(X),N=nrow(X))$nFactors[1] # take Bartlett method
  if(method == "bentler")
    nTot = nFactors::nBentler(data.frame(X),N=nrow(X))$nFactors
  if(method == "cng")
    nTot = nFactors::nCng(data.frame(X),model="factors")$nFactors# cattell, nelson, gorsuch
  if(method == "mreg")
    nTot = nFactors::nMreg(data.frame(X),model="factors")$nFactors[1] # take b
  if(method == "paran")
    nTot= paran(data.frame(X),cfa=T)$Retained
  if(method == "scree")
    nTot = nFactors::nScree(data.frame(X),model="factors")$Components$noc # take the optimal coordinates results (can implement others later)
  if(method == "seScree")
    nTot = nFactors::nSeScree(data.frame(X),model="factors")$nFactors[2] # take the R2: se was unreasonably high

  return(nTot)
}

#' Estimates the number of shared and study-specific factors present in multi-study data.
#'
#' @param X_s List of length \eqn{S}{S}, corresponding to number of different studies considered.
#' Each element of the list contains a data matrix, with the same number of columns \eqn{P}{P} for all the studies.
#' No standardization is carried out by the function.
#' @param method Method for estimating the number of factors in the shared data, from the \code{nFactors} package. Default is "cng".
#' Options include "bartlett","bentler","cng","mreg","paran","scree","seScree". See \code{nFactors} documentation for details.
#' @return A list  containing  \code{bics}, a list of the BIC values for each number of factors considered,
#' \code{pooledBIC}, the BIC for a standard factor analysis model that does not split the studies,
#' \code{nTot}, the total number of factors estimated for each study according to minimum BIC,
#' \code{k}, the number of shared factors estimated according to minimum BIC,
#' \code{j_s}, the number of study-specific factors in each study estimated according to minimum BIC.
#' @import nFactors
#' @references Raiche, Gilles, David Magis, and Maintainer Gilles Raiche. "Package ‘nFactors’." Repository CRAN (2020): 1-58.
#' @export
get_factor_count = function(X_s, method = "cng")
{
  print(paste("[get_factor_count] METHOD:",method))
  nTot = list()
  S = length(X_s)
  p = ncol(X_s[[1]])

  # compare with a baseline model that does not include study-specific factors, via AIC/BIC
  nTotPooled = NA
  X_s_pooled = do.call(rbind,X_s)
  nTotPooled = get_n_factors_vanilla(X_s_pooled,method = method)
  modelPooled = factanal(X_s_pooled,factors = nTotPooled)

  # the covariance to test is PhiPhi^T + sigma
  modelCov = modelPooled$loadings %*% t(modelPooled$loadings) + diag(modelPooled$uniquenesses)
  pooledLL = loglik_vanilla(Sig_s1 = solve(modelCov),  ds_s = det(modelCov), n_s = nrow(X_s_pooled), cov_s = cov(X_s_pooled))

  # number of parameters: p for Psi and k*(p-(k-1))/2 for Phi?
  nparPooled <- p + nTotPooled * (p - (nTotPooled - 1) / 2)
  pooledAIC = -2 * pooledLL + nparPooled * 2
  pooledBIC = -2 * pooledLL + nparPooled * log(nrow(X_s_pooled))

  print(paste("Pooled likelihood:",pooledLL))
  print(paste("Pooled number of parameters:",nparPooled))
  print("ICs: lower is better")
  print(paste("factanal AIC:", pooledAIC))
  print(paste("factanal BIC:", pooledBIC))

  for(s in 1:S)
  {
    nTot[[s]] = get_n_factors_vanilla(X_s[[s]], method = method)
  }

  # check to see if there are too many factors, as in factanal
  p = ncol(X_s[[1]]) # assume same number of predictors for every study
  #dof <- 0.5 * ((p - factors)^2 - p - factors) from factanal
  #dof = 0.5*(p^2 - 2*p*factors + factors^2 - p - factors)
  #dof = 0.5*(factors^2 - (2*p+1)*factors + p^2-p)
  roots = c((2*p+1 - sqrt((2*p+1)^2 - 4*(p^2-p)))/2,(2*p+1 + sqrt((2*p+1)^2 - 4*(p^2-p)))/2)
  # problem: what if nTot is bigger than roots[1]
  for(s in 1:S)
  {
    if(nTot[[s]] > floor(roots[1]))
    {
      nTot[[s]] = floor(roots[1]) # can't handle more factors than this
    }
  }

  maxShared = min(unlist(nTot))-1
  # for every possible value of maxShared
  bics = rep(NA,maxShared)
  for(k in maxShared:1)
  {
    j_s = unlist(nTot) - k
    print(paste("Testing k = ",k))
    if(check_constraint(p,k,j_s,length(j_s)))
    {

      start_k = tryCatch(expr = {start_msfa(X_s,k=k,j_s=unlist(j_s),method="fa")},
                         error = {function(e) NULL})

      if(!is.null(start_k))
      {
        start_k = start_msfa(X_s,k=k,j_s=unlist(j_s),method="fa")
        start_k$Phi = as.matrix(start_k$Phi)
        start_k$Lambda_s = lapply(start_k$Lambda_s,as.matrix)
        mod_k = ecm_msfa(X_s,start=start_k,extend=T,tol=1e-3) #,constraint="block_lower1")
      }

      else
      {
        print(paste("[get_factor_count] Not able to start fa with k=",k))
      }
    }
    bics[k] =mod_k$BIC
  }

  return(list("bics"=bics,
              "pooledBIC"=pooledBIC,
              "nTot"=nTot,
              "k"=which.min(bics),
              "j_s"=unlist(nTot)-which.min(bics)))
}

#' @export
get_factor_count_bayes = function(X_s, method = "cng")
{
  # run factanal on each dataset to determine k* and j* (j_s* = tot_s - k*)

  print(paste("[get_factor_count] METHOD:",method))
  nTot = list()
  S = length(X_s)
  p = ncol(X_s[[1]])

  # compare with a baseline model that does not include study-specific factors, via AIC/BIC
  nTotPooled = NA
  X_s_pooled = do.call(rbind,X_s)
  nTotPooled = get_n_factors_vanilla(X_s_pooled,method = method)
  modelPooled = factanal(X_s_pooled,factors = nTotPooled)

  # the covariance to test is PhiPhi^T + sigma
  modelCov = modelPooled$loadings %*% t(modelPooled$loadings) + diag(modelPooled$uniquenesses)
  pooledLL = loglik_vanilla(Sig_s1 = solve(modelCov),  ds_s = det(modelCov), n_s = nrow(X_s_pooled), cov_s = cov(X_s_pooled))

  # number of parameters: p for Psi and k*(p-(k-1))/2 for Phi?
  nparPooled <- p + nTotPooled * (p - (nTotPooled - 1) / 2)
  pooledAIC = -2 * pooledLL + nparPooled * 2
  pooledBIC = -2 * pooledLL + nparPooled * log(nrow(X_s_pooled))

  print(paste("Pooled likelihood:",pooledLL))
  print(paste("Pooled number of parameters:",nparPooled))
  print("ICs: lower is better")
  print(paste("factanal AIC:", pooledAIC))
  print(paste("factanal BIC:", pooledBIC))

  for(s in 1:S)
  {
    nTot[[s]] = get_n_factors_vanilla(X_s[[s]], method = method)
  }

  print(nTot)

  # check to see if there are too many factors, as in factanal
  p = ncol(X_s[[1]]) # assume same number of predictors for every study
  #dof <- 0.5 * ((p - factors)^2 - p - factors) from factanal
  #dof = 0.5*(p^2 - 2*p*factors + factors^2 - p - factors)
  #dof = 0.5*(factors^2 - (2*p+1)*factors + p^2-p)
  roots = c((2*p+1 - sqrt((2*p+1)^2 - 4*(p^2-p)))/2,(2*p+1 + sqrt((2*p+1)^2 - 4*(p^2-p)))/2)
  # problem: what if nTot is bigger than roots[1]
  for(s in 1:S)
  {
    if(nTot[[s]] > floor(roots[1]))
    {
      nTot[[s]] = floor(roots[1]) # can't handle more factors than this
    }
  }

  maxShared = min(unlist(nTot))-1

  # fit MSFA with the k* and j*
  k = maxShared
  mod_k = NULL # in case it is not able to fit, detect the problem
  minShared = 1

  j_s = list()
  for(s in 1:S)
    j_s[[s]]=nTot[[s]]-minShared # intentional? This gives more factors than maybe possible, makes sense in Bayes case perhaps

  while(k > 0)
  {
    print(paste("Testing k = ",k))
    if(check_constraint(p,k,j_s,length(j_s)))
    {

      start_k = tryCatch(expr = {start_msfa(X_s,k=k,j_s=unlist(j_s),method="fa")},
                         error = {function(e) NULL})

      if(!is.null(start_k))
      {
        start_k = start_msfa(X_s,k=k,j_s=unlist(j_s),method="fa")
        start_k$Phi = as.matrix(start_k$Phi)
        start_k$Lambda_s = lapply(start_k$Lambda_s,as.matrix)
        mod_k = ecm_msfa(X_s,start=start_k,extend=T,tol=1e-3) #,constraint="block_lower1")
        k = 0
      }

      else
      {
        print(paste("[get_factor_count] Not able to start fa with k=",k))
        k = k-1
      }
    }
    else
    {k = k-1}
  }

  # get eigenvalues that explain > 0.05 of total sum of eigenvalues
  if(is.null(mod_k))
  {
    print("[get_factor_count] Unable to optimize.")
    return(NA)
  }

  sigma_phi = mod_k$Phi %*% t(mod_k$Phi)
  val_eigen = eigen(sigma_phi)$values
  prop_var = val_eigen / sum(val_eigen)

  nShared = sum(prop_var > 0.05)

  j_s = list()
  for(s in 1:S)
  {
    # do same eigenvalue approach for the study specific matrices
    sigma_lambda_s = mod_k$Lambda_s[[s]] %*% t(mod_k$Lambda_s[[s]])
    val_eigen = eigen(sigma_lambda_s)$values
    #scree.plot(val_eigen, title=paste("Screeplot Study",s))
    prop_var = val_eigen / sum(val_eigen)
    j_s[[s]]= sum(prop_var > 0.05)

  }

  k = nShared

  msfa_AIC = mod_k$AIC
  msfa_BIC = mod_k$BIC
  print(paste("MSFA likelihood:",mod_k$loglik))
  tot_s = k + unlist(j_s)
  npar <- p * S + sum(tot_s * (p - (tot_s - 1) / 2))
  print(paste("MSFA number of parameters:",npar))
  print(paste("MSFA AIC:",mod_k$AIC))
  print(paste("MSFA BIC:",mod_k$BIC))
  return(list("k"=k, "j_s" = j_s))
}

#' Immune System Data
#
#'
#'
#' A data set used by De Vito et al. (2019) as an illustrative example.
#'
#' @format A list of two matrices, each with 63 columns and 285 and 140 rows, respectively.
#' @examples
#' \dontrun{
#' The commands below show how the dataset was obtained from libraries on the Bioconductor repository.
#' source("http://bioconductor.org/biocLite.R")
#' biocLite(c("limma", "curatedOvarianData", "RTCGAToolbox"), suppressUpdates=TRUE)
#' library(curatedOvarianData)
#' library(RTCGAToolbox)
#' data(package="curatedOvarianData")
#' data(GSE20565_eset)
#' data(GSE9891_eset)
#' im_response <- c(
#'  "TP53",
#'  "TUBA3C","PRKACG","FGF6","FGF23","FGF22","FGF20","ASB4","TUBB1","LAT","ULBP1","NCR1",
#'  "SIGLEC5","CD160","KLRD1","NCR3","TRIM9","FGF18","ICOSLG",
#'  "MYH2","C9","MBL2",
#'  "GRIN2B","POLR2F","CSF2","IL5","CRP","C8A","SPTA1","GRIN2A","CCR6","FGA","LBP",
#'  "DUSP9","FCN2","PRKCG","ADCY8","IL5RA","GRIN1","C8B",
#'  "GH2","TNFSF18","GRIN2D","FGB","PRL","SPTBN5","CD70","FGG","RASGRF1","IFNG",
#'  "SPTBN4","TRIM10","ACTN2","LTA","TNFSF11","GRIN2C","CAMK2B","SPTB","IL1A","TNFRSF13B",
#'  "ITGA2B","CAMK2A","TRIM31","EREG")
#' GSE20_eset <- t(as.matrix(GSE20565_eset))
#' GSE98_eset <- t(as.matrix(GSE9891_eset))
#' leD <- length(im_response)
#' mat1 <- matrix(sapply(1:leD, function(i) which(colnames(GSE20_eset)==im_response[i])))
#' mat2 <- matrix(sapply(1:leD, function(i) which(colnames(GSE98_eset)==im_response[i])))
#' GSE20 <- matrix(c(GSE20_eset[,mat1]), 140, leD)
#' GSE98 <-matrix(c(GSE98_eset[,mat2]), 285, leD)
#' data_immune <- list(GSE98, GSE20)}
"data_immune"


