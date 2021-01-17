###======================================###
###============== LORSEN ================###
###
### Contributor: Cheng Gao
### Email: chenggao@mtu.edu
### Date: 1/16/2021
### Ref: Gao, C., et al, LORSEN: fast and efficient eQTL mapping with low rank penalized regression

#' asso_res_ns: this function is used to calculate the precision given the number of associations 
#' for joint modeling without SNP screening in simulation studies.
#' @param mat: beta matrix yielded from joint modeling
#' @param m: the number of associations
#' @param true.beta: a TRUE-FALSE beta index matrix. In simulations, if \beta_{ij} is nonzero, then 
#' the corresponding entry of true.beta is TRUE; otherwise, FALSE. 

asso_res_ns <- function(mat, m, true.beta){
  
  nonzero.num <- sum((mat != 0))
  
  if(nonzero.num < m){
    
    PPV <- NA
    
  }else{
    
    crit.val <- sort(abs(mat),decreasing = TRUE)[m]
    
    if(sum((abs(mat) == crit.val))>1){
      
      PPV <- NA
      
    }else{
      
      mat.bin <- (abs(mat) >= crit.val)
      
      TP <- sum(mat.bin * true.beta)
      PPV <- TP/m ## precision
    }
  }
  
  return(PPV)
}

#' asso_res_ws: this function is used to calculate the precision given the number of associations 
#' for joint modeling with SNP screening in simulation studies.
#' @param mat: beta matrix yielded from joint modeling
#' @param m: the number of associations
#' @param true.beta: a TRUE-FALSE beta index matrix. In simulations, if \beta_{ij} is nonzero, then 
#' the corresponding entry of true.beta is TRUE; otherwise, FALSE. 
#' @param select.snp.index: the index vector of selected SNPs from SNP screening 
#' @param true.snp.index: the index vector of causal SNPs in simulations
 
asso_res_ws <- function(mat, m, select.snp.index, true.snp.index, true.beta){
  
  nonzero.num <- sum((mat != 0))
  true.sel.inter.num <- length(intersect(select.snp.index,true.snp.index))
  
  if(nonzero.num < m | true.sel.inter.num == 0){
    
    PPV <- NA
    
  }else{
    
    crit.val <- sort(abs(mat),decreasing = TRUE)[m]
    
    if(sum((abs(mat) == crit.val))>1){
      
      PPV <- NA
      
    }else{
      
      mat.bin <- (abs(mat) >= crit.val)
      
      TP <- sum(mat.bin * true.beta)
      PPV <- TP/m ## precision
    }
  }
  
  return(PPV)
}


#' roc: this function returns the TPR sequence and FPR sequence given a sequence of thresholds of 
#' beta's for drawing ROC curve. [Adapted from original code in LORS.]
#' @param B: the beta matrix yielded from joint modeling 
#' @param thr.seq: sequence of thresholds of beta's
#' @param true.beta: a TRUE-FALSE beta index matrix. In simulations, if \beta_{ij} is nonzero, then 
#' the corresponding entry of true.beta is TRUE; otherwise, FALSE.

roc <- function(B,thr.seq,true.beta){
  
  num.thr <- length(thr.seq)
  TPR.seq <- FPR.seq <- rep(NA,num.thr)
  
  for(i in 1:num.thr){
    
    S <- abs(B) > thr.seq[i]
    TPR.seq[i] <- sum(S*true.beta)/sum(true.beta)
    FPR.seq[i] <- sum(S*(!true.beta))/sum(!true.beta)
  }
  
  AUC <- sum((FPR.seq[1:(num.thr-1)]-FPR.seq[2:num.thr])*(TPR.seq[1:(num.thr-1)]+TPR.seq[2:num.thr])/2)
  
  return(list(TPR.seq = TPR.seq, FPR.seq = FPR.seq, AUC = AUC))
}

#' LORS_eval: perform LORS. min 1/2 ||Y - XB - L -1\mu||^{2}_{F} + \rho ||B||_{1} + \lambda ||L||_{*}   
#' @param Y: n * q gene expression matrix, n is sample size, q is the number of genes
#' @param X: n * p SNP matrix, n is sample size, p is the number of SNPs 
#' @param CV: cross-validation, please set as FALSE, left for future development
#' @param nfold: once CV is FALSE, ignore this parameter
#' @param nrho: the number of \rho's to be searched over, 20 is good
#' @param quickRho: TRUE or FALSE
#' @Omega: TRUE-FALSE index matrix. If some entry is observed, TRUE; otherwise, FALSE.
#' @param eps: 2.2204e-16 
#' @param tol: 1e-4 

LORS_eval <- function(Y, X, CV, nfold, nrho, quickRho, Omega, eps, tol){
  
  start <- proc.time()
  start_param <- proc.time()
  
  LORSPara <- LORS_BestPara(Y, X, CV, nfold, nrho, quickRho, Omega, maxiter = 5000)
  lambda <- LORSPara$lambda
  rho <- LORSPara$rho 
  
  end_param <- proc.time()
  param_time <- end_param[3] - start_param[3]
  
  start_model <- proc.time()
  
  LORS_Obj <- LORS0(Y, X, rho, lambda, maxiter = 5000, eps, tol, verbose = TRUE)
  
  end_model <- proc.time()
  model_time <- end_model[3] - start_model[3]
  
  end <- proc.time()
  total_time <- end[3] - start[3]
  
  return(list("LORS_Obj" = LORS_Obj, "param_time" = param_time,
              "model_time" = model_time, "total_time" = total_time, "rho" = rho, "lambda" = lambda))
}

#' FastLORS_eval: perform FastLORS. min 1/2 ||Y - XB - L -1\mu||^{2}_{F} + \rho ||B||_{1} + \lambda ||L||_{*}   
#' @param Y: n * q gene expression matrix, n is sample size, q is the number of genes
#' @param X: n * p SNP matrix, n is sample size, p is the number of SNPs 
#' @param CV: cross-validation, please set as FALSE, left for future development
#' @param nfold: once CV is FALSE, ignore this parameter
#' @param nrho: the number of \rho's to be searched over, 20 is good
#' @param quickRho: TRUE or FALSE
#' @Omega: TRUE-FALSE index matrix. If some entry is observed, TRUE; otherwise, FALSE.
#' @param eps: 2.2204e-16 
#' @param tol: 1e-4 

FastLORS_eval <- function(Y, X, CV, nfold, nrho, quickRho, Omega, eps, tol){
  
  start <- proc.time()
  start_param <- proc.time()
  
  FastLORSPara <- FastLORS_BestPara(Y, X, CV, nfold, nrho, quickRho, Omega, maxiter = 5000)
  lambda <- FastLORSPara$lambda
  rho <- FastLORSPara$rho
  
  end_param <- proc.time()
  param_time <- end_param[3] - start_param[3]
  
  start_model <- proc.time()
  
  Fast_LORS_Obj <- Fast_LORS(Y, X, rho, lambda, maxiter = 5000, eps, tol, verbose = FALSE, omega_SOR = omega_SOR)
  
  end_model <- proc.time()
  model_time <- end_model[3] - start_model[3]
  
  end <- proc.time()
  total_time <- end[3] - start[3]
  
  return(list("Fast_LORS_Obj" = Fast_LORS_Obj, "param_time" = param_time,
              "model_time" = model_time, "total_time" = total_time, "rho" = rho, "lambda" = lambda))
}

#' LORSEN_eval: perform LORSEN with parameter tuning. min 1/2 ||Y - XB - L -1\mu||^{2}_{F} + \lambda(\alpha ||B||_{1} + 
#' (1 - \alpha)/2 ||B||^{2}_{F}) + \rho ||L||_{*}   
#' @param Y: n * q gene expression matrix, n is sample size, q is the number of genes
#' @param X: n * p SNP matrix, n is sample size, p is the number of SNPs 
#' @param alpha_EN_seq: sequence of alpha values to be searched over, for example, (0.2,0.4,0.6,0.8,0.9)
#' @param nlambda: the number of alpha values, for example, 50
#' @param CV: cross-validation, please set as FALSE, left for future development
#' @param nfold: once CV is FALSE, ignore this parameter
#' @param nrho: the number of \rho's to be searched over, 20 is good
#' @param quickRho: TRUE or FALSE
#' @Omega: TRUE-FALSE index matrix. If some entry is observed, TRUE; otherwise, FALSE.
#' @param eps: 2.2204e-16 
#' @param tol: 1e-4 
#' @param Alg: please set as "FISTA", left for future development
#' @param stepsize: "constant" or "backtracking"
#' @param maxiter: the maximum number of iterations, for example, 5000 
#' @Omega: TRUE-FALSE index matrix. If some entry is observed, TRUE; otherwise, FALSE.

LORSEN_eval <- function(Y, X, alpha_EN_seq, nlambda, maxiter = 5000, Omega, eps = 2.2204e-16, tol = 1e-4, verbose = TRUE, Alg, stepsize){
  
  start <- proc.time()
  
  stopp <- TRUE
  
  while(stopp){
    
    start_param <- proc.time()
    LORSENPara <- LORSEN_BestPara(Y, X, alpha_EN_seq, nlambda, Omega, maxiter, Alg, stepsize)
    
    if(length(LORSENPara$alpha_EN.best) > 1){
      
      break
      
    }else{
      
      stopp <- FALSE
    }
    
    lambda.best <- LORSENPara$lambda.best
    alpha_EN.best <- LORSENPara$alpha_EN.best
    lambda_EN.best <- LORSENPara$lambda_EN.best
  
    end_param <- proc.time()
    param_time <- end_param[3] - start_param[3]
  
    start_model <- proc.time()
  
    LORSEN_Obj <- LORSEN(Y, X, lambda=lambda.best, lambda_EN=lambda_EN.best, alpha_EN=alpha_EN.best, maxiter = 5000, eps = 2.2204e-16, tol = 1e-4, verbose = TRUE, Alg, stepsize)
  
    end_model <- proc.time()
    model_time <- end_model[3] - start_model[3]
  }
  end <- proc.time()
  total_time <- end[3] - start[3]
  
  return(list("LORSEN_Obj" = LORSEN_Obj, "param_time" = param_time,
              "model_time" = model_time, "total_time" = total_time, "lambda.best" = lambda.best, "alpha_EN.best" = alpha_EN.best, "lambda_EN.best" = lambda_EN.best))
}

#' EN_tune: this function is used to tune parameter \lambda_EN
#' Please note the difference in parameters used in paper and code here, left for further improvemnet
#' paper: min 1/2 ||Y - XB - L -1\mu||^{2}_{F} + \lambda(\alpha ||B||_{1} + (1 - \alpha)/2 ||B||^{2}_{F}) + \rho ||L||_{*} 
#' code: min 1/2 ||Y - XB - L -1\mu||^{2}_{F} + \lambda_EN(\alpha_EN ||B||_{1} + (1 - \alpha_EN)/2 ||B||^{2}_{F}) + \lambda ||L||_{*}
#' @param Training: index matrix of observations used in training, as \Omega_{1} in paper
#' @param Validation: index matrix of observations used in validation, as \Omega_{2} in paper 
#' @param Alg: please set as "FISTA", left for future development
#' @param stepsize: "constant" or "backtracking"
#' @param maxiter: the maximum number of iterations, for example, 5000  

EN_tune <- function(Y,X,L,lambda,alpha_EN,nlambda,Training,Validation,maxiter,Alg,stepsize){

    ENErr <- matrix(0, nrow = nlambda, ncol = 1)
    B <- matrix(0, nrow = ncol(X), ncol = ncol(Y))
    mu <- matrix(0, nrow = 1, ncol = ncol(Y))
	
    lambda.seq <- lambda_seq(X,Y,L,nlambda,Training,alpha_EN)
	  
	for(ilambda in 1:nlambda){
	    
	   lambda_EN <- lambda.seq[ilambda]
	   out <- LORSEN_Tuning(Y, X, lambda, lambda_EN, alpha_EN, Training, Validation, maxiter = 5000, eps = 2.2204e-16, tol = 1e-4, verbose = TRUE, L, Alg, stepsize)
       Err <- out[["Err"]]
	   ENErr[ilambda,1] <- Err	   
	}	
	
	bestind_lambda <- ENErr==min(ENErr)
    lambda_EN <- mean(lambda.seq[bestind_lambda]) 
    return(list(lambda_EN = lambda_EN, minErr = min(ENErr)))	
}	

#' LORSEN_BestPara returns the optimal parameters for LORSEN.

LORSEN_BestPara <- function(Y, X, alpha_EN_seq, nlambda, Omega, maxiter, Alg, stepsize){

    nfold <- 2
	mask.ls <- mask(nfold, Y)
	   
	params <- ParamTune(Training = Omega & mask.ls[[1]], Validation = Omega & mask.ls[[2]], tune_method)
    lambda.best <- params[["lambda"]]
    L <- params[["L"]]
	   
	nalpha <- length(alpha_EN_seq)
	lambda_EN_seq <- minErr_seq <- numeric(nalpha)
	   
	for(ialpha in 1:nalpha){
	   
	   alpha_EN <- alpha_EN_seq[ialpha]
	   EN.info <- EN_tune(Y,X,L,lambda = lambda.best,alpha_EN,nlambda,Training = Omega & mask.ls[[1]], Validation = Omega & mask.ls[[2]],maxiter,Alg,stepsize)
	   lambda_EN_seq[ialpha] <- EN.info$lambda_EN
	   minErr_seq[ialpha] <- EN.info$minErr
	}
	
    print(paste0("lambda_EN_seq: ", lambda_EN_seq))
    print(paste0("minErr_seq: ", minErr_seq))	
	best.ind <- which(minErr_seq == min(minErr_seq))
	print(paste0("best.ind: ", best.ind))	
	alpha_EN.best <- alpha_EN_seq[best.ind]
	lambda_EN.best <- lambda_EN_seq[best.ind]
    return(list(lambda.best = lambda.best, alpha_EN.best = alpha_EN.best, lambda_EN.best = lambda_EN.best))
}

#' GetMaxLambda calculates the smallest \lambda such that beta matrix B = 0 based on coordinate descent algorithm.
#' min 1/2 ||Y - XB - L -1\mu||^{2}_{F} + \lambda(\alpha ||B||_{1} + (1 - \alpha)/2 ||B||^{2}_{F}) + \rho ||L||_{*} 

GetMaxLambda <- function(X, Y, L, Omega,alpha_EN){
 
  maxlambda = matrix(0, nrow = q, ncol = 1)
  
  for(i in 1:q){
    maxlambda[i,1] = max(abs(t(X) %*% ((matrix(Y[,i], nrow = nrow(Y), ncol = 1) - matrix(L[,i], nrow = nrow(L), ncol = 1)) * matrix(Omega[,i], nrow = nrow(Omega), ncol = 1))/alpha_EN))
  }
  
  MaxLambda = max(maxlambda)
  return(MaxLambda)
}

#' lambda_seq generates the sequence of \lambda_EN values 
#' paper: min 1/2 ||Y - XB - L -1\mu||^{2}_{F} + \lambda(\alpha ||B||_{1} + (1 - \alpha)/2 ||B||^{2}_{F}) + \rho ||L||_{*} 
#' code: min 1/2 ||Y - XB - L -1\mu||^{2}_{F} + \lambda_EN(\alpha_EN ||B||_{1} + (1 - \alpha_EN)/2 ||B||^{2}_{F}) + \lambda ||L||_{*}

lambda_seq <- function(X,Y,L,nlambda,Omega,alpha_EN){

    ## Set a sequence of rho
	
    MaxLambda <- GetMaxLambda(X, Y, L, Omega, alpha_EN)  
    
    lambdaseq <- logspace(log10(MaxLambda),log10(MaxLambda*.02),nlambda)
    
	return(lambdaseq)
}

#' LORSEN: perform LORSEN with fixed tuning parameter [PGD-CD algorithm]
#' paper: min 1/2 ||Y - XB - L -1\mu||^{2}_{F} + \lambda(\alpha ||B||_{1} + (1 - \alpha)/2 ||B||^{2}_{F}) + \rho ||L||_{*} 
#' code: min 1/2 ||Y - XB - L -1\mu||^{2}_{F} + \lambda_EN(\alpha_EN ||B||_{1} + (1 - \alpha_EN)/2 ||B||^{2}_{F}) + \lambda ||L||_{*}

LORSEN <- function(Y, X, lambda, lambda_EN, alpha_EN, maxiter = 5000, eps = 2.2204e-16, tol = 1e-4, verbose = TRUE, Alg, stepsize){
  
  ### Initial Setup
  
  n <- nrow(Y)
  q <- ncol(Y)
  p <- ncol(X)
  
  ones <- matrix(1, nrow = n, ncol = 1)
  
  B <- matrix(0, nrow = p, ncol = q)
  mu <- matrix(0, nrow = 1, ncol = q)
  L <- matrix(0, nrow = n, ncol = q)
  
  fval_old <- fval <- 0
  f_val_vec <- c()
  nor_vec <- c()
  
  t_L <- 1
  t_B <- 1/(max(svd(X)$d)^2)
  t_mu <- 1/max(nrow(Y), ncol(Y)) ## should be 1/max(n,q), not n, especially when n < q
  
  ## initialize Lx0,Ly,t0,Bx0,By,mux0,muy
  
  Lx0 <- Lx <- L
  Bx0 <- Bx <- B
  mux0 <- mux <- mu
  t0 <- t <- 1
  
  #### Update B, mu, and L
  
  ##=====================================
  
  if(Alg == "FISTA"){
     	
    
    if(stepsize == "constant"){
      
      iter_LORS <- 1
      stopp <- FALSE	  
	  
      for(iter in iter_LORS:maxiter){
        
        print("Beginning LORSEN--PGD Updates") 
		
        t0 <- t            
		
		if((iter > iter_LORS)&(iter > 1)){
          fval_old <- fval   
        }
		
        t <- (1 + sqrt(1 + 4 * t0^2))/2 
        Ly <- Lx + (t0 - 1)/t * (Lx - Lx0) 
        By <- Bx + (t0 - 1)/t * (Bx - Bx0)
        muy <- mux + (t0 - 1)/t * (mux - mux0)
		
        Lx0 <- Lx  
        Bx0 <- Bx
        mux0 <- mux
		
        Lx <- prox_tr(X,By,muy,Ly,Y,lambda,t_L)   
        Bx <- prox_EN(X,By,muy,Lx,Y,lambda_EN,alpha_EN,t_B)
        mux <- prox_null(X,Bx,muy,Lx,Y,t_mu)   
		
        #### Check Convergence    
        
        fval <- FF_EN(X,Bx,Lx,mux,Y,lambda,alpha_EN,lambda_EN) 
        
        nor = abs(fval-fval_old)/abs(fval_old+eps)  
        
        f_val_vec <- c(f_val_vec, fval) 
        nor_vec <- c(nor_vec, nor) 
        
        print(paste('Iter ', iter, 'fval', fval, 'nor', nor))
        
        if (nor < tol){
          break
        }
        
        if(iter > iter_LORS){
          
          if(f_val_vec[iter] > f_val_vec[iter-1]){
            
			print("Beginning LORSEN--glmnet Updates")
			
            Lx1 <- Lx	
			Bx1 <- Bx
			mux1 <- mux
			
            fval1 <- fval 
            
            L0 <- L1 <- Lx0  
            B0 <- B1 <- Bx0
			mu0 <- mu1 <- mux0
            
            for(iter_LORS in iter:maxiter){
                
              #### Compute L
              
			  svd_st_res <- svd_st(Y - X%*%B1 - ones %*% mu1,lambda)  
              Lx <- svd_st_res$L   
			  
			  for(j in 1:q){
                fit <- glmnet(X, matrix(Y[,j]) - matrix(Lx[,j]), family = "gaussian", alpha = alpha_EN, lambda = lambda_EN/n, standardize = FALSE)
                a0 <- fit[["a0"]]
                old_beta <- fit[["beta"]]
                my_beta <- as(old_beta, "matrix")
                Bx[,j] <- my_beta
                mux[,j] <- a0
              }
			  
              #### Check Convergence    
              
              fval <- FF_EN(X,B,L,mu,Y,lambda,alpha_EN,lambda_EN)  
              nor <- abs(fval-fval_old)/abs(fval_old+eps) 
              
              if(iter_LORS == iter){
                
				print(paste0("CD-fval: ",fval," ---- PGD-fval: ",fval1))
				print(paste0("fval < fval1: ",fval < fval1))
				
                if(fval < fval1){     
                  
                  f_val_vec <- f_val_vec[-iter] 
                  nor_vec <- nor_vec[-iter]  
				  
				  L1 <- Lx
				  B1 <- Bx
				  mu1 <- mux
				  fval_old <- fval 
                  
                }else{
                  
                  iter_LORS <- iter + 1   
                  Lx <- Lx1      
				  Bx <- Bx1
				  mux <- mux1
				  fval_old <- fval1
                  break
                }			
              }
              
              f_val_vec <- c(f_val_vec, fval) 
              nor_vec <- c(nor_vec, nor) 
              
              print(paste('iter_LORS ', iter_LORS, 'fval', fval, 'nor', nor))
              
              if (nor < tol){
                stopp <- TRUE
                break
              }else{
                stopp <- FALSE
              }
              
              if((f_val_vec[iter_LORS] > f_val_vec[iter_LORS-1])&(iter_LORS > iter)){
                
				f_val_vec <- f_val_vec[-iter_LORS] 
                nor_vec <- nor_vec[-iter_LORS] 
				
                Lx <- L1 
                Lx0 <- L0 
                mux <- mu1
				mux0 <- mu0
				Bx <- B1
				Bx0 <- B0
			
                break
                
              }else if((f_val_vec[iter_LORS] <= f_val_vec[iter_LORS-1])&(iter_LORS > iter)){
			  
                L0 <- L1 
                L1 <- Lx 
				B0 <- B1
				B1 <- Bx
				mu0 <- mu1
				mu1 <- mux
                fval_old <- fval 
              }      
            }
          }
        }
		
        if(stopp){
          break
        }
      }
    }       
    
    if(stepsize == "backtracking"){
      
      eta_L <- eta_B <- eta_mu <- 0.5
      iter_LORS <- 1
      stopp <- FALSE	  
      for(iter in iter_LORS:maxiter){
        
        print("Beginning BerhuLORS--PGD Updates") 
		
        t0 <- t            
		
		if((iter > iter_LORS)&(iter > 1)){
          fval_old <- fval   
        }
		
        t <- (1 + sqrt(1 + 4 * t0^2))/2 
        Ly <- Lx + (t0 - 1)/t * (Lx - Lx0) 
        By <- Bx + (t0 - 1)/t * (Bx - Bx0)
        muy <- mux + (t0 - 1)/t * (mux - mux0)
		
        Lx0 <- Lx  
        Bx0 <- Bx
        mux0 <- mux
		
        Lx <- Update_L(X,By,muy,Ly,Y,lambda,t_L,eta_L)  
        Bx <- Update_B_EN(X,By,muy,Lx,Y,lambda_EN,alpha_EN,t_B,eta_B)
        mux <- Update_mu(X,Bx,muy,Lx,Y,t_mu,eta_mu)  
		
        #### Check Convergence    
        
        fval <- FF_EN(X,Bx,Lx,mux,Y,lambda,alpha_EN,lambda_EN) 
        
        nor = abs(fval-fval_old)/abs(fval_old+eps) # new-to-old ratio  
        
        f_val_vec <- c(f_val_vec, fval) 
        nor_vec <- c(nor_vec, nor) 
        
        print(paste('Iter ', iter, 'fval', fval, 'nor', nor))
        
        if (nor < tol){
          break
        }
        
        if(iter > iter_LORS){
          
          if(f_val_vec[iter] > f_val_vec[iter-1]){
            
			print("Beginning LORSEN--glmnet Updates")
			
            Lx1 <- Lx	
			Bx1 <- Bx
			mux1 <- mux
			
            fval1 <- fval 
            
            L0 <- L1 <- Lx0  
            B0 <- B1 <- Bx0
			mu0 <- mu1 <- mux0
            
            for(iter_LORS in iter:maxiter){
                
              #### Compute L
              
			  svd_st_res <- svd_st(Y - X%*%B1 - ones %*% mu1,lambda)  
              Lx <- svd_st_res$L   
			  
			  for(j in 1:q){
                fit <- glmnet(X, matrix(Y[,j]) - matrix(Lx[,j]), family = "gaussian", alpha = alpha_EN, lambda = lambda_EN/n, standardize = FALSE)
                a0 <- fit[["a0"]]
                old_beta <- fit[["beta"]]
                my_beta <- as(old_beta, "matrix")
                Bx[,j] <- my_beta
                mux[,j] <- a0
              }
			  
              #### Check Convergence    
              
              fval <- FF_EN(X,B,L,mu,Y,lambda,alpha_EN,lambda_EN)  
              nor <- abs(fval-fval_old)/abs(fval_old+eps) 
              
              if(iter_LORS == iter){
                
                if(fval < fval1){     
                  
                  f_val_vec <- f_val_vec[-iter] 
                  nor_vec <- nor_vec[-iter]  
				  
				  L1 <- Lx
				  B1 <- Bx
				  mu1 <- mux
				  fval_old <- fval 
                  
                }else{
                  
                  iter_LORS <- iter + 1   
                  Lx <- Lx1      
				  Bx <- Bx1
				  mux <- mux1
				  fval_old <- fval1
                  break
                }			
              }
              
              f_val_vec <- c(f_val_vec, fval) 
              nor_vec <- c(nor_vec, nor) 
              
              print(paste('iter_LORS ', iter_LORS, 'fval', fval, 'nor', nor))
              
              if (nor < tol){
                stopp <- TRUE
                break
              }else{
                stopp <- FALSE
              }
              
              if((f_val_vec[iter_LORS] > f_val_vec[iter_LORS-1])&(iter_LORS > iter)){
                
				f_val_vec <- f_val_vec[-iter_LORS] 
                nor_vec <- nor_vec[-iter_LORS] 
				
                Lx <- L1 
                Lx0 <- L0 
                mux <- mu1
				mux0 <- mu0
				Bx <- B1
				Bx0 <- B0
			
                break
                
              }else if((f_val_vec[iter_LORS] <= f_val_vec[iter_LORS-1])&(iter_LORS > iter)){
			  
                L0 <- L1 
                L1 <- Lx 
				B0 <- B1
				B1 <- Bx
				mu0 <- mu1
				mu1 <- mux
                fval_old <- fval 
              }      
            }
          }
        }
		
        if(stopp){
          break
        }
      }
    }   			
  }
  
  if(iter >= iter_LORS){
    
    stepp <- iter
    
  }else{
    
    stepp <- iter_LORS
  }
  ### Return the B, L, mu, objective function values, residual values, and total iterates
  return(list("B" = Bx, "L" = Lx, "mu" = mux, "f_val_vec" = f_val_vec, "nor_vec" = nor_vec, "iter" = length(f_val_vec), "step" = stepp))
}

#' LORSEN_Tuning: this function is used in parameter tuning of LORSEN

LORSEN_Tuning <- function(Y, X, lambda, lambda_EN, alpha_EN, Training, Validation, maxiter = 5000, eps = 2.2204e-16, tol = 1e-4, verbose = TRUE, L, Alg, stepsize){
  
  ### Initial Setup
  
  n <- nrow(Y)
  q <- ncol(Y)
  p <- ncol(X)
  
  ones <- matrix(1, nrow = n, ncol = 1)
  B <- matrix(0, nrow = p, ncol = q)
  mu <- matrix(0, nrow = 1, ncol = q)
  
  if(is.null(L) == TRUE){
    L <- matrix(0, nrow = n, ncol = q)
  }
  
  fval_old <- fval <- 0
  f_val_vec <- c()
  nor_vec <- c()
  
  t_L <- 1
  t_B <- 1/(max(svd(X)$d)^2)
  t_mu <- 1/nrow(Y)
  
  ## initialize Lx0,Ly,t0,Bx0,By,mux0,muy
  
  Lx0 <- Lx <- L
  Bx0 <- Bx <- B
  mux0 <- mux <- mu
  t0 <- t <- 1
  
  #### Update B, mu, and L
  
  ##=====================================
  
  if(Alg == "FISTA"){
    
    if(stepsize == "constant"){
      
	  iter_LORS <- 1
	  stopp <- FALSE
      for(iter in iter_LORS:maxiter){
	  
        print("Beginning LORSEN--PGD Updates")
		
		t0 <- t            
		
		if((iter > iter_LORS)&(iter > 1)){
          fval_old <- fval   
        }
		
        t <- (1 + sqrt(1 + 4 * t0^2))/2 
        Ly <- Lx + (t0 - 1)/t * (Lx - Lx0) 
        By <- Bx + (t0 - 1)/t * (Bx - Bx0)
        muy <- mux + (t0 - 1)/t * (mux - mux0)
		
        Lx0 <- Lx  
        Bx0 <- Bx
        mux0 <- mux
		
        Lx <- prox_tr_Tuning(X,By,muy,Ly,Y,lambda,t_L,Training)   
        Bx <- prox_EN_Tuning(X,By,muy,Lx,Y,lambda_EN,alpha_EN,t_B,Training)
        mux <- prox_null_Tuning(X,Bx,muy,Lx,Y,t_mu,Training)  
        
        #### Check Convergence    
        
        fval <- FF_EN_Tuning(X,Bx,Lx,mux,Y,lambda,alpha_EN,lambda_EN,Training)
		
        nor <- abs(fval-fval_old)/abs(fval_old+eps) # new-to-old ratio
        
        f_val_vec <- c(f_val_vec, fval)
        nor_vec <- c(nor_vec, nor)
        
        print(paste('Iter ', iter, 'fval', fval, 'nor', nor))
        
        if (nor < tol){
          break
        }
        
        if(iter > iter_LORS){
          if(f_val_vec[iter] > f_val_vec[iter-1]){
        
		    print("Beginning LORSEN--CD Updates") 
			
            Lx1 <- Lx
            Bx1 <- Bx
            mux1 <- mux			
			fval1 <- fval
			
            L0 <- L1 <- Lx0
            B0 <- B1 <- Bx0
            mu0 <- mu1 <- mux0
            
			for (iter_LORS in iter:maxiter){
         
		        energy = Inf
				
				Z = Y - X %*% B1 - ones %*% mu1
				LL <- L1
                for (innerIts in 1:50){
                    
                    C = Z*Training + LL*(!Training)
            
                    svd_st_res <- svd_st(C,lambda)
                    Lx <- svd_st_res$L
                    W <- svd_st_res$W
            
                    Lnorm = sum(abs(diag(W)))
                    energy_old = energy
            
                    cL = c(Lx)
                    cZ = c(Z)
                    newvec = matrix(cL[c(Training == 1)]-cZ[c(Training == 1)], ncol = 1)
                    energy = lambda * Lnorm + norm(newvec,'F')**2 / 2
            
                    mydiff = abs(energy - energy_old) / abs(energy_old)
            
                    if (!is.na(mydiff)){
                       if (abs(energy - energy_old) / abs(energy_old) < tol){
                          break
                        }
                    }
					
					LL <- Lx
                }
			     
				for(j in 1:q){
                   fit <- glmnet(X[Training[,j],], matrix(Y[Training[,j],j]) - matrix(Lx[Training[,j],j]), family = "gaussian", alpha = alpha_EN, lambda = lambda_EN/sum(Training[,j]), standardize = FALSE)
                   a0 <- fit[["a0"]]
                   old_beta <- fit[["beta"]]
                   my_beta <- as(old_beta, "matrix")
                   Bx[,j] <- my_beta
                   mux[,j] <- a0
                }
          
                # Convergence
          
                fval <- FF_EN_Tuning(X,B,L,mu,Y,lambda,alpha_EN,lambda_EN,Training)
                nor <- abs(fval-fval_old)/abs(fval_old+eps) # new-to-old ratio
          
		        if(iter_LORS == iter){
                
                  if(fval < fval1){     
                  
                    f_val_vec <- f_val_vec[-iter] 
                    nor_vec <- nor_vec[-iter]  
				  
				    L1 <- Lx
				    B1 <- Bx
				    mu1 <- mux
				    fval_old <- fval 
                  
                  }else{
                  
                    iter_LORS <- iter + 1   
                    Lx <- Lx1      
				    Bx <- Bx1
				    mux <- mux1
				    fval_old <- fval1
                    break
                  }			
                }
		  
                f_val_vec <- c(f_val_vec, fval)
                nor_vec <- c(nor_vec, nor)
          
			    print(paste('iter_LORS ', iter_LORS, 'fval', fval, 'nor', nor))
          
                if (nor < tol){
				   stopp <- TRUE
                   break
                }else{
				   stopp <- FALSE
				}
			
			    if((f_val_vec[iter_LORS] > f_val_vec[iter_LORS-1])&(iter_LORS > iter)){
                
				  f_val_vec <- f_val_vec[-iter_LORS] 
                  nor_vec <- nor_vec[-iter_LORS] 
				
                  Lx <- L1 
                  Lx0 <- L0 
                  mux <- mu1
				  mux0 <- mu0
				  Bx <- B1
				  Bx0 <- B0
			
                  break
                
                }else if((f_val_vec[iter_LORS] <= f_val_vec[iter_LORS-1])&(iter_LORS > iter)){
			  
                  L0 <- L1 
                  L1 <- Lx 
				  B0 <- B1
				  B1 <- Bx
				  mu0 <- mu1
				  mu1 <- mux
                  fval_old <- fval 
                }  
            }
          }       
        } 
         
        if(stopp){
		   break
		}		 
      }
    }  
   
        
    if(stepsize == "backtracking"){
	
      eta_L <- eta_B <- eta_mu <- 0.5
	  iter_LORS <- 1
	  stopp <- FALSE
      for(iter in iter_LORS:maxiter){
	  
        print("Beginning LORSEN--PGD Updates")
		
		t0 <- t            
		
		if((iter > iter_LORS)&(iter > 1)){
          fval_old <- fval   
        }
		
        t <- (1 + sqrt(1 + 4 * t0^2))/2 
        Ly <- Lx + (t0 - 1)/t * (Lx - Lx0) 
        By <- Bx + (t0 - 1)/t * (Bx - Bx0)
        muy <- mux + (t0 - 1)/t * (mux - mux0)
		
        Lx0 <- Lx  
        Bx0 <- Bx
        mux0 <- mux
		
        Lx <- Update_L_Tuning(X,By,muy,Ly,Y,lambda,t_L,eta_L,Training)  
        Bx <- Update_B_EN_Tuning(X,By,muy,Lx,Y,lambda_EN,alpha_EN,t_B,eta_B,Training)
        mux <- Update_mu_Tuning(X,Bx,muy,Lx,Y,t_mu,eta_mu,Training)
        
        #### Check Convergence    
        
        fval <- FF_EN_Tuning(X,Bx,Lx,mux,Y,lambda,alpha_EN,lambda_EN,Training)
		
        nor <- abs(fval-fval_old)/abs(fval_old+eps) # new-to-old ratio
        
        f_val_vec <- c(f_val_vec, fval)
        nor_vec <- c(nor_vec, nor)
        
        print(paste('Iter ', iter, 'fval', fval, 'nor', nor))
        
        if (nor < tol){
          break
        }
        
        if(iter > iter_LORS){
          if(f_val_vec[iter] > f_val_vec[iter-1]){
        
		    print("Beginning LORSEN--CD Updates") 
			
            Lx1 <- Lx
            Bx1 <- Bx
            mux1 <- mux			
			fval1 <- fval
			
            L0 <- L1 <- Lx0
            B0 <- B1 <- Bx0
            mu0 <- mu1 <- mux0
            
			for (iter_LORS in iter:maxiter){
         
		        energy = Inf
				
				Z = Y - X %*% B1 - ones %*% mu1
				LL <- L1
                for (innerIts in 1:50){
                    
                    C = Z*Training + LL*(!Training)
            
                    svd_st_res <- svd_st(C,lambda)
                    Lx <- svd_st_res$L
                    W <- svd_st_res$W
            
                    Lnorm = sum(abs(diag(W)))
                    energy_old = energy
            
                    cL = c(Lx)
                    cZ = c(Z)
                    newvec = matrix(cL[c(Training == 1)]-cZ[c(Training == 1)], ncol = 1)
                    energy = lambda * Lnorm + norm(newvec,'F')**2 / 2
            
                    mydiff = abs(energy - energy_old) / abs(energy_old)
            
                    if (!is.na(mydiff)){
                       if (abs(energy - energy_old) / abs(energy_old) < tol){
                          break
                        }
                    }
					
					LL <- Lx
                }
			     
				for(j in 1:q){
                   fit <- glmnet(X[Training[,j],], matrix(Y[Training[,j],j]) - matrix(Lx[Training[,j],j]), family = "gaussian", alpha = alpha_EN, lambda = lambda_EN/sum(Training[,j]), standardize = FALSE)
                   a0 <- fit[["a0"]]
                   old_beta <- fit[["beta"]]
                   my_beta <- as(old_beta, "matrix")
                   Bx[,j] <- my_beta
                   mux[,j] <- a0
                }
          
                # Convergence
          
                fval <- FF_EN_Tuning(X,B,L,mu,Y,lambda,alpha_EN,lambda_EN,Training)
                nor <- abs(fval-fval_old)/abs(fval_old+eps) # new-to-old ratio
          
		        if(iter_LORS == iter){
                
                  if(fval < fval1){     
                  
                    f_val_vec <- f_val_vec[-iter] 
                    nor_vec <- nor_vec[-iter]  
				  
				    L1 <- Lx
				    B1 <- Bx
				    mu1 <- mux
				    fval_old <- fval 
                  
                  }else{
                  
                    iter_LORS <- iter + 1   
                    Lx <- Lx1      
				    Bx <- Bx1
				    mux <- mux1
				    fval_old <- fval1
                    break
                  }			
                }
		  
                f_val_vec <- c(f_val_vec, fval)
                nor_vec <- c(nor_vec, nor)
          
			    print(paste('iter_LORS ', iter_LORS, 'fval', fval, 'nor', nor))
          
                if (nor < tol){
				   stopp <- TRUE
                   break
                }else{
				   stopp <- FALSE
				}
			
			    if((f_val_vec[iter_LORS] > f_val_vec[iter_LORS-1])&(iter_LORS > iter)){
                
				  f_val_vec <- f_val_vec[-iter_LORS] 
                  nor_vec <- nor_vec[-iter_LORS] 
				
                  Lx <- L1 
                  Lx0 <- L0 
                  mux <- mu1
				  mux0 <- mu0
				  Bx <- B1
				  Bx0 <- B0
			
                  break
                
                }else if((f_val_vec[iter_LORS] <= f_val_vec[iter_LORS-1])&(iter_LORS > iter)){
			  
                  L0 <- L1 
                  L1 <- Lx 
				  B0 <- B1
				  B1 <- Bx
				  mu0 <- mu1
				  mu1 <- mux
                  fval_old <- fval 
                }  
            }
          }       
        } 
         
        if(stopp){
		   break
		}		 
      }
    }   
            
  }  
  
  residual = Y - X%*%Bx - ones %*% mux - Lx
  err = residual*Validation
  Err = t(c(err)) %*% c(err) / sum(Validation)
  
  if(iter >= iter_LORS){
  
    stepp <- iter
	
  }else{
  
    stepp <- iter_LORS
  }
  
  ### Return the B, L, mu, objective function values, residual values, and total iterates
  return(list("B" = Bx, "L" = Lx, "mu" = mux, "Err" = Err, "f_val_vec" = f_val_vec, "nor_vec" = nor_vec, "iter" = length(f_val_vec), "step" = stepp))
}

### Update B

Update_B_EN <- function(X,B,mu,L,Y,lambda_EN,alpha_EN,t_B,eta_B){

      i_k <- 1
	  arg1 <- prox_EN(X,B,mu,L,Y,lambda_EN,alpha_EN,t_B)
	  arg2 <- B
	  
	  F_B_val <- F(X,arg1,L,mu,Y)
	  Q_B_val <- Q_B(arg1,arg2,X,L,mu,Y,t_B)
	  
      while(F_B_val > Q_B_val){
	  
	       i_k <- i_k + 1
	       t_B <- eta_B^{i_k} * t_B
		   arg1 <- prox_EN(X,B,mu,L,Y,lambda_EN,alpha_EN,t_B)
	       
	       F_B_val <- F(X,arg1,L,mu,Y)  
	       Q_B_val <- Q_B(arg1,arg2,X,L,mu,Y,t_B)   
	  }
	  
	  return(arg1)
}


### Update B

Update_B_EN_Tuning <- function(X,B,mu,L,Y,lambda_EN,alpha_EN,t_B,eta_B,Training){
      
      i_k <- 1
	  arg1 <- prox_EN_Tuning(X,B,mu,L,Y,lambda_EN,alpha_EN,t_B,Training)
	  arg2 <- B
	  
	  F_B_val <- F_Tuning(X,arg1,L,mu,Y,Training)
	  
	  Q_B_val <- Q_B_Tuning(arg1,arg2,X,L,mu,Y,t_B,Training)
	  
      while(F_B_val > Q_B_val){
	  
	       i_k <- i_k + 1
	       t_B <- eta_B^{i_k} * t_B
		   arg1 <- prox_EN_Tuning(X,B,mu,L,Y,lambda_EN,alpha_EN,t_B,Training)
	       
	       F_B_val <- F_Tuning(X,arg1,L,mu,Y,Training)  
	       Q_B_val <- Q_B_Tuning(arg1,arg2,X,L,mu,Y,t_B,Training)
	 
	  }
	  return(arg1)

}


FF_EN_Tuning <- function(X,B,L,mu,Y,lambda,alpha_EN,lambda_EN,Training){

     lambda1 <- lambda_EN * alpha_EN
	 lambda2 <- (1 - alpha_EN) * lambda_EN
	 
     res <- F_Tuning(X,B,L,mu,Y,Training) + lambda * tr_norm(L) + lambda1 * l1_norm(B) + 0.5 * lambda2 * l2_norm(B)
	 return(res)
}

FF_EN <- function(X,B,L,mu,Y,lambda,alpha_EN,lambda_EN){

     lambda1 <- lambda_EN * alpha_EN
	 lambda2 <- (1 - alpha_EN) * lambda_EN
	 
     res <- F(X,B,L,mu,Y) + lambda * tr_norm(L) + lambda1 * l1_norm(B) + 0.5 * lambda2 * l2_norm(B)
	 return(res)
}

prox_EN_Tuning <- function(X,B,mu,L,Y,lambda_EN,alpha_EN,t_B,Training){

    lambda1 <- lambda_EN * alpha_EN
	lambda2 <- (1 - alpha_EN) * lambda_EN
	
    L1_res <- prox_1(B - t_B * t(X) %*% (Training * (X %*% B + ones %*% mu + L - Y)), lambda1 * t_B)

    L2_res <- c()

    for(i in 1:ncol(Y)){

       L2 <- prox_L2(L1_res[,i], lambda2 * t_B)
       L2_res <- cbind(L2_res,L2)
    }
	return(L2_res)
}

prox_EN <- function(X,B,mu,L,Y,lambda_EN,alpha_EN,t_B){

    lambda1 <- lambda_EN * alpha_EN
	lambda2 <- (1 - alpha_EN) * lambda_EN
	
    L1_res <- prox_1(B - t_B * t(X) %*% (X %*% B + ones %*% mu + L - Y), lambda1 * t_B)

    L2_res <- c()

    for(i in 1:ncol(Y)){

       L2 <- prox_L2(L1_res[,i], lambda2 * t_B)
       L2_res <- cbind(L2_res,L2)
    }
	return(L2_res)
}


prox_L2 <- function(v, lambda2){

    res <- (1 - lambda2/max(norm(v,"2"),lambda2)) * v
    return(res)	
}

## soft thresholding operator

STO <- function(a,b){

    if(a>0 & b < abs(a)){
	
	   res <- a - b
	}
	
	if(a < 0 & b < abs(a)){
	
	   res <- a + b
	}
	
	if(b >= abs(a)){
	
	   res <- 0
	}
	
	return(res)
}

## calculate L2 norm

calnorm2 <- function(vec){

   return(norm(vec, type = "2")**2)
}

### softImpute，used to obtain the best lambda [soft imputing algorithm]

softImpute <- function(X,Z,Omega,Omega1,Omega2,alpha0,maxRank){
  
  X_0 = X*Omega
  if (is.null(Z) == TRUE) {
    Z = X_0
  }
  
  if (is.null(alpha0) == TRUE){
  
    my_svd = svd(X_0)
    DD = diag(my_svd[["d"]])
    alpha0 = DD[2,2]
  }
  
  if (is.null(maxRank) == TRUE){
    maxRank = -1
  }
  
  # parameters
  eta = 0.9
  epsilon = 1e-4
  maxInnerIts = 50
  
  ## soft-impute
  
  # 1. initialize
  alpha = alpha0
  
  Err = c()
  rank_alpha = c()
  znorm = c()
  Alpha = c()
  
  while (TRUE){
    energy = Inf
    for (innerIts in 1:maxInnerIts){
      # (a)i
      C = X*Omega1 + Z*(!Omega1)
      C.Obj <- svd_st(C,alpha)
      Z <- C.Obj$L
	  W <- C.Obj$W
	  
      Znorm = sum(abs(diag(W)))
      energy_old = energy
      cZ = c(Z)
      cX = c(X)
      newvec = matrix(cZ[c(Omega1 == 1)]-cX[c(Omega1 == 1)], ncol = 1)
      energy = alpha*Znorm + norm(newvec, 'F')**2 / 2
      
      mydiff = abs(energy - energy_old) / energy_old
      if (!is.na(mydiff)){
        if (abs(energy - energy_old) / energy_old < epsilon){
          break
        }
      }
    }
    
    e = X * Omega2 - Z * Omega2
    err2 = t(c(e)) %*%  c(e)
    Err = cbind(Err, err2)
    znorm = cbind(znorm, Znorm)
    k = length(diag(W))
    rank_alpha = cbind(rank_alpha, k)
    Alpha = cbind(Alpha, alpha)
    
    if (k <= maxRank && alpha > 1e-3){
      alpha = alpha * eta
    }else{
      break
    }
    
  }
  return(list("Z" = Z, "Err" = Err, "rank_alpha" = rank_alpha, "znorm" = znorm, "Alpha" = Alpha))
}


### svd_st: soft thresholding SVD

svd_st <- function(X, lambda){
    mysvd <- svd(X)
    U <- mysvd[["u"]]
    V <- mysvd[["v"]]
    VT <- t(V)
    D <- diag(mysvd[["d"]])
    d <- mysvd[["d"]]
    index_list <- which(d >= lambda)
    if(length(index_list) > 1){
      W <- diag(c(d[index_list] - lambda))
      L <- U[,index_list] %*% W %*% VT[index_list,]
    }
    if(length(index_list) == 1){
      W <- matrix(d[index_list] - lambda, nrow = 1, ncol = 1)
      L <- U[,index_list] %*% W %*% VT[index_list,]
    }
    if(length(index_list) == 0){
      L <- matrix(0, nrow(X), ncol(X))
	  W <- 0
    }
    
    return(list(L = L, W = W))
}

## logspace is a Matlab built-in function

logspace <- function (x1, x2, n = 50) {
  if (x2 == pi){
    x2 <- log10(x2)
  }
  10^seq(x1, x2, length.out = n)
}

### LORSscreen, used to estimate \beta_{i} for each SNP, used in LORS-Screening

LORSscreen <- function(Y, X, lambda, tol){
  
  eps = 2.2204e-16
  
  n = nrow(Y)
  q = ncol(Y)
  
  B = matrix(0, nrow = 1, ncol=q)
  mu = matrix(0, nrow=1, ncol=q)
  
  maxIter = 100
  
  fval_old = 0
  
  for (iter in 1:maxIter){
    
    #### Compute L
    svd_st_res <- svd_st(Y - X%*%B - ones %*% mu, lambda)
	L <- svd_st_res$L
	W <- svd_st_res$W
    
    #### Solve the Least Squares Problem
    
    Z <- Y - L
    A <- cbind(X,1)
    myqr <- qr(A)
    Q <- qr.Q(myqr)
    R <- qr.R(myqr)
    B_QR <- ginv(t(R) %*% R) %*% t(A) %*% Z
    
    B <- matrix(B_QR[1,], nrow = 1, ncol = q)
    mu <- matrix(B_QR[2,], nrow = 1, ncol = q)
    
    dum = Y - X%*%B - ones%*%mu - L
    dum = c(dum)
    fval = 0.5 * t(dum) %*% dum + lambda*sum(abs(diag(W)))
    
    res = abs(fval-fval_old)/abs(fval_old+eps)
    
    if (res < tol){
      break
    }
    
    fval_old <- fval
  }
  return (list("B" = B, "L" = L, "mu" = mu))
}


## InitialEst, obtain initial estimate of B

InitialEst <- function(Y, X, lambda = NULL){
    
  n <- nrow(Y)
  q <- ncol(Y)
  p <- ncol(X) 
	
  if(is.null(lambda) == TRUE){
    # first use soft-impute to select a reasonable lambda
    
    mask = matrix(runif(n*q) > 0.5, nrow = n, ncol = q)
    Omega1 = Omega & mask
    Omega2 = Omega & !mask
    
    maxRank = min(n,q)/2
    mySI= softImpute(Y,NULL,Omega, Omega1, Omega2, NULL,maxRank)
    Err = mySI[["Err"]]
    Alpha = mySI[["Alpha"]]
    bestind_lam = min(Err)==Err;
    
    lambda = min(Alpha[bestind_lam])
  }
  
  myB <- matrix(NA, p, q)
  print("Building Initial Estimate")
  for(SNP_col in 1:ncol(X)){
    print(paste("Building initial estimate: On column", SNP_col,"of",ncol(X)))
    X1 <- matrix(X[,SNP_col], ncol = 1)
    LS <- LORSscreen(Y, X1, lambda, 0.01)
    myB[SNP_col,] <- LS$B
  }
  
  return(list("B"  = myB))
}


### standardizeBhat，standardization of B

standardizeBhat <- function(Y, X, Bhat){
  Bhat_standard <- Bhat
  print("Standardizing Initial Estimate")
  for(i in 1:ncol(X)){
    for(j in 1:ncol(Y)){
      var_B_ij <- (t(X[,i]) %*% X[,i])^(-1) * var(Y[,j] - X[,i] * Bhat[i,j])
      Bhat_standard[i,j] <- Bhat[i,j]/sqrt(var_B_ij)
    }
  }
  return(Bhat_standard)
}


### LORS2, used in parameter tuning in LORS

LORS2 <- function(Y, X, L, Omega1, Omega2, rho, lambda, tol, maxIter = 5000){
  
  #p <- ncol(X)
  B = matrix(0, nrow = ncol(X), ncol = q)
  mu = matrix(0, nrow=1, ncol=q)
  fval_old = 0
  
  maxInnerIts = 50
  
  #### initialization
  
  if (is.null(L) == 1){
    L = Y
  }
  
  energy = Inf
  
  for (iter in 1:maxIter){
    
    for (innerIts in 1:maxInnerIts){
      
      Z = Y - X %*% B - ones %*% mu;
      C = Z*Omega1 + L*(!Omega1)
      
      svd_st_res <- svd_st(C,lambda)
	  L <- svd_st_res$L
	  W <- svd_st_res$W
      
      Lnorm = sum(abs(diag(W)))
      energy_old = energy
      
      cL = c(L)
      cZ = c(Z)
      newvec = matrix(cL[c(Omega1 == 1)]-cZ[c(Omega1 == 1)], ncol = 1)
      energy = lambda * Lnorm + norm(newvec,'F')**2 / 2
      
      mydiff = abs(energy - energy_old) / energy_old
      
      if (!is.na(mydiff)){
        if (abs(energy - energy_old) / energy_old < tol){
          break
        }
      }
    }
    
    ### Compute B using glmnet
    
    for (j in 1:q){
      fit <- glmnet(X[Omega1[,j],], matrix(Y[Omega1[,j],j]) - matrix(L[Omega1[,j],j]), family = "gaussian", lambda = rho/sum(Omega1[,j]), standardize = FALSE)
      a0 <- fit[["a0"]]
      old_beta <- fit[["beta"]]
      my_beta <- as(old_beta, "matrix")
      B[,j] <- my_beta
      mu[,j] <- a0
    }
    
    # Convergence
    residual = Y - X%*%B - ones %*% mu - L
    dum = residual*Omega1
    dum = c(dum)
    fval = 0.5 * t(dum) %*% dum + rho * sum(abs(c(B))) + lambda*sum(abs(diag(W)))
    res = abs(fval-fval_old)/abs(fval_old+eps)
    
    if (res < tol){
      break
    }
    
    fval_old <- fval
    
  }
  
  err = residual*Omega2
  Err = t(c(err)) %*% c(err) / sum(Omega2)
  
  return (list("B" = B, "mu" = mu, "L" = L, "Err" = Err))
}


#' GetMaxRho
#' min 1/2 ||Y - XB - L -1\mu||^{2}_{F} + \rho ||B||_{1} + \lambda ||L||_{*} 

GetMaxRho <- function(X, Y, L, Omega){
 
  maxrho = matrix(0, nrow = q, ncol = 1)
  
  for(i in 1:q){
    maxrho[i,1] = max(abs(t(X) %*% ((matrix(Y[,i], nrow = nrow(Y), ncol = 1) - matrix(L[,i], nrow = nrow(L), ncol = 1)) * matrix(Omega[,i], nrow = nrow(Omega), ncol = 1))))
  }
  
  MaxRho = max(maxrho)
  return(MaxRho)
}


### prox_1, Soft Threshold Function for L1 norm

prox_1 <- function(b, tau){
    prox_b <- sign(b)*(sapply(abs(b) - tau, FUN=function(x) {max(x,0)}))
    return(prox_b)
}


### Fast_LORS_Tuning

Fast_LORS_Tuning <- function(Y, X, rho, lambda, Training, Validation, maxiter = 5000, eps = 2.2204e-16, tol = 1e-4, L, omega_SOR = 1.999) {
  
  #p <- ncol(X)
  B <- matrix(0, nrow = ncol(X), ncol = q)
  mu <- matrix(0, nrow = 1, ncol = q)
  
  if(is.null(L) == TRUE){
    L <- matrix(0, nrow = n, ncol = q)
  }
 
  fval_old <- 0
  f_val_vec <- c()
  res_vec <- c()
  
  t_L <- 1
  t_B <- 1/(max(svd(X)$d)^2)
  t_mu <- 1/nrow(Y)
  
  for(iter in 1:maxiter){
    
    #### Update B, mu, and L
    L <- svd_st(L - t_L * (Training * (X %*% B + ones %*% mu + L - Y)), t_L * lambda)$L
    B <- prox_1(B - t_B * t(X) %*% (Training * (X %*% B + ones %*% mu + L - Y)), t_B * rho)
    mu <- mu - t_mu * t(ones) %*% (Training * (X %*% B + ones %*% mu + L - Y))
    
    ## Update via SOR
    
    if(iter > 1){
      L_update <- (1-omega_SOR) * L_update + omega_SOR * L
      B_update <- (1-omega_SOR) * B_update + omega_SOR * B
      mu_update <- (1-omega_SOR) * mu_update + omega_SOR * mu
      
      L <- L_update
      B <- B_update
      mu <- mu_update
    }
    
    if(iter == 1){
      L_update <- L
      B_update <- B
      mu_update <- mu
      
      L <- L_update
      B <- B_update
      mu <- mu_update
    }
    
    #### Check Convergence
    
    dum <- c(Training * (Y - X %*% B - ones%*%mu - L))
    fval <- 0.5 * norm(dum, type = "2")^2 + rho * sum(abs(B)) + lambda*sum(svd(L)[["d"]])
	
    residual = Y - X%*%B - ones %*% mu - L
	
    res = abs(fval-fval_old)/abs(fval_old+eps)
    
    fval_old <- fval
    f_val_vec <- c(f_val_vec, fval)
    res_vec <- c(res_vec, res)
    
    print(paste('Iter', iter, 'fval', fval, 'res', res))
    
    if (res < tol){
      break
    }
    
    if(iter > 1){
      if(f_val_vec[iter] > f_val_vec[iter-1]){
        break
      }
    }
 
  }
  
  if(iter < maxiter & res >= tol){
    print("Beginning FastLORS Updates")
    energy = Inf
    
    for (iter_LORS in (iter + 1):maxiter){
      
      for (innerIts in 1:50){
        
        Z = Y - X %*% B - ones %*% mu;
        C = Z*Training + L*(!Training)
        
        svd_st_res <- svd_st(C,lambda)
        L <- svd_st_res$L
		W <- svd_st_res$W
		
        Lnorm = sum(diag(W))
        energy_old = energy
        
        cL = c(L)
        cZ = c(Z)
        newvec = matrix(cL[c(Training == 1)]-cZ[c(Training == 1)], ncol = 1)
        energy = lambda * Lnorm + norm(newvec,'F')**2 /2
        
        mydiff = abs(energy - energy_old) / energy_old
        
        if (!is.na(mydiff)){
          if (abs(energy - energy_old) / energy_old < tol){
            break
          }
        }
      }
      
      ### Compute B using glmnet
      
      for (j in 1:q){
        fit <- glmnet(X[Training[,j],], matrix(Y[Training[,j],j]) - matrix(L[Training[,j],j]), family = "gaussian", lambda = rho/sum(Training[,j]), standardize = FALSE)
        a0 <- fit[["a0"]]
        old_beta <- fit[["beta"]]
        my_beta <- as(old_beta, "matrix")
        B[,j] <- my_beta
        mu[,j] <- a0
      }
      
      # Convergence
      residual = Y - X%*%B - ones %*% mu - L
      dum = residual*Training
      dum = c(dum)
      fval = 0.5 * t(dum) %*% dum + rho * sum(abs(c(B))) + lambda*sum(abs(diag(W)))
      res = abs(fval-fval_old)/abs(fval_old+eps)
      
      fval_old <- fval
      f_val_vec <- c(f_val_vec, fval)
      res_vec <- c(res_vec, res)
      
      print(paste('Iter ', iter_LORS, 'fval', fval, 'res', res))
      
      if (res < tol){
        break
      }
      
    }
     
  }
  
  err = residual*Validation
  Err = t(c(err)) %*% c(err) / sum(Validation)
  
  ### Return the B, L, mu, objective function values, residual values, and total iterates
  return(list("B" = B, "L" = L, "mu" = mu, "Err" = Err, "f_val_vec" = f_val_vec, "res_vec" = res_vec, "iter" = iter))
}



### Fast_LORS

Fast_LORS <- function(Y, X, rho, lambda, maxiter = 5000, eps = 2.2204e-16, tol = 1e-4, verbose = TRUE, omega_SOR = 1.999) {
  
  ### Initial Setup
  #p <- ncol(X)
  B <- matrix(0, nrow = ncol(X), ncol = q)
  mu <- matrix(0, nrow = 1, ncol = q)
  L <- matrix(0, nrow = n, ncol = q)
  
  
  fval_old <- 0
  f_val_vec <- c()
  res_vec <- c()
  
  t_L <- 1
  t_B <- 1/(max(svd(X)$d)^2)
  t_mu <- 1/nrow(Y)
  
  ones <- matrix(1, nrow = n, ncol = 1)
  
  if(verbose == TRUE){
    for(iter in 1:maxiter){
      
      #### Update B, mu, and L
      L <- svd_st(L - t_L * (X %*% B + ones %*% mu + L - Y), t_L * lambda)$L
      B <- prox_1(B - t_B * t(X) %*% (X %*% B + ones %*% mu + L - Y), t_B * rho)
      mu <- mu - t_mu * t(ones) %*% (X %*% B + ones %*% mu + L - Y)
      
      ## Update via SOR
      
      if(iter > 1){
        L_update <- (1-omega_SOR) * L_update + omega_SOR * L
        B_update <- (1-omega_SOR) * B_update + omega_SOR * B
        mu_update <- (1-omega_SOR) * mu_update + omega_SOR * mu
        
        L <- L_update
        B <- B_update
        mu <- mu_update
      }
      
      if(iter == 1){
        L_update <- L
        B_update <- B
        mu_update <- mu
        
        L <- L_update
        B <- B_update
        mu <- mu_update
      }
      
      #### Check Convergence
      
      dum <- c(Y - X %*% B - ones%*%mu - L)
      fval <- 0.5 * norm(dum, type = "2")^2 + rho * sum(abs(B)) + lambda*sum(svd(L)[["d"]])
      
      res = abs(fval-fval_old)/abs(fval_old+eps)
      
      print(paste('Iter ', iter, 'fval', fval, 'res', res))
      
      fval_old <- fval
      f_val_vec <- c(f_val_vec, fval)
      res_vec <- c(res_vec, res)
      
      if (res < tol){
        break
      }
      
      if(iter > 1){
        if(f_val_vec[iter] > f_val_vec[iter-1]){
          break
        }
      }
      
    }
    
    if(iter < maxiter & res >= tol){
      print("Beginning LORS Updates")
      for (iter_LORS in (iter+1):maxiter){
        
        #### Compute L
        svd_st_res <- svd_st(Y - X%*%B - ones %*% mu,lambda)
        W <- svd_st_res$W
        L <- svd_st_res$L
        
        for (j in 1:q){
          fit <- glmnet(X, Y[,j]-L[,j], family = "gaussian", lambda = rho/n, standardize = FALSE)
          a0 <- fit[["a0"]]
          old_beta <- fit[["beta"]]
          my_beta <- as(old_beta, "matrix")
          B[,j] <- my_beta
          mu[,j] <- a0
        }
        
        dum <- c(Y - X %*% B - ones %*% mu - L)
        fval <- 0.5 * t(dum) %*% dum + rho * sum(abs(B)) + lambda*sum(abs(diag(W)))
        
        res <- abs(fval-fval_old)/abs(fval_old+eps)
        
        print(paste('Iter ', iter_LORS, 'fval', fval, 'res', res))
        
        fval_old <- fval
        f_val_vec <- c(f_val_vec, fval)
        res_vec <- c(res_vec, res)
        
        if (res < tol){
          break
        }
        
      }
    }
  }
  
  if(verbose == FALSE){
    for(iter in 1:maxiter){
      
      #### Update B, mu, and L
      L <- svd_st(L - t_L * (X %*% B + ones %*% mu + L - Y), t_L * lambda)$L
      B <- prox_1(B - t_B * t(X) %*% (X %*% B + ones %*% mu + L - Y), t_B * rho)
      mu <- mu - t_mu * t(ones) %*% (X %*% B + ones %*% mu + L - Y)
      
      ## Update via SOR
      
      if(iter > 1){
        L_update <- (1-omega_SOR) * L_update + omega_SOR * L
        B_update <- (1-omega_SOR) * B_update + omega_SOR * B
        mu_update <- (1-omega_SOR) * mu_update + omega_SOR * mu
        
        L <- L_update
        B <- B_update
        mu <- mu_update
      }
      
      if(iter == 1){
        L_update <- L
        B_update <- B
        mu_update <- mu
        
        L <- L_update
        B <- B_update
        mu <- mu_update
      }
      
      #### Check Convergence
      
      dum <- c(Y - X %*% B - ones%*%mu - L)
      fval <- 0.5 * norm(dum, type = "2")^2 + rho * sum(abs(B)) + lambda*sum(svd(L)[["d"]])
      
      res = abs(fval-fval_old)/abs(fval_old+eps)
      
      fval_old <- fval
      f_val_vec <- c(f_val_vec, fval)
      res_vec <- c(res_vec, res)
      
      if (res < tol){
        break
      }
      
      if(iter > 1){
        if(f_val_vec[iter] > f_val_vec[iter-1]){
          break
        }
      }
      
    }
    
    if(iter < maxiter & res >= tol){
      for (iter_LORS in (iter+1):maxiter){
        
        #### Compute L
        svd_st_res <- svd_st(Y - X%*%B - ones %*% mu,lambda)
        W <- svd_st_res$W
        L <- svd_st_res$L
        
        for (j in 1:q){
          fit <- glmnet(X, Y[,j]-L[,j], family = "gaussian", lambda = rho/n, standardize = FALSE)
          a0 <- fit[["a0"]]
          old_beta <- fit[["beta"]]
          my_beta <- as(old_beta, "matrix")
          B[,j] <- my_beta
          mu[,j] <- a0
        }
        
        dum <- c(Y - X%*%B - ones%*%mu - L)
        fval <- 0.5 * t(dum) %*% dum + rho * sum(abs(B)) + lambda*sum(abs(diag(W)))
        
        res <- abs(fval-fval_old)/abs(fval_old+eps)
        
        fval_old <- fval
        f_val_vec <- c(f_val_vec, fval)
        res_vec <- c(res_vec, res)
        
        if (res < tol){
          break
        }
        
        
      }
    }
  }
  
  ### Return the B, L, mu, objective function values, residual values, and total iterates
  return(list("B" = B, "L" = L, "mu" = mu, "f_vals" = f_val_vec, "res_vec" = res_vec, "iter" = length(f_val_vec)))
}


### LORS0, used to solve the LORS optimization problem

LORS0 <- function(Y, X, rho, lambda, maxiter = 5000, eps = 2.2204e-16, tol = 1e-4, verbose = FALSE){
  
  B <- matrix(0, nrow = ncol(X), ncol=q)
  mu <- matrix(0, nrow=1, ncol=q)
  
  fval_old <- 0
  f_val_vec <- c()
  res_vec <- c()
  
  if(verbose == TRUE){
    for (iter in 1:maxiter){
      
      #### Compute L
      svd_st_res <- svd_st(Y - X %*% B - ones %*% mu,lambda)
      L <- svd_st_res$L
	  W <- svd_st_res$W
      
      for (j in 1:q){
        fit <- glmnet(X, Y[,j]-L[,j], family = "gaussian", lambda = rho/n, standardize = FALSE)
        a0 <- fit[["a0"]]
        old_beta <- fit[["beta"]]
        my_beta <- as(old_beta, "matrix")
        B[,j] <- my_beta
        mu[,j] <- a0
      }
      
      dum <- c(Y - X%*%B - ones%*%mu - L)
      fval <- 0.5 * t(dum) %*% dum + rho * sum(abs(B)) + lambda*sum(abs(diag(W)))
      
      res <- abs(fval-fval_old)/abs(fval_old+eps)
      
      print(paste('Iter ', iter, 'fval', fval, 'res', res))
      
      if (res < tol){
        break
      }
      
      fval_old <- fval
      f_val_vec <- c(f_val_vec, fval)
      res_vec <- c(res_vec, res)
    }
  }
  
  if(verbose == FALSE){
    for (iter in 1:maxiter){
      
      #### Compute L
      svd_st_res <- svd_st(Y - X%*%B - ones %*% mu,lambda)
      L <- svd_st_res$L
	  W <- svd_st_res$W
      
      for (j in 1:q){
        fit <- glmnet(X, Y[,j]-L[,j], family = "gaussian", lambda = rho/n, standardize = FALSE)
        a0 <- fit[["a0"]]
        old_beta <- fit[["beta"]]
        my_beta <- as(old_beta, "matrix")
        B[,j] <- my_beta
        mu[,j] <- a0
      }
      
      dum <- c(Y - X%*%B - ones%*%mu - L)
      fval <- 0.5 * t(dum) %*% dum + rho * sum(abs(B)) + lambda*sum(abs(diag(W)))
      
      res <- abs(fval-fval_old)/abs(fval_old+eps)
      
      if (res < tol){
        break
      }
      
      fval_old <- fval
      f_val_vec <- c(f_val_vec, fval)
      res_vec <- c(res_vec, res)
    }
  }
  return (list("B" = B, "L" = L, "mu" = mu, "f_val_vec" = f_val_vec, "res_vec" = res_vec, "iter" = length(f_val_vec)))
}


### Run_LORS_Screening

Run_LORS_Screening <- function(myB){
  
  index_list = c()
  for (i in 1:ncol(myB)){
    cands = abs(myB[,i])
    sorted_cands <- sort(cands, index.return=TRUE, decreasing=TRUE)$ix[1:nrow(X)]
    index_list <- c(index_list, sorted_cands)
  }
  selectedSNPs <- sort(unique(index_list))
  return(selectedSNPs)
}


###===========================================###
###============== norm functions =============###

### trace norm

tr_norm <- function(mat){

   res <- sum(svd(mat)[["d"]])
   return(res)  
}


### L1 norm

l1_norm <- function(mat){

   res <- sum(abs(mat))
   return(res)
}

### L2 norm

l2_norm <- function(mat){

   res <- norm(mat,"F")**2
   return(res)
}


###==================================================###
###================ proximal mapping ================###

### B, L1 penalty

prox_l1 <- function(X,B,mu,L,Y,rho,t_B){

    res <- prox_1(B - t_B * t(X) %*% (X %*% B + ones %*% mu + L - Y),rho * t_B)
	return(res)
}


### mu, no penalty

prox_null <- function(X,B,mu,L,Y,t_mu){

    res <- mu - t_mu * t(ones) %*% (X %*% B + ones %*% mu + L - Y)
	return(res)
}


### L, nuclear penalty

prox_tr <- function(X,B,mu,L,Y,lambda,t_L){

    res <- svd_st(L - t_L * (X %*% B + ones %*% mu + L - Y),lambda * t_L)$L
	return(res)
}


###=================================================###
###================= gradient functions ============###

### gradf_B

gradf_B <- function(X,B,mu,L,Y){

    res <- t(X) %*% (X %*% B + ones %*% mu + L - Y)
	return(res)
}


### gradf_L

gradf_L <- function(X,B,mu,L,Y){

    res <- X %*% B + ones %*% mu + L - Y
    return(res)
}


### gradf_mu

gradf_mu <- function(X,B,mu,L,Y){

    res <- t(ones) %*% (X %*% B + ones %*% mu + L - Y)
	return(res)
}


###==================================================###
###================== F、Q functions =================###

F <- function(X,B,L,mu,Y){

    res <- 0.5 * norm(Y - X %*% B - ones %*% mu - L, "F")**2 
	return(res)
}

## objective function
FF <- function(X,B,L,mu,Y,lambda,rho,M){

     res <- F(X,B,L,mu,Y) + lambda * tr_norm(L) + rho * Bh_norm(B,M)
	 return(res)
}

## reference Q(z,x)
Q_L <- function(arg1,arg2,X,B,mu,Y,t_L){

    res <- F(X,B,arg2,mu,Y) + sum(diag(t(arg1 - arg2) %*% gradf_L(X,B,mu,arg2,Y))) + (norm(arg1 - arg2, "F")**2) / (2 * t_L)
	return(res)
}

Q_B <- function(arg1,arg2,X,L,mu,Y,t_B){

    res <- F(X,arg2,L,mu,Y) + sum(diag(t(arg1 - arg2) %*% gradf_B(X,arg2,mu,L,Y))) + (norm(arg1 - arg2, "F")**2) / (2 * t_B)
    return(res)
}


Q_mu <- function(arg1,arg2,X,B,L,Y,t_mu){

     res <- F(X,B,L,arg2,Y) + sum(diag(t(arg1 - arg2) %*% gradf_mu(X,B,arg2,L,Y))) + (norm(arg1 - arg2, "F")**2) / (2 * t_mu)
     return(res)
}


###=================================================###
###============= used to update L、B、mu =============###

### Update_L

Update_L <- function(X,B,mu,L,Y,lambda,t_L,eta_L){
 
      i_k <- 1
	  arg1 <- prox_tr(X,B,mu,L,Y,lambda,t_L)
	  arg2 <- L
	  
	  F_L_val <- F(X,B,arg1,mu,Y)
	  Q_L_val <- Q_L(arg1,arg2,X,B,mu,Y,t_L)
	  
      while(F_L_val > Q_L_val){
	  
	       i_k <- i_k + 1
	       t_L <- eta_L^{i_k} * t_L
		   arg1 <- prox_tr(X,B,mu,L,Y,lambda,t_L)
	      
	       F_L_val <- F(X,B,arg1,mu,Y)    
	       Q_L_val <- Q_L(arg1,arg2,X,B,mu,Y,t_L)
	  
	  }
	  
	  return(arg1)
}


### Update_B

Update_B <- function(X,B,mu,L,Y,rho,M,t_B,eta_B){

      i_k <- 1
	  arg1 <- prox_Bh(X,B,mu,L,Y,rho,M,t_B)
	  arg2 <- B
	  
	  F_B_val <- F(X,arg1,L,mu,Y)
	  Q_B_val <- Q_B(arg1,arg2,X,L,mu,Y,t_B)
	  
      while(F_B_val > Q_B_val){
	  
	       i_k <- i_k + 1
	       t_B <- eta_B^{i_k} * t_B
		   arg1 <- prox_Bh(X,B,mu,L,Y,rho,M,t_B)
	       
	       F_B_val <- F(X,arg1,L,mu,Y)  
	       Q_B_val <- Q_B(arg1,arg2,X,L,mu,Y,t_B)
		   
	  }
	  
	  return(arg1)

}


### Update_mu

Update_mu <- function(X,B,mu,L,Y,t_mu,eta_mu){

      i_k <- 1
	  arg1 <- prox_null(X,B,mu,L,Y,t_mu)
	  arg2 <- mu
	  
	  F_mu_val <- F(X,B,L,arg1,Y)
	  Q_mu_val <- Q_mu(arg1,arg2,X,B,L,Y,t_mu)
	  
      while(F_mu_val > Q_mu_val){
	  
	       i_k <- i_k + 1
	       t_mu <- eta_mu^{i_k} * t_mu
		   arg1 <- prox_null(X,B,mu,L,Y,t_mu)
	       
	       F_mu_val <- F(X,B,L,arg1,Y)  
	       Q_mu_val <- Q_mu(arg1,arg2,X,B,L,Y,t_mu)
	  
	  }
	  
	  return(arg1)

}

###=================================================###

### mask generates \Omega1, \Omega2,...

mask <- function(nfold, Y){
    
	Omega <- !(is.na(Y))
	n <- nrow(Y)
	q <- ncol(Y)
    Omega.ls <- list()
	randn <- matrix(runif(n*q),nrow = n, ncol = q)
	
	for(i in 1:nfold){
	
	   Omega.ls[[i]] <- (randn <= i/nfold) & (randn > (i-1)/nfold)
	}
	
	return(Omega.ls)
}


### B, L1 penalty

prox_l1_Tuning <- function(X,B,mu,L,Y,rho,t_B,Training){

    res <- prox_1(B - t_B * t(X) %*% (Training * (X %*% B + ones %*% mu + L - Y)),rho * t_B)
	return(res)
}


### mu, no penalty

prox_null_Tuning <- function(X,B,mu,L,Y,t_mu,Training){

    res <- mu - t_mu * t(ones) %*% (Training * (X %*% B + ones %*% mu + L - Y))
	return(res)
}


### L, nuclear penalty

prox_tr_Tuning <- function(X,B,mu,L,Y,lambda,t_L,Training){

    res <- svd_st(L - t_L * (Training * (X %*% B + ones %*% mu + L - Y)),lambda * t_L)$L
	return(res)
}


###=============================================================###
###================= gradient functions with \Omega1 ===========###

### gradf_B

gradf_B_Tuning <- function(X,B,mu,L,Y,Training){

    res <- t(X) %*% (Training * (X %*% B + ones %*% mu + L - Y))
	return(res)
}


### gradf_L

gradf_L_Tuning <- function(X,B,mu,L,Y,Training){

    res <- Training * (X %*% B + ones %*% mu + L - Y)
    return(res)
}


### gradf_mu

gradf_mu_Tuning <- function(X,B,mu,L,Y,Training){

    res <- t(ones) %*% (Training * (X %*% B + ones %*% mu + L - Y))
	return(res)
}


###===============================================================###
###================== F、Q functions with \Omega1 =================###

F_Tuning <- function(X,B,L,mu,Y,Training){

    res <- 0.5 * norm(Training * (Y - X %*% B - ones %*% mu - L), "F")**2 
	return(res)
}

FF_Tuning <- function(X,B,L,mu,Y,lambda,rho,M,Training){

     res <- F_Tuning(X,B,L,mu,Y,Training) + lambda * tr_norm(L) + rho * Bh_norm(B,M)
	 return(res)
}



Q_L_Tuning <- function(arg1,arg2,X,B,mu,Y,t_L,Training){

    res <- F_Tuning(X,B,arg2,mu,Y,Training) + sum(diag(t(arg1 - arg2) %*% gradf_L_Tuning(X,B,mu,arg2,Y,Training))) + (norm(arg1 - arg2, "F")**2) / (2 * t_L)
	return(res)
}

Q_B_Tuning <- function(arg1,arg2,X,L,mu,Y,t_B,Training){

    res <- F_Tuning(X,arg2,L,mu,Y,Training) + sum(diag(t(arg1 - arg2) %*% gradf_B_Tuning(X,arg2,mu,L,Y,Training))) + (norm(arg1 - arg2, "F")**2) / (2 * t_B)
    return(res)
}


Q_mu_Tuning <- function(arg1,arg2,X,B,L,Y,t_mu,Training){

     res <- F_Tuning(X,B,L,arg2,Y,Training) + sum(diag(t(arg1 - arg2) %*% gradf_mu_Tuning(X,B,arg2,L,Y,Training))) + (norm(arg1 - arg2, "F")**2) / (2 * t_mu)
     return(res)
}


###=================================================###
###============= used to update L、B、mu =============###

### Update_L

Update_L_Tuning <- function(X,B,mu,L,Y,lambda,t_L,eta_L,Training){
 
      i_k <- 1
	  arg1 <- prox_tr_Tuning(X,B,mu,L,Y,lambda,t_L,Training)
	  arg2 <- L
	  
	  F_L_val <- F_Tuning(X,B,arg1,mu,Y,Training)
	  Q_L_val <- Q_L_Tuning(arg1,arg2,X,B,mu,Y,t_L,Training)
	  
      while(F_L_val > Q_L_val){
	  
	       i_k <- i_k + 1
	       t_L <- eta_L^{i_k} * t_L
		   arg1 <- prox_tr_Tuning(X,B,mu,L,Y,lambda,t_L,Training)
	      
	       F_L_val <- F_Tuning(X,B,arg1,mu,Y,Training)    
	       Q_L_val <- Q_L_Tuning(arg1,arg2,X,B,mu,Y,t_L,Training)
	  
	  }
	  
	  return(arg1)
}


### Update_B

Update_B_Tuning <- function(X,B,mu,L,Y,rho,M,t_B,eta_B,Training){
      
      i_k <- 1
	  arg1 <- prox_Bh_Tuning(X,B,mu,L,Y,rho,M,t_B,Training)
	  arg2 <- B
	  
	  F_B_val <- F_Tuning(X,arg1,L,mu,Y,Training)
	  
	  Q_B_val <- Q_B_Tuning(arg1,arg2,X,L,mu,Y,t_B,Training)
	  
      while(F_B_val > Q_B_val){
	  
	       i_k <- i_k + 1
	       t_B <- eta_B^{i_k} * t_B
		   arg1 <- prox_Bh_Tuning(X,B,mu,L,Y,rho,M,t_B,Training)
	       
	       F_B_val <- F_Tuning(X,arg1,L,mu,Y,Training)  
	       Q_B_val <- Q_B_Tuning(arg1,arg2,X,L,mu,Y,t_B,Training)
	 
	  }
	  
	  return(arg1)

}


### Update_mu

Update_mu_Tuning <- function(X,B,mu,L,Y,t_mu,eta_mu,Training){

      i_k <- 1
	  arg1 <- prox_null_Tuning(X,B,mu,L,Y,t_mu,Training)
	  arg2 <- mu
	  
	  F_mu_val <- F_Tuning(X,B,L,arg1,Y,Training)
	  Q_mu_val <- Q_mu_Tuning(arg1,arg2,X,B,L,Y,t_mu,Training)
	  
      while(F_mu_val > Q_mu_val){
	  
	       i_k <- i_k + 1
	       t_mu <- eta_mu^{i_k} * t_mu
		   arg1 <- prox_null_Tuning(X,B,mu,L,Y,t_mu,Training)
	       
	       F_mu_val <- F_Tuning(X,B,L,arg1,Y,Training)  
	       Q_mu_val <- Q_mu_Tuning(arg1,arg2,X,B,L,Y,t_mu,Training)
	  
	  }
	  
	  return(arg1)

}


###==================================================###

### rho_seq, generate a sequence of rho values

rho_seq <- function(X,Y,L,nrho,Omega,quickRho){

    ## Set a sequence of rho
	
    MaxRho <- GetMaxRho(X, Y, L, Omega)  
    
    rhoseq <- logspace(log10(MaxRho),log10(MaxRho*.05),nrho)
    
    if(quickRho == TRUE){
      rhoseq2 <- rhoseq[8:13]
      rhoseq <- rhoseq2
    }
    nrho <- length(rhoseq)
	
	return(list(rhoseq = rhoseq, nrho = nrho))
}


### ParamTune, modify ParamTune in R package FastLORS to output L and lambda

ParamTune <- function(Training, Validation, tune_method){
 
    maxRank = min(n,q)/2
    mySI <- softImpute(Y,NULL,Omega, Training, Validation, NULL,maxRank)
    Err <- mySI[["Err"]]
    Alpha <- mySI[["Alpha"]]
    bestind_lam <- min(Err)==Err;
    
    lambda <- mean(Alpha[bestind_lam])
    
    ## Get a good initialization of L
	
    L <- svd_st(Y,lambda)$L               
    
	return(list(lambda = lambda, L = L))
}


### rho_tune: to obtain the best \rho

rho_tune <- function(tune_method,Y,X,L,lambda,rhoseq,nrho,quickRho,Training,Validation,maxiter,Alg,stepsize){

    rhoErr <- matrix(0, nrow = nrho, ncol = 1)
    B <- matrix(0, nrow = ncol(X), ncol = q)
    mu <- matrix(0, nrow = 1, ncol = q)

    for(irho in 1:length(rhoseq)){
	
      print(paste("On rho",irho,"of", length(rhoseq)))
      rho <- rhoseq[irho]
	  
	  if(tune_method == "FastLORS"){
        out <- Fast_LORS_Tuning(Y = Y, X = X, rho = rho, lambda = lambda, Training = Training, Validation = Validation, tol = tol, maxiter = 5000, L = L, omega_SOR = omega_SOR)
      }
      if(tune_method == "LORS"){
        out <- LORS2(Y = Y, X = X, L = L, Omega1 = Training, Omega2 = Validation, rho = rho, lambda = lambda, tol = tol, maxIter = maxiter)
      }
	  
      Err <- out[["Err"]]
      rhoErr[irho,1] <- Err	  	  

    }
	
	bestind_rho <- rhoErr==min(rhoErr)
    rho <- mean(rhoseq[bestind_rho]) 
    return(rho)	
}	


### LORS_BestPara

LORS_BestPara <- function(Y, X, CV, nfold, nrho, quickRho, Omega, maxiter){

    if(CV == "FALSE"){
	
	   nfold <- 2
	   mask.ls <- mask(nfold, Y)
	   
	   params <- ParamTune(Training = Omega & mask.ls[[1]], Validation = Omega & mask.ls[[2]], tune_method)
       lambda.best <- params[["lambda"]]
       L <- params[["L"]]
	   
	   rho_info <- rho_seq(X,Y,L,nrho,Omega,quickRho)
       nrho <- rho_info$nrho
       rhoseq <- rho_info$rhoseq
	   
	   rho.best <- rho_tune(tune_method = "LORS", Y = Y, X = X, L = L, lambda = lambda.best, rhoseq = rhoseq, nrho = nrho, quickRho, Training = Omega & mask.ls[[1]], Validation = Omega & mask.ls[[2]], maxiter = maxiter,Alg,stepsize)
	   
	}else{


    lambda.best.vec <- rho.best.vec <- nrho.vec <- numeric(nfold)
	L.list <- rhoseq.list <- list(nfold)
	
    mask.ls <- mask(nfold, Y)
   
    for(cvn in 1:nfold){
   
        params <- ParamTune(Training = Omega & !mask.ls[[cvn]], Validation = Omega & mask.ls[[cvn]], tune_method)
        lambda.best.vec[cvn] <- params[["lambda"]]
		L.list[[cvn]] <- params[["L"]]
	   
	    rho_info <- rho_seq(X,Y,L.list[[cvn]],nrho,Omega,quickRho)
        nrho.vec[cvn] <- rho_info$nrho
        rhoseq.list[[cvn]] <- rho_info$rhoseq
	   
	    rho.best.vec[cvn] <- rho_tune("LORS", Y, X, L.list[[cvn]], lambda.best.vec[cvn], rhoseq.list[[cvn]], nrho.vec[cvn], quickRho, Omega & !mask.ls[[cvn]], Omega & mask.ls[[cvn]], maxiter,Alg,stepsize)
    }
 
    lambda.best <- mean(lambda.best.vec)
    rho.best <- mean(rho.best.vec)
	
	}
    return(list(lambda = lambda.best, rho = rho.best))
}


### FastLORS_BestPara

FastLORS_BestPara <- function(Y, X, CV, nfold, nrho, quickRho, Omega, maxiter){

    if(CV == "FALSE"){
	
	   nfold <- 2
	   mask.ls <- mask(nfold, Y)
	   
	   params <- ParamTune(Training = Omega & mask.ls[[1]], Validation = Omega & mask.ls[[2]], tune_method)
       lambda.best <- params[["lambda"]]
       L <- params[["L"]]
	   
	   rho_info <- rho_seq(X,Y,L,nrho,Omega,quickRho)
       nrho <- rho_info$nrho
       rhoseq <- rho_info$rhoseq
	   
	   rho.best <- rho_tune(tune_method = "FastLORS", Y = Y, X = X, L = L, lambda = lambda.best, rhoseq = rhoseq, nrho = nrho, quickRho, Training = Omega & mask.ls[[1]], Validation = Omega & mask.ls[[2]], maxiter = maxiter,Alg,stepsize)
	   
	}else{

 
    lambda.best.vec <- rho.best.vec <- nrho.vec <- numeric(nfold)
	L.list <- rhoseq.list <- list(nfold)
	
    mask.ls <- mask(nfold, Y)
   
    for(cvn in 1:nfold){
   
        params <- ParamTune(Training = Omega & !mask.ls[[cvn]], Validation = Omega & mask.ls[[cvn]], tune_method)
        lambda.best.vec[cvn] <- params[["lambda"]]
		L.list[[cvn]] <- params[["L"]]
	   
	    rho_info <- rho_seq(X,Y,L.list[[cvn]],nrho,Omega,quickRho)
        nrho.vec[cvn] <- rho_info$nrho
        rhoseq.list[[cvn]] <- rho_info$rhoseq
	   
	    rho.best.vec[cvn] <- rho_tune("FastLORS", Y, X, L.list[[cvn]], lambda.best.vec[cvn], rhoseq.list[[cvn]], nrho.vec[cvn], quickRho, Omega & !mask.ls[[cvn]], Omega & mask.ls[[cvn]], maxiter, Alg, stepsize)
    }
 
    lambda.best <- mean(lambda.best.vec)
    rho.best <- mean(rho.best.vec)
	}
	
    return(list(lambda = lambda.best, rho = rho.best))
}


### rankHC2 is copied from FastLORS

rankHC2 <- function(Bhat_standardized){
  HC_vec <- c()
  q <- ncol(Bhat_standardized)
  for (j in 1:nrow(Bhat_standardized)){
    if(length(which(Bhat_standardized[j,] == 0)) != ncol(Bhat_standardized)){  #### for very sparse matrices, some rows may be all 0.  Don't want to calculate HC for these
      t_vec <- sort(abs(Bhat_standardized[j,]), decreasing=FALSE)
      #t_vec <- t_vec[1:(length(t_vec) - 2)] ### values at the end of t_vec can be extreme
      t_vec <- t_vec[1:(length(t_vec) - 1)] ### values at the end of t_vec can be extreme
      
      HC = max(sqrt(q)*(S(t_vec, Bhat_standardized[j,])/q - 2 * survival(t_vec)) / sqrt(2 * survival(t_vec) * (1 - 2 * survival(t_vec))))
      HC_vec <- c(HC_vec, HC)
    }
    
    if(length(which(Bhat_standardized[j,] == 0)) == ncol(Bhat_standardized)){  #### for very sparse matrices, some rows may be all 0.  HC = 0 for these
      HC_vec <- c(HC_vec, 0)
    }
    
  }
  
  sorted_HC <- sort(HC_vec, index.return=TRUE, decreasing=TRUE)
  
  return(list("index" = sorted_HC$ix, "HC_vec" = HC_vec))
}

### S is copied from FastLORS

S <- function(t,my_matrix){
  length(my_matrix)*(1 - ecdf(abs(my_matrix))(t))
}

### survival is copied from FastLORS

survival <- function(t){
  s <- 1 - pnorm(t, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  return (s)
}

