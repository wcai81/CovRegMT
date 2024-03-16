#####################################################################################################
library(Rcpp)
library(RcppArmadillo)

sourceCpp("Functions.cpp")




is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol




reg_kFWER <- function(data, center=TRUE, k=1, step="SS", B=100, alpha=0.05){
	
	Tsize <- dim(data)[1]
	Nsize <- dim(data)[2]
	M <- (Nsize^2 - Nsize)/2
		
	if (center==TRUE){
		### Center returns by subtracting means 

		mu <- apply(data, 2, mean)
		Y <- data - matrix(mu, Tsize, Nsize, byrow=TRUE) 		
	} else{
		Y <- data
	}		

	# Sample covariance and correlation matrices

	S <- matrix(0, Nsize, Nsize)
	ss <- apply(Y^2, 2, sum)		
	diag(S) <- sqrt(ss)	
	cormat <- solve(S) %*%  t(Y) %*% Y  %*% solve(S)
	STDs <- S/sqrt(Tsize)
	covmat <- STDs %*% cormat %*% STDs				
			
	if (is.wholenumber(alpha*B)==FALSE){		

		print("Warning: The value of B should be chosen so that alpha*B is an integer.")
		cat("\n")			

	} else{
		
		if (step=="SS") pvalues <- SSpvalues_Cpp(Y=data, k, B)

		if (step=="SD") pvalues <- SDpvalues_Cpp(Y=data, k, B)
		
		idx <- 1*( pvalues <= alpha)
		cormatShrink <- Shrink_Cpp(cormat, cormat*idx, Tsize)
		regcovmat <- STDs %*% cormatShrink %*% STDs				

		return(list("covmat"=covmat, "cormat"=cormat, "pvalues"=pvalues, "regcovmat"= regcovmat))														
	}
}




reg_FDP <- function(data, center=TRUE, gamma=0.1, step="SS", B=100, alpha=0.05){
	
	Tsize <- dim(data)[1]
	Nsize <- dim(data)[2]
	M <- (Nsize^2 - Nsize)/2
						
	if (Nsize >= 500){
		print("Warning: Previous testing with N=500 took approximately 6 minutes to execute. Expect longer times for larger N values.")
		cat("\n")
	}	
		
	if (center==TRUE){
		### Center returns by subtracting means 

		mu <- apply(data, 2, mean)
		Y <- data - matrix(mu, Tsize, Nsize, byrow=TRUE) 		
	} else{
		Y <- data
	}		
		
	# Sample covariance and correlation matrices

	S <- matrix(0, Nsize, Nsize)
	ss <- apply(Y^2, 2, sum)		
	diag(S) <- sqrt(ss)	
	cormat <- solve(S) %*%  t(Y) %*% Y  %*% solve(S)
	STDs <- S/sqrt(Tsize)
	covmat <- STDs %*% cormat %*% STDs				
		
	if (is.wholenumber(alpha*B)==FALSE){
		
		print("Warning: The value of B should be chosen so that alpha*B is an integer.")
		cat("\n")		
	
	} else{

		kl <- 1
		kr <- M

		while ((kr-kl) > 1){	

			km <- floor( (kl + kr)/2 )
			k <- km
			set.seed(123456)

			if (step=="SS") pvaluesFDP <- SSpvalues_Cpp(Y=data, k, B)
			if (step=="SD") pvaluesFDP <- SDpvalues_Cpp(Y=data, k, B)
						
			NumRejsm <- CheckFDP_Cpp(pvaluesFDP, alpha)
	
			if (km <= gamma*(NumRejsm + 1)){ 
				kl <- km				
			}else{
				kr <- km
			}	
		}
		
		k <- kl
		set.seed(123456)
		
		if (step=="SS") pvaluesFDP <- SSpvalues_Cpp(Y=data, k, B)
		if (step=="SD") pvaluesFDP <- SDpvalues_Cpp(Y=data, k, B)
				
		NumRejs <- CheckFDP_Cpp(pvaluesFDP, alpha)

		if (k > gamma*(NumRejs + 1) ){   
		
			if (gamma == 0){	

				pvaluesFDP.final <- pvaluesFDP			

			}else if (gamma > 0){		

				pvaluesFDP.final <- "FDP-adjusted p-values cannot be produced."		
			}		
		} 

		while (k <= gamma*(NumRejs + 1) ){
								
			pvaluesFDP.final <- pvaluesFDP
			k <- k + 1			
			set.seed(123456)

			if (step=="SS") pvaluesFDP <- SSpvalues_Cpp(Y=data, k, B)
			if (step=="SD") pvaluesFDP <- SDpvalues_Cpp(Y=data, k, B)

			NumRejs <- CheckFDP_Cpp(pvaluesFDP, alpha)							
		}

		pvalues <- pvaluesFDP.final
		
		idx <- 1*( pvalues <= alpha)
		cormatShrink <- Shrink_Cpp(cormat, cormat*idx, Tsize)
		regcovmat <- STDs %*% cormatShrink %*% STDs				
		
		return(list("covmat"=covmat, "cormat"=cormat, "pvalues"=pvalues, "regcovmat"= regcovmat))							
	}	
}







