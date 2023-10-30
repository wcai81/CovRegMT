#####################################################################################################
library(Rcpp)
library(RcppArmadillo)

sourceCpp("Functions.cpp")



#####################################################################################################
#### Illustrative example with simulated returns data

Rets <- read.table("Returns.txt")
Rets <- as.matrix(Rets)
	
Tsize <- dim(Rets)[1]
Nsize <- dim(Rets)[2]
M <- (Nsize^2 - Nsize)/2

	
### Need to center returns by subtracting means 

mu <- apply(Rets, 2, mean)
Y <- Rets - matrix(mu, Tsize, Nsize, byrow=TRUE) 


# Sample covariance and correlation matrices

S <- matrix(0, Nsize, Nsize)
ss <- apply(Y^2, 2, sum)		
diag(S) <- sqrt(ss)	
cormat <- solve(S) %*%  t(Y) %*% Y  %*% solve(S)
STDs <- S/sqrt(Tsize)

covmat <- STDs %*% cormat %*% STDs				

options(digits = 3)

print(covmat)
print(cormat)



#####################################################################################################
#### Procedures for k-FWER control 


# k is the value that defines the k-FWER error rate measure (1 <= k <= M)
k <- 1

# alpha is the desired k-FWER (0 < alpha < 1)	
alpha <- 0.05	
	
# totsim should be chosen so that alpha*totsim is an integer (preferably, totsim >= 100)	
totsim <- 100	 
	
	
### The next function call returns the single-step (SS) k-FWER-adjusted Monte Carlo p-values

set.seed(123456)	
pvaluesSS <- SSpvalues_Cpp(Y, k, totsim)
print(pvaluesSS)


## Shrinkage step 

idx <- 1*( pvaluesSS <= alpha)
cormatShrink <- Shrink_Cpp(cormat, cormat*idx, Tsize)
covmatReg.SS <- STDs %*% cormatShrink %*% STDs				

# Regularized covariance matrix
print(covmatReg.SS)


### The next function call returns the step-down (SD) k-FWER-adjusted Monte Carlo p-values

set.seed(123456)	
pvaluesSD <- SDpvalues_Cpp(Y, k, totsim)
print(pvaluesSD)


## Shrinkage step 

idx <- 1*( pvaluesSD <= alpha)
cormatShrink <- Shrink_Cpp(cormat, cormat*idx, Tsize)
covmatReg.SD <- STDs %*% cormatShrink %*% STDs				

# Regularized covariance matrix
print(covmatReg.SD)



#####################################################################################################
#### Procedures for FDP control 


# gamma is the lower bound FDP value (0 <= gamma < 1)
gamma <- 0.10

# alpha is the desired significance level (0 < alpha < 1)	
alpha <- 0.05
	
# totsim should be chosen so that alpha*totsim is an integer (preferably, totsim >= 100)	
totsim <- 100	 


### The following code chunk computes the FDP-adjusted Monte Carlo p-values using the SS procedure

kl <- 1
kr <- M

while ((kr-kl) > 1){	

	km <- floor( (kl + kr)/2 )
	k <- km
	set.seed(123456)
	pvaluesFDP.SS <- SSpvalues_Cpp(Y, k, totsim)
	NumRejsSSm <- CheckFDP_Cpp(pvaluesFDP.SS, alpha)
	
	if (km <= gamma*(NumRejsSSm + 1)){ 
		kl <- km				
	}else{
		kr <- km
	}	
}

k <- kl
set.seed(123456)
pvaluesFDP.SS <- SSpvalues_Cpp(Y, k, totsim)
NumRejsSS <- CheckFDP_Cpp(pvaluesFDP.SS, alpha)

if (k > gamma*(NumRejsSS + 1) ){   
		
	if (gamma == 0){	

		pvaluesFDP.SS.final <- pvaluesFDP.SS			

	}else if (gamma > 0){		

		pvaluesFDP.SS.final <- "FDP-adjusted p-values cannot be produced"		
	}	
	
} 

while (k <= gamma*(NumRejsSS + 1) ){
		
	pvaluesFDP.SS.final <- pvaluesFDP.SS
	k <- k + 1			
	set.seed(123456)
	pvaluesFDP.SS <- SSpvalues_Cpp(Y, k, totsim)			
	NumRejsSS <- CheckFDP_Cpp(pvaluesFDP.SS, alpha)	
}


# FDP-adjusted p-values using SS procedure
print(pvaluesFDP.SS.final)


## Shrinkage step 

idx <- 1*( pvaluesFDP.SS.final <= alpha)
cormatShrink <- Shrink_Cpp(cormat, cormat*idx, Tsize)
covmatReg.FDP.SS <- STDs %*% cormatShrink %*% STDs				

# Regularized covariance matrix
print(covmatReg.FDP.SS)




#####
### The next code chunk is the same as above except that now it uses the SD procedure

kl <- 1
kr <- M

while ((kr-kl) > 1){	

	km <- floor( (kl + kr)/2 )
	k <- km
	set.seed(123456)
	pvaluesFDP.SD <- SDpvalues_Cpp(Y, k, totsim)
	NumRejsSDm <- CheckFDP_Cpp(pvaluesFDP.SD, alpha)
	
	if (km <= gamma*(NumRejsSDm + 1)){ 
		kl <- km				
	}else{
		kr <- km
	}	
}

k <- kl
set.seed(123456)
pvaluesFDP.SD <- SDpvalues_Cpp(Y, k, totsim)
NumRejsSD <- CheckFDP_Cpp(pvaluesFDP.SD, alpha)

if (k > gamma*(NumRejsSD + 1) ){   
		
	if (gamma == 0){	

		pvaluesFDP.SD.final <- pvaluesFDP.SD			

	}else if (gamma > 0){		

		pvaluesFDP.SD.final <- "FDP-adjusted p-values cannot be produced"		
	}	
	
} 

while (k <= gamma*(NumRejsSD + 1) ){
		
	pvaluesFDP.SD.final <- pvaluesFDP.SD
	k <- k + 1			
	set.seed(123456)
	pvaluesFDP.SD <- SDpvalues_Cpp(Y, k, totsim)			
	NumRejsSD <- CheckFDP_Cpp(pvaluesFDP.SD, alpha)	
}


# FDP-adjusted p-values using SD procedure
print(pvaluesFDP.SD.final)


## Shrinkage step 

idx <- 1*( pvaluesFDP.SD.final <= alpha)
cormatShrink <- Shrink_Cpp(cormat, cormat*idx, Tsize)
covmatReg.FDP.SD <- STDs %*% cormatShrink %*% STDs				

# Regularized covariance matrix
print(covmatReg.FDP.SD)




	
