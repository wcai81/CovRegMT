#include <RcppArmadillo.h>
// [[Rcpp::plugins(cpp11)]]
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

   


//[[Rcpp::export]]
double ranklex_Cpp(vec X, vec U){

    int n = X.n_rows;
    vec r(n, fill::zeros);
    
    for (int i =0; i <= (n-2); i++){

        if ( X(n-1) > X(i) ) r(i) = 1;
        
        if ( ( X(n-1) == X(i) ) && ( U(n-1) > U(i) ) ) r(i) = 1;
    }
    
    return(sum(r)+1);
}




//[[Rcpp::export]]
vec vechs_Cpp(mat cormat){
    
    int Nsize = cormat.n_rows, k;
    vec corvechs( Nsize*(Nsize-1)/2 );
    
    k = 0;
    for (int j = 1; j <= Nsize-1; j++){
        for (int i = j+1; i <= Nsize; i++){
            
            corvechs(k) = cormat(i-1,j-1);
            k = k + 1;
        }
    }
    
    return(corvechs);
}




//[[Rcpp::export]]
mat unvechs_Cpp(vec corvechs, int Nsize){
        
    mat cormat(Nsize, Nsize, fill::ones);
    
    int k = 0;
    for (int j = 1; j <= Nsize-1; j++){
        for (int i = j+1; i <= Nsize; i++){
            
            cormat(i-1,j-1) = corvechs(k);
            cormat(j-1,i-1) = cormat(i-1,j-1);
            k = k + 1;
        }
    }
    
    return(cormat);
}




//[[Rcpp::export]]
double kmax_Cpp(vec X, int k){

    vec sorted_X = sort(X, "descend");
                        
    return(sorted_X(k-1));
}




//[[Rcpp::export]]
mat SSpvalues_Cpp(mat &YY, int k, int totsim){
    
    int Tsize = YY.n_rows, Nsize = YY.n_cols;
    
    int totcors = Nsize*(Nsize-1)/2;
    
    vec abscors(totcors), pvaluesSS(totcors, fill::zeros);
        
    mat pvaluesmat(Nsize, Nsize, fill::zeros);
                       
    vec ss(Nsize);
    
    mat D(Nsize, Nsize, fill::zeros), cormat(Nsize, Nsize, fill::zeros), cormatDATA(Nsize, Nsize, fill::zeros);
    
    mat statsSS(totsim, totcors, fill::zeros);
    
    mat YYtil(Tsize, Nsize, fill::zeros), Stil(Tsize, Nsize, fill::zeros);
    
    vec U(totsim, fill::randu);
    
    double kmaxabscors;
                    
    for (int j = 0; j <= Nsize-1; j++){
        ss(j) = sum( YY(span(0, Tsize-1), j) % YY(span(0, Tsize-1), j) );
        D(j,j) = 1/sqrt(ss(j));
    }
    
    cormat = D * trans(YY) * YY * D;
    cormatDATA = cormat;
        
    abscors = vechs_Cpp(abs(cormat));
            
    statsSS(totsim-1, span(0, totcors-1)) = trans(abscors);
               
    for (int isim = 0; isim <= totsim-2; isim++){
                
        Stil = sign(randn(Tsize, Nsize));
        YYtil = Stil % YY;

        cormat = D * trans(YYtil) * YYtil * D;

        abscors = vechs_Cpp(abs(cormat));
            
        kmaxabscors = kmax_Cpp(abscors, k);
                                
        for (int icol=0; icol<= totcors-1; icol++) statsSS(isim, icol)= kmaxabscors;
  
    }
        
    for (int j = 0; j <= totcors-1; j++) pvaluesSS(j) = (totsim - ranklex_Cpp( statsSS( span(0,totsim-1),j ), U ) + 1)/totsim;
        
    pvaluesmat = unvechs_Cpp(pvaluesSS(span(0, totcors-1), 0), Nsize);
    
    pvaluesmat.diag() = zeros(Nsize);
        
    return(pvaluesmat);
}




//[[Rcpp::export]]
mat SDpvalues_Cpp(mat &YY, int k, int totsim){
    
    int Tsize = YY.n_rows, Nsize = YY.n_cols;
    
    int totcors = Nsize*(Nsize-1)/2;
    
    vec abscors(totcors), pvaluesSD(totcors), mmaxs(totcors), orderedpvalues(totcors);
    
    vec umaxs(totcors);
    
    mat pvaluesmat(Nsize, Nsize, fill::zeros);
        
    vec temp(2);
       
    uvec index(totcors);
        
    vec ss(Nsize);
    
    mat D(Nsize, Nsize, fill::zeros), cormat(Nsize, Nsize, fill::zeros), cormatDATA(Nsize, Nsize, fill::zeros);
    
    mat statsSD(totsim, totcors, fill::zeros);
    
    mat YYtil(Tsize, Nsize, fill::zeros), Stil(Tsize, Nsize, fill::zeros);
    
    vec U(totsim, fill::randu);
    
    double kmaxabscors;
                    
    for (int j = 0; j <= Nsize-1; j++){
        ss(j) = sum( YY(span(0, Tsize-1), j) % YY(span(0, Tsize-1), j) );
        D(j,j) = 1/sqrt(ss(j));
    }
    
    cormat = D * trans(YY) * YY * D;
    cormatDATA = cormat;
        
    abscors = vechs_Cpp(abs(cormat));
    
    vec sorted_abscors = sort(abscors, "descend");
    
    index = sort_index(abscors, "descend");
                
    statsSD(totsim-1, span(0, totcors-1)) = trans(sorted_abscors);
               
    for (int isim = 0; isim <= totsim-2; isim++){
                
        Stil = sign(randn(Tsize, Nsize));
        YYtil = Stil % YY;

        cormat = D * trans(YYtil) * YYtil * D;

        abscors = vechs_Cpp(abs(cormat));
            
        kmaxabscors = kmax_Cpp(abscors, k);
                                        
        umaxs(totcors-1) = abscors(index(totcors-1));
        
        for (int i=(totcors-2); i>= 1; i--){
                        
            temp(0)  = umaxs(i+1);
            temp(1)  = abscors(index(i));
            umaxs(i) = max(temp);
        }

        kmaxabscors = kmax_Cpp(abscors, k);
        mmaxs(0) = kmaxabscors;
        
        for (int i=1; i<= (totcors-1); i++){
                                    
            temp(0)  = umaxs(i);
            temp(1)  = mmaxs(i-1);
            mmaxs(i) = min(temp);
        }
                
        statsSD(isim, span(0, totcors-1)) = trans(mmaxs);
    }
        
    for (int j = 0; j <= totcors-1; j++) orderedpvalues(j) = (totsim - ranklex_Cpp( statsSD( span(0,totsim-1),j ), U ) + 1)/totsim;
    
    for (int j = 1; j <= totcors-1; j++){
        
        temp(0) = orderedpvalues(j-1);
        temp(1) = orderedpvalues(j);
        orderedpvalues(j) = max(temp);
    }
        
    for (int j = 0; j <= totcors-1; j++) pvaluesSD(index(j)) = orderedpvalues(j);
        
    pvaluesmat = unvechs_Cpp(pvaluesSD(span(0, totcors-1), 0), Nsize);
    
    pvaluesmat.diag() = zeros(Nsize);
    
    return(pvaluesmat);
}




//[[Rcpp::export]]
int CheckFDP_Cpp(mat &pvalues, double alpha){
    
    int Nsize = pvalues.n_cols, NumRej = 0;
                
    for (int i=(1-1); i<=(Nsize-1-1); i++ ){
        
        for (int j=(i+1); j<= (Nsize-1); j++){
            
            if (pvalues(i,j) <= alpha) NumRej = NumRej + 1;
        }
    }
    
    return(NumRej);
}




//[[Rcpp::export]]
mat Shrink_Cpp(mat &rho, mat &rhotil, int Tsize){
    
    int Nsize = rho.n_cols;
        
    mat RtilSMT(Nsize, Nsize, fill::zeros), A(Nsize, Nsize, fill::zeros);;
            
    vec eigval;
    mat eigvec;
    
    double lambdamin, tol, inc, xi0, xi, xistar;
    
    eig_sym(eigval, eigvec, rhotil);
    
    RtilSMT = rhotil;
    
    lambdamin = min(eigval);
    
    tol = 0.01; inc= tol/2;
    
    xi0 = 0;
    if (lambdamin <= tol) xi0 = (tol-lambdamin)/(1-lambdamin);
    
    double num = 0;
    double den1 = 0;
    double den2 = 0;
    double temp1 = 0;
    double temp2 = 0;
    
    for (int i =0; i <= Nsize-1; i++){
        for (int j = 0; j <= Nsize-1; j++ ){
            
            if (i != j){

                temp1 = 1 - pow(rho(i,j), 2);

                temp2 = rho(i,j) - rho(i,j) * temp1/(2*Tsize);

                num = num + rho(i,j) * temp2;
                                
                den1 = den1 + pow(temp1, 2);
                den2 = den2 + pow(temp2, 2);
            }
        }
    }
    
    den1 = den1/Tsize;
    double thetastar = 1 - num/(den1 + den2);
    if (thetastar < 0) thetastar = 0;
    if (thetastar > 1) thetastar = 1;
            
    mat I = eye(Nsize, Nsize);
    mat R0 = thetastar * I + (1-thetastar) * rho;
    mat R0inv = inv(R0);
    
    int numxi = 1 + (1-xi0)/inc;
    mat Results(numxi, 2, fill::zeros);
        
    for (int i=0; i <= numxi-1; i++){
        
        xi = xi0 + i*inc;
        
        Results(i,0) = xi;
        
        RtilSMT = xi * I + (1-xi) * rhotil;
        
        eig_sym(eigval, eigvec, RtilSMT);
        
        lambdamin = min(eigval);
        
        A = R0inv - inv(RtilSMT);
                        
        Results(i,1) = sum(diagvec(trans(A)*A));
    }
    
    uword idx = Results(span(0, numxi-1),1).index_min();
            
    xistar = Results(idx, 0);
    
    RtilSMT = xistar * I + (1-xistar) * rhotil;
                            
    return(RtilSMT);
}

