#include <RcppArmadillo.h>
// [[Rcpp::plugins(cpp11)]]
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

   


//[[Rcpp::export]]
double ranklex_Cpp(vec X, vec U){

    int n = X.n_rows;
    vec r(n, fill::zeros);
    
    for (int i = 0; i <= n-2; i++){

        if ( X(n-1) > X(i) ) r(i) = 1;
        
        if ( ( X(n-1) == X(i) ) && ( U(n-1) > U(i) ) ) r(i) = 1;
    }
    
    return(sum(r)+1);
}


