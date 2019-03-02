#include <Rcpp.h>
using namespace Rcpp;

//' rcpp_kronDBS
//' 
//' @param A Vector
//' @param dA Vector
//' @param B Vector
//' @param p Vector
//' @return kronDBS calculation
//' @export
//' @examples
//' rcpp_gkronDBS(c(1,1), c(0,0), c(.75), c(1,1))
// [[Rcpp::export]]
void rcpp_fastmatclcr(NumericMatrix I, NumericVector w, NumericMatrix MSEmat, NumericVector S, int maxlevel){
  
  int ns = S.size();
  int nb = w.size();
  int d = I.ncol();
  
  NumericVector RR(ns);
  
  for(int bindex=0; bindex<=(nb-1); bindex++){
    if(abs(w[bindex])>0.5){
    RR = MSEmat( _ ,I(bindex,0)-1);
    for(int dim=1;dim<=(d-1); dim++){
      RR = RR*MSEmat( _ ,maxlevel*dim+I(bindex,dim)-1);
    }
    RR = -w[bindex]*RR;
     
     S += RR;
    }
  }
  return;
}


#include <Rcpp.h>
using namespace Rcpp;

//' rcpp_kronDBS
//' 
//' @param A Vector
//' @param dA Vector
//' @param B Vector
//' @param p Vector
//' @return kronDBS calculation
//' @export
//' @examples
//' rcpp_gkronDBS(c(1,1), c(0,0), c(.75), c(1,1))
// [[Rcpp::export]]
void rcpp_fastmatclcranddclcr(NumericMatrix I, NumericVector w, NumericMatrix MSEmat, NumericMatrix dMSEmat, NumericVector S, NumericMatrix dS, int maxlevel, int numpara){
  
  int ns = S.size();
  int nb = w.size();
  int d = I.ncol();
  
  NumericVector RR(ns);
  
  for(int bindex=0; bindex<=(nb-1); bindex++){
    if(abs(w[bindex])>0.5){
    RR = MSEmat( _ ,I(bindex,0)-1);
    for(int dim=1;dim<=(d-1); dim++){
      RR = RR*MSEmat( _ ,maxlevel*dim+I(bindex,dim)-1);
    }
    RR = -w[bindex]*RR;
    
    S += RR;
    
    for(int dim=0;dim<=(d-1); dim++){
      for(int para=0;para<=(numpara-1); para++){
        dS( _ , numpara*dim+para) = dS( _ , numpara*dim+para)+RR*dMSEmat( _ ,maxlevel*numpara*dim+para*maxlevel+I(bindex,dim)-1);
      }
    }
    }
  }
  return;
}
