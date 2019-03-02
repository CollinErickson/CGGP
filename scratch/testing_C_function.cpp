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
void rcpp_fastmatclcr(IntegerMatrix I, IntegerMatrix w, NumericMatrix MSEmat, NumericVector S){
  
  int ns = S.size();
  int nb = w.size();
  int d = I.size();
  
  NumericVector RR(ns);
  
  for(int bindex=0; bindex<=(nb-1); bindex++){
     RR = (-w[bindex])*MSEmat( _ ,I(bindex,0));
     for(int dim=1;dim<=(d-1); dim++){
       RR = RR*MSEmat( _ ,I(bindex,dim));
     }
     
     S += RR;
  }
  void;
}
