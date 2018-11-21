#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//' rcpp_kronDBS
//' 
//' @param A Vector
//' @param B Vector
//' @param p Vector
//' @param Al Length of A
//' @param Bl Length of B
//' @param d Dimension
//' @return kronDBS calculation
//' @export
// [[Rcpp::export]]
NumericVector rcpp_kronDBS(NumericVector A, NumericVector B, NumericVector p, int Al, int Bl,  int d){
  
  int sv = 0;
  
  NumericVector x(Bl);
  NumericVector y(Bl);
  int p0 = 0;
  int n0 = 0;
  
  sv = Al;
  for(int dim = d-1; dim>=0;dim--){
    sv = sv-p(dim)*p(dim);
    p0 = p(dim);
    n0 = Bl/p0;
    // if(p0 > 1.5){
    for(int h = 0; h < n0; h++)
    {
      for ( int i = 0; i < p0; i++)
      {
        x(h*p0+i) = B(h*p0+i);
        for ( int k = i-1; k >=0; k-- ) x(h*p0+i) -= A(sv+i*p0+k) * x(h*p0+k);
        x(h*p0+i) /= A(sv+i*p0+i);
      }
      
      for ( int i = p0-1; i >=0; i--)
      {
        y(h*p0+i) = x(h*p0+i);
        for ( int k=i+1; k < p0; k++) y(h*p0+i) -= A(sv+k*p0+i) *y(h*p0+k);
        y(h*p0+i) /= A(sv+i*p0+i);
      }
    }
    
    int c = 0;
    for ( int i = 0; i < p0; i++)
    {
      for(int h = 0; h < n0; h++)
      {
        B(c) = y(h*p0+i);
        c++;
      }
    }
    // }else{
    //   B = B/(A(sv)*A(sv));
    // }
  }
  
  return B;
}


// [[Rcpp::export]]
NumericMatrix rcpp_gkronDBS(NumericVector A,NumericVector dA, NumericVector B, NumericVector p, int Bl,  int d){
  
  int sv = 0;
  int sv2 = 0;
  
  NumericVector x(Bl);
  NumericVector y(Bl);
  NumericVector B2(Bl);
  NumericVector x2(Bl);
  NumericVector y2(Bl);
  
  
  NumericMatrix gz(Bl,d);
  int c = 0;
  int p0 = 0;
  int n0 = 0;
  
  for(int dim = d-1; dim>=0;dim--) sv = sv+p(dim)*p(dim);
  
  for(int dim = d-1; dim>=0;dim--){
    sv = sv-p(dim)*p(dim);
    p0 = p(dim);
    n0 = Bl/p0;
    
    // if(p0 > 1.5){
    for(int h = 0; h < n0; h++)
    {
      for ( int i = 0; i < p0; i++)
      {
        x(h*p0+i) = B(h*p0+i);
        for ( int k = i-1; k >=0; k-- ) x(h*p0+i) -= A(sv+i*p0+k) * x(h*p0+k);
        x(h*p0+i) /= A(sv+i*p0+i);
      }
      for ( int i = p0-1; i >=0; i--)
      {
        y(h*p0+i) = x(h*p0+i);
        for ( int k=i+1; k < p0; k++) y(h*p0+i) -= A(sv+k*p0+i)*y(h*p0+k);
        y(h*p0+i) /= A(sv+i*p0+i);
      }
    }
    
    c = 0;
    for ( int i = 0; i < p0; i++)
    {
      for(int h = 0; h < n0; h++)
      {
        B(c) = y(h*p0+i);
        c++;
      }
    }
    
    for(int h = 0; h < n0; h++)
    {
      for ( int i = 0; i < p0; i++)
      {
        B2(h*p0+i) = 0;
        for ( int k = 0; k <p0; k++) B2(h*p0+i) += dA(sv+k*p0+i)*y(h*p0+k);
      }
    }
    
    for(int h = 0; h < n0; h++)
    {
      for ( int i = 0; i < p0; i++)
      {
        x(h*p0+i) = B2(h*p0+i);
        for ( int k = i-1; k >=0; k-- ) x(h*p0+i) -= A(sv+i*p0+k) * x(h*p0+k);
        x(h*p0+i) /= A(sv+i*p0+i);
      }
      for ( int i = p0-1; i >=0; i--)
      {
        y(h*p0+i) = x(h*p0+i);
        for ( int k=i+1; k < p0; k++) y(h*p0+i) -= A(sv+k*p0+i)*y(h*p0+k);
        y(h*p0+i) /= A(sv+i*p0+i);
      }
    }
    
    c = 0;
    for ( int i = 0; i < p0; i++)
    {
      for(int h = 0; h < n0; h++)
      {
        B2(c) = y(h*p0+i);
        c++;
      }
    }
    // }else{
    //   B = B/(A(sv)*A(sv));
    //   B2 = dA(sv)*y/(A(sv)*A(sv)*A(sv)*A(sv));
    // }
    
    
    //if(dA(sv)>0.0000000001){
    sv2 = sv;
    for(int dim2 = dim-1; dim2>=0;dim2--){
      sv2 = sv2-p(dim2)*p(dim2);
      p0 = p(dim2);
      n0 = Bl/p0;
      
      // if(p0 > 1.5){
      for(int h = 0; h < n0; h++)
      {
        for ( int i = 0; i < p0; i++)
        {
          x2(h*p0+i) = B2(h*p0+i);
          for ( int k = i-1; k >=0; k-- ) x2(h*p0+i) -= A(sv2+i*p0+k) * x2(h*p0+k);
          x2(h*p0+i) /= A(sv2+i*p0+i);
        }
        
        for ( int i = p0-1; i >=0; i--)
        {
          y2(h*p0+i) = x2(h*p0+i);
          for (int k=i+1; k < p0; k++) y2(h*p0+i) -= A(sv2+k*p0+i) *y2(h*p0+k);
          y2(h*p0+i) /= A(sv2+i*p0+i);
        }
      }
      c = 0;
      for ( int i = 0; i < p0; i++)
      {
        for(int h = 0; h < n0; h++)
        {
          B2(c) = y2(h*p0+i);
          c++;
        }
      }
      // }else{
      //   B2 = B2/(A(sv)*A(sv));
      // }
    }
    //}
    
    
    for(int h2 = 0; h2 < Bl; h2++){
      gz(h2,dim) = B2(h2); 
    }
  }
  
  return gz;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
