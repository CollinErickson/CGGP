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
//' @return kronDBS calculation
//' @export
// [[Rcpp::export]]
void rcpp_kronDBS(NumericVector A, NumericVector B, NumericVector p){
  int sv=A.size(); int Bl=B.size(); int d=p.size();
  NumericVector x(Bl); NumericVector y(Bl);
  int p0; int sv0; int c; int i; int k; int h;
  
  for(int dim=d-1; dim>=0;dim--){ //loop over the all demensions
    if(p[dim]>1.5){ //do not need to do much when kron with respect to one thing
      p0 = p[dim]; //look at our one demension
      sv = sv-p0*p0;
      for(h=0; h<Bl; h+=p0)  //loop over the leftover demensions
      {
        x[h]=B[h]/A[sv];  //do first outside loop
        for (i=1; i<p0; i++)  //backsolve with respect to dim+1 the first time
        {
          x[h+i]=B[h+i];
          sv0 = sv+i*p0;   //speed up index reference
          for(k=i-1; k>=0; k--) x[h+i]-=A[sv0+k]*x[h+k];
          x[h+i]/=A[sv0+i];
        }
        
        y[h+p0-1]=x[h+p0-1]/A[sv+p0*p0-1];   //do first outside loop
        for (i=p0-2; i>=0; i--) //backsolve with respect to dim+1 the second time
        {
          y[h+i]=x[h+i];
          sv0=sv+i*p0;   //speed up index reference
          for(k=i+1; k<p0; k++) y[h+i]-=A[sv0+k]*y[h+k];
          y[h+i]/=A[sv0+i];
        }
      }
      c=0;
      for(i=0; i<p0; i++) for(h=0; h<Bl; h+=p0){B[c]=y[h+i]; c++;} //spinning the vector to right orientation
    }else{sv--; B=B/(A[sv]*A[sv]);} // kron with respect to one thing shortcut
  }
  
  return;}


//' rcpp_kronDBS
//' 
//' @param A Vector
//' @param B Vector
//' @param p Vector
//' @return kronDBS calculation
//' @export
//' @examples
//' rcpp_gkronDBS(c(1,1), c(0,0), c(.75), c(1,1))
// [[Rcpp::export]]
NumericMatrix rcpp_gkronDBS(NumericVector A,NumericVector dA, NumericVector B, NumericVector p){
  int sv=A.size(); int Bl=B.size(); int d=p.size();
  NumericVector x(Bl); NumericVector y(Bl); NumericVector B1(Bl); NumericVector B2(Bl); 
  NumericVector dBn(Bl); NumericVector dBn2(Bl); NumericMatrix dB(d,Bl);
  int p0; int sv0; int c; int i; int k; int h; int dim2; int p1;
  
  for(int dim=d-1; dim>=0;dim--){ //loop over the all demensions
    if(p[dim]>1.5){ //do not need to do much when kron with respect to one thing
      p0 = p[dim]; //look at our one demension
      sv = sv-p0*p0;
      for(h=0; h<Bl; h+=p0)  //loop over the leftover demensions
      {
        x[h]=B[h]/A[sv];  //do first outside loop
        for (i=1; i<p0; i++)  //backsolve with respect to dim+1 the first time
        {
          x[h+i]=B[h+i];
          sv0 = sv+i*p0;   //speed up index reference
          for(k=i-1; k>=0; k--) x[h+i]-=A[sv0+k]*x[h+k];
          x[h+i]/=A[sv0+i];
        }
        
        y[h+p0-1]=x[h+p0-1]/A[sv+p0*p0-1];   //do first outside loop
        for (i=p0-2; i>=0; i--) //backsolve with respect to dim+1 the second time
        {
          y[h+i]=x[h+i];
          sv0=sv+i*p0;   //speed up index reference
          for(k=i+1; k<p0; k++) y[h+i]-=A[sv0+k]*y[h+k];
          y[h+i]/=A[sv0+i];
        }
      }
      c=0;
      for(i=0; i<p0; i++) for(h=0; h<Bl; h+=p0){B[c]=y[h+i]; c++;} //spinning the vector to right orientation
    }else{sv--; B=B/(A[sv]*A[sv]);} // kron with respect to one thing shortcut
  }
  
  B1=clone(B); //save this value, it will spin
  sv=A.size();
  
  for(int dim=d-1; dim>=0;dim--){ //loop over the all demensions
    if(p[dim]>1.5){  //do kron with respect to more than one thing
      p0 = p[dim]; //look at our one demension
      sv = sv-p0*p0;
      for(h = 0; h < Bl; h+=p0){
        for(i=0; i<p0; i++) {
          sv0=sv+i*p0;
          dBn(h+i) = 0;
          for (k = 0; k <p0; k++) dBn[h+i] += dA[sv+i*p0+k]*B1[h+k];
        }
      }
      
      for(h=0; h<Bl; h+=p0)  //loop over the leftover demensions
      {
        x[h]=dBn[h]/A[sv];  //do first outside loop
        for (i=1; i<p0; i++)  //backsolve with respect to dim+1 the first time
        {
          x[h+i]=dBn[h+i];
          sv0 = sv+i*p0;   //speed up index reference
          for(k=i-1; k>=0; k--) x[h+i]-=A[sv0+k]*x[h+k];
          x[h+i]/=A[sv0+i];
        }
        
        y[h+p0-1]=x[h+p0-1]/A[sv+p0*p0-1];   //do first outside loop
        for (i=p0-2; i>=0; i--) //backsolve with respect to dim+1 the second time
        {
          y[h+i]=x[h+i];
          sv0=sv+i*p0;   //speed up index reference
          for(k=i+1; k<p0; k++) y[h+i]-=A[sv0+k]*y[h+k];
          y[h+i]/=A[sv0+i];
        }
      }
      
      B2 = clone(B1);
      c=0;
      for(i=0; i<p0; i++) for(h=0; h<Bl; h+=p0){B1[c]=B2[h+i]; c++;} //spinning the vector to next orientation
      
      c=0;
      for(i=0; i<p0; i++) for(h=0; h<Bl; h+=p0){dBn[c]=y[h+i]; c++;} //spinning the vector to right orientation
      for(dim2=dim-1;dim2>=0;dim2--){
        p1 = p(dim2); //look at our one demension
       if(p(dim2)>1.5){  //spin if we have something to spin over
          dBn2 = clone(dBn);
          c=0;
          for(i=0; i<p1; i++) for(h=0; h<Bl; h+=p1){dBn[c]=dBn2[h+i]; c++;} //spinning the vector to right orientation
        }
      }
      for(i = 0; i < Bl; i++) dB(dim,i) = B[i]*(dA[sv]/(A[sv]*A[sv]));
    }else{sv--; for(i = 0; i < Bl; i++) dB(dim,i) = B[i]*(dA[sv]/(A[sv]*A[sv]));} // kron with respect to one thing shortcut
  }
  return dB;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
