
library(Rcpp)

cppFunction("
            NumericMatrix DouBS(NumericMatrix A, NumericVector B, int n,  int m){
            
            NumericMatrix x(n,m);
            NumericMatrix y(n,m);
            
            for(int h = 0; h < m; h++)
            {
            for ( int i = 0; i < n; i++)
            {
            x(i,h) = B(i,h);
            for ( int k = i-1; k >=0; k-- ) x(i,h) -= A(k,i) * x(k,h);
            x(i,h) /= A(i,i);
            }
            
            
            for ( int i = n - 1; i >= 0; i-- )
            {
            y(i,h) = x(i,h);
            for ( int k = i + 1; k < n; k++ ) y(i,h) -= A(i,k) * y(k,h);
            y(i,h) /= A(i,i);
            }
            }
            
            return y;
            }")

cppFunction("
            NumericVector kronDBS(NumericVector A, NumericVector B, NumericVector p, int Al, int Bl,  int d){
            
            int sv = 0;

            NumericVector x(Bl);
            NumericVector y(Bl);
          
            sv = Al;
           for(int dim = d-1; dim>=0;dim--){
            sv = sv-p(dim)*p(dim);
            int p0 = p(dim);
            int n0 = Bl/p0;
          
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
            c+=1;
            }
            }
            }

            return B;
            } ")




A1 = t(cbind(c(1,3,2),c(0,1,4),c(0,0,6)))
#A1 =  t(cbind(c(1,0,0),c(0,1,0),c(0,0,1)))
A2 = t(cbind(c(1,2),c(0,4)))
A = kronecker(A1,A2)
b = c(4, 5, 10,6,7,14)
backsolve(A,backsolve(A,as.matrix(b),6,transpose = TRUE))

Av = c(as.vector(A1),as.vector(A2))
kronDBS(Av, b, c(3,2), length(Av), length(b),  2)
  