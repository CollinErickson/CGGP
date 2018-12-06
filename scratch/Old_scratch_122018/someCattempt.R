install.packages("inline")
  
library(Rcpp)
cppFunction('int add(int x, int y, int z) {
  int sum = x + y + z;
            return sum;
            }')
# add works like a regular R function
add(1, 2, 3)



cppFunction("
NumericVector exfun(NumericVector x, int i){
x = x*i;
return x;
}")
exfun(1:5, 3)

cppFunction("
            NumericMatrix exfun(NumericMatrix x, NumericVector k){
            x(k(0),k(1)) = x(k(0),k(1))*3;
            return x;
            }")


cppFunction("
            NumericMatrix exfun(NumericMatrix x, NumericVector k){
            
            return val;
            }")
exfun(cbind(1:5,2:6), as.integer(c(1,1)))
exfun(cbind(1:5,2:6), as.integer(c(1,1)))

R = array(c(5,2,2))
cppFunction("
            NumericVector MEEAGA(NumericMatrix Fmat, NumericVector I, NumericVector w, NumericVector S, int uc, int L, int d){
            double V = 1;
          
            for (int k = 0; k <= L; k++) {
            for (int i = 0; i <= uc; i++) {
                V=1;
                for (int j = 0; j <= d; j++){
                  V = V*(1-Fmat(I(i,j),j));
                }
                S(i) += w(i)*V;
            } 
            }
               return S;
            }")
exfun(cbind(1:5,2:6), as.integer(c(1,1)))

