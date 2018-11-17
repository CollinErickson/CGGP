
cppFunction("
            NumericVector MEEAGA(NumericMatrix Fmat, NumericMatrix I, NumericVector w, NumericVector S, int uc, int L, int d){
            double V = 1;
            
            for (int k = 0; k <= (L-1); k++) {
            for (int i = 0; i <= (uc-1); i++) {
            V=1;
            for (int j = 0; j <= (d-1); j++){
            V = V*(1-Fmat(I(i,j)-1,k));
            }
            S(k) += -w(i)*V;
            } 
            }
            return S;
            }")


