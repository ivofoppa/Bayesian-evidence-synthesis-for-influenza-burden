library(Rcpp)
# 
# fluhospls[seas] ~ dpois(rfluhospls[seas]*FSNpopls[seas]) ## The # of flu hosp in FSN with not fatal outcome; unobserved
# FSNfluhospls[seas] ~ dbin(ptarr[seas,1],fluhospls[seas])

cppFunction('int binomialN(unsigned int x){
            unsigned long long res=1;
            for (int i=1; i <= x; i++){
                res*=i;
            }
            return res;
}')


cppFunction('int bincoeff(int n, int k){
            long int sol=1;
if ( k > (n - k))
k = n - k;
            for (int i=0; i < k; i++){
                sol*=(n-i);
                sol/=(1+i);
            }
            return sol;
            
}')

cppFunction('double binlik(int n, int x, double p){
            int c = a - b;
            return c;
            }')

cppFunction('NumericVector cs(NumericVector x){
            int n=x.size();
            NumericVector y=x;
            for (int i = 1; i < n; i++){
            y[i] = y[i-1] + y[i];
            }
            return y;
            }'