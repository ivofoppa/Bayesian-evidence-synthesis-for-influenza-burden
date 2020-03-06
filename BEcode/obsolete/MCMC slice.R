library(Rcpp)

cppFunction('double choosecpp(int n, int k)
            {
            int res ;
            if ( n < k ) 
            {
            return(0) ; 
            }
            else 
            {
            if ( k > ( n - k ) ){
              k=n-k ;
            }
            res = 1 ;
            for (int i =  1  ; i <= k   ; i = i+1) {
            res = res * (n - i + 1)/i ;
            }
            return(res) ;
            } 
            }
            
            ')


cppFunction('double binompdf(int x, int n, double p) {
    double f ;
    if ( (p < 0) || (p > 1) || (x > n) ) 
    {
    return(0) ; 
    }
    else 
    {
    f = choosecpp(n,x)*pow(p,x)*pow((1-p),(n-x)) ;
    return(f) ;
    }
  }
  ')

 

cppFunction('int factorial(unsigned int x){
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
            