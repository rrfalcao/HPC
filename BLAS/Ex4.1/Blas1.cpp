#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
using namespace std;

#define F77NAME(x) x##_

extern "C" {
    double F77NAME(ddot) (const int& n,
        const double *x, const int& incx,
        const double *y, const int& incy);
    int F77NAME(idamax)(const int& n,
        const double *x, const int& incx);
    double F77NAME(dasum)(const int& n,
        const double *x, const int& incx);
}
int main(){
    const int N=8;
    double* x= new double[N];
    
    double d = 0.0;
    double s = 0.0;
    int m = 0.0;

    srand(time(0));
    for (int i=0;i<N;++i){
        x[i]=(double)(rand())/ RAND_MAX;
        
    }
    

    d=F77NAME(ddot)(N,x,1,x,1);
    d=sqrt(d);
    s=F77NAME(dasum)(N,x,1);
    m=F77NAME(idamax)(N,x,1);
    
    std::cout<<"Quantity 1: "<< d <<endl;
    std::cout<<"Quantity 2: "<< s <<endl;
    std::cout<<"Quantity 3: index - "<< m <<" Value - "<< x[m-1]<<endl;
    for (int j=0;j<N;++j){
        std::cout<<x[j]<<endl;

    }
}

