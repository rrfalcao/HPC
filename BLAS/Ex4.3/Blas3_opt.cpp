#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
using namespace std;

#define F77NAME(x) x##_

extern "C" {

    void F77NAME(dcopy) (const int& n,
        const double *x, const int& incx,
        double *y, const int& incy);

    double F77NAME(ddot) (const int& n,
        const double *x, const int& incx,
        const double *y, const int& incy);
    
    void F77NAME(daxpy) (const int& n,
        const double& alpha,
        const double *x, const int& incx,
        double *y, const int& incy);

    void F77NAME(dscal) (const int& n,
        const double& alpha,
        double *x, const int& incx);

    
    void F77NAME(dspmv) (const char& uplo, 
        const int& n, 
        const double& alpha, 
        const double* AP,  // Packed storage matrix
        const double* x, const int& incx, 
        const double& beta, 
        double* y, const int& incy);
}

void conj_grad_packed(double* A,double* x, double *b,int n,double tol, int max_iter){
    double alpha, beta, rcur, rnext;
    double* r = new double[n];
    double* p = new double[n];
    double* Ap = new double[n];

    F77NAME(dspmv)('U',n,1.0,A,x,1,0.0,Ap,1); 
    F77NAME(dcopy) (n,b,1,r,1); 
    F77NAME(daxpy) (n,-1.0,Ap,1,r,1); // r0 = b - Ax0
    F77NAME(dcopy) (n,r,1,p,1); // p0 = r0

    for (int i=0; i<max_iter;i++){

        rcur=F77NAME(ddot) (n,r,1,r,1); 
        F77NAME(dspmv)('U',n,1.0,A,p,1,0.0,Ap,1);
        alpha=rcur/(F77NAME(ddot) (n,p,1,Ap,1)); // alpha = ((rk)^T*(rk))/((pk)^T*A(pk))

        F77NAME(daxpy) (n,alpha,p,1,x,1); // xk+1 = xk + alphak*pk
        F77NAME(daxpy) (n,-alpha,Ap,1,r,1); // rk+1 = rk - alphak*Apk

        rnext=F77NAME(ddot) (n,r,1,r,1); //if |rk+1|<tol, break
        if (sqrt(rnext)<tol){
            std::cout<< "Converged in "<< i+1<< " iterations"<<endl;
            break;
        }

        beta=rnext/rcur;
        F77NAME(dscal) (n,beta,p,1);
        F77NAME(daxpy) (n,1.0,r,1,p,1);
    }
    delete[] r;
    delete[] p;
    delete[] Ap;

}



void writemat (double* X,int N){
    for (int i=0;i<N;++i){
        for (int j=0;j<N;++j){
            std::cout<<X[i*N+j]<<" ";
        }
        std::cout<<endl;
    }

    std::cout<<endl;
    
}
void writevec (double* X,int N){
    for (int i=0;i<N;++i){
        
        std::cout<<X[i]<<endl;
           
    }

    std::cout<<endl;
    
}


int main(){
    const int N=61;
    double* A=new double[N*(N+1)/2];
    double* x=new double[N];
    double* b=new double[N];
    double* y=new double[N];

    srand(time(0));
    for(int i=0;i<N;++i){
        for(int j=0; j<=i; ++j){
            // Creating column-major matrix, good for BLAS because it's random it could be either upper or lower
            A[i*(i+1)/2 +j] =double(rand())/RAND_MAX; 
              
        }
        x[i]=double(rand())/RAND_MAX-1; 
        y[i]=0.0;
    } 
    

    F77NAME(dspmv)('U',N,1.0,A,x,1,0.0,b,1);

    conj_grad_packed(A,y,b,N,0.000001,100000);
    std::cout << "Expected x:\n";
    writevec(x,N);

    std::cout << "Computed x:\n";
    writevec(y,N);

    delete[] A;
    delete[] x;
    delete[] b;
    delete[] y;
}