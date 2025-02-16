#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
using namespace std;

#define F77NAME(x) x##_

extern "C" {
    void F77NAME(dgemm) (const char& transA,const char& transB, 
        const int& m, const int& n, const int& k, 
        const double& alpha,
        const double* A, const int& lda,
        const double* B, const int& ldB,
        const double& beta, 
        double* C, const int& ldC);

    void F77NAME(dcopy) (const int& n,
        const double *x, const int& incx,
        double *y, const int& incy);

    double F77NAME(ddot) (const int& n,
        const double *x, const int& incx,
        const double *y, const int& incy);
    
    double F77NAME(daxpy) (const int& n,
        const double& alpha,
        const double *x, const int& incx,
        double *y, const int& incy);

    void F77NAME(dscal) (const int& n,
        const double& alpha,
        const double *x, const int& incx);

    void F77NAME(dgemv) (const char& trans, 
        const int& m, const int& n, 
        const double& alpha,
        const double* A, const int& lda,
        const double* x, const int& incx,
        const double& beta, double* y,
        const int& incy);
    
}

void conj_grad(double* A,double* x, double *b,int n,double tol, int max_iter){
    double alpha, beta, rcur, rnext;
    double* r = new double[n];
    double* p = new double[n];
    double* Ap = new double[n];

    F77NAME(dgemv)('N',n,n,1.0,A,n,x,1,0.0,Ap,1); 
    F77NAME(dcopy) (n,b,1,r,1); 
    F77NAME(daxpy) (n,-1.0,Ap,1,r,1); // r0 = b - Ax0
    F77NAME(dcopy) (n,r,1,p,1); // p0 = r0

    for (int i=0; i<max_iter;i++){

        rcur=F77NAME(ddot) (n,r,1,r,1); 
        F77NAME(dgemv)('N',n,n,1.0,A,n,p,1,0.0,Ap,1);
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
    const int N=4;
    double* M=new double[N*N];
    double* A=new double[N*N];
    double* x=new double[N];
    double* b=new double[N];
    double* y=new double[N];
    srand(time(0));
    for(int i=0;i<N;++i){
        for(int j=0; j<N; ++j){
            // Creating column-major matrix, good for BLAS
            M[j*N+i] =double(rand())/RAND_MAX-1;    
        }
        x[i]=double(rand())/RAND_MAX-1; 
        y[i]=0.0;
    } 

    
    F77NAME(dgemm)('T','N',N,N,N,1.0,M,N,M,N,0.0,A,N);
    F77NAME(dgemv)('N',N,N,1.0,A,N,x,1,0.0,b,1);

    conj_grad(A,y,b,N,0.000001,10000);
    writevec(x,N);
    writevec(y,N);
}