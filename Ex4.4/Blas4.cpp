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

    
    void F77NAME(dsymv) (const char& uplo, 
    const int& n, 
    const double& alpha, 
    const double* A, const int& lda,
    const double* x, const int& incx, 
    const double& beta, double* y, const int& incy);
}
void conj_grad(double* A,double* x, double *b,int n,double tol, int max_iter){
    double alpha, beta, rcur, rnext;
    double* r = new double[n];
    double* p = new double[n];
    double* Ap = new double[n];

    F77NAME(dsymv)('U',n,1.0,A,n,x,1,0.0,Ap,1); 
    F77NAME(dcopy) (n,b,1,r,1); 
    F77NAME(daxpy) (n,-1.0,Ap,1,r,1); // r0 = b - Ax0
    F77NAME(dcopy) (n,r,1,p,1); // p0 = r0

    for (int i=0; i<max_iter;i++){

        rcur=F77NAME(ddot) (n,r,1,r,1); 
        F77NAME(dsymv)('U',n,1.0,A,n,p,1,0.0,Ap,1);
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

void print_upper_tri(const double* A, int N) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i <= j) {
                std::cout  << A[j * (j + 1) / 2 + i] << " "; // Upper triangular stored
            } else {
                std::cout << A[i * (i + 1) / 2 + j] << " "; // Use symmetry (A_ji)
            }
        }
        std::cout << std::endl;
    }
}
void generate_reduced_helmholtz_system(double* A,double* f,double* u_an,int n,double a,double b,double lambda,double u0,double un1){
    double dx = (b-a)/(n-1);
    double beta = 1.0/(dx*dx);
    double alpha = -(2.0*beta + lambda);
     for (int i=0; i<n-2;++i){
        for(int j=0;j<=i;++j){
            if (j==i){
                A[i*(i+1)/2 +j] =alpha;
            }
            else if (j==i-1) 
            {
                A[i*(i+1)/2 +j] =beta;
            }
            else{
                A[i*(i+1)/2 +j] =0.0;
            }
              
        }
        for (int k=0;k<n-2;++k){
            double x_i=a+i*dx;
            f[i]=-(lambda+(M_PI*M_PI))*sin(M_PI*x_i);//Forcing function
            u_an[i]=sin(M_PI*x_i);
        }
        
        
     }
        f[0]-=u0*beta;
        f[n-3]-=un1*beta;
        

}
void writevec (double* X,int N){
    for (int i=0;i<N;++i){
        
        std::cout<<X[i]<<endl;
           
    }

    std::cout<<endl;
    
}
int main(){
    int N=21;
    double a = 0.0;
    double b=2.0;
    double lambda =1.0;
    double u0=0.0;  
    double un1=0.0;
    double* A= new double[(N-1)*(N-2)/2];
    double* f= new double[N-2];
    double* uanalytic=new double[N-2];
    double* u=new double[N-2];
    generate_reduced_helmholtz_system(A,f,uanalytic,N,a,b,lambda,u0,un1);
    std::cout << "Upper triangular matrix (stored in A):\n";
    print_upper_tri(A, N-2);
    // writevec(f,N-2);
    conj_grad(A,u,f,N-2,0.000001,10000);
    std::cout << "Expected u:\n";
    writevec(uanalytic,N-2);
    std::cout << "Computed u:\n";
    writevec(u,N-2);
}