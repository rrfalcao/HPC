#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
using namespace std;

#define F77NAME(x) x##_

extern "C"
{
    double F77NAME(dnrm2) (const int& n,
        const double *x, const int& incx);

    void F77NAME(dgbmv) (const char& trans, 
    const int& n, const int& m, const int& kl, const int& ku,   
                         const double& alpha, 
                         const double* A, const int& lda, // Banded matrix (leading dimension)
                         const double* x, const int& incx, 
                         const double& beta, 
                         double* y, const int& incy);
}

void prep_f(double *f, double dx, int Nx)
    {
        double x_i = 0.0;
        for (int i = 0; i < Nx; ++i)
        {

            x_i = i * dx;
            f[i] = sin(M_PI * x_i);
            
        }
        
    }
void prep_upper_banded_D(double *A,double Nx, int k)
    {
    

    for (int j = 0; j < Nx; ++j) {  // Column-major storage (traverse columns)
        if (j >= 3) A[(k - 3) + (k + 1) * j] = -1; // 3rd upper diagonal (theta)
        if (j >= 2) A[(k - 2) + (k + 1) * j] = 4;  // 2nd upper diagonal (gamma)
        if (j >= 1) A[(k - 1) + (k + 1) * j] = -5; // 1st upper diagonal (beta)
        A[k + (k + 1) * j] = 2; // Main diagonal (alpha)
    }
    }

void prep_lower_banded_D(double *A, int Nx, int k) {
    for (int j = 0; j < Nx; ++j) {  // Column-major storage (traverse columns)
        A[(k - 0) + (k + 1) * j] = 2; // Main diagonal (alpha)
        if (j < Nx - 1) A[(k - 1) + (k + 1) * j] = -5; // 1st lower diagonal (beta)
        if (j < Nx - 2) A[(k - 2) + (k + 1) * j] = 4;  // 2nd lower diagonal (gamma)
        if (j < Nx - 3) A[(k - 3) + (k + 1) * j] = -1; // 3rd lower diagonal (theta)
    }
}
void gen_B(double *B,double dx,double a, double b){
    // I want to fill the first 997 values of ~B with 0s and specify the last three indices
    double L=b-a;
    for (int i = 0; i < 998; ++i)
    {
        B[i] = 0;
    }
    B[998] = -sin(M_PI * (L+dx));
    B[999] = 4.0*sin(M_PI * (L+dx))-sin(M_PI * (L+2.0*dx));
    B[1000] = 5.0*sin(M_PI * (L+dx))+4.0*sin(M_PI * (L+2.0*dx))-sin(M_PI * (L+3.0*dx));

}
void gen_fdashdash(double *fdashdash, double dx, int Nx){
    for (int i = 0; i < Nx; ++i)
    {
        fdashdash[i] = -M_PI*M_PI*sin(M_PI*i*dx);
    }
}
void difference(double *A, double *B, double *C, int Nx){
    for (int i = 0; i < Nx; ++i)
    {
        C[i] = A[i] - B[i];
    }
}

int main(){
    int N=1001;
    double* D=new double[N*(N+1)/2];
    double* f=new double[N];
    double* B=new double[N];
    double* dif=new double[N];
    int k =3;
    double a=0.0;
    double b=1.0;
    double dx=(b-a)/(N-1);
    double* fdashdash=new double[N];
    
    prep_f(f,dx,N);
    prep_upper_banded_D(D,N,k);
    gen_B(B,dx,a,b);
    
    F77NAME(dgbmv)('N', N, N, 0,3, 1.0,D, k + 1, f, 1, 1.0, B, 1);
    gen_fdashdash(fdashdash,dx,N);
    difference(B,fdashdash,dif,N);
    
    double norm = F77NAME(dnrm2)(N,dif,1);
    std::cout<<"Norm of the difference is: "<<norm<<endl;

    delete[] D;
    delete[] f;
    delete[] B;
    delete[] dif;
    delete[] fdashdash;

}
