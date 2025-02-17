
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
using namespace std;

#define F77NAME(x) x##_
extern "C"{
    double F77NAME(dznrm2) (const int& n,
        std::complex<double>*x, const int& incx);

    void F77NAME(zgemv) (const char& trans, 
        const int& m, const int& n, 
        const double& alpha,
        std::complex<double>* A, const int& lda,
        std::complex<double>* x, const int& incx,
        const double& beta, std::complex<double>* y,
        const int& incy);
    
}
void gen_x(std::complex<double>* x,int N, double dx,double start) {//Generate vector x
    double x_i = start;
    // Example complex function
    for (int i = 0; i < N; ++i) {
        x_i = start + i * dx;
        x[i]=cos(2*M_PI*x_i) +cos(6*M_PI*x_i)+1;
        
    }
}
void difference(std::complex<double> *A, std::complex<double> *B, std::complex<double> *C, int Nx){//find the difference between two matrices
    for (int i = 0; i < Nx; ++i)
    {
        C[i] = A[i] - B[i];
    }
}
void gen_W(std::complex<double>* W,int N){ //Generate matrix W
    for (int j = 0; j < N; ++j) {
        for (int k = 0; k <N; ++k) {
            W[j*N +k] = std::complex<double>((1/sqrt(N))*sin(-2*j*k*M_PI/N), (1/sqrt(N))*cos(-2*j*k*M_PI/N));
        }
    }
}

void writevec (std::complex<double>* X,int N){ //Write a vector to terminal
    for (int i=0;i<N;++i){
        
        std::cout<<X[i]<<endl;
           
    }

    std::cout<<endl;
    
}
void write_halfvec (std::complex<double>* X,int N){ //Write half a vector to terminal
    for (int i=0;i<N/2;++i){
        
        std::cout<<X[i].real() <<" + "<<X[i].imag()<<"i"<<endl;
           
    }

    std::cout<<endl;
    
}
void writemat (std::complex<double>* X,int N){ // Write a matrix to terminal
    for (int i=0;i<N;++i){
        for (int j=0;j<N;++j){
            std::cout<<X[i*N+j]<<" ";
        }
        std::cout<<endl;
    }

    std::cout<<endl;
    
}

int main() {
    const int N = 16;
    const double start = 0.0;
    const double end = 1.0;
    const double dx = (end - start) / (N - 1);
    
    std::complex<double>*  W=new std::complex<double>[N*N];
    std::complex<double>*  X=new std::complex<double>[N];
    std::complex<double>* x= new std::complex<double>[N];
    std::complex<double>* dif= new std::complex<double>[N];
    
    gen_x(x, N, dx, start);
    gen_W(W,N);
    
    F77NAME(zgemv)('N', N, N, 1.0, W, N, x, 1, 0.0, X, 1);
    write_halfvec(X,N);
    F77NAME(zgemv)('N', N, N, 1.0, W, N, X, 1, 0.0, X, 1);
    F77NAME(zgemv)('N', N, N, 1.0, W, N, X, 1, 0.0, X, 1);
    F77NAME(zgemv)('N', N, N, 1.0, W, N, X, 1, 0.0, X, 1);
    

    difference(X,x,dif,N);
    
    double norm = F77NAME(dznrm2)(N,dif,1);
    std::cout<<"Norm of the difference between x and y is: "<<norm<<endl;

delete[] W;
delete[] x;
delete[] X;
delete[] dif;

    

}