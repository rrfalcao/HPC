#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
using namespace std;

#define F77NAME(x) x##_
extern "C"
{
    

  void F77NAME(zgemm) (const char& transA, const char& transB,
                         const int& m, const int& n, const int& k,
                         const std::complex<double>& alpha, 
                         const std::complex<double>* A, const int& lda,
                         const std::complex<double>* B, const int& ldb,
                         const std::complex<double>& beta, 
                         std::complex<double>* C, const int& ldc);}

void complex_gen(std::complex<double>* A, int N){
    srand(time(0));
    double im=0.0;
    double re=0.0;
    for(int i=0;i<N;++i){
        for(int j=0; j<N; ++j){

            re=double(rand())/RAND_MAX-1;
            im=double(rand())/RAND_MAX-1;
            A[i*N+j] =std::complex<double>(re,im);
        }
        
    } 
}
void writemat (std::complex<double>* X,int N){
    for (int i=0;i<N;++i){
        for (int j=0;j<N;++j){
            std::cout<<X[i*N+j]<<" ";
        }
        std::cout<<endl;
    }

    std::cout<<endl;
    
}
bool is_hermitian(std::complex<double>* C, int N, double tol = 1e-10) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < i; ++j) { // Only check upper triangle
            if (std::abs(C[i * N + j] - std::conj(C[j * N + i])) > tol) {
                std::cout << "Element (" << i << "," << j << ") is not equal to its conjugate transpose" << std::endl;
                std::cout << C[i * N + j] << " != " << std::conj(C[j * N + i]) << std::endl;
                return false;
            }
        }
    }
    return true;
}
void conj_transpose(std::complex<double>* A,std::complex<double>* B,int N){
    for (int i=0;i<N;++i){
        for (int j=0;j<N;++j){
            B[i*N+j]=std::conj(A[j*N+i]);
        }
    }
}
int main(){
    int N=4;
    std::complex<double>*  A=new std::complex<double>[N*N];
    std::complex<double>*  B=new std::complex<double>[N*N];
    std::complex<double>*  C=new std::complex<double>[N*N];

    complex_gen(A,N);
    writemat(A,N);

    conj_transpose(A,B,N);
    writemat(B,N);

    const std::complex<double> alpha(1.0, 0.0);
    const std::complex<double> beta(0.0, 0.0);
    F77NAME(zgemm)('N','N',N,N,N,alpha,A,N,B,N,beta,C,N);
    writemat(C,N);

    std::cout<< "Hermitian C: "<<is_hermitian(C,N)<<endl;
}