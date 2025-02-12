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
        const double* C, const int& ldC);
    
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
int main(){
    const int n=1;
    const int N = pow(2,n);

    double* A=new double[N*N];
    double* B=new double[N*N];
    double* C=new double[N*N];
    double* D=new double[N*N];

    srand(time(0));
    for(int i=0;i<N;++i){
        for(int j=0; j<N; ++j){
            // Creating row-major matrices, good for c++
            A[i*N+j] =double(rand())/RAND_MAX;
            B[i*N+j] =double(rand())/RAND_MAX;
        }
    }
    // A and B matrix are transposed to make them column-major, which are compatible with BLAS
    F77NAME(dgemm)('T','T',N,N,N,1.0,A,N,B,N,0.0,C,N);
    
    for (int i=0;i<N;++i){
        for (int j=0;j<N;++j){
            D[i*N+j]=0;
            for (int k=0;k<N;++k){
                D[i*N+j]+=A[i*N+k]*B[k*N+j];
            }
        }
    }

    double dif=0.0;

    for (int i=0;i<N;++i){
        for (int j=0;j<N;++j){
            dif+=abs(C[i*N+j]-D[i*N+j]);
        }
    }


    writemat(A,N);
    writemat(B,N);
    writemat(C,N);
    writemat(D,N);

    std::cout<<dif;
}