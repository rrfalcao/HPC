#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <boost/timer/timer.hpp>
using namespace std;

#define F77NAME(x) x##_

extern "C"
{
    
        // LAPACK routine for banded matrix eigenvalue calculation
    void F77NAME(dsbev)(const char& jobz, const char& ul, const int& n,
                            const int& kd,
                            double* ab, const int& ldab, double* w,
                            double* z,  const int& ldz,
                            double* work, int* info );
        // LAPACK routine for banded matrix factorisation
     void F77NAME(dgbtrf)(const int& M, const int& N, const int& KL,
                     const int& KU, double* A, const int& lda,
                     int* ipiv, int* info);
        // LAPACK routine for banded matrix solve
     void F77NAME(dgbtrs)(const char& TRANS, const int& N, const int& KL,
                     const int& KU, const int& NRHS,
                     const double* A, const int& lda,
                     int* ipiv, double* B, const int& ldb, int* info);
}
    void prep_banded_A(double *A, double nu, double Nx, int lda)
    {
        
    double altnu = 1 - 2 * nu;

    for (int i = 0; i < Nx; ++i) {
        A[i*lda+1] = nu;
        A[i*lda+2] = altnu;
        A[i*lda+3] = nu;
    }}

    void prep_u(double *u, double dx, int Nx, double nu)
    {
        double x_i = 0.0;
        for (int i = 0; i < Nx; ++i)
        {

            x_i = (i + 1) * dx;
            u[i] = sin(M_PI * x_i);
        }
        
    }

    
    void writevec(double *X, int N)
    {
        for (int i = 0; i < N; ++i)
        {

            std::cout << X[i] << endl;
        }

        std::cout << endl;
    }
    void eig_calc(double *A, int Nx, int ldab)
    {
        char jobz = 'N';  // Eigenvalues only
        char uplo = 'L';  // Upper triangular
        int info;
        double* w = new double[Nx];  // Array for eigenvalues
        double* work = new double[3 * Nx];
        double* z=0;
        F77NAME(dsbev)(jobz, uplo, Nx, 1, A+2, ldab, w, z, 1, work, &info);
        if (info != 0)
        {
            std::cerr << "Error computing eigenvalues, info = " << info << std::endl;
        }
        else
        {
            std::cout << "Eigenvalues: ";
            for (int i = 0; i < Nx - 2; ++i)
            {
                std::cout << w[i] << " ";
            }
            std::cout << std::endl;
        }
        delete[] w;
        delete[] work;
        delete[] z;
    }
    void PACK_integrate(double* A,double* u,int n,int kl,int ku, int lda, int Nt){
        int info;
        boost::timer::cpu_timer timer;
        int* ipiv = new int[n]();
        F77NAME(dgbtrf)(n, n, kl, ku, A, lda, ipiv, &info);
        if (info) {
            cout << "Error in DGBTRF: " << info << endl;
           
        }

        // Time integrate
        for (int t = 0; t < Nt; t++) {
            F77NAME(dgbtrs)('N', n, kl, ku, 1, A, lda, ipiv, u, n, &info);
            if (info) {
                cout << "Error in DGBTRS: " << info << endl;
            }
        }
        timer.stop();
        std::cout << "Time taken for time integration: " << timer.format() << std::endl;
        delete[] ipiv;


    }
    

    int main()
    {
        const int Nx=20;
        const int kl = 1;
        const int ku = 1;
        const int lda = 2*kl + ku + 1;
        double* A = new double[lda*Nx]();
        double* w = new double[Nx]();
        double* z = 0;
        double* work = new double[3*Nx - 2]();
        int info = 0;
        double nu = 0.4;
        double Nt = 100;
        const double dx = 1.0/(Nx-1);
        double* u = new double[Nx]();
        
            
        prep_banded_A(A, nu, Nx, lda);
        prep_u(u, dx, Nx, nu);
        eig_calc(A, Nx,lda);
        
        // Modify A to represent the implicit operator (note change of signs)
        for (int i = 0; i < Nx; ++i) {
            A[i*lda+1] = -nu;
            A[i*lda+2] = 1.0 + 2.0*nu;
            A[i*lda+3] = -nu;
        }

        // To impose Dirichlet boundary conditions, we zero the off-diagonal entries
        A[lda+1] = 0.0;
        A[(Nx-2)*lda+3] = 0.0;
        
        
        PACK_integrate(A,u,Nx,kl,ku,lda,Nt);
        writevec(u,Nx);
        
        
        
    }
    