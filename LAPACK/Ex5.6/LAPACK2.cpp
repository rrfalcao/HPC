#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
using namespace std;

#define F77NAME(x) x##_

extern "C"
{
    void F77NAME(dcopy) (const int& n,
        const double *x, const int& incx,
        double *y, const int& incy);

    void F77NAME(dgbmv) (const char& trans, 
                         const int& m,
                         const int& n, 
                         const int& kl,  
                         const int& ku, // Bandwidth
                         const double& alpha, 
                         const double* A, const int& lda, // Banded matrix (leading dimension)
                         const double* x, const int& incx, 
                         const double& beta, 
                         double* y, const int& incy);
    void F77NAME(dsbev)(char* jobz, char* uplo, int* n, 
        int kd, double* ab, int* ldab, double* w, 
        double* z, int ldz, double* work, int* info);
    void F77NAME(dpbtrf)(char* uplo, int* n, int* kd, double* ab, int* ldab, int* info);
}
    void prep_banded_A(double *A, double nu, double Nx, int kl, int ku)
    {
        int lda = kl + ku + 1;
    double altnu = 1 - 2 * nu;

    for (int i = 0; i < Nx - 2; ++i)
    {   
        A[i*lda] = nu;
        A[i*lda + 1] = altnu;
        A[i*lda + 2] = nu;
        
    }}

    void prep_u(double *u, double dx, int Nx, double nu, double u0, double u1)
    {
        double x_i = 0.0;
        for (int i = 0; i < Nx - 2; ++i)
        {

            x_i = (i + 1) * dx;
            u[i] = sin(M_PI * x_i);
        }
        u[0] -= u0 * nu;
        u[Nx - 3] -= u1 * nu;
    }

    void time_integrate(double *u, double *A, double Nt, int n,int kl,int ku)
    {
        double* u_new = new double[n];
        
        for (int i = 0; i < Nt; ++i)
        {
            
            F77NAME(dgbmv)('N', n,n, kl,ku, 1.0, A, kl+ku+1, u, 1, 0.0, u_new, 1);
            F77NAME(dcopy)(n, u_new, 1, u, 1);
        }

        delete[] u_new;
    }
    void writevec(double *X, int N)
    {
        for (int i = 0; i < N; ++i)
        {

            std::cout << X[i] << endl;
        }

        std::cout << endl;
    }
    void eig_calc(double *A, double *u, int Nx,char uplo, int ldab)
    {
        char jobz = 'N';  // Eigenvalues only
        int info;
        double* w = new double[Nx];  // Array for eigenvalues
        double* work = new double[3 * (Nx)];
        F77NAME(dsbev)(&jobz, &uplo, &Nx, 1, A, &ldab, w, u, 1, work, &info);
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
    }
    void PACK_integrate(double* A,int n,int kd, int ldab, char uplo){
        int info;
        F77NAME(dpbtrf)(&uplo, &n, &kd, A, &ldab, &info);
        if (info != 0)
        {
            std::cerr << "Error computing eigenvalues, info = " << info << std::endl;
        }

    }
    void print_banded_matrix(double* A, int n, int kl, int ku)
    {
        int lda = kl + ku + 1; // Number of rows in the banded matrix
        std::cout << "Banded Matrix (Column-major order):" << std::endl;

        for (int i = 0; i < lda; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                // Print the banded matrix in memory layout
                std::cout << A[j * lda + i] << "\t";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    int main()
    {
        int Nx=20;
        int Nt=100;
        double dt=0.001;
        int kl = 1;
        int ku = 1;
        double gamma0 = 0.0;
        double gamma1 = 0.0;

        

        double *A = new double[(kl +ku + 1) * (Nx - 2)];
        double *u = new double[Nx - 2];
        double dx = 1.0 / (Nx - 1);
        double nu = dt / (dx * dx);
        char uplo = 'U';
        // int ldab = k + 1;
        
      
        
        prep_banded_A(A, nu, Nx, kl,ku);
        prep_u(u, dx, Nx, nu, gamma0, gamma1);

        // eig_calc(A, u, Nx-2, uplo,kl+ku+1);
        // PACK_integrate(A,Nx-2,k,ldab,uplo);
        print_banded_matrix(A, Nx - 2, kl, ku);
        time_integrate(u, A, Nt, Nx - 2,kl,ku);
        writevec(u,Nx-2);
    }
    