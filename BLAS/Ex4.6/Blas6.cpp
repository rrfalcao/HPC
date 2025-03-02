
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
    void F77NAME(dcopy) (const int& n,
        const double *x, const int& incx,
        double *y, const int& incy);

    void F77NAME(dsbmv) (const char& uplo, 
                         const int& n, 
                         const int& k,   // Bandwidth
                         const double& alpha, 
                         const double* A, const int& lda, // Banded matrix (leading dimension)
                         const double* x, const int& incx, 
                         const double& beta, 
                         double* y, const int& incy);
}
    void prep_banded_A(double *A, double nu, double Nx, int k)
    {
    double altnu = 1 - 2 * nu;

    for (int i = 0; i < Nx - 2; ++i)
    {
        A[((k + 1) * i) + 1] = altnu;
        if (i > 0)
        {
            A[(k + 1) * i] = nu;
        }
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

    void time_integrate(double *u, double *A, double Nt, int n,int k)
    {
        double* u_new = new double[n];
        boost::timer::cpu_timer timer;
        for (int i = 0; i < Nt; ++i)
        {
            
            F77NAME(dsbmv)('U', n, k, 1.0, A, k + 1, u, 1, 0.0, u_new, 1);
            F77NAME(dcopy)(n, u_new, 1, u, 1);
        }
        timer.stop();
        std::cout << "Time taken for time integration: " << timer.format() << std::endl;

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
    int main()
    {
        int Nx, Nt;
        double dt;
        int k = 1;
        double gamma0 = 0.0;
        double gamma1 = 0.0;

        std::cout << "Enter number of grid points (Nx): ";
        std::cin >> Nx;

        std::cout << "Enter time-step size (dt): ";
        std::cin >> dt;

        std::cout << "Enter number of time steps (Nt): ";
        std::cin >> Nt;

        std::cout << "Nx: " << Nx << ", dt: " << dt << ", Nt: " << Nt << std::endl;

        double *A = new double[(k + 1) * (Nx - 2)];
        double *u = new double[Nx - 2];
        double dx = 1.0 / (Nx - 1);
        double nu = dt / (dx * dx);

        prep_banded_A(A, nu, Nx, k);
        prep_u(u, dx, Nx, nu, gamma0, gamma1);
        time_integrate(u, A, Nt, Nx - 2,k);
        writevec(u,Nx-2);
    }