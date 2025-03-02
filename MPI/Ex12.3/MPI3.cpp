#include <iostream>
#include <mpi.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>

using namespace std;


#define F77NAME(x) x##_

extern "C" {
    double F77NAME(ddot) (const int& n,
        const double *x, const int& incx,
        const double *y, const int& incy);
    double F77NAME(dnrm2) (const int& n,
        const double *x, const int& incx);
    }
void generate_vector(double* x,double* y,int n,int rank){
    srand(time(0)+rank);
    for (int i = 0; i < n; ++i)
    {
        x[i]=(double)rand()/RAND_MAX;
        y[i]=(double)rand()/RAND_MAX;
    }
}

int main(int argc, char **argv) {
    int err = MPI_Init(&argc, &argv);
    if (err != MPI_SUCCESS) {
        cout << "Failed to initialise MPI" << endl;
        return -1;
    }
    int rank = 0;
    int size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n = 1024 / size;
    double* x = new double[n];
    double* y = new double[n];
    generate_vector(x,y,n,rank);

    double dot = F77NAME(ddot)(n,x,1,y,1);
    double nrmx = F77NAME(dnrm2)(n,x,1);
    double nrmy = F77NAME(dnrm2)(n,y,1);

    double dot_r;
    double nrmx_r;
    double nrmy_r;

    MPI_Reduce(&dot,  &dot_r,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&nrmx, &nrmx_r, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&nrmy, &nrmy_r, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        cout << "Dot product: " << dot_r << endl;
        cout << "Norm of x: " << nrmx_r << endl;
        cout << "Norm of y: " << nrmy_r << endl;
    }

    delete[] x;
    delete[] y;

    // Finalise MPI.
    MPI_Finalize();
}