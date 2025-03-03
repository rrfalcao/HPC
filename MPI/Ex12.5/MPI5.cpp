#include <iostream>
#include <mpi.h>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>

using namespace std;
#define F77NAME(x) x##_

extern "C" {
    void F77NAME(dgemv) (const char
    &trans, const int& m, const int& n, const double& alpha, 
    const double* A, const int& lda, const double* x, 
    const int& incx, const double& beta, double* y, const int& incy);
    }
void gen_matrix(double* x, int rows,int cols,int rank){
    srand(time(0)+rank);
    for(int i=0;i<rows;++i){
        for(int j=0;j<cols;++j){
            x[i*cols+j]=(double)rand()/RAND_MAX;
        }
    }
}
void gen_vector(double* x, int n){
    srand(time(0));
    for (int i = 0; i < n; ++i)
    {
        x[i]=(double)rand()/RAND_MAX;
    }
}

void writevec (double* X,int N){
    for (int i=0;i<N;++i){     
        std::cout<<X[i]<<endl;       
    }
    std::cout<<endl;
 
}
void writematrix(double* mat, int rows, int cols, int rank) {
    std::cout << "Process " << rank << " matrix:" << endl;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << mat[i * cols + j] << " ";
        }
        std::cout << endl;
    }
    std::cout << endl;
}

int main(int argc, char **argv){

    int err = MPI_Init(&argc, &argv);
    if (err != MPI_SUCCESS) {
        cout << "Failed to initialise MPI" << endl;
        return -1;
    }

    int rank = 0;
    int size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int n=pow(2,4);
    if (n % size != 0) {
        if (rank == 0) {
            cerr << "Matrix rows must be divisible by the number of processes." << endl;
        }
        MPI_Finalize();
        return -1;
    }

   
    int set = n/size;
    double* vec = new double[n];
    if (rank==0){
        
        gen_vector(vec,n);
        writevec(vec,n);
    }
    MPI_Bcast( vec , n , MPI_DOUBLE , 0 , MPI_COMM_WORLD);

    double* mat = new double[set*n];
    gen_matrix(mat,set,n,rank);
    
    
    double* result = new double[set]();
    F77NAME(dgemv)('N',set,n,1.0,mat,set,vec,1,0.0,result,1);

    // Debugging: print the local result on each process
    cout << "Process " << rank << " local result: ";
    for (int i = 0; i < set; ++i) {
        cout << result[i] << " ";
    }
    cout << endl;


    double* final_result = nullptr;
    if (rank == 0) {
        final_result = new double[n];
    }
    
    MPI_Gather(result, set, MPI_DOUBLE, final_result, set, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Debugging: print the gathered result on the root process
    if (rank == 0) {
        cout << "Gathered final result: ";
        for (int i = 0; i < n; ++i) {
            cout << final_result[i] << " ";
        }
        cout << endl;
    }

    if (rank == 0) {
        cout << "Matrix-vector multiplication result: " << endl;
        writevec(final_result, n);
        delete[] final_result;
    }

    // Clean up
    delete[] vec;
    delete[] mat;
    delete[] result;

    MPI_Finalize();
    return 0;
}