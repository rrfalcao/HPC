#include <iostream>
#include <mpi.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
using namespace std;

double f(double x){
    return 1.0/(1.0+(x*x));
}

double local_trap(double a, double b, int n) {
    if (n==0) return 0.0;
    double h = (b - a) / n;
    double integral = (f(a) + f(b)) / 2.0; // First and last points
    for (int i = 1; i < n; i++) {
        integral += f(a + i * h); // Middle points
    }
    return integral * h;
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

    int intervals = 0;
    if(rank==0){
        string input;
        std::cout << "Enter the number of intervals: ";
        std::cin >> input;
        intervals=stoi(input);
    }
    
    
    MPI_Bcast(&intervals, 1, MPI_INT, 0, MPI_COMM_WORLD);
    double dx = 1.0/intervals;

    int base_n = intervals / size;
    int remainder = intervals % size;
    int n = base_n + (rank < remainder ? 1 : 0);
    double a=0.0;
    double b=0.0;

    // Calculate interval to integrate for each rank
    if (rank < remainder) {
        a = rank * (n) * dx;
        b = a + (n)*dx;
    }
    else {
        a = n * remainder * dx + (rank - remainder) * n * dx;
        b = a + n * dx;
    }
    
    
    // Calculate local contribution to integral
    double local_integral = local_trap(a, b, n);
    double full_integral=0.0;

    MPI_Reduce(&local_integral, &full_integral, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        cout << "Approx pi: " << 4.0*full_integral << endl;
    }

    MPI_Finalize();
}