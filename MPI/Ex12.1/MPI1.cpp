#include <iostream>
#include <mpi.h>
using namespace std;
int main(int argc, char * argv[]){
    MPI_Init(&argc,&argv);
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    int n;
    if (rank==0){
        
        std::cout << "Enter a number: ";
        std::cin >> n;
    }
    MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
    std::cout << "This is rank: " << rank <<" The number is: "<<n<<std::endl;
    MPI_Finalize(); 
    return 0;
    
}