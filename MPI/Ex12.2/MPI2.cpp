#include <iostream>
#include <mpi.h>
using namespace std;
int main(int argc, char * argv[]){
    MPI_Init(&argc,&argv);
    
    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    int n;
    if (rank==0){
        
        std::cout << "Enter a number: ";
        std::cin >> n;
        MPI_Ssend(&n,1,MPI_INT,1,0,MPI_COMM_WORLD);
    }
    if (rank!=0){
        MPI_Recv(&n,1,MPI_INT,rank-1,rank-1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        std::cout << "This is rank: " << rank <<" The number is: "<<n<<std::endl;
        if (rank!=size-1){
            MPI_Ssend(&n,1,MPI_INT,rank+1,rank,MPI_COMM_WORLD);
        }
        else{
            MPI_Ssend(&n,1,MPI_INT,0,rank,MPI_COMM_WORLD);
        }

    }
    
    if (rank==0){ 
        MPI_Recv( &n,1,MPI_INT,size-1,size-1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        std::cout << "This is rank: " << rank <<" The loop is over!"<<std::endl;
        
        
        }
    MPI_Finalize();   
    return 0;
    }
    

    000000000000000000000000000
    