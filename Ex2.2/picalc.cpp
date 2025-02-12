#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;

const double realpi = 3.1415926536;
const double tol = 5e-6;

double form1 (int n){
    double tot =0;
    for (int i = 0; i<=n;i++){
        tot += (pow(-1.0,i))/((2.0*i)+1.0);
    }
    return 4.0*tot;
}
double form2(int n){
    double tot=1;
    for (int i =1; i<=n;i++){
        tot*= 4.0*i*i/((4.0*i*i)-1.0);
    }
    return 2.0*tot;
}
int main(){
    int n=1;
    bool trial = true;
    double dif=0.0;
    while(trial){
        dif = abs(form1(n)-realpi);
        if (dif<tol){
            trial=false;
        }
        else{n++;}
    }
    std::cout<<"Formula 1 converged at n1="<<n<<endl;

    n=1;
    trial = true;
    while(trial){
        dif = abs(form2(n)-realpi);
        if (dif<tol){
            trial=false;
        }
        else{n++;}
    }
    std::cout<<"Formula 2 converged at n2="<<n;
}