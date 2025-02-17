#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;


const double tol = 5e-8;

double form1 (int i){
    float sign=1.0;
    if(i%2==0){
        sign=1.0;
    }
    else{
        sign=-1.0;
    }
    return 4.0*(sign)/((2.0*i)+1.0);
}
double form2(int i){
    return (4.0*i*i/((4.0*i*i)-1.0));
}
int main(){
    int n=0;
    int m=1;
    double S1=0.0;
    double S2=1.0;
    bool conv1=false;
    bool conv2=false;
    while(!(conv1 && conv2)){
        
        if(!conv1){
            S1+=form1(n);
            conv1 = (abs(S1-M_PI)<=tol);
            if(!conv1) ++n;
        }
        if(!conv2){
            S2*=form2(m);
            conv2 = (abs((2.0*S2)-M_PI)<=tol);
            if(!conv2) ++m;
        }
        
        
        
    }
    std::cout<<"Formula 1 converged at n1="<<n<<endl;
    std::cout<<"Formula 2 converged at n2="<<m<<endl;
}