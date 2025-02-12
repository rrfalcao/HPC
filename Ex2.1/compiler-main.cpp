#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;

#include "cblas.h"
#include "compiler-norm.h"


/*
 * The main routine
 */
int main() {
    // Declare variables
    int y = 0;
#if 0
    int n;
    cin >> n
#else
    int n = 100000000;
#endif

    if (n == 0) {
        cout << "N must be > 0" << endl;
    }

    vector<double> x(n, 0.0);

    // Seed the random number generate with the current epoch time
    srand(time(0));

    // Generate random numbers and print them to the screen
    generate(x.begin(), x.end(), [] { return (double)rand()/RAND_MAX; });

    // Print out the values computed by BLAS and manually
    cout << "BLAS Norm:   " << cblas_dnrm2(n, &x[0], 1) << endl;
    cout << "Manual norm: " << ComputeNorm(x) << endl;

    return 0;
}
