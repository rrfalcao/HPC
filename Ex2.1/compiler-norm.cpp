#include <vector>
#include <cmath>
using namespace std;

double ComputeNorm(const vector<double>& x) {
    double result = 0.0;
    for (auto a : x) {
         result = result + a*a;
    }
    return sqrt(result);
}
