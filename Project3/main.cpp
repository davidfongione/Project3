#include <iostream>
#include "protoplanet.h"

using namespace std;

int main()
{
    int n;
    double final_time;
    cout << "Number of integration points : ";
    cin >> n;
    cout << "Final time : ";
    cin >> final_time;

    protoplanet Earth(0.000003,1.,0.0,0.0,0.0,6.28,0.);

    Earth.EulerBinary(n,final_time);
    cout << "Euler ok" << endl;
    Earth.VerletBinary(n, final_time);
    cout << "Verlet ok" << endl;
    return 0;
}
