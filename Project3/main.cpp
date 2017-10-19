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

    protoplanet Earth(0.000003,1.,0.0,0.0,0.0,6.3,0.);


    ofstream output_one("Euler");
    ofstream output_two("Verlet");

    Earth.EulerBinary(n,final_time, output_one);
    cout << "Euler ok" << endl;
    Earth.VerletBinary(n, final_time, output_two);
    cout << "Verlet ok" << endl;
    return 0;
}
