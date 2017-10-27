#include <iostream>
#include "planetsimplified.h"


using namespace std;

int main()
{
    int n;
    double final_time, vx, vy;
    cout << "Number of integration points : ";
    cin >> n;
    cout << "Final time : ";
    cin >> final_time;
    cout << "Initial value for Vx : ";
    cin >> vx;
    cout << "Initial value for Vy : ";
    cin >> vy;

    planetsimplified Earth(0.000003,1.,0,0.,vx, vy ,0.);

    for (double b=2.0;b<=3.0;b += 0.1)
    {
        cout << "beta : " << b <<endl;
        string terminaison = to_string(b) + ".txt";
        string filename_euler = "Euler beta = " + terminaison;
        string filename_verlet = "Verlet beta = " + terminaison;

        Earth.EulerBinary(n,final_time,b, filename_euler);
        cout << "Euler ok" << endl;
        Earth.VerletBinary(n, final_time,b, filename_verlet);
        cout << "Verlet ok" << endl;
}

    return 0;
}
