#include <iostream>
#include "planet.h"
#include "System.h"
#include <string>

using namespace std;

int main()
{
beginning:
    int n;
    double final_time, b;
    cout << "Number of integration points : ";
    cin >> n;
    cout << "Final time : ";
    cin >> final_time;
    cout << "Value for beta : ";
    cin >> b;


    planet Earth(0.000003,9.776084084140449E-01, 2.252978885335488E-01,-1.428938088034023E-04,-4.062408436217812E-03, 1.671862543878177E-02 ,-1.806432427648686E-07);
    planet Sun(1.,0.,0.,0.,0.,0.,0.);
    planet Jupiter(0.00095,-4.627833494978254,-2.854717085491127,1.153513368222256e-01,-4.627833494978254,-2.854717085491127, 1.153513368222256E-01);
    planet Mars (0.00000033,-1.510556500968014, 7.016843752438682e-01, 5.158146872545729e-02,-5.328917141208186E-03,-1.151197416056758E-02,-1.105814511135091E-04);
    planet Venus (0.00000245,-5.122691445263290E-01, 5.057193230515091E-01, 3.642042956185625E-02,-1.417569102004634E-02, -1.460561845394588E-02, 6.174282365295180E-04);
    planet Saturn (0.000275,-4.106878033627468E-01,-1.004669994513129E+01, 1.910286401217743E-01,5.268654062679371E-03 ,-2.458299042190405E-04 ,-2.054537524214767E-04);
    planet Mercury (0.000000165,-3.896475834228487E-01,-2.816359732979998E-02, 3.305924612088227E-02,-3.438017523837078E-03,-2.681547759811458E-02,-1.876572045563162E-03);
    planet Uranus (0.000044,1.787905976395623E+01, 8.773022260100626,-1.990427195843487E-01,-1.761322338307447E-03, 3.347557466700209E-03, 3.516031928207541E-05);
    planet Neptune (0.0000515,2.860387845165146E+01,-8.854538501256135,-4.768631616180988E-01,9.072716295025177E-04, 3.016975531099665E-03,-8.333936157206056E-05);
    planet Pluto(0.00000000655,1.051353424271748E+01,-3.171623355522280E+01, 3.526949474225156E-01,3.046469278872717E-03, 3.317169711406130E-04,-9.167186249993658E-04);

    System three_body_problem((100.0));
    three_body_problem.add(Sun);
    three_body_problem.add(Earth);
    three_body_problem.add(Jupiter);

    System solar((three_body_problem));
    solar.add(Mars);
    solar.add(Mars);
    solar.add(Venus);
    solar.add(Saturn);
    solar.add(Mercury);
    solar.add(Uranus);
    solar.add(Neptune);
    solar.add(Pluto);
    three_body_problem.VelocityVerlet(3, n, final_time,3,b, 10e-6);
    solar.VelocityVerlet(3, n, final_time,3,b, 10e-6);

    string interact;
    cout << "Repeat ? ";
    cin >> interact;
    if (interact == "yes" || interact == "YES" || interact == "Yes" || interact == "Y" || interact == "y") goto beginning;
    return 0;
}
