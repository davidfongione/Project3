#ifndef PLANETSIMPLIFIED_H
#define PLANETSIMPLIFIED_H
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
using std::vector;
using std::string;

class planetsimplified
{
public:
    // Properties
    double mass;
    double position[3];
    double velocity[3];
    double potential;
    double kinetic;
    double angular[3];

    // Initializers
    planetsimplified();
    planetsimplified(double M,double x,double y,double z,double vx, double vy,double vz);

    // Functions
    double r();
    double distance(planetsimplified otherPlanet);
    double GravitationalForce(planetsimplified otherPlanet);
    double Acceleration(planetsimplified otherPlanet);
    double KineticEnergy();
    double PotentialEnergy(planetsimplified &otherPlanet, double epsilon);
    double * LinearMomentum();
    double * AngularMomentum();
    void EulerBinary(int integration_points, double final_time, double beta, string filename);
    void VerletBinary(int integration_points, double final_time, double beta, string filename);
    bool Bound();
};

#endif // PLANETSIMPLIFIED_H
