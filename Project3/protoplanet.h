#ifndef PROTOPLANET_H
#define PROTOPLANET_H
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
using std::vector;


class protoplanet
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
    protoplanet();
    protoplanet(double M,double x,double y,double z,double vx, double vy,double vz);

    // Functions
    double r();
    double distance(protoplanet otherPlanet);
    double GravitationalForce(protoplanet otherPlanet);
    double Acceleration(protoplanet otherPlanet);
    double KineticEnergy();
    double PotentialEnergy(protoplanet &otherPlanet, double epsilon);
    double * LinearMomentum();
    double * AngularMomentum();
    void EulerBinary(int integration_points, double final_time);
    void VerletBinary(int integration_points, double final_time);
    bool Bound();
};

#endif // PROTOPLANET_H
