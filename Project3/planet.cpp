#include "planet.h"

planet::planet()
{
    mass = 1.;
    position[0] = 1.;
    position[1] = 0.;
    position[2] = 0.;
    velocity[0] = 0.;
    velocity[1] = 0.;
    velocity[2] = 0.;
    potential = 0.;
    kinetic = 0.;
    angular[0] = 0.;
    angular[1] = 0.;
    angular[2] = 0.;
}

planet::planet(double M, double x, double y, double z, double vx, double vy, double vz)
{
    mass = M;
    position[0] = x;
    position[1] = y;
    position[2] = z;
    velocity[0] = 365.256363*vx; //nasa data is in Au/day, whereas we're working on Au/year
    velocity[1] = 365.256363*vy;
    velocity[2] = 365.256363*vz;
    potential = 0.;
    kinetic = 0.;
    angular[0] = 0.;
    angular[1] = 0.;
    angular[2] = 0.;
}



double planet::r()
{
    double x, y, z;

    x = this->position[0];
    y = this->position[1];
    z = this->position[2];

    return sqrt(x*x + y*y + z*z);
}

double planet::distance(planet otherPlanet)
{
    double x1,y1,z1,x2,y2,z2,xx,yy,zz;

    x1 = this->position[0];
    y1 = this->position[1];
    z1 = this->position[2];

    x2 = otherPlanet.position[0];
    y2 = otherPlanet.position[1];
    z2 = otherPlanet.position[2];

    xx = x1-x2;
    yy = y1-y2;
    zz = z1-z2;

    return sqrt(xx*xx + yy*yy + zz*zz);
 }

double planet::GravitationalForce(planet otherPlanet,double Gconst)
{
    double r = this->distance(otherPlanet);
    if(r!=0) return Gconst*this->mass*otherPlanet.mass/(r*r);
    else return 0;
}

double planet::Acceleration(planet otherPlanet, double Gconst)
{
    double r = this->distance(otherPlanet);
    if(r!=0) return this->GravitationalForce(otherPlanet,Gconst)/(this->mass*r);
    else return 0;
}

double planet::KineticEnergy()
{
    double velocity2 = (this->velocity[0]*this->velocity[0]) + (this->velocity[1]*this->velocity[1]) + (this->velocity[2]*this->velocity[2]);
    return 0.5*this->mass*velocity2;
}

double planet::PotentialEnergy(planet &otherPlanet, double Gconst, double epsilon)
{
    if(epsilon==0.0) return -Gconst*this->mass*otherPlanet.mass/this->distance(otherPlanet);
    else return (Gconst*this->mass*otherPlanet.mass/epsilon)*(atan(this->distance(otherPlanet)/epsilon) - (0.5*M_PI));
}

double * planet::LinearMomentum()
{
    double * linear_momentum;
    linear_momentum = new double[3];

    for (int i = 0;i<3;i++) linear_momentum[i] = this->mass * this->velocity[i];

    return linear_momentum;
}

double * planet::AngularMomentum()
{
    double * angmom;
    angmom = new double[3];

    double * linmom = this->LinearMomentum();

    for (int i=0;i<3;i++)
        angmom[i] = this->position[(i+1)%3] * linmom[(i+2)%3] - this->position[(i+2)%3] * linmom[(i+1)%3];
    return angmom;
}
