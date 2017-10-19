#include "protoplanet.h"

protoplanet::protoplanet()
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
}

protoplanet::protoplanet(double M, double x, double y, double z, double vx, double vy, double vz)
{
    mass = M;
    position[0] = x;
    position[1] = y;
    position[2] = z;
    velocity[0] = vx;
    velocity[1] = vy;
    velocity[2] = vz;
    potential = 0.;
    kinetic = 0.;
}

double protoplanet::r()
{
    double x, y, z;

    x = this->position[0];
    y = this->position[1];
    z = this->position[2];

    return sqrt(x*x + y*y + z*z);
}

double protoplanet::distance(protoplanet otherPlanet)
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

double protoplanet::GravitationalForce(protoplanet otherPlanet,double Gconst)
{
    double r = this->distance(otherPlanet);
    if(r!=0) return Gconst*this->mass*otherPlanet.mass/(r*r);
    else return 0;
}

double protoplanet::Acceleration(protoplanet otherPlanet, double Gconst)
{
    double r = this->distance(otherPlanet);
    if(r!=0) return this->GravitationalForce(otherPlanet,Gconst)/(this->mass*r);
    else return 0;
}

double protoplanet::KineticEnergy()
{
    double velocity2 = (this->velocity[0]*this->velocity[0]) + (this->velocity[1]*this->velocity[1]) + (this->velocity[2]*this->velocity[2]);
    return 0.5*this->mass*velocity2;
}

double protoplanet::PotentialEnergy(protoplanet &otherPlanet, double Gconst, double epsilon)
{
    if(epsilon==0.0) return -Gconst*this->mass*otherPlanet.mass/this->distance(otherPlanet);
    else return (Gconst*this->mass*otherPlanet.mass/epsilon)*(atan(this->distance(otherPlanet)/epsilon) - (0.5*M_PI));
}

void protoplanet::EulerBinary(int integration_points, double final_time, std::ofstream &output)
{
    double h = final_time/((double) integration_points);
    double time = 0.0;
    time += h;
    double FourPi2 = 4 *M_PI*M_PI;

    double r3 = this->r() * this->r() * this->r();
    output << std::setw(5) << "time" << std::setw(15) << "x" << std::setw(15) << "y" << std::setw(15) << "z" << std::setw(15) << "vx" << std::setw(15) << "vy" << std::setw(15) << "vz" << std::endl;
    while (time <= final_time)
    {

        output << std::setw(5) << std::setprecision(3) << time;
        for(int j=0;j<3;j++) output << std::setw(15) << std::setprecision(8) << this->position[j];
        for(int j=0;j<3;j++) output << std::setw(15) << std::setprecision(8) << this->velocity[j];
        output << std::endl;

        this->position[0] += h * this->velocity[0];
        this->position[1] += h * this->velocity[1];
        this->position[2] += h * this->velocity[2];

        this->velocity[0] += -h * FourPi2 * this->position[0]/r3;
        this->velocity[1] += -h * FourPi2 * this->position[1]/r3;
        this->velocity[2] += -h * FourPi2 * this->position[2]/r3;

        r3 = this->r()*this->r()*this->r();

        time += h;

    }
    output << time ;
    for(int j=0;j<3;j++) output << std::setw(15) << std::setprecision(8) << this->position[j];
    for(int j=0;j<3;j++) output << std::setw(15) << std::setprecision(8) << this->velocity[j];
    output << std::endl;
}

void protoplanet::VerletBinary(int integration_points, double final_time, std::ofstream &output)
{
    double h = final_time/((double) integration_points);
    double h2 = h*h;
    double time = 0.0;
    time += h;
    double FourPi2 = 4 *M_PI*M_PI;

    double r3 = this->r() * this->r() * this->r();
    double * x_prev, * x_plus;
    x_prev = new double [3];
    x_plus = new double [3];
    output << std::setw(5) << "time" << std::setw(15) << "x" << std::setw(15) << "y" << std::setw(15) << "z" << std::setw(15) << "vx" << std::setw(15) << "vy" << std::setw(15) << "vz" << std::endl;
    while (time <= final_time)
    {
        for(int j=0;j<3;j++) { x_prev[j]=this->position[j]; this->position[j]=x_plus[j];}
        r3 = this->r()*this->r()*this->r();

        for(int j=0;j<3;j++) x_plus[j] = (2 - h2*FourPi2/r3) * this->position[j] - x_prev[j];
        for(int j=0;j<3;j++) this->velocity[j] = (x_plus[j] - x_prev[j])/(2 * h);

        output << std::setw(5) << std::setprecision(3) << time;
        for(int j=0;j<3;j++) output << std::setw(15) << std::setprecision(8) << this->position[j];
        for(int j=0;j<3;j++) output << std::setw(15) << std::setprecision(8) << this->velocity[j];
        output << std::endl;

        time += h;

    }
}
