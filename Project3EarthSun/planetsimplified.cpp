#include "planetsimplified.h"
#include "time.h"

planetsimplified::planetsimplified()
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

planetsimplified::planetsimplified(double M, double x, double y, double z, double vx, double vy, double vz)
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
    angular[0] = 0.;
    angular[1] = 0.;
    angular[2] = 0.;
}

double planetsimplified::r()
{
    double x, y, z;

    x = this->position[0];
    y = this->position[1];
    z = this->position[2];

    return sqrt(x*x + y*y + z*z);
}

double planetsimplified::distance(planetsimplified otherPlanet)
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

double planetsimplified::GravitationalForce(planetsimplified otherPlanet)
{
    double r = this->distance(otherPlanet);
    if(r!=0) return (4 * M_PI * M_PI)*this->mass*otherPlanet.mass/(r*r);
    else return 0;
}

double planetsimplified::Acceleration(planetsimplified otherPlanet)
{
    double r = this->distance(otherPlanet);
    if(r!=0) return this->GravitationalForce(otherPlanet)/(this->mass*r);
    else return 0;
}

double planetsimplified::KineticEnergy()
{
    double velocity2 = (this->velocity[0]*this->velocity[0]) + (this->velocity[1]*this->velocity[1]) + (this->velocity[2]*this->velocity[2]);
    return 0.5*this->mass*velocity2;
}

double planetsimplified::PotentialEnergy(planetsimplified &otherPlanet, double epsilon)
{
    double G = 4 * M_PI *M_PI;
    if(epsilon==0.0) return -G*this->mass*otherPlanet.mass/this->distance(otherPlanet);
    else return (G*this->mass*otherPlanet.mass/epsilon)*(atan(this->distance(otherPlanet)/epsilon) - (0.5*M_PI));
}

double * planetsimplified::LinearMomentum()
{
    double * linear_momentum;
    linear_momentum = new double[3];

    for (int i = 0;i<3;i++) linear_momentum[i] = this->mass * this->velocity[i];

    return linear_momentum;
}

double * planetsimplified::AngularMomentum()
{
    double * angmom;
    angmom = new double[3];

    double * linmom = this->LinearMomentum();

    for (int i=0;i<3;i++)
        angmom[i] = this->position[(i+1)%3] * linmom[(i+2)%3] - this->position[(i+2)%3] * linmom[(i+1)%3];
    return angmom;
}

void planetsimplified::EulerBinary(int integration_points, double final_time, double beta, string filename)
{
    double h = final_time/((double) integration_points);
    double time = 0.0;
    double FourPi2 = 4 *M_PI*M_PI;
    planetsimplified center(1.,0.,0.,0.,0.,0.,0.);

    string filename_energy = "Energy_" + filename;
    std::ofstream output(filename);
    std::ofstream output_energy(filename_energy);

    double denominator = pow(this->r(), beta + 1.0);
    output << std::setw(6) << "time" << std::setw(15) << "x" << std::setw(15) << "y" << std::setw(15) << "z" << std::setw(15) << "vx" << std::setw(15) << "vy" << std::setw(15) << "vz" << std::endl;
    output_energy << std::setw(6) << "time" << std::setw(15) << "kinetic" << std::setw(15) << "potential" << std::setw(15) << "total" << std::setw(15) << "angular_x" << std::setw(15) << "angular_y" << std::setw(15) << "angular_z" << std::endl;

    clock_t start, finish;
    start = clock();

    while (time < final_time)
    {
        this->kinetic = this->KineticEnergy();
        this->potential = this->PotentialEnergy(center, 1.0e-6);
        for(int j=0;j<3;j++) this->angular[j] = this->AngularMomentum()[j];

        output_energy << std::setw(6) << std::setprecision(3) << time;
        output_energy << std::setw(15) << std::setprecision(8) << this->kinetic;
        output_energy << std::setw(15) << std::setprecision(8) << this->potential;
        output_energy << std::setw(15) << std::setprecision(8) << this->potential + this->kinetic;
        for(int j=0;j<3;j++) output_energy << std::setw(15) << std::setprecision(8) << this->angular[j];
        output_energy << std::endl;

        output << std::setw(6) << std::setprecision(3) << time;
        for(int j=0;j<3;j++) output << std::setw(15) << std::setprecision(8) << this->position[j];
        for(int j=0;j<3;j++) output << std::setw(15) << std::setprecision(8) << this->velocity[j];
        output << std::endl;

        this->position[0] += h * this->velocity[0];
        this->position[1] += h * this->velocity[1];
        this->position[2] += h * this->velocity[2];

        this->velocity[0] += -h * FourPi2 * this->position[0]/denominator;
        this->velocity[1] += -h * FourPi2 * this->position[1]/denominator;
        this->velocity[2] += -h * FourPi2 * this->position[2]/denominator;

        denominator = pow(this->r(), beta + 1.0);

        time += h;

    }
    finish = clock();
    std::cout << "elapsed time :" << (double) (finish-start)/((double) CLOCKS_PER_SEC) << std::endl;
    std::cout << "bound : " << this->Bound() << std::endl;
}

void planetsimplified::VerletBinary(int integration_points, double final_time, double beta, string filename)
{
    double h = final_time/((double) integration_points);
    double h2 = h*h;
    double time = 0.0;
    double TwoPi2 = 2 *M_PI*M_PI;
    planetsimplified center(1.,0.,0.,0.,0.,0.,0.);

    string filename_energy = "Energy_" + filename;
    std::ofstream output(filename);
    std::ofstream output_energy(filename_energy);

    double denominator = pow(this->r(), beta + 1.0);
    double denominator_prev;
    double * x_prev;
    x_prev = new double [3];
    output << std::setw(6) << "time" << std::setw(15) << "x" << std::setw(15) << "y" << std::setw(15) << "z" << std::setw(15) << "vx" << std::setw(15) << "vy" << std::setw(15) << "vz" << std::endl;
    output_energy << std::setw(6) << "time" << std::setw(15) << "kinetic" << std::setw(15) << "potential" << std::setw(15) << "total" << std::setw(15) << "angular_x" << std::setw(15) << "angular_y" << std::setw(15) << "angular_z";

    clock_t start, finish;
    start = clock();

    while (time < final_time)
    {
        this->kinetic = this->KineticEnergy();
        this->potential = this->PotentialEnergy(center, 1.0e-6);
        for(int j=0;j<3;j++) this->angular[j] = this->AngularMomentum()[j];

        output_energy << std::setw(6) << std::setprecision(3) << time;
        output_energy << std::setw(15) << std::setprecision(8) << this->kinetic;
        output_energy << std::setw(15) << std::setprecision(8) << this->potential;
        output_energy << std::setw(15) << std::setprecision(8) << this->potential + this->kinetic;

        for(int j=0;j<3;j++) output_energy << std::setw(15) << std::setprecision(8) << this->angular[j];
        output_energy << std::endl;

        output << std::endl;
        output << std::setw(6) << std::setprecision(3) << time;
        for(int j=0;j<3;j++) output << std::setw(15) << std::setprecision(8) << this->position[j];
        for(int j=0;j<3;j++) output << std::setw(15) << std::setprecision(8) << this->velocity[j];

        for(int j=0;j<3;j++) x_prev[j]=this->position[j];
        denominator_prev = denominator;

        for(int j=0;j<3;j++)
            this->position[j] = (1 - h2*TwoPi2/denominator) * this->position[j] + h * this->velocity[j];
        denominator = pow(this->r(), beta + 1.0);
        for(int j=0;j<3;j++)
            this->velocity[j] += - h * TwoPi2 * ( (x_prev[j]/denominator_prev) + this->position[j]/denominator );

        time += h;
    }
    finish = clock();
    std::cout << "elapsed time :" << (double) (finish-start)/((double) CLOCKS_PER_SEC) << std::endl;
    std::cout << "bound : " << this->Bound() << std::endl;

}

bool planetsimplified::Bound()
{
    return ((this->kinetic + this->potential) < 0.0);
}
