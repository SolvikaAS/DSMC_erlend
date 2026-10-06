#include <iostream>
#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include <fstream>

struct Particle
{
	double vx;
	double vy;
	double vz;
};


void collide(Particle& p1, Particle& p2)
{
	double Vx = 0.5*(p1.vx+p2.vx);
	double Vy = 0.5*(p1.vy+p2.vy);
	double Vz = 0.5*(p1.vz+p2.vz);

	double gx = p1.vx -p2.vx;
	double gy = p1.vy -p2.vy;
	double gz = p1.vz -p2.vz;

	double g = std::sqrt(gx*gx+gy*gy+gz*gz);

	static std::mt19937 generator (12345);
	static std::uniform_real_distribution<double> distribution(-1,1);

	double nx, ny, nz;
	double n2;

	do
	{
		nx = distribution(generator);
		ny = distribution(generator);
		nz = distribution(generator);

		n2 = nx*nx +ny*ny +nz*nz;
	}
	while (n2 ==0.0 || n2 >1.0);

	double n = std::sqrt(n2);

	nx /= n;
	ny /= n;
	nz /= n;

	double gx_new = g*nx;
	double gy_new = g*ny;
	double gz_new = g*nz;

	p1.vx = Vx + 0.5*gx_new;
	p1.vy = Vy + 0.5*gy_new;
	p1.vz = Vz + 0.5*gz_new;

	p2.vx = Vx - 0.5*gx_new;
        p2.vy = Vy - 0.5*gy_new;
        p2.vz = Vz - 0.5*gz_new;
}


int main()
{
	//simulation parameters
	const int N = 100000;
	const double T = 300;

	//constants
	const double kB = 1.380649e-23;
	const double m = 6.63e-26;

	const double sigma = std::sqrt(kB*T/m);

	std::cout << "DSMC 0D test\n";
	std::cout << "Particles: " <<N<<"\n";
	std::cout << "Temperature: " <<T<<" K\n";
	std::cout << "Velocity sigma: " << sigma <<"m/s\n";

	//random number generator
	std::mt19937 generator(42);

	//gaussian distribution
	std::normal_distribution<double> velocity_distribution(0.0,sigma);

	//create particles
	std::vector<Particle> particles;
	particles.reserve(N);

	for (int i=0;i<N;++i)
	{
		Particle p;
		p.vx = velocity_distribution(generator);
		p.vy = velocity_distribution(generator);
		p.vz = velocity_distribution(generator);

		particles.push_back(p);
	}

	double sum_v2 = 0.0;
	for (const Particle& p :particles)
	{
		double v2 = p.vx*p.vx
			  + p.vy*p.vy
			  + p.vz*p.vz;
		sum_v2 +=v2;
	}

	double mean_v2 = sum_v2/N;

	double measured_temperature= m*mean_v2/(3*kB);

	std::ofstream file("speed.dat");
	for (const Particle& p :particles)
	{
		double speed = std::sqrt(p.vx*p.vx
                          + p.vy*p.vy
                          + p.vz*p.vz);

		file << speed <<"\n";
	}

	file.close();

	std::cout << "\nResults:\n";
    	std::cout << "Mean v^2: " << mean_v2 << " m^2/s^2\n";
    	std::cout << "Measured temperature: " << measured_temperature << " K\n";

	Particle p1{100,0,0};
	Particle p2{-200,0,0};

	double px_before = m*(p1.vx + p2.vx);
	double E_before = 0.5*m*(p1.vx*p1.vx+ p1.vy*p1.vy+ p1.vz*p1.vz
			       + p2.vx*p2.vx+ p2.vy*p2.vy+ p2.vz*p2.vz);


	collide(p1, p2);

	std::cout << "p1.vx: " << p1.vx << "\n";
	std::cout << "p1.vy: " << p1.vy << "\n";
	std::cout << "p1.vz: " << p1.vz << "\n";

	double px_after = m*(p1.vx + p2.vx);
	double E_after = 0.5*m*(p1.vx*p1.vx+ p1.vy*p1.vy+ p1.vz*p1.vz
                               + p2.vx*p2.vx+ p2.vy*p2.vy+ p2.vz*p2.vz);

	std::cout << "Momentum before: " << px_before << "\n";
	std::cout << "Momentum after:  " << px_after << "\n";

	std::cout << "Energy before:   " << E_before << "\n";
	std::cout << "Energy after:    " << E_after << "\n";

	std::cout << "p2.vx: " << p2.vx << "\n";
        std::cout << "p2.vy: " << p2.vy << "\n";
        std::cout << "p2.vz: " << p2.vz << "\n";

	return 0;
}
