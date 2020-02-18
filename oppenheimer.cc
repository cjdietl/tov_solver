// This code integrates the Oppenheimer-Volkov equation
// The calculation follows http://prola.aps.org/abstract/PR/v55/i4/p374_1
#include <iostream>
#include <fstream>
#include "math.h"

using namespace std;

// Maximum amount of steps
#define STEPS 150000

// Increase Core Density if we have more steps than this
#define MAX_STEPS 20000

// // Dimensionless mass unit -> EarthMass
// #define MASS_C 0.302481 // EarthMass
// // Earth Mass
// #define M_EARTH 5.9742E24 // Kilogram
// // Dimensionless mass unit -> Kilogram
// #define MASS_C_KG (MASS_C*M_EARTH) // Kilogram
// //Dimensionless length unit -> cm
// #define RADIUS_C 0.134196  // Centi Meter
// //Dimensionless length unit -> Meter
// #define RADIUS_C_M RADIUS_C/100.0 // Meter

// ****************************************
// Mass in eV
//#define ParticleMass 3E12
double ParticleMass;
// ****************************************

// Earth Mass
#define M_EARTH 5.9742E24 // Kilogram
//Gravitational Constant
#define G 6.67428E-11 // m^3/(kg * s^2)
// Gravitational Acceleration Earth
#define G_EARTH 9.80665
//Speed of Light
#define C 2.99792458E8
// Mass conversion eV -> gram
#define EV_TO_GRAM 1.78266E-33
// Mass in gram
//#define ParticleMassGram (ParticleMass*EV_TO_GRAM)
double ParticleMassGram;
// Length conversion dimensionless -> cm
//#define DLtoCM (3.83815E-42/(ParticleMassGram*ParticleMassGram))
double DLtoCM;
// Length conversion dimensionless -> m
//#define DLtoM ((3.83815E-42/(ParticleMassGram*ParticleMassGram))*(1E-2))
double DLtoM;
// Mass Conversion dimensionless -> kg
//#define DLtoKG (5.16843E-17/(ParticleMassGram*ParticleMassGram))
double DLtoKG;
// Earth Mass
#define M_EARTH 5.9742E24 // Kilogram
// Solar Mass
#define M_SOL 1.988435E30 // Kilogram
//Mass is considered constant if error is below this
#define MASS_EPS 1E-11
//Density Conversion g/cm^3 -> dimensionless
//#define DensitytoDL (1.09397E-111/(ParticleMassGram*ParticleMassGram*ParticleMassGram*ParticleMassGram))
double DensitytoDL;
//Polytrope Exponential (gamma)
#define POLY_EXP 5.0/3.0
// Scaling Rho
#define RHO_ZERO 1E28
// If t/tzero below T_EPS and mass constant, it is assumed that the surface is reached
#define T_EPS 5E-2
//Degeneracy Factor
#define DEG 2

#define BISECTION_ERR 1E-6

// Constants for utility function get Sign
#define NEG -1
#define POS +1

#define DLTOKGCONST 5.16843E-17
#define DLTOCMCONST 3.83815E-42

void recalcConst()
{
	ParticleMassGram = ParticleMass*EV_TO_GRAM;
	DLtoCM = (DLTOCMCONST/(ParticleMassGram*ParticleMassGram));
	DLtoM = ((DLTOCMCONST/(ParticleMassGram*ParticleMassGram))*(1E-2));
	DLtoKG = (DLTOKGCONST/(ParticleMassGram*ParticleMassGram));
	DensitytoDL = (1.09397E-111/(ParticleMassGram*ParticleMassGram*ParticleMassGram*ParticleMassGram));
}
	
short getSign(double x)
{
	if (x < 0.0)
		return NEG;
	else
		return POS;
}

// Provides the equation for the bisection function
double edensity_bi(double t, double lhs)
{	
	return ((1.0/(4.0*M_PI))*(sinh(t)-t)) - DensitytoDL*lhs;
}

//Dimensionless energy density
double edensity(double t)
{
	return (1.0/(4.0*M_PI))*(sinh(t)-t);
}

//Dimensionless pressure
double p(double t)
{
	return (1.0/(12.0*M_PI))*(sinh(t)-8.0*sinh(0.5*t)+3.0*t);
}

//Solves epsilon_bi with the bisection method
double bisection(double a, double b, double epsilon, double lhs ,double (*func)(double,double))
{	
	double left = a;
	double right = b;
	
	double mid;
	
	int i = 0;
	
	while ((right-left) > epsilon)
	{	
		mid = (left+right)*0.5;
		
		double midVal = func(mid,lhs);
		double leftVal = func(left,lhs);
		double rightVal = func(right,lhs);
		
		if(rightVal == 0)
			return right;
		else if (leftVal == 0)
			return left;
		else if (midVal == 0)
			return mid;
		
		if (getSign(midVal) != getSign(rightVal))
		{
			left = mid;
		}
		else if (getSign(midVal) != getSign(leftVal))
		{
			right = mid;
		}
		i++;
	}
	
	return (left+right)/(2.0);
}

//Returns the derivatives in dydr
void f(double r, const double y[], double dydr[])
{				
	dydr[0] = r*r*0.5*DEG*(sinh(y[1])-y[1]);
	
	//For easier debugging, the factors have been split
	double a =-(4.0/(r*(r-2.0*y[0])));
	
	double b = (sinh(y[1])-(2.0*sinh(0.5*y[1])));
	if (y[1] <= 1E-3)
	{
		//~ cout << "b before division: " << b << endl;
		//~ cout << "0.03125*y[1]*y[1]*y[1]*y[1] " << 0.03125*y[1]*y[1]*y[1]*y[1] << endl;
		b/=(0.03125*y[1]*y[1]*y[1]*y[1]);
		//~ cout << "b after division: " << b << endl;
	}
	else
		b/=(cosh(y[1])-4.0*cosh(0.5*y[1])+3.0);
		
	//double b =(sinh(y[1])-(2.0*sinh(0.5*y[1])))/(cosh(y[1])-4.0*cosh(0.5*y[1])+3.0);
	double c = ((1.0/3.0)*r*r*r*DEG*0.5*(sinh(y[1])-8.0*sinh(0.5*y[1])+3.0*y[1])+y[0]);
	
	double test = 1+0.5*y[1]*y[1];
	
	//~ cout << test-1 << endl;
	//~ cout << cosh(y[1])-1 << endl;
	//~ cout << 4*cosh(y[1]*0.5)-4 << endl;
	//~ cout << "y[1] " << y[1] << endl;
	//~ cout << "cosh terme " << (cosh(y[1])-4.0*cosh(0.5*y[1])+3.0) << endl;
	//~ cout << "a " << a << endl;
	//~ cout << "b " << b << endl;
	//~ cout << "c " << c << endl;
	
	if (r == 0.0)
		dydr[1] = 0.0;
	else
		dydr[1] = a*b*c;
}

//Runge-Kutta fourth order integrator
double rk4(double r, double h, const double y_arg[], double y_sol[], void (*f)(double r, const double[], double dydr[]))
{
	double dydr[4][2];
	double tmp[2];
	
	f(r,y_arg,dydr[0]);
	
	//~ cout << "dydr[0][0] " << dydr[0][0] << endl;
	//~ cout << "dydr[0][1] " << dydr[0][1] << endl;
	//~ cin.get();
	
	for(int i=0;i<2;i++)
		tmp[i] = y_arg[i] + 0.5*h*dydr[0][i];
		
	f(r+0.5*h,tmp,dydr[1]);
	
	//~ cout << "dydr[1][0] " << dydr[1][0] << endl;
	//~ cout << "dydr[1][1] " << dydr[1][1] << endl;
	//~ cin.get();
	
	for(int i=0;i<2;i++)
		tmp[i] = y_arg[i] + 0.5*h*dydr[1][i];
		
	f(r+0.5*h,tmp,dydr[2]);
	
	//~ cout << "dydr[2][0] " << dydr[2][0] << endl;
	//~ cout << "dydr[2][1] " << dydr[2][1] << endl;
	//~ cin.get();
	
	for(int i=0;i<2;i++)
		tmp[i] = y_arg[i] + h*dydr[2][i];
		
	f(r+h,tmp,dydr[3]);
	
	//~ cout << "dydr[3][0] " << dydr[3][0] << endl;
	//~ cout << "dydr[3][1] " << dydr[3][1] << endl;
	//~ cin.get();
	
	for(int i=0;i<2;i++)
		y_sol[i] = y_arg[i] + h*(1.0/6.0)*(dydr[0][i]+2.0*(dydr[1][i]+dydr[2][i])+dydr[3][i]);
}	

//Integrates the OV-Equation with rk4 and a simple step-adaption		
int integrate(double RhoCore, double h, ofstream& massradius, ofstream& gacc, ofstream& massdefect)
{
	double y[3][2];
	double r = 0.0;
	double step = h;
	int i = 0;
	
	char filename[64];
	
	//Get dimensionless pressure (t) with given central energy density in g/cm^3
	double tzero = bisection(0,50,BISECTION_ERR,RhoCore,&edensity_bi);
	
	y[0][0] = 0.0; // dimensionless mass (u)
	y[0][1] = tzero; // dimensionless pressure (t)
	
//	cout << "T zero:" << tzero << endl;
	
	double tmpy[2];
	
	//desired precision for step adaption
	double const eps = 1E-7;
	
	while(1)
	{
		// Integrate one step and compare with two steps with half the step size
		// Bisect the step size until error <= eps
		while(1)
		{
			rk4(r,step,y[0],y[1],f);
			
			tmpy[0] = y[1][0];
			tmpy[1] = y[1][1];
					
			rk4(r,step*0.5,y[0],y[1],f);
			rk4(r+step*0.5,step*0.5,y[1],y[2],f);
			
			//~ cout << "y[2][1]:" << y[2][1] << endl;
			//~ cout << "tmpy[1]:" << tmpy[1] << endl;
			//~ cout << "deviation:" << fabs(y[2][1] - tmpy[1]) << endl;
			//~ 
			//~ cout << "deviation q:" << y[2][1]/tmpy[1] << endl;
			//~ if (isnan(tmpy[1]) || isinf(tmpy[1]) || tmpy[1]<0)
			//~ {
				//~ cout << "tmpy nan, inf or neg" << endl;
				
				//~ cin.get();
			//~ }
			
			if((fabs(y[2][1] - tmpy[1]) <= eps) && !(isnan(tmpy[1]) || isinf(tmpy[1])) && (tmpy[0]>0 && tmpy[1]>0))
			{
				y[0][0] = tmpy[0];
				y[0][1] = tmpy[1];
				break;
			}
			else
			{
				//~ cout << "step size:" << step << endl;
				//~ cin.get();
				
				step*=0.5;
			}
		}
		

		
		// This should not happen. It indicates that the step size is too large near t~0
		if (isnan(y[0][1]) || y[0][1] < 0.0 )
		{
			cout << y[0][1]/tzero << " mass: " << y[0][0] << endl;
			cout << "t is NaN or negative" << endl;
			cin.get();
			return -2;
		}
		// Integration has to be stopped since array reached maximum size. This could be solved more elegantly,
		// but the integration does usually not need more than a few thousand steps
	/*	else if (i >= MAX_STEPS)
		{
			cout << "Surface could not be reached within the maximum steps" << endl;
			return -1;
		}*/
		// Surface reached if t/tzero small
		else if (y[0][1]/tzero < T_EPS)
		{
			double R = r;
			double M = y[0][0];
			
			char str[64];
			
			cout << "Surface reached, steps: " << i << " t: " << y[0][1]/tzero << " M: " << (M*DLtoKG)/M_SOL << " M_SOL R: " << R*DLtoM << " m" << endl;
			
			sprintf(str,"%.2E",RhoCore);
			massradius  << R*DLtoCM << "\t" << (M*DLtoKG/M_SOL) <<  "\t" << str << endl;
			massdefect << RhoCore << "\t" << (2.0*M*DLtoKG*G)/(C*C*R*DLtoM) << endl;
			// gravitational acceleration at surface in g_earth
			// gaccfactor << ParticleMass << "\t" << ((G*M*DLtoKG/(R*R*DLtoM*DLtoM))/9.81)/pow(RhoCore/RHO_ZERO,POLY_EXP*0.5) << endl;
			gacc << RhoCore << "\t" << (G*M*DLtoKG/(R*R*DLtoM*DLtoM))/G_EARTH << endl;
			break;
		}
		r += step;

		//~ cout << "radius: " << r << " ratio: " << y[0][1]/tzero << " mass: " << y[0][0] << " step: " << step;
		//~ cout << "\r";
		//~ cin.get();
		step = h;
		i++;	
	}
	return 0;
}

#define START_RHO 	1E23
#define END_RHO 	1E30

int main(int argc, char **argv)
{
	/*
	double pmasseV = PARTICLE_MASS;
	double pmass = PARTICLE_MASS_G;
	double density = DENSITY_TO_DL;
	*/
	/*
	cout << "pmass eV " << pmasseV << endl;
	cout << "pmass g" << pmass << endl;
	cout << "density factor" << density << endl;
	*/
	char filename[64];
	
	/*sprintf(filename, "gravityfactor%.2E.dat", PARTICLE_MASS);
	ofstream gravityfactor(filename);*/
		
	/*
	gravityfactor << "# particle mass " << PARTICLE_MASS << " eV" << endl;
	gravityfactor << "# rho from " << START_RHO << " to " << END_RHO << " g/cm^3" << endl;
	*/
	
	// initial step size for runge-kutta
	double h = 1E-2;
	double rhoinc = 1.3;
	bool write_file;
	
	// double masses[5] = {5E9,1E12,3E12,5E12,10E12};
	double masses[3] = {5E9,1E12,3E12};
	double rhorange[5][2] = {{1E10,1E31},{1E10,1E31},{1E10,1E33},{1E10,1E34},{1E10,1E35}};
	
	for (int a =0;a<3;a++)
	{
		ParticleMass = masses[a];
		cout << "Particle Mass: " << ParticleMass << endl;
		recalcConst();
		
		sprintf(filename, "massradius%.2E.dat", ParticleMass);
		ofstream fmassradius(filename);
		sprintf(filename, "acceleration%.2E.dat", ParticleMass);
		ofstream gravity(filename);
		sprintf(filename, "massdefect%.2E.dat", ParticleMass);
		ofstream massdefect(filename);
		
		fmassradius << "# particle mass " << ParticleMass << " eV" << endl;
		fmassradius << "# rho from " << rhorange[a][0] << " to " << rhorange[a][1] << " g/cm^3" << endl;
		fmassradius << "# degeneracy " << DEG << endl;
		
		gravity << "# particle mass " << ParticleMass << " eV" << endl;
		gravity << "# rho from " << rhorange[a][0] << " to " << rhorange[a][1] << " g/cm^3" << endl;
		
		massdefect << "# particle mass " << ParticleMass << " eV" << endl;
		massdefect << "# rho from " << rhorange[a][0] << " to " << rhorange[a][1] << " g/cm^3" << endl;
		
		for(double rho=rhorange[a][0];rho<=rhorange[a][1];rho*=rhoinc)
		{			
			integrate(rho,h,fmassradius,gravity,massdefect);
		}
					
		fmassradius.close();
		gravity.close();
		massdefect.close();
	}
	
	return 0;
}
