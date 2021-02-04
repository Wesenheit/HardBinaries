#include <iostream>
#include <fstream>
#include <math.h>
#include <array>
using namespace std;


array<double,4> add(array<double,4> w1,array<double,4> w2) //adding arrays
{
array<double,4> wyn0;
for (int i=0; i<4; i++)
	{
	wyn0[i]=w1[i]+w2[i];
	}
return wyn0;
}

array<double,4> mul(array<double,4> w, double c) //multiplicating array by scalar
{
array<double,4> wyn1;
for (int i=0; i<4; i++)
	{
	wyn1[i]=w[i]*c;
	}
return wyn1;
}

array<double,2> force (array<double,4> wekt,double k) //calculating force based on position
{
double r=sqrt(wekt[0]*wekt[0]+wekt[1]*wekt[1]);
double fx=-k*wekt[0]/(r*r*r);
double fy=-k*wekt[1]/(r*r*r);
array<double,2> wyn2;
wyn2[0]=fx;
wyn2[1]=fy;
return wyn2;
}

array<double,4> diff(array<double,4> wekt,double k) //calculating derivation of vecotor
{
array<double,2> temp=force(wekt,k);
array<double,4> wyn3;
wyn3[0]=wekt[2];	//vx
wyn3[1]=wekt[3];	//vy
wyn3[2]=temp[0];	//fx
wyn3[3]=temp[1];	//fy
return wyn3;
}

array<double,4> rungekutta(array<double,4> wekt,double k, double dt) //runge kutta methode, 4th order
{
array<double,4> wektp;
array<double,4> k1=diff(wekt,k);
array<double,4> k2=diff(add(wekt,mul(k1,dt/2)),k);
array<double,4> k3=diff(add(wekt,mul(k2,dt/2)),k);
array<double,4> k4=diff(add(wekt,mul(k3,dt)),k);
for (int i=0;i<4;i++)
	{
	wektp[i]=dt/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i])+wekt[i];
	}
return wektp;
}

int main()
{
double k=1000;       //force constant
double dt=0.0001;	//time step
array<double,4> wekt={10.,0.,0.,12.}; //inital position
ofstream fout;
fout.open("file.txt");	//output file
double T=100;		//Time of simulation
double t=0;		//initial time
while (t<T)		//start of simulation
	{
	wekt=rungekutta(wekt,k,dt);
	fout<<wekt[0]<<" "<<wekt[1]<<" "<<wekt[2]<<" "<<wekt[3]<<endl;
	t=t+dt;
	}
fout.close();
return 0;
}
