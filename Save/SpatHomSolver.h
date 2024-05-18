#pragma once
#include"Vec3.h"
#include<iostream>
#include<vector>
#include<list>
#include<math.h>

using std::vector;
using std::list;
using std::cout;
using std::endl;
struct Particle
{
	double x;
	Vec3 v;
	double mass;
};
struct DataInfo
{
	int colCount = 0;
	double lambdaMax = 0;
	double lambdaAvg = 0;
	double lMax = 0;
	double lMin = 1e16;
	double lAvg = 0;
	double Echange = 0;
	double mass = 0;
	double Energy = 0;
	Vec3 p = { 0,0,0 };

	void Add(DataInfo& info)
	{
		mass += info.mass;
		Energy += info.Energy;
		p = p + info.p;
	}

	void Average(int N)
	{
		mass /= N;
		Energy /= N;
		p = p / N;
	}

	void Print()
	{
		cout << "Mass change = " << mass << endl;
		cout << "Energy change = " << Energy << endl;
		cout << "Impuls change = (" << p.x << ", " << p.y << ", " << p.z << endl;
	}
};
struct SpatHomSolver
{
	int N;							//число частиц
	int initN;						//начальное число частиц в €чейке
	double Kn = 1.0;						//„исло  нудсена
	double dt;						//шаг по времени
	double sdt;						//корень из шага по времени
	double n;						//числова€ плотность
	double D = 0.001;				//радиус частицы
	double dw;						//объемный угол
	double left, right;				//граница области
	vector<Particle> prev_layer;	//предыдущий временный слой
	vector<Particle> layer;			//текущий слой по времени
	SpatHomSolver(int N, double dt, double Kn, double left, double right);
	void InitialDistribution(double m, double u, double T);
	void InitialDistribution(double m, double u, double E, int n);
	void InitialSmoothedDistribution(double rho1, double u1, double E1, double rho2, double u2, double E2, int n, double l, double r);
	void DistribUniformly();
	DataInfo CalculateStep();
	DataInfo CalculateLangevinStep(double E);
	vector<Particle>::iterator PrepareToChangeCell();
	void PlayParticle(Particle& p, double u, double E);
	friend double Energy(vector<Particle>& v);
};

