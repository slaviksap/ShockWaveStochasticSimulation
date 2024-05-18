#pragma once
#include"Vec3.h"
#include<vector>
#include"SpatHomSolver.h"
#include"Histogram.h"
using std::vector;

struct ShockWaveParams
{
	double E1 = -1, E2 = -1;
	double u1 = -1, u2 = -1;
	double rho1 = -1, rho2 = -1;
	double D = 0;
};

struct ShockWaveSolver
{
	int k;						//число шагов по времени
	int M;						//число испытаний
	double tmax;				//интервал по времени
	double dt;
	vector<SpatHomSolver> cells;
	vector<double> cellsBounds = {-24, -22,-20,-18,-16,-14,-12, -10,-8,-6, -5,-4,-3,-2,-1.5,-1.2,-0.9,-0.7,-0.5,-0.3,-0.2,-0.1,0,
		0.1,0.2,0.3,0.5,0.7,0.9,1.2,1.5,2,3,4,5,6,8,10,12,14,16,18,20,22,24 };
	ShockWaveParams params;
	Histogram initHistoX;
	Histogram lastHistoX;
	Histogram initHistoV;
	Histogram lastHistoV;
	double visuaDt;
	vector<Histogram> rhoX;
	vector<Histogram> vX;
	ShockWaveSolver(double t, int k, int M, int histSect, double histMin, double histMax, double visua_dt);
	void CheckParams();
	void InitialDistribution();
	void Solve();
	void UpdateCells();
	void ReturnInDomain();
	void Refill(vector<Particle>& movingParticles);
	double CalcD();
	SpatHomSolver& front() { return cells.front(); }
	SpatHomSolver& back() { return cells.back(); }
};