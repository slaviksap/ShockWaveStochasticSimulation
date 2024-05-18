#pragma once
#include"Vec3.h"
#include<vector>
#include"SpatHomSolver.h"
#include"Histogram.h"
#include"ShockWaveParams.h"
using std::vector;

struct ShockWaveSolver
{
	int k;						//число шагов по времени
	int M;						//число испытаний
	double tmax;				//интервал по времени
	double dt;
	int Ninit = -1;
	bool useSmooth;
	double leftSmooth, rightSmooth;
	bool useDiffusion;
	double leftDiffusion, rightDiffusion;
	vector<SpatHomSolver> cells;
	vector<double> cellsBounds = {-24, -22,-20,-18,-16,-14,-12, -10,-8,-6, -5,-4,-3,-2,-1.5,-1.2,-0.9,-0.7,-0.5,-0.3,-0.2,-0.15,-0.1,-0.05,-0.03, -0.01,0, 0.01, 0.03, 0.05,
		0.1,0.15, 0.2,0.3,0.5,0.7,0.9,1.2,1.5,2,3,4,5,6,8,10,12,14,16,18,20,22,24 };
	ShockWaveParams params;
	Histogram initHistoX;
	Histogram initHistoV;
	Histogram initHistoT;
	double visuaDt;
	vector<Histogram> rhoX;
	vector<Histogram> vX;
	vector<Histogram> TX;
	vector<Histogram> distr;
	ShockWaveSolver(double t, int k, int M, int N, int histSect, double histMin, double histMax, double visua_dt);
	void InitialDistribution();
	void Solve();
	void UpdateCells();
	void ReturnInDomain();
	void Refill(vector<Particle>& movingParticles);
	SpatHomSolver& front() { return cells.front(); }
	SpatHomSolver& back() { return cells.back(); }
};