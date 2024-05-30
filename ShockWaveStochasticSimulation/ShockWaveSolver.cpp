#include "ShockWaveSolver.h"
#include <fstream>
#include <iostream>
#include "Randomizer.h"
#include "string"
#include <algorithm>
#include <omp.h>
using namespace std;

ShockWaveSolver::ShockWaveSolver(double t, int k, int M, int N, int histSect, double histMin, double histMax, double visuaDt) :tmax(t),
k(k), M(M), Ninit(N), initHistoX(histSect, histMin, histMax), initHistoV(histSect, histMin, histMax), initHistoT(histSect, histMin, histMax), visuaDt(visuaDt)
{
	double dt = t / k;
	this->dt = dt;
	//cellsBounds = {-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,10,11,12,13,14,15};
	cellsBounds = {};
	for (double p = -15; p <= 15; p += 0.1)
		cellsBounds.push_back(p);
	/*cellsBounds = {-14, -12, -10, -8, -6, -5, -4, -3, -2, -1.5, -1.2, -0.9, -0.7, -0.5, -0.3, -0.2, -0.15, -0.1, -0.05, -0.03, -0.01, 0, 0.01, 0.03, 0.05,
		0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 0.9, 1.2, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 14};*/
	for (int i = 0; i < cellsBounds.size() - 1; ++i)
	{
		const double left = cellsBounds[i];
		const double right = cellsBounds[i + 1];
		const double Kn = 1;
		if (Ninit <= 0)
			throw("Wrong initial N!\n\a");
		cells.emplace_back(SpatHomSolver(Ninit, dt, Kn, left, right));
	}
	Randomizer::InitOmegas(20, 20);

	int visuaTimeSteps = int(t / visuaDt);
	rhoX.assign(visuaTimeSteps, Histogram(histSect, histMin, histMax));
	vX.assign(visuaTimeSteps, Histogram(histSect, histMin, histMax));
	TX.assign(visuaTimeSteps, Histogram(histSect, histMin, histMax));
	distr.assign(visuaTimeSteps, Histogram(300, -5, 20));
	//omp_set_num_threads(8);
}
double ShockWaveParams::sigma(double T)
{
	return sqrt(T);
}
void ShockWaveParams::CalcParams()
{
	const double Ma2 = Ma * Ma;
	rho2 = 1;
	rho1 = rho2 * (gamma + 1) * Ma2 / (2 + (gamma - 1) * Ma2);
	D = Ma;
	u1 = D - Ma * (2 + (gamma - 1) * Ma2) / ((gamma + 1) * Ma2);
	u2 = 0;
	p2 = 1 / gamma;
	p1 = p2 * (2 * gamma * Ma2 - gamma + 1) / (gamma + 1);
	T1 = gamma * p1 / rho1;
	T2 = gamma * p2 / rho2;
	sigma1 = sigma(T1);
	sigma2 = sigma(T2);
	pprintf("rho1 = %g\tu1 = %g\tp1 = %g\tT1 = %g\n", rho1, u1, p1, T1);
	pprintf("rho2 = %g\tu2 = %g\tp2 = %g\tT2 = %g\n", rho2, u2, p2, T2);
	pprintf("D = %g\n", D);
}

void ShockWaveSolver::InitialDistribution()
{
	for (int i = 0; i < cells.size() / 2; ++i)
	{
		if (useSmooth && cells[i].left >= leftSmooth)
		{
			cells[i].InitialSmoothedDistribution(params.rho1, params.u1, params.T1, params.rho2, params.u2, params.T2, cells[i].initN, leftSmooth, rightSmooth);
		}
		else
		{
			double m1 = (cells[i].right - cells[i].left) * params.rho1 / cells[i].initN;
			cells[i].InitialDistribution(m1, params.u1, params.sigma1);
		}
	}
	for (int i = cells.size() / 2; i < cells.size(); ++i)
	{
		if (useSmooth && cells[i].right <= rightSmooth)
		{
			cells[i].InitialSmoothedDistribution(params.rho1, params.u1, params.T1, params.rho2, params.u2, params.T2, cells[i].initN, leftSmooth, rightSmooth);
		}
		else
		{
			double m2 = (cells[i].right - cells[i].left)* params.rho2 / cells[i].initN;
			cells[i].InitialDistribution(m2, params.u2, params.sigma2);
		}
	}
	double minMass = 1000000;
	for (auto& cell : cells)
		if (cell.mInit < minMass)
			minMass = cell.mInit;
	for (auto& cell : cells)
		cell.mAmb = 1.0 / minMass;
}
double Energy(vector<Particle>& v)
{
	double E = 0;
	for (const auto& particle : v)
	{
		const double U = magnitide(particle.v);
		E += 0.5 * U * U * particle.mass;
	}
	return E;
}
double Energy(vector<SpatHomSolver>& cells)
{
	double E = 0;
	for (auto& cell : cells)
	{
		E += Energy(cell.layer);
	}
	return E;
}
Vec3 Impulse(vector<Particle>& v)
{
	Vec3 p = { 0,0,0 };
	for (auto& particle : v)
	{
		p = p + particle.v * particle.mass;
	}
	return p;
}
Vec3 Impulse(vector<SpatHomSolver>& cells)
{
	Vec3 p = { 0,0,0 };
	for (auto& cell : cells)
	{
		p = p + Impulse(cell.layer);
	}
	return p;
}
double Mass(vector<Particle>& v)
{
	double m = 0;
	for (const auto& particle : v)
	{
		m += particle.mass;
	}
	return m;
}
double Mass(vector<SpatHomSolver>& cells)
{
	double m = 0;
	for (auto& cell : cells)
	{
		m += Mass(cell.layer);
	}
	return m;
}
void ShockWaveSolver::Solve()
{
	DataInfo allInfo;
	for (int m = 1; m <= M; ++m)
	{
		InitialDistribution();
		initHistoX.ReplenishX(cells);
		initHistoV.ReplenishV(cells);
		initHistoT.ReplenishT(cells);
		int curVisua = 0;
		double t = 0;
		DataInfo info;
		double StartEnergy = Energy(cells);
		Vec3 StartImpuls = Impulse(cells);
		double StartMass = Mass(cells);
		for (int timeStep = 1; timeStep <= k; ++timeStep)
		{
			#pragma omp parallel for schedule(dynamic)
			for (int i = 0; i < cells.size(); ++i)
			{
				if (!useDiffusion || cells[i].left >= leftDiffusion && cells[i].right <= rightDiffusion)
					cells[i].CalculateStep();
				else
					cells[i].CalculateLangevinStep();
			}
			UpdateCells();
			t = timeStep * (tmax / k);
			if (t > (curVisua + 1) * visuaDt)
			{
				rhoX[curVisua].ReplenishX(cells);
				vX[curVisua].ReplenishV(cells);
				TX[curVisua].ReplenishT(cells);
				distr[curVisua].ReplenishDistr(cells[25]);
				++curVisua;
			}
		}
		info.Energy = StartEnergy / Energy(cells);
		info.p.x = StartImpuls.x / Impulse(cells).x;
		info.p.y = StartImpuls.y / Impulse(cells).y;
		info.p.z = StartImpuls.z / Impulse(cells).z;
		info.mass = StartMass / Mass(cells);
		allInfo.Add(info);
		pprintf("%d\n", m);
	}
	allInfo.Average(M);
	allInfo.Print();

	initHistoX.Average(M);
	initHistoX.Write("init.txt");
	initHistoV.Average(M);
	initHistoV.Write("initV.txt");
	initHistoT.Average(M);
	initHistoT.Write("initT.txt");
	for (int i = 0; i < rhoX.size(); ++i)
	{
		rhoX[i].Average(M);
		rhoX[i].Write("rhoX_" + to_string((i + 1) * visuaDt) + ".txt");

		vX[i].Average(M);
		vX[i].Write("vX_" + to_string((i + 1) * visuaDt) + ".txt");

		TX[i].Average(M);
		TX[i].Write("TX_" + to_string((i + 1) * visuaDt) + ".txt");

		distr[i].Average(M);
		distr[i].Write("distr_" + to_string((i + 1) * visuaDt) + ".txt");
	}
}

void ShockWaveSolver::UpdateCells()
{
	for (auto& cell : cells)
	{
		for (auto& particle : cell.layer)
		{
			particle.x += (particle.v.x - params.D) * cell.dt;
		}
	}
	vector<Particle> movingParticles;
	for (auto& cell : cells)
	{
		auto forMoves = cell.PrepareToChangeCell();
		movingParticles.insert(movingParticles.end(), forMoves, cell.layer.end());
		cell.layer.erase(forMoves, cell.layer.end());
	}
	for (auto& particle : movingParticles)
	{
		double x = particle.x;
		if (x < cellsBounds.front())
		{
			cells.front().layer.push_back(particle);
		}
		else if (x > cellsBounds.back())
		{
			cells.back().layer.push_back(particle);
		}
		else
		{
			auto comp = [=](SpatHomSolver& cell) {return x > cell.left && x <= cell.right; };
			std::find_if(cells.begin(), cells.end(), comp)->layer.push_back(particle);
		}
	}
	ReturnInDomain();
}
void ShockWaveSolver::ReturnInDomain()
{
	double m1 = (front().right - front().left) * params.rho1 / front().initN;
	front().InitialDistribution(m1, params.u1, params.sigma1);
	double m2 = (back().right - back().left) * params.rho2 / back().initN;
	cells.back().InitialDistribution(m2, params.u2, params.sigma2);
}
void ShockWaveSolver::Refill(vector<Particle>& movingParticles)
{
	for (auto& particle : movingParticles)
	{
		double x = particle.x;
		if (x > cellsBounds.front() && x < cellsBounds.back())
		{
			auto comp = [=](SpatHomSolver& cell) {return x > cell.left && x <= cell.right; };
			std::find_if(cells.begin(), cells.end(), comp)->layer.push_back(particle);
		}
	}
	double l = dt * params.D;
	double dx = l / back().initN;
	int Ncreate = back().initN * l / (back().right - back().left);
	double x = back().right - dx;
	for (int i = 0; i < Ncreate; ++i)
	{
		Particle p;
		normal_distribution<> norm(0, params.sigma2);
		uniform_real_distribution<> uni(back().left, back().right);
		p.x = uni(Randomizer::gen());
		p.v.x = norm(Randomizer::gen()) + params.u2;
		p.v.y = norm(Randomizer::gen());
		p.v.z = norm(Randomizer::gen());
		p.mass = back().mInit;
		back().layer.push_back(p);
	}
}