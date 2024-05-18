#include "ShockWaveSolver.h"
#include <fstream>
#include <iostream>
#include "Randomizer.h"
#include "string"
#include <omp.h>
using namespace std;

ShockWaveSolver::ShockWaveSolver(double t, int k, int M, int histSect, double histMin, double histMax, double visuaDt) :tmax(t),
k(k), M(M), initHistoX(histSect, histMin, histMax), lastHistoX(histSect, histMin, histMax),
initHistoV(histSect, histMin, histMax), lastHistoV(histSect, histMin, histMax), visuaDt(visuaDt)
{
	double dt = t / k;
	this->dt = dt;
	
	for (int i = 0; i < cellsBounds.size() - 1; ++i)
	{
		const double left = cellsBounds[i];
		const double right = cellsBounds[i + 1];
		const double Kn = 1;
		int N = 300;
		if (left >= -3 && right <= 3)
			N = 300;
		cells.emplace_back(SpatHomSolver(N, dt, Kn, left, right));
	}
	Randomizer::InitOmegas(20, 20);

	int visuaTimeSteps = int(t / visuaDt);
	rhoX.assign(visuaTimeSteps, Histogram(histSect, histMin, histMax));
	//omp_set_num_threads(8);
}

void ShockWaveSolver::CheckParams()
{
	constexpr double gamma = 1.6;
	constexpr double gam_inv = 1.0 / (gamma - 1);
	const double r0 = params.rho2;
	const double r1 = params.rho1;
	//const double u = params.u1;
	//const double D = r1 * u / (r1 - r0);
	//const double E1 = -r0 * D * (D * u * gam_inv + 0.5 * u * u) / (r0 * D - r1 * D - (gamma - 1) * r1 * u);
	//const double E0 = r1 / r0 * E1 - D * u * gam_inv;

	const double D = r1 * params.u1 / (r1 - r0);
	const double u1 = -(D - params.u1);
	const double u0 = -D;
	const double E0 = (u0 * u0 / (2 * gamma) + u1 * u1 * gam_inv - r0 * u0 * u0 * gam_inv / r1) / (r0 / r1 - 1);
	const double E1 = E0 + (u0 * u0 - u1 * u1) / (2 * gamma);
	cout << "->->->->->D = " << D << endl;
	cout << "->->->->->E1 = " << E1 << endl;
	cout << "->->->->->E0 = " << E0 << endl;
	params.D = D;
	params.E1 = E1;
	params.E2 = E0;
}

void ShockWaveSolver::InitialDistribution()
{
	double left = -4, right = 4;
	for (int i = 0; i < cells.size() / 2; ++i)
	{
		if (cells[i].left >= left)
		{
			cells[i].InitialSmoothedDistribution(params.rho1, params.u1, params.E1, params.rho2, params.u2, params.E2, cells[i].initN, left, right);
		}
		else
		{
			double m1 = (cells[i].right - cells[i].left) * params.rho1 / cells[i].initN;
			cells[i].InitialDistribution(m1, params.u1, params.E1);
		}
	}
	for (int i = cells.size() / 2; i < cells.size(); ++i)
	{
		if (cells[i].right <= right)
		{
			cells[i].InitialSmoothedDistribution(params.rho1, params.u1, params.E1, params.rho2, params.u2, params.E2, cells[i].initN, left, right);
		}
		else
		{
			double m2 = (cells[i].right - cells[i].left) * params.rho2 / cells[i].initN;
			cells[i].InitialDistribution(m2, params.u2, params.E2);
		}
	}
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
		int curVisua = 0;
		double t = 0;
		DataInfo info;
		double StartEnergy = Energy(cells);
		Vec3 StartImpuls = Impulse(cells);
		double StartMass = Mass(cells);
		for (int timeStep = 1; timeStep <= k; ++timeStep)
		{
			//#pragma omp parallel for
			for (int i = 0; i < cells.size(); ++i)
			{
				if (cells[i].left >= -4 && cells[i].right <= 4)
					cells[i].CalculateStep();
				else
					cells[i].CalculateLangevinStep(cells[i].left < -4 ? params.E1 : params.E2);
			}
			UpdateCells();
			t = timeStep * (tmax / k);
			if (t > (curVisua + 1) * visuaDt)
			{
				rhoX[curVisua].ReplenishX(cells);
				++curVisua;
			}
		}
		info.Energy = StartEnergy / Energy(cells);
		info.p.x = StartImpuls.x / Impulse(cells).x;
		info.p.y = StartImpuls.y / Impulse(cells).y;
		info.p.z = StartImpuls.z / Impulse(cells).z;
		info.mass = StartMass / Mass(cells);
		allInfo.Add(info);
		//lastHistoX.ReplenishX(cells);
		std::cout << m << "\n";
	}
	allInfo.Average(M);
	allInfo.Print();
	//initHistoX.Normalize();
	ofstream file("init.txt");
	for (auto& cell : initHistoX.hist)
	{
		cell /= M;
	}
	initHistoX.Write(file);
	file.close();
	for (int i = 0; i < rhoX.size(); ++i)
	{
		//rhoX[i].Normalize();
		for (auto& cell : rhoX[i].hist)
		{
			cell /= M;
		}
		ofstream file("rhoX_" + to_string((i + 1) * visuaDt) + ".txt");
		rhoX[i].Write(file);
		file.close();
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
	//Refill(movingParticles);
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
	for (auto& cell : cells)
	{
		for (auto& particle : cell.layer)
		{
			if (particle.x < cellsBounds.front())
			{
				//particle.x += abs(particle.x - cellsBounds.front()) * 2;
				//particle.v = particle.v * (-1.0);
				cells.front().PlayParticle(particle, params.u1, params.E1);
			}
			if (particle.x > cellsBounds.back())
			{
				//particle.x -= abs(particle.x - cellsBounds.back()) * 2;
				//particle.v = particle.v * (-1.0);
				cells.back().PlayParticle(particle, params.u2, params.E2);

			}
		}
	}

	double m1 = (front().right - front().left) * params.rho1 / front().layer.size();
	//double m12 = (cells[1].right - cells[1].left) * params.rho1 / cells[1].layer.size();
	double m2 = (back().right - back().left) * params.rho2 / back().layer.size();
	//double m22 = (cells[cells.size() - 2].right - cells[cells.size() - 2].left) * params.rho2 / cells[cells.size() - 2].layer.size();
	cells.front().InitialDistribution(m1, params.u1, params.E1, front().layer.size());
	//cells[1].InitialDistribution(m12, params.u1, params.E1, cells[1].layer.size());
	cells.back().InitialDistribution(m2, params.u2, params.E2, back().layer.size());
	//cells[cells.size() - 2].InitialDistribution(m22, params.u2, params.E2, cells[cells.size() - 2].layer.size());
	//cells.front().DistribUniformly();
	//cells.back().DistribUniformly();
}
double ShockWaveSolver::CalcD()
{
	const double gamma = 1.6;
	const double p1 = (gamma - 1) * params.rho1 * params.E1;
	const double p2 = (gamma - 1) * params.rho2 * params.E2;
	return params.u1 * sqrt((p2 - p1) / (params.u1 - params.u2));
}
void ShockWaveSolver::Refill(vector<Particle>& movingParticles)
{
	//for (auto& particle : movingParticles)
	//{
	//	double x = particle.x;
	//	if (x < cellsBounds.front() || x > cellsBounds.back())
	//		continue;
	//	else
	//	{
	//		auto comp = [=](SpatHomSolver& cell) {return x > cell.left && x <= cell.right; };
	//		std::find_if(cells.begin(), cells.end(), comp)->layer.push_back(particle);
	//	}
	//}
	//double dt = tmax / k;
	//double dx = (params.D - params.u1) * dt / 2;
	//double V1 = front().right - front().left;
	//int Ncreate1 = cells.front().initN * dx / V1;
	//double m1 = params.rho1 * dx / Ncreate1;
	//for (int i = 0; i < Ncreate1; ++i)
	//{
	//	Particle p;
	//	normal_distribution<> norm(0, sqrt(params.E1));
	//	uniform_real_distribution<> uni(front().left, front().left + dx);
	//	//double guess = -1;
	//	//while (guess < 0)
	//	//{
	//	//	guess = norm(Randomizer::gen());
	//	//}
	//	//p.v.x = guess + params.u1;
	//	p.v.x = norm(Randomizer::gen()) + params.u1;
	//	p.v.y = norm(Randomizer::gen());
	//	p.v.z = norm(Randomizer::gen());

	//	p.x = uni(Randomizer::gen());
	//	p.mass = m1;
	//	front().layer.push_back(p);
	//}
	//dx = params.D*dt;
	//double V2 = back().right - back().left;
	//int Ncreate2 = back().initN * dx / V2;
	//double m2 = params.rho2 * dx / Ncreate2;
	//for (int i = 0; i < Ncreate2; ++i)
	//{
	//	Particle p;
	//	normal_distribution<> norm(0, sqrt(params.E2));
	//	uniform_real_distribution<> uni(back().right - dx, back().right);
	//	double guess = 1;
	//	while (guess > 0)
	//	{
	//		guess = norm(Randomizer::gen());
	//	}
	//	p.v.x = guess + params.u2;
	//	p.v.y = norm(Randomizer::gen());
	//	p.v.z = norm(Randomizer::gen());

	//	p.x = uni(Randomizer::gen());
	//	p.mass = m2;
	//	back().layer.push_back(p);
	//}
	for (auto& particle : movingParticles)
	{
		double x = particle.x;
		if (x < cellsBounds.front())
		{
			continue;
			particle.x = cells.back().right + (particle.v.x - params.D) * dt;
			/*particle.v = particle.v * (-1);*/
			normal_distribution<> norm(0, sqrt(params.E2));
			particle.v.x = norm(Randomizer::gen()) + params.u2;
			particle.v.y = norm(Randomizer::gen());
			particle.v.z = norm(Randomizer::gen());
			back().layer.push_back(particle);
		}
		else if (x > cellsBounds.back())
		{
			particle.x -= (particle.v.x - params.D) * dt;
			//particle.v = particle.v * (-1);
			normal_distribution<> norm(0, sqrt(params.E2));
			particle.v.x = norm(Randomizer::gen()) + params.u2;
			particle.v.y = norm(Randomizer::gen());
			particle.v.z = norm(Randomizer::gen());
			back().layer.push_back(particle);
		}
		else
		{
			auto comp = [=](SpatHomSolver& cell) {return x > cell.left && x <= cell.right; };
			std::find_if(cells.begin(), cells.end(), comp)->layer.push_back(particle);
		}
	}
}