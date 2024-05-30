#include "SpatHomSolver.h"
#include "Randomizer.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include "Histogram.h"
#include"ShockWaveParams.h"
using namespace std;
constexpr double C_2_3 = 2.0 / 3.0;
SpatHomSolver::SpatHomSolver(int N, double dt, double Kn, double left, double right)
	: N(N), dt(dt), Kn(Kn), prev_layer(N), layer(N), n(1.0), left(left), right(right), initN(N)
{
	sdt = sqrt(dt);
	dw = 3.14159265358979 / 4.0;
}

void SpatHomSolver::InitialDistribution(double m, double u, double sigma)
{
	N = initN;
	mInit = m;
	layer.assign(initN, Particle());
	normal_distribution<> normx(u, sigma);
	normal_distribution<> norm(0, sigma);
	uniform_real_distribution<> uni(left, right);
	for (auto& particle : layer)
	{
		particle.v.x = normx(Randomizer::gen());
		particle.v.y = norm(Randomizer::gen());
		particle.v.z = norm(Randomizer::gen());

		particle.mass = m;
	}
	DistribUniformly();
}
void SpatHomSolver::DistribUniformly()
{
	double dx = (right - left) / layer.size();
	double cur_x = left + dx / 2;
	for (int i = 0; i < layer.size(); ++i)
	{
		auto& particle = layer[i];
		particle.x = cur_x;
		cur_x += dx;
	}
}
DataInfo SpatHomSolver::CalculateStep()
{
	N = layer.size();
	//random_shuffle(layer.begin(), layer.end());
	DataInfo info;
	Histogram3D F(100, -50.0, 50.0, 0);
	F.Build(layer);
	double M = 0;
	for (auto& particle : layer)
		M += particle.mass* mAmb;
	assert(M >= N, "wrong M value");
	double T = Temperature();
	Vec3 V = MacroVelocity();
	double rho = Rho();
	for (int i = 0; i < N; ++i)
	{
		for (int j = i + 1; j < N; ++j)
		{
			auto w = Randomizer::sampleOmega();
			double prod = dot(w, layer[i].v - layer[j].v);
			double l = abs(prod) * dt;
			double lambda = l / M / Kn / 2 * dw * rho * F.Get(layer[j].v) * mAmb;

			if (Randomizer::sampleUni() < lambda)
			{
				Vec3 dv = w * prod;
				layer[i].v = layer[i].v - dv;
				layer[j].v = layer[j].v + dv;
				++info.colCount;
			}
		}

	}

	return info;
}
DataInfo SpatHomSolver::CalculateLangevinStep()
{
	N = layer.size();
	double avgV = MacroVelocity().x;
	double T = Temperature();
	const double a = sqrt(T);
	const double sigma = sqrt(2 * T * a);
	DataInfo info;
	for (int i = 0; i < N; ++i)
	{
		Vec3 v = layer[i].v;
		Vec3 dv = -(v - Vec3(avgV,0,0)) * a * dt + Randomizer::sampleNormVec3() * sdt * sigma;
		layer[i].v += dv;
	}

	return info;
}
vector<Particle>::iterator SpatHomSolver::PrepareToChangeCell()
{
	auto end_iter = layer.end();
	for (auto iter = layer.begin(); iter != end_iter;)
	{
		if (iter->x >= left && iter->x < right)
		{
			++iter;
		}
		else
		{
			--end_iter;
			swap(*iter, *end_iter);
		}
	}
	return end_iter;
}
Vec3 SpatHomSolver::MacroVelocity()
{
	Vec3 avgVel = Vec3(0,0,0);
	double massCount = 0;
	for (Particle& particle : layer)
	{
		avgVel += particle.v * particle.mass;
		massCount += particle.mass;
	}
	return avgVel / massCount;
}
double SpatHomSolver::Temperature()
{
	double avgVel = 0;
	double avgVelSquared = 0;
	double massCount = 0;
	for (const Particle& particle : layer)
	{
		const double x = particle.x;
		avgVel += particle.v.x * particle.mass;
		avgVelSquared += particle.v.x * particle.v.x * particle.mass;
		massCount += particle.mass;
	}
	avgVel /= massCount;
	avgVelSquared /= massCount;
	return avgVelSquared - avgVel*avgVel;
}

double SpatHomSolver::Rho()
{
	double mass = 0;
	for (const Particle& particle : layer)
	{
		mass += particle.mass;
	}
	return mass / (right - left);
}

double SmoothFunction(double x, double left, double right, double leftValue, double rightValue)
{
	constexpr double C23 = 2.0 / 3.0;
	constexpr double C43 = 4.0 / 3.0;
	const double M = 2 / (right - left);
	const double L = left + (right - left) / 2;
	const double scale = leftValue - rightValue;
	return (((-x * M + L) - pow(-x * M + L, 3) / 3 + C23) / C43) * scale + rightValue;
}

void SpatHomSolver::InitialSmoothedDistribution(double rho1, double u1, double T1, double rho2, double u2, double T2, int n, double l, double r)
{
	N = n;
	layer.assign(n, Particle());
	double dx = (right - left) / n;
	double x = left + dx / 2;
	for (int i = 0; i < n; ++i)
	{
		double rho = SmoothFunction(x, l, r, rho1, rho2);
		layer[i].mass = rho * dx;
		layer[i].x = x;
		double u = SmoothFunction(x, l, r, u1, u2);
		double T = SmoothFunction(x, l, r, T1, T2);
		normal_distribution<> norm(0, ShockWaveParams::sigma(T));
		layer[i].v.x = norm(Randomizer::gen()) + u;
		layer[i].v.y = norm(Randomizer::gen());
		layer[i].v.z = norm(Randomizer::gen());
		x += dx;
	}
}