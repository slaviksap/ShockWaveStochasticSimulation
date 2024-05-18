#include "SpatHomSolver.h"
#include "Randomizer.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include "Histogram.h"
using namespace std;
constexpr double C_2_3 = 2.0 / 3.0;
SpatHomSolver::SpatHomSolver(int N, double dt, double Kn, double left, double right)
	: N(N), dt(dt), Kn(Kn), prev_layer(N), layer(N), n(1.0), left(left), right(right), initN(N)
{
	sdt = sqrt(dt);
	dw = 3.14159265358979 / 4.0;
}

void SpatHomSolver::InitialDistribution(double m, double u, double E)
{
	N = initN;
	layer.assign(initN, Particle());
	normal_distribution<> normx(u, sqrt(C_2_3 * E));
	normal_distribution<> norm(0, sqrt(C_2_3 * E));
	uniform_real_distribution<> uni(left, right);
	for (auto& particle : layer)
	{
		particle.v.x = normx(Randomizer::gen());
		particle.v.y = norm(Randomizer::gen());
		particle.v.z = norm(Randomizer::gen());

		//particle.x = uni(Randomizer::gen());

		particle.mass = m;
	}
	DistribUniformly();
}
void SpatHomSolver::InitialDistribution(double m, double u, double E, int n)
{
	N = n;
	layer.assign(n, Particle());
	normal_distribution<> normx(u, sqrt(C_2_3 * E));
	normal_distribution<> norm(0, sqrt(C_2_3 * E));
	uniform_real_distribution<> uni(left, right);
	for (auto& particle : layer)
	{
		particle.v.x = normx(Randomizer::gen());
		particle.v.y = norm(Randomizer::gen());
		particle.v.z = norm(Randomizer::gen());

		//particle.x = uni(Randomizer::gen());

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
	//prev_layer = layer;
	N = layer.size();
	DataInfo info;
	Histogram3D F(5, -5.0, 5.0, 0);
	F.Build(layer);
	for (int i = 0; i < N; ++i)
	{
		for (int j = i + 1; j < N; ++j)
		{
			auto w = Randomizer::sampleOmega();
			double prod = dot(w, layer[i].v - layer[j].v);
			double l = abs(prod) * dt;
			double lambda = l / N / Kn / 2 * dw * layer[j].mass * F.Get(layer[j].v);

			if (Randomizer::sampleUni() < lambda * 2)
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
DataInfo SpatHomSolver::CalculateLangevinStep(double E)
{
	N = layer.size();
	double avgV = 0;
	for (auto& particle : layer)
	{
		avgV += particle.v.x;
	}
	avgV /= N;
	const double a = 3;// abs(avgV);
	const double sigma = sqrt(3);//sqrt(4.0 / 3.0 * E * a);
	DataInfo info;
	for (int i = 0; i < N; ++i)
	{
		Vec3 v = layer[i].v;
		layer[i].v = (v - (v - Vec3(avgV, 0, 0)) * a * dt + Randomizer::sampleNormVec3() * sdt * sigma);
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
void SpatHomSolver::PlayParticle(Particle& p, double u, double E)
{
	normal_distribution<> norm(0, sqrt(C_2_3 * E));
	uniform_real_distribution<> uni(left, right);
	p.v.x = norm(Randomizer::gen()) + u;
	p.v.y = norm(Randomizer::gen());
	p.v.z = norm(Randomizer::gen());

	p.x = uni(Randomizer::gen());
}

double SmoothFunction(double x, double left, double right, double leftValue, double rightValue)
{
	constexpr double C23 = 2.0 / 3.0;
	constexpr double C43 = 4.0 / 3.0;
	const double M = 2 / (right - left);
	const double L = 0;// left + (right - left) / 2;
	const double scale = leftValue - rightValue;
	return (((-x * M + L) - pow(-x * M + L, 3) / 3 + C23) / C43) * scale + rightValue;
}

void SpatHomSolver::InitialSmoothedDistribution(double rho1, double u1, double E1, double rho2, double u2, double E2, int n, double l, double r)
{
	ofstream file("smooth.txt");
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
		double E = SmoothFunction(x, l, r, E1, E2);
		normal_distribution<> norm(0, sqrt(C_2_3 * E));
		layer[i].v.x = norm(Randomizer::gen()) + u;
		layer[i].v.y = norm(Randomizer::gen());
		layer[i].v.z = norm(Randomizer::gen());
		x += dx;
	}
}