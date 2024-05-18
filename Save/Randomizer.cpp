#include "Randomizer.h"
#include<iostream>
double constexpr _Pi = 3.14159265358979;
void Randomizer::init(int threads_num, int seed, int jump)
{
	//gens.assign(threads_num, mt19937(seed));
	for (int i = 0; i < threads_num; ++i)
		gens[i].seed(seed);
	for (int i = 0; i < threads_num; ++i)
		gens[i].discard(jump * i);
}

void Randomizer::InitOmegas(int cPhi, int cThette)
{
	double dPhi = 2 * _Pi / cPhi;
	double dThette = _Pi / cThette;
	for (double phi = 0; phi < 2 * _Pi; phi += dPhi)
	{
		for (double thetta = 0; thetta < _Pi; thetta += dThette)
		{
			omegas.emplace_back(Vec3{ sin(thetta) * cos(phi),sin(thetta) * sin(phi),cos(thetta) });
		}
	}
}

Vec3 Randomizer::sampleOmega()
{
	//std::cout << "Size = " << omegas.size() - 1 << "\n";
	static std::uniform_int_distribution<>uni_omega(0, omegas.size() - 1);
	return omegas[uni_omega(gen())];
}

double Randomizer::sampleNorm()
{
	return distr(gens[omp_get_thread_num()]);
}

double Randomizer::sampleUni()
{
	return uni(gens[omp_get_thread_num()]);
}

Vec3 Randomizer::sampleNormVec3()
{
	Vec3 v;
	int threadnum = omp_get_thread_num();
	v.x = distr(gens[threadnum]);
	v.y = distr(gens[threadnum]);
	v.z = distr(gens[threadnum]);
	return v;
}

std::mt19937& Randomizer::gen()
{
	return gens[omp_get_thread_num()];
}
