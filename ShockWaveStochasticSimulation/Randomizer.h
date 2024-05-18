#pragma once
#include<vector>
#include<random>
#include<cmath>
#include"Vec3.h"
#include<omp.h>
using std::vector;
using std::mt19937;
using std::mt19937_64;

typedef mt19937_64 generator;
namespace Randomizer
{
	static vector<generator> gens(8, generator(23232));
	static std::normal_distribution<> distr(0, 1);
	static std::uniform_real_distribution<> uni(0, 1);

	void init(int threads_num, int seed, int jump);
	void MPI_Init(int size, int rank, int seed, unsigned long long jump);
	static vector<Vec3> omegas;
	void InitOmegas(int cPhi, int cThette);
	Vec3 sampleOmega();
	double sampleNorm();
	double sampleUni();
	Vec3 sampleNormVec3();

	generator& gen();
};

void pprintf(char const* const _Format, ...);