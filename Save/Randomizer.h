#pragma once
#include<vector>
#include<random>
#include<cmath>
#include"Vec3.h"
#include<omp.h>
using std::vector;
using std::mt19937;

namespace Randomizer
{
	static vector<mt19937> gens(8, mt19937(23232));
	static std::normal_distribution<> distr(0, 1);
	static std::uniform_real_distribution<> uni(0, 1);

	void init(int threads_num, int seed, int jump);
	static vector<Vec3> omegas;
	void InitOmegas(int cPhi, int cThette);
	Vec3 sampleOmega();
	double sampleNorm();
	double sampleUni();
	Vec3 sampleNormVec3();

	mt19937& gen();
};