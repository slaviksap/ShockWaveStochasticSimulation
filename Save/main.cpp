#include"ShockWaveSolver.h"
#include<iostream>
int main()
{
	double t = 1;
	int k = 100;
	int M = 100;
	int sect = 200;
	double visuaDt = 0.1;
	ShockWaveSolver solver(t, k, M, sect, -20, 20, visuaDt);
	ShockWaveParams params;
	params.rho1 = 2, params.rho2 = 1;
	params.u1 = 2, params.u2 = 0;
	solver.params = params;
	solver.CheckParams();
	solver.Solve();
	std::cout << "Ready!\n\a";
	//std::cin.get();
	return 0;
}