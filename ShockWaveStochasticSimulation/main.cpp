#include"ShockWaveSolver.h"
#include<iostream>
#include<mpi.h>
#include "Randomizer.h"

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);
	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double t = 1;
	int k = 400;
	int M = 40;
	int N = 400;
	int sect = 200;
	double visuaDt = 0.01;
	int Mlocal = M / size;
	unsigned long long jump = unsigned long long(Mlocal) * unsigned long long(N) * unsigned long long(N) * unsigned long long(k + 3);
	Randomizer::MPI_Init(1, rank, 12345, jump);
	ShockWaveSolver solver(t, k, Mlocal, N, sect, -15, 15, visuaDt);
	ShockWaveParams params;
	params.gamma = 5.0 / 3.0;
	params.Ma = 9;
	params.CalcParams();
	solver.params = params;
	solver.useSmooth = true;
	solver.leftSmooth = -3;
	solver.rightSmooth = 3;
	solver.useDiffusion = true;
	solver.leftDiffusion = -5;
	solver.rightDiffusion = 5;
	MPI_Barrier(MPI_COMM_WORLD);
	double start = MPI_Wtime();
	solver.Solve();
	pprintf("Ready!\n\a");
	pprintf("Calculation time %g\n", MPI_Wtime() - start);
	MPI_Finalize();
	return 0;
}