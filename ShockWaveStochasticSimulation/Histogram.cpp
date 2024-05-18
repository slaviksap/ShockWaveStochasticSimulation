#include "Histogram.h"
#include <fstream>
#include <mpi.h>
#define MCW MPI_COMM_WORLD
using namespace std;


void Histogram::ReplenishX(vector<SpatHomSolver>& cells)
{
	vector<double> temp(sect, 0);
	const double dx = (max - min) / sect;
	for (const auto& cell : cells)
	{
		for (const Particle& particle : cell.layer)
		{
			const double x = particle.x;
			if (min < x && x <= max)
			{
				temp[size_t((x - min) / dx)] += particle.mass;
			}
		}
	}
	for (double& x : temp)
		x /= dx;
	for (int i = 0; i < hist.size(); ++i)
		hist[i] += temp[i];
}
void Histogram::ReplenishV(vector<SpatHomSolver>& cells)
{
	vector<double> temp(sect, 0);
	vector<double> counter(sect, 0);
	const double dx = (max - min) / sect;
	for (const auto& cell : cells)
	{
		for (const Particle& particle : cell.layer)
		{
			const double x = particle.x;
			if (min < x && x <= max)
			{
				size_t idx = size_t((x - min) / dx);
				temp[idx] += particle.v.x * particle.mass;
				counter[idx] += particle.mass;
			}
		}
	}
	for (int i = 0; i < temp.size(); ++i)
		temp[i] /= counter[i];
	for (int i = 0; i < hist.size(); ++i)
		hist[i] += temp[i];
}
void Histogram::ReplenishT(vector<SpatHomSolver>& cells)
{
	vector<double> avgV(sect, 0);
	vector<double> avgVV(sect, 0);
	vector<double> counter(sect, 0);
	const double dx = (max - min) / sect;
	for (const auto& cell : cells)
	{
		for (const Particle& particle : cell.layer)
		{
			const double x = particle.x;
			if (min < x && x <= max)
			{
				size_t idx = size_t((x - min) / dx);
				avgV[idx] += particle.v.x * particle.mass;
				avgVV[idx] += particle.v.x * particle.v.x * particle.mass;
				counter[idx] += particle.mass;
			}
		}
	}
	for (int i = 0; i < avgV.size(); ++i)
	{
		avgV[i] /= counter[i];
		avgVV[i] /= counter[i];
	}
	for (int i = 0; i < hist.size(); ++i)
		hist[i] += avgVV[i] - avgV[i] * avgV[i];
}
void Histogram::ReplenishDistr(SpatHomSolver& cell)
{
	double integral = 0;
	const double dx = (max - min) / sect;
	for (int i = 0; i < cell.layer.size(); ++i)
	{
		Particle& p = cell.layer[i];
		double v = p.v.x;
		if (v > min && v < max) 
		{
			hist[int((v - min) / dx)] +=p.mass;
			integral += p.mass / dx;
		}
	}
	for (int i = 0; i < hist.size(); ++i)
	{
		hist[i] /= integral;
	}
}
void Histogram::Average(int M)
{
	int size, rank;
	MPI_Comm_size(MCW, &size);
	MPI_Comm_rank(MCW, &rank);
	MPI_Barrier(MCW);
	int globalM = M;
	if (size > 1)
	{
		vector<double> rcvBuf;
		if (rank == 0)
			rcvBuf.assign(hist.size(), 0);
		MPI_Reduce(&hist[0], &rcvBuf[0], hist.size(), MPI_DOUBLE, MPI_SUM, 0, MCW);
		hist = rcvBuf;
		MPI_Reduce(&M, &globalM, 1, MPI_INT, MPI_SUM, 0, MCW);
	}

	for (auto& cell : hist)
	{
		cell /= globalM;
	}
}
void Histogram::Write(std::string fileName)
{
	int size, rank;
	MPI_Comm_size(MCW, &size);
	MPI_Comm_rank(MCW, &rank);
	if (rank == 0)
	{
		ofstream file(fileName);
		const double dx = (max - min) / sect;
		for (int i = 0; i < hist.size(); ++i)
		{
			double x = min + (i + 0.5) * dx;
			file << x << "\t" << hist[i] << "\n";
		}
		file.close();
	}
}

void Histogram3D::Build(vector<Particle>& particles)
{
	for (int i = 0; i < hist.size(); ++i)
		for (int j = 0; j < hist[i].size(); ++j)
			for (int k = 0; k < hist[i][j].size(); ++k)
				hist[i][j][k] = 0;
	double dx = (maxX - minX) / sect;
	double dy = (maxY - minY) / sect;
	double dz = (maxZ - minZ) / sect;
	N = particles.size();
	for (const auto& particle : particles)
	{
		int xIdx = int((particle.v.x - minX) / dx);
		int yIdx = int((particle.v.y - minY) / dy);
		int zIdx = int((particle.v.z - minZ) / dz);
		if (xIdx < 0 || xIdx >= sect) continue;
		if (yIdx < 0 || yIdx >= sect) continue;
		if (zIdx < 0 || zIdx >= sect) continue;
		hist[xIdx][yIdx][zIdx] += particle.mass / volume;
	}
}

void Histogram3D::Write(std::ostream& out)
{
}

double Histogram3D::Get(Vec3& v)
{
	double dx = (maxX - minX) / sect;
	double dy = (maxY - minY) / sect;
	double dz = (maxZ - minZ) / sect;
	int xIdx = int((v.x - minX) / dx);
	int yIdx = int((v.y - minY) / dy);
	int zIdx = int((v.z - minZ) / dz);
	if (xIdx < 0 || xIdx >= sect) return 1.0;
	if (yIdx < 0 || yIdx >= sect) return 1.0;
	if (zIdx < 0 || zIdx >= sect) return 1.0;
	return hist[xIdx][yIdx][zIdx];
}

double Histogram3D::Get(Vec3& v, Vec3& V, double T)
{
	const double dv = 3;
	const double vv = (v.x - V.x) * (v.x - V.x) + (v.y - V.y) * (v.y - V.y) + (v.z - V.z) * (v.z - V.z);
	return pow(1.0 / (sqrt(2 * T * 3.14159265)), 1.5) * exp(-vv / (T)) * dv;
}
