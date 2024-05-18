#include "Histogram.h"
using namespace std;
//void Histogram::ReplenishX(vector<Particle>& particles)
//{
//	vector<double> temp(sect, 0);
//	const double dx = (max - min) / sect;
//	double M = 0;
//	for (const auto& particle : particles)
//	{
//		const double x = particle.x;
//		if (min < x && x <= max)
//		{
//			temp[size_t((x - min) / dx)] += particle.mass;
//			M == particle.mass;
//		}
//	}
//	for (double& x : temp)
//		x /= M * dx;
//	for (int i = 0; i < hist.size(); ++i)
//		hist[i] += temp[i];
//}

void Histogram::ReplenishX(vector<SpatHomSolver>& cells)
{
	double M = 0;
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
				M += particle.mass;
			}
		}
	}
	for (double& x : temp)
		x /= dx;
	for (int i = 0; i < hist.size(); ++i)
		hist[i] += temp[i];
}
void Histogram::ReplenishV(vector<Particle>& particles)
{
	vector<double> temp(sect, 0);
	const double dx = (max - min) / sect;
	for (const auto& particle : particles)
	{
		const double vx = particle.v.x;
		if (min < vx && vx <= max)
			temp[size_t((vx - min) / dx)] += 1;
	}
	for (double& x : temp)
		x /= particles.size() * dx;
	for (int i = 0; i < hist.size(); ++i)
		hist[i] += temp[i];
}
void Histogram::Normalize()
{
	double integral = 0;
	const double dx = (max - min) / sect;
	for (int i = 0; i < hist.size(); ++i)
		integral += hist[i] * dx;
	for (auto& cell : hist)
		cell /= integral;
}

void Histogram::Write(std::ostream& out)
{
	const double dx = (max - min) / sect;
	for (int i = 0; i < hist.size(); ++i)
	{
		double x = min + (i + 0.5) * dx;
		out << x << "\t" << hist[i] << "\n";
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
		hist[xIdx][yIdx][zIdx] += 1.0;
	}
}

void Histogram3D::Normalize()
{
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
