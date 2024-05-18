#pragma once 
#include<vector>
#include<ostream>
#include"SpatHomSolver.h"

using std::vector;

struct Histogram
{
	vector<double> hist;
	int sect;
	double min, max;

	Histogram(int sect, double min, double max) : sect(sect), min(min), max(max), hist(sect,0) {}
	//void ReplenishX(vector<Particle>& particles);
	void ReplenishX(vector<SpatHomSolver>& cells);
	void ReplenishV(vector<Particle>& particles);
	void Normalize();
	void Write(std::ostream& out);
};

struct Histogram3D
{
	vector<vector<vector<double>>> hist;
	int sect;				//размер гистограммы ПО ОДНОМУ НАПРАВЛЕНИЮ
	double minX, maxX, minY, maxY, minZ, maxZ;		//по каджому направлению
	int N;
	Histogram3D(int sect, double min, double max, double Ex) :
		sect(sect), minY(min), maxY(max), minZ(min), maxZ(max), minX(0), maxX(0)
		, hist(sect, vector<vector<double>>(sect, vector<double>(sect, 0)))
	{
		double avg = (max - min) / 2;
		minX = Ex - avg;
		maxX = Ex + avg;
	}
	void Build(vector<Particle>& particles);
	void Normalize();
	void Write(std::ostream& out);
	double Get(Vec3& v);
};