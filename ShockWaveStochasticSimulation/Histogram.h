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
	void ReplenishV(vector<SpatHomSolver>& cells);
	void ReplenishT(vector<SpatHomSolver>& cells);
	void ReplenishDistr(SpatHomSolver& cell);
	void Average(int M);
	void Write(std::string fileName);
};

struct Histogram3D
{
	vector<vector<vector<double>>> hist;
	int sect;				//размер гистограммы ПО ОДНОМУ НАПРАВЛЕНИЮ
	double minX, maxX, minY, maxY, minZ, maxZ;		//по каджому направлению
	int N;
	double volume;
	Histogram3D(int sect, double min, double max, double Ex) :
		sect(sect), minY(min), maxY(max), minZ(min), maxZ(max), minX(0), maxX(0)
		, hist(sect, vector<vector<double>>(sect, vector<double>(sect, 0)))
	{
		double avg = (max - min) / 2;
		minX = Ex - avg;
		maxX = Ex + avg;
		volume = pow((max - min) / sect, 3);
	}
	void Build(vector<Particle>& particles);
	void Write(std::ostream& out);
	double Get(Vec3& v);
	double Get(Vec3& v, Vec3& V, double T);
};