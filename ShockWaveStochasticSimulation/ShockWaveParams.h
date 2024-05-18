#pragma once
struct ShockWaveParams
{
	double sigma1 = -1, sigma2 = -1;
	double u1 = -1, u2 = -1;
	double p1 = -1, p2 = -1;
	double rho1 = -1, rho2 = -1;
	double T1 = -1, T2 = -1;
	double gamma = -1;
	double Ma = -1;
	double D = 0;

	static double sigma(double T);

	void CalcParams();
};