#pragma once
#include<cassert>
class Vec3
{
public:
	double x;
	double y;
	double z;

	Vec3(double x, double y, double z) : x(x), y(y), z(z) {}
	Vec3() : x(0), y(0), z(0) {}
	Vec3 operator+(const Vec3& v)
	{
		return Vec3(x + v.x, y + v.y, z + v.z);
	}
	Vec3 operator-(const Vec3& v)
	{
		return Vec3(x - v.x, y - v.y, z - v.z);
	}
	Vec3 operator-()
	{
		x = -x;
		y = -y;
		z = -z;
		return *this;
	}
	Vec3 operator*(const double c)
	{
		return Vec3(x * c, y * c, z * c);
	}
	Vec3 operator/(const double c)
	{
		//if (c < 1e-16)
		//	throw "Division by zero";
		assert(c > 1e-16, "Division by zero");
		return Vec3(x / c, y / c, z / c);
	}
	Vec3& operator+=(const Vec3& v)
	{
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}
	double friend dot(const Vec3& left, const Vec3& right)
	{
		return left.x * right.x + left.y * right.y + left.z * right.z;
	}
	double friend magnitide(const Vec3& v)
	{
		return dot(v, v);
	}
	Vec3 Unit()
	{
		double l = magnitide(*this);
		return *this / l;
	}
};

