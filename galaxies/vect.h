#ifndef __VECT_H__
#define __VECT_H__

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;
class Vect
{
public:
	double x, y, z;

	Vect()
	{
		this->x = 0;
		this->y = 0;
		this->z = 0;
	}

	Vect(double _x, double _y, double _z)
	{
		this->x = _x;
		this->y = _y;
		this->z = _z;
	}

	Vect operator-(const Vect &v)
	{
		Vect r;
		r.x = this->x - v.x;
		r.y = this->y - v.y;
		r.z = this->z - v.z;
		return r;
	}

	Vect operator+(const Vect &v)
	{
		Vect r;
		r.x = this->x + v.x;
		r.y = this->y + v.y;
		r.z = this->z + v.z;
		return r;
	}

	Vect operator+=(const Vect &v)
	{
		this->x += v.x;
		this->y += v.y;
		this->z += v.z;
		return *this;
	}

	Vect operator*=(const double &scalar)
	{
		this->x *= scalar;
		this->y *= scalar;
		this->z *= scalar;
		return *this;
	}

	Vect operator/=(const double &scalar)
	{
		this->x /= scalar;
		this->y /= scalar;
		this->z /= scalar;
		return *this;
	}

	Vect operator*(const double &scalar)
	{
		Vect re;
		re.x = this->x * scalar;
		re.y = this->y * scalar;
		re.z = this->z * scalar;
		return re;
	}

	Vect operator/(const double &scalar)
	{
		Vect re;
		re.x = this->x / scalar;
		re.y = this->y / scalar;
		re.z = this->z / scalar;
		return re;
	}

	double Int()
	{
		return sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
	}

	double SquareInt()
	{
		return this->x * this->x + this->y * this->y + this->z * this->z;
	}

	double CubeInt()
	{
		double square = this->x * this->x + this->y * this->y + this->z * this->z;
		return sqrt(square) * square;
	}

	Vect Rotate(double theta)
	{
		Vect ret;
		ret.z = this->z * cos(theta) - this->y * sin(theta);
		ret.y = this->z * sin(theta) + this->y * cos(theta);
		ret.x = this->x;
		return ret;
	}

	Vect RotateY(double theta)
	{
		Vect ret;
		ret.z = this->z * cos(theta) - this->x * sin(theta);
		ret.x = this->z * sin(theta) + this->x * cos(theta);
		ret.y = this->y;

		return ret;
	}
};

#endif
