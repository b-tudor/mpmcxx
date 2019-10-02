#include <iomanip>
#include <chrono>
#include <sstream>

#include "Vector3D.h"
#include "Rando.h"



Vector3D::Vector3D()
{
	X = 0;
	Y = 0;
	Z = 0;
}
Vector3D::Vector3D(double x, double y, double z)
{
	X = x;
	Y = y;
	Z = z;
}

Vector3D::~Vector3D(void)
{
}

void Vector3D::set(double x, double y, double z)
{
	X = x;
	Y = y;
	Z = z;
}

Vector3D & Vector3D::operator=(const Vector3D &rhs) {
	if (this != &rhs) {
		this->X = rhs.x();
		this->Y = rhs.y();
		this->Z = rhs.z();
	}

	return *this;
}

bool Vector3D::equals(const Vector3D &rhs, const double threshold) {

	if (squaredDistance(*this, rhs) <= (threshold*threshold))
		return true;
	return false;
}

Vector3D & Vector3D::operator=(const double rhs) {

	this->X = rhs;
	this->Y = rhs;
	this->Z = rhs;

	return *this;
}

Vector3D & Vector3D::operator+=(const Vector3D &rhs) {
	this->X = this->X + rhs.x();
	this->Y = this->Y + rhs.y();
	this->Z = this->Z + rhs.z();

	return *this;
}

Vector3D & Vector3D::operator-=(const Vector3D &rhs) {
	this->X = this->X - rhs.x();
	this->Y = this->Y - rhs.y();
	this->Z = this->Z - rhs.z();

	return *this;
}

Vector3D & Vector3D::operator*=(const double scalar) {
	this->X *= scalar;
	this->Y *= scalar;
	this->Z *= scalar;
	return (*this);
}

Vector3D & Vector3D::operator/=(const double scalar) {
	this->X /= scalar;
	this->Y /= scalar;
	this->Z /= scalar;
	return (*this); 
}

std::string Vector3D::xyzOut(std::string symbol) {
	using namespace std;
	stringstream out;
	out << left << setw(5) << symbol << right;
	out << fixed << setprecision(5);
	out << setw(10) << X << " "
		<< setw(10) << Y << " "
		<< setw(10) << Z;

	return out.str();
}

Vector3D & Vector3D::normalize() {
	double mag = norm();
	if (mag != 0) {
		X = X / mag;
		Y = Y / mag;
		Z = Z / mag;
	}
	else
		X = Y = Z = 0;
	return *this;
}

void Vector3D::randomize() {
	X = Rando::rand_normal();
	Y = Rando::rand_normal();
	Z = Rando::rand_normal();
	normalize();
}

// Non-member functions
// Constant * vector
Vector3D operator*(const double d, const Vector3D rhs) {
	return Vector3D(d*rhs.x(), d*rhs.y(), d*rhs.z());
}

Vector3D Vector3D::operator-() const {
	return Vector3D(-X, -Y, -Z);
}
