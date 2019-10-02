#ifndef VECTOR3D_3D_H
#define	VECTOR3D_3D_H

#include "Rando.h"
/*
	Vector3D();                            - creates the vector (0,0,0)
	Vector3D(double x, double y, double z) - creates the vector (x,y,z)
	~Vector3D(void);

	GETTERS/SETTERS
		x() - returns X coord
		y() - returns Y coord
		z() - returns Z coord

		set(double x, double y, double z) - sets x, y and z coords
		setX(double x)                    - sets x coord
		setY(double y)                    - sets y coord
		setZ(double z)                    - sets z coord


	SUPPORTED OPERATORS/OPERATIONS
		Vector3D & operator=(const Vector3D &rhs)            - vector assignment
		Vector3D & operator=(const double rhs)               - assigns all vector components a value of rhs
		bool operator==(const Vector3D &rhs)                 - boolean check for exact equality
		bool equals(const Vector3D &rhs, double threshold)   - boolean equality check such that the magnitude of the difference vector is < threshold

	ADDITION/SUBTRACTION OPERATORS
		Vector3D & operator+=(const Vector3D &rhs)
		Vector3D operator+(const Vector3D &rhs)
		Vector3D & operator-=(const Vector3D &rhs)
		Vector3D operator-(const Vector3D &rhs) const {

	MULTIPLICATIVE OPERATORS/OPERATIONS
		operator*(const Vector3D &rhs)                      - dot (scalar) product
		cross(const Vector3D &rhs)                          - cross (vector) product
		Vector3D operator*(double  rhs)                     - uniformly scales vector by rhs
		friend Vector3D operator*(double lhs, Vector3D rhs) - uniformly scales vector by d
		Vector3D operator/(const double  rhs)               - uniformly scales vector by 1/rhs

	MISCELLANEOUS OPERATIONS
		normalize()                                                      - (void)     maintains vector direction, but scales it to magnitude 1
		randomize()                                                      - (void)     changes vector into a randomly oriented unit vector, selected from a uniform distribution across the surface of a unit sphere
		squaredDistance(Vector3D A, Vector3D B)                          - (double)   returns the distance between points A & B, squared (i.e. distance * distance)
		squaredDistancePBC(Vector3D A, Vector3D B, double BoxLength)     - (double)   returns the distance between points A & B, squared (i.e. distance * distance), accounting for periodic boundary conditions in a cubic cell
		Vector3D differencePBC(Vector3D A, Vector3D B, double BoxLength) - (Vector3D) returns the distance between points A & B, accounting for periodic boundary conditions in a cubic cell
		magnitude()                                                      - (double)   returns magnitude of the vector (absolute distance of head from origin)

	OUTPUT
		xyzOut(std::string symbol) - (std::string) returns a string suitable as a line in an XYZ format text file, using symbol as the atomic symbol, e.g.  xyzOut("H") produces something like "H   1.0    2.0   -3.0"
*/



class Vector3D
{

private:
	double X, Y, Z;


public:
	Vector3D();
	Vector3D(double x, double y, double z);
	~Vector3D(void);

	void set(double x, double y, double z);
	void setX(double x) { X = x; }
	void setY(double y) { Y = y; }
	void setZ(double z) { Z = z; }
	inline double x() const { return X; }
	inline double y() const { return Y; }
	inline double z() const { return Z; }

	

	Vector3D & operator=(const Vector3D &rhs);
	Vector3D & operator=(const double rhs);

	bool operator==(const Vector3D &rhs) {
		return  (this->x() == rhs.x()) && (this->y() == rhs.y()) && (this->z() == rhs.z());
	}

	bool equals(const Vector3D &rhs, const double threshold);



	// Addition
	Vector3D & operator+=(const Vector3D &rhs);
	inline const Vector3D operator+(const Vector3D &rhs) const {
		return Vector3D(*this) += rhs;
	}

	// Subtraction
	Vector3D & operator-=(const Vector3D &rhs);
	const Vector3D operator-(const Vector3D &rhs) const {
		return Vector3D(*this) -= rhs;
	}
	
	
	// Dot (scalar) product
	inline double operator*(const Vector3D &rhs) const {
		return this->x()*rhs.x() + this->y()*rhs.y() + this->z()*rhs.z();
	}

	// Cross (vector) product
	inline Vector3D cross(const Vector3D &rhs) const {
		return Vector3D(
			this->y()*rhs.z() - this->z()*rhs.y(),
			this->z()*rhs.x() - this->x()*rhs.z(),
			this->x()*rhs.y() - this->y()*rhs.x()
		);
	}

	// Scalar * vector
	Vector3D operator-() const;
	Vector3D & operator*=(const double scalar);
	inline Vector3D operator*(const double  rhs) const {
		return Vector3D(this->x()*rhs, this->y()*rhs, this->z()*rhs);
	}
	friend Vector3D operator*(const double lhs, Vector3D rhs);

	// vector / scalar
	Vector3D & operator/=(const double scalar);
	inline Vector3D operator/(const double  rhs) const {
		return Vector3D(this->x() / rhs, this->y() / rhs, this->z() / rhs);
	}



	inline static double squaredDistance(Vector3D A, Vector3D B) { return (A - B)*(A - B); }
	inline static double squaredDistancePBC(Vector3D A, Vector3D B, double BoxLength) {
		Vector3D dr = A - B;
		dr.set(
			dr.x() - BoxLength * nint(dr.x() / BoxLength),
			dr.y() - BoxLength * nint(dr.y() / BoxLength),
			dr.z() - BoxLength * nint(dr.z() / BoxLength)
		);
		return dr * dr;
	}
	inline static double angle(Vector3D A, Vector3D B) { return acos(   (A*B)  /  (A.norm()*B.norm())   ); }

	inline static Vector3D differencePBC(Vector3D A, Vector3D B, double BoxLength) {
		Vector3D dr = A - B;
		dr.set(dr.x() - BoxLength * nint(dr.x() / BoxLength),
		       dr.y() - BoxLength * nint(dr.y() / BoxLength),
		       dr.z() - BoxLength * nint(dr.z() / BoxLength));
		return dr;
	}

	Vector3D & normalize();
	inline double norm2() const { return (*this)*(*this); }
	inline double norm () const { return sqrt(norm2()); }

	// changes Vector3D to be a random unit Vector3D, selected from a uniform distribution across the surface of a unit sphere
	void randomize();



	std::string xyzOut(std::string symbol);





private:
	inline double static nint(double x) { return floor(x + 0.5); }
};

#endif	// VECTOR3D_H