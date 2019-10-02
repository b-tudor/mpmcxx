// ©2013 Adam Hogan
// Space Research Group
// Department of Chemistry
// University of South Florida

// see http://www.cprogramming.com/tutorial/3d/quaternions.html
// and http://content.gpwiki.org/index.php/OpenGL:Tutorials:Using_Quaternions_to_represent_rotation

#pragma once

#include "Vector3D.h"


class Quaternion
{
public:
	Quaternion( double x, double y, double z, double w, int mode );
	Quaternion(Vector3D v);
	~Quaternion();

	static const int XYZW = 1;
	static const int AXIS_ANGLE_RADIAN = 2;
	static const int AXIS_ANGLE_DEGREE = 3;

	void normalize();
	Quaternion  operator*( const Quaternion &rhs ) const;
	Quaternion& operator=( const Quaternion &rhs );
	Quaternion  conjugate() const;
	Vector3D rotate( const Vector3D &rotatee ) const;

	inline double x() { return X; }
	inline double y() { return Y; }
	inline double z() { return Z; }
	inline double w() { return W; }

private:
	double X, Y, Z, W;
  
};

