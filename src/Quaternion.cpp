// ©2013 Adam Hogan
// Space Research Group
// Department of Chemistry
// University of South Florida

// see http://www.cprogramming.com/tutorial/3d/quaternions.html
// and http://content.gpwiki.org/index.php/OpenGL:Tutorials:Using_Quaternions_to_represent_rotation

#include <cmath>

#include "constants.h"
#include "Quaternion.h"



Quaternion::Quaternion( double x, double y, double z, double w, int construction_mode )
{
	double angle     = w;
	double sinAngle  = 0,
	       magnitude = 0;

	switch( construction_mode ) {

		case XYZW:
			// Construct quaternion from components
			this->X = x;
			this->Y = y;
			this->Z = z;
			this->W = w;
			break;

		case AXIS_ANGLE_DEGREE:
			// Convert angle from degrees to radians, prior to construction of quaternion
			angle /= 57.2957795; // fall-through here is intentional...

		case AXIS_ANGLE_RADIAN:

			// Construct quaternion from an axis and an angle (in radians)
			// Normalizes the axis vector

			magnitude = sqrt( x*x + y*y + z*z );
			if( magnitude == 0.0 ) // edge case, if the axis to rotate around doesn't exist just return a quaternion with no rotation
			{
				this->X = 0;
				this->Y = 0;
				this->Z = 0;
				this->W = 1;
				return;
			}

			x = x/magnitude;
			y = y/magnitude;
			z = z/magnitude;

			sinAngle = sin(angle/2.0);
		
			this->X = x*sinAngle;
			this->Y = y*sinAngle;
			this->Z = z*sinAngle;
			this->W = cos(angle/2.0);
			break;
	
		default:
			throw invalid_quaternion_mode;
			break;
	}
}


Quaternion::Quaternion(Vector3D v) {
	X = v.x();
	Y = v.y();
	Z = v.z();
	W = 0;
}


Quaternion::~Quaternion() {}


// Normalize quaternion
void Quaternion::normalize()
{
	double magnitude = sqrt( X*X + Y*Y + Z*Z + W*W );
	X = X/magnitude;
	Y = Y/magnitude;
	Z = Z/magnitude;
	W = W/magnitude;
}

Quaternion& Quaternion::operator=(const Quaternion &rhs)
{
	if( &rhs != this ) {  // Avoids self assignment
		X = rhs.X;
		Y = rhs.Y;
		Z = rhs.Z;
		W = rhs.W;
	}
	return *this; // allows chaining, i.e. x = y = z;
}


// QuaternionStore = Q1 * Q2
// Order matters!
Quaternion Quaternion::operator*( const Quaternion &rhs ) const
{
  double result_w = W*rhs.W - X*rhs.X - Y*rhs.Y - Z*rhs.Z;
  double result_x = W*rhs.X + X*rhs.W + Y*rhs.Z - Z*rhs.Y;
  double result_y = W*rhs.Y - X*rhs.Z + Y*rhs.W + Z*rhs.X;
  double result_z = W*rhs.Z + X*rhs.Y - Y*rhs.X + Z*rhs.W;
  
  return Quaternion( result_x, result_y, result_z, result_w, XYZW );
}


// A conjugate quaternion performs the opposite rotation
Quaternion Quaternion::conjugate() const
{
	return Quaternion( -X, -Y, -Z, W, XYZW );
}



Vector3D Quaternion::rotate( const Vector3D &rotatee ) const {
	Quaternion newR = (*this) * rotatee * this->conjugate();
	return Vector3D(newR.x(), newR.y(), newR.z());
}