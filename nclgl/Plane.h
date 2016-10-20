/******************************************************************************
Class:Plane
Implements:
Author:Rich Davison	<richard.davison4@newcastle.ac.uk>
Description:VERY simple Plane class. Students are encouraged to modify 
this as necessary!

-_-_-_-_-_-_-_,------,   
_-_-_-_-_-_-_-|   /\_/\   NYANYANYAN
-_-_-_-_-_-_-~|__( ^ .^) /
_-_-_-_-_-_-_-""  ""   

*//////////////////////////////////////////////////////////////////////////////

#pragma once
#include "vector3.h"

class Plane	{
public:
	Plane(void){};
	Plane(const Vector3 &normal, float distance, bool normalise = false);
	~Plane(void){};

	//Sets the planes normal, which should be UNIT LENGTH!!!
	void	SetNormal(const Vector3 &normal) {this->normal = normal;}
	//Gets the planes normal.
	Vector3 GetNormal() const				 {return normal;}
	//Sets the planes distance from the origin
	void	SetDistance(float dist)	{distance = dist;}
	//Gets the planes distance from the origin
	float	GetDistance() const		{return distance;}
	//Performs a simple sphere / plane test
	bool SphereInPlane(const Vector3 &position, float radius) const; 
	//Performs a simple sphere / point test
	bool PointInPlane(const Vector3 &position) const;

	float DistFromPlane(const Vector3& P) const
	{
		// if N is not normalized this is *not* really the distance, 
		// but the computations work just the same.
		return Vector3::Dot(normal, P) + distance;
	}

	bool GetSegmentPlaneIntersection(const Vector3& P1, const Vector3& P2, Vector3& outP) const 
	{
		float d1 = DistFromPlane(P1),
			d2 = DistFromPlane(P2);

		if (d1*d2 > 0)  // points on the same side of plane
			return false;

		float t = d1 / (d1 - d2); // 'time' of intersection point on the segment
		outP = P1 + (P2 - P1) * t;

		return true;
	}

protected:
	//Unit-length plane normal
	Vector3 normal;
	//Distance from origin
	float	distance;
};