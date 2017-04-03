#pragma once

#include <cmath>
#include <iostream>

class Vector3	{
public:

	//ctors
	Vector3(void) {
		ToZero();
	}

	Vector3(const float x, const float y, const float z)
	{
		vec_array[0] = x;
		vec_array[1] = y;
		vec_array[2] = z;
	}

	~Vector3(void){}



	//Default States
	inline void ToZero() { x = y = z = float(0.0); }



	//Accessors
	inline float  operator[](const int index) const { return vec_array[index]; }
	inline float& operator[](const int index) { return vec_array[index]; }



	//Common Functionality (non-static)
	inline float Length() const
	{
		return sqrt((x*x) + (y*y) + (z*z));
	}
	inline float LengthSquared() const
	{
		return ((x*x + y*y) + (z*z));
	}
	inline float Dot(const Vector3 &b) const
	{
		return (x*b.x) + (y*b.y) + (z*b.z);
	}
	inline Vector3  Cross(const Vector3 &b) const
	{
		return Vector3((y*b.z) - (z*b.y), (z*b.x) - (x*b.z), (x*b.y) - (y*b.x));
	}
	static Vector3 ClampVars(const Vector3& in, float lwr, float upr)
	{
		Vector3 out = in;
		for (int i = 0; i < 3; ++i)
		{ 
			if (out[i] > lwr)
			{
				if (out[i] > upr)
					out[i] = upr;
			}
			else
				out[i] = lwr;
		};

		return out;
	}

	Vector3& Normalise() {
		float length = Length();

		if (length != float(0.0))	{
			length = float(1.0) / length;
			x = x * length;
			y = y * length;
			z = z * length;
		}
		else
		{
			x = float(0.0);
			y = float(0.0);
			z = float(0.0);
		}
		return *this;
	}




	//Common Functionality (Static call)
	static Vector3 Normalise(const Vector3& rhs) {
		float length = rhs.Length();

		Vector3 out;
		if (length != float(0.0))	{
			length = float(1.0) / length;
			out.x = rhs.x * length;
			out.y = rhs.y * length;
			out.z = rhs.z * length;
		}
		else
		{
			out.x = float(1.0);
			out.y = float(0.0);
			out.z = float(0.0);
		}
		return out;
	}
	static float Dot(const Vector3 &a, const Vector3 &b)
	{
		return (a.x*b.x) + (a.y*b.y) + (a.z*b.z);
	}
	static Vector3	Cross(const Vector3 &a, const Vector3 &b)
	{
		return Vector3((a.y*b.z) - (a.z*b.y),
			(a.z*b.x) - (a.x*b.z),
			(a.x*b.y) - (a.y*b.x));
	}
	static Vector3 InterpolateLinear(const Vector3& pStart, const Vector3& pEnd, float pFactor)
	{
		return (pStart * (float(1.0) - pFactor) + pEnd * pFactor);
	}






	//Operator Overloading
	inline Vector3  operator+(const Vector3& a) const{
		return Vector3(x + a.x, y + a.y, z + a.z);
	}

	inline Vector3  operator-(const Vector3& a) const{
		return Vector3(x - a.x, y - a.y, z - a.z);
	}

	inline Vector3  operator-() const{
		return Vector3(-x, -y, -z);
	}
	inline Vector3& operator+=(const Vector3& a){
		x += a.x;
		y += a.y;
		z += a.z;
		return *this;
	}

	inline Vector3& operator-=(const Vector3& a){
		x -= a.x;
		y -= a.y;
		z -= a.z;
		return *this;
	}

	inline Vector3  operator+(const float a) const{
		return Vector3(x + a, y + a, z + a);
	}

	inline Vector3  operator-(const float a) const{
		return Vector3(x - a, y - a, z - a);
	}

	inline Vector3  operator*(const float a) const{
		return Vector3(x * a, y * a, z * a);
	}

	inline Vector3& operator*=(const float& a) {
		x *= a;
		y *= a;
		z *= a;
		return *this;
	}

	inline Vector3  operator*(const Vector3& a) const{
		return Vector3(x * a.x, y * a.y, z * a.z);
	}

	inline Vector3  operator/(const Vector3& v) const{
		return Vector3(x / v.x, y / v.y, z / v.z);
	};

	inline Vector3  operator/(const float v) const{
		return Vector3(x / v, y / v, z / v);
	};

	inline bool	operator==(const Vector3 &A) const {
		return (A.x == x && A.y == y && A.z == z) ? true : false;
	};

	inline bool	operator!=(const Vector3 &A)const {
		return (A.x == x && A.y == y && A.z == z) ? false : true;
	};


	inline friend std::ostream& operator<<(std::ostream& o, const Vector3& v) {
		o << "Vector3(" << v.x << "," << v.y << "," << v.z << ")" << std::endl;
		return o;
	}



	//Union Array Representation
	union
	{
		float vec_array[3];
		struct {
			float x;
			float y;
			float z;
		};
	};
};



