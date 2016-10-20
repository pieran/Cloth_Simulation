#pragma once
#include <cmath>
#include <iostream>

class Vector2	{
public:

	//ctors
	Vector2(void) : x(0.0f), y(0.0f) {}
	Vector2(const float nx, const float ny) : x(nx), y(ny) {}
	~Vector2(void){}



	//Default States
	inline void ToZero() { x = 0.0f; y = 0.0f; }



	//Accessors
	inline const float&  operator[](const int index) const	{ return vec_array[index]; }
	inline float&		 operator[](const int index)		{ return vec_array[index]; }



	//Common Functionality (non-static)
	inline float Length() const
	{
		return sqrt((x*x) + (y*y));
	}
	inline float LengthSquared() const
	{
		return ((x*x) + (y*y));
	}
	inline float Dot(const Vector2& b) const
	{
		return (x*b.x) + (y*b.y);
	}
	inline float  Cross(const Vector2 &b) const
	{
		return x * b.y - y * b.x;
	}

	Vector2& Normalise() {
		float length = Length();

		if (length > 0.0f)	{
			length = 1.0f / length;
			x = x * length;
			y = y * length;
		}

		return *this;
	}




	//Common Functionality (Static call)
	static Vector2 Normalise(const Vector2& rhs) {
		float length = rhs.Length();

		Vector2 out;
		if (length != float(0.0))	{
			length = float(1.0) / length;
			out.x = rhs.x * length;
			out.y = rhs.y * length;
		}

		return out;
	}

	static inline float Dot(const Vector2& a, const Vector2& b)
	{
		return (a.x*b.x) + (a.y*b.y);
	}

	static inline float  Cross(const Vector2& a, const Vector2 &b)
	{
		return a.x * b.y - a.y * b.x;
	}

	static Vector2 InterpolateLinear(const Vector2& pStart, const Vector2& pEnd, float pFactor)
	{
		return (pStart * (1.0f - pFactor) + pEnd * pFactor);
	}






	//Operator Overloading
	inline Vector2  operator+(const Vector2& a) const{
		return Vector2(x + a.x, y + a.y);
	}

	inline Vector2  operator-(const Vector2& a) const{
		return Vector2(x - a.x, y - a.y);
	}

	inline Vector2  operator-() const{
		return Vector2(-x, -y);
	}

	inline Vector2& operator+=(const Vector2& a){
		x += a.x;
		y += a.y;
		return *this;
	}

	inline Vector2& operator-=(const Vector2& a){
		x -= a.x;
		y -= a.y;
		return *this;
	}

	inline Vector2  operator+(const float a) const{
		return Vector2(x + a, y + a);
	}

	inline Vector2  operator-(const float a) const{
		return Vector2(x - a, y - a);
	}

	inline Vector2  operator*(const float a) const{
		return Vector2(x * a, y * a);
	}

	inline Vector2  operator/(const float v) const{
		return Vector2(x / v, y / v);
	};

	inline bool	operator==(const Vector2 &A) const {
		return (A.x == x && A.y == y) ? true : false;
	};

	inline bool	operator!=(const Vector2 &A)const {
		return (A.x == x && A.y == y) ? false : true;
	};



	//Union Array Representation
	union
	{
		float vec_array[2];
		struct {
			float x;
			float y;
		};
	};
};


