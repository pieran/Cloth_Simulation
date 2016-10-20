
#pragma once
#ifndef MAT33_H
#define MAT33_H

#include "Vector3.h"

class Matrix4;

class Matrix3
{
public:

	static const Matrix3 Identity;
	static const Matrix3 ZeroMatrix;

	//ctor
	Matrix3();

	Matrix3(float elements[9]);

	Matrix3(const Vector3& c1, const Vector3& c2, const Vector3& c3);

	Matrix3(float a1, float a2, float a3,
			float b1, float b2, float b3,
			float c1, float c2, float c3);

	Matrix3(const Matrix4& mat44);

	~Matrix3(void);




	//Default States
	void	ToZero();
	void	ToIdentity();



	//Default Accessors
	inline float   operator[](int index) const        { return mat_array[index]; }
	inline float&  operator[](int index)              { return mat_array[index]; }
	inline float   operator()(int row, int col) const { return mat_array[row + col * 3]; }
	inline float&  operator()(int row, int col)       { return mat_array[row + col * 3]; }

	inline const Vector3&	GetCol(int idx) const				{ return *((Vector3*)&mat_array[idx * 3]); }
	inline void				SetCol(int idx, const Vector3& row)	{ memcpy(&mat_array[idx * 3], &row.x, 3 * sizeof(float)); }

	inline Vector3			GetRow(int idx)	const				{ return Vector3(mat_array[idx], mat_array[3 + idx], mat_array[6 + idx]); }
	inline void				SetRow(int idx, const Vector3& col)	{ mat_array[idx] = col.x; mat_array[3 + idx] = col.y; mat_array[6 + idx] = col.z; }



	//Common Matrix Properties
	inline Vector3			GetScalingVector() const			{ return Vector3(_11, _22, _33); }
	inline void				SetScalingVector(const Vector3& in)	{ _11 = in.x, _22 = in.y, _33 = in.z; }



	//Transformation Matrix
	static Matrix3 Rotation(float degrees, const Vector3 &axis);
	static Matrix3 Scale(const Vector3 &scale);



	// Standard Matrix Functionality
	static Matrix3 Inverse(const Matrix3& rhs);
	static Matrix3 Transpose(const Matrix3& rhs);
	static Matrix3 Adjugate(const Matrix3& m);
	static Matrix3 OrthoNormalise(const Matrix3& m);
	static Matrix3 OrthoNormalise(const Vector3 rows[3], const int primary_row = 0);
	static Matrix3 OrthoNormalise(const Vector3& row1, const Vector3& row2, const Vector3& row3, const int primary_row = 0);
	static Matrix3 OuterProduct(const Vector3& a, const Vector3& b);
	static Matrix3 TransposeMirror(const Matrix3& rhs);


	// Additional Functionality
	float Trace() const;
	float Determinant() const;



	//Other representation of data.
	union
	{
		float	mat_array[9];
		struct {
			float _11, _21, _31;
			float _12, _22, _32;
			float _13, _23, _33; 
		};
	};
};
/*
inline Matrix3& operator+=(Matrix3& a, const Matrix3& rhs);
inline Matrix3& operator-=(Matrix3& a, const Matrix3& rhs);

inline Matrix3 operator+(const Matrix3& a, const Matrix3& rhs);
inline Matrix3 operator-(const Matrix3& a, const Matrix3& rhs);
inline Matrix3 operator*(const Matrix3& a, const Matrix3& rhs);


inline Matrix3& operator+=(Matrix3& a, const float b);
inline Matrix3& operator-=(Matrix3& a, const float b);
inline Matrix3& operator*=(Matrix3& a, const float b);
inline Matrix3& operator/=(Matrix3& a, const float b);

inline Matrix3 operator+(const Matrix3& a, const float b);
inline Matrix3 operator-(const Matrix3& a, const float b);
inline Matrix3 operator*(const Matrix3& a, const float b);
inline Matrix3 operator/(const Matrix3& a, const float b);

inline Vector3 operator*(const Matrix3& a, const Vector3& b);

#include "Matrix3.inl"*/


inline Matrix3& operator+=(Matrix3& a, const Matrix3& b)
{
	for (unsigned int i = 0; i < 9; ++i)
		a.mat_array[i] += b.mat_array[i];
	return a;
}

inline Matrix3& operator-=(Matrix3& a, const Matrix3& b)
{
	for (unsigned int i = 0; i < 9; ++i)
		a.mat_array[i] -= b.mat_array[i];
	return a;
}

inline Matrix3 operator+(const Matrix3& a, const Matrix3& b)
{
	Matrix3 m;
	for (unsigned int i = 0; i < 9; ++i)
		m.mat_array[i] = a.mat_array[i] + b.mat_array[i];
	return m;
}

inline Matrix3 operator-(const Matrix3& a, const Matrix3& b)
{
	Matrix3 m;
	for (unsigned int i = 0; i < 9; ++i)
		m.mat_array[i] = a.mat_array[i] - b.mat_array[i];
	return m;
}


inline Matrix3& operator+=(Matrix3& a, const float b)
{
	for (unsigned int i = 0; i < 9; ++i)
		a.mat_array[i] += b;
	return a;
}

inline Matrix3& operator-=(Matrix3& a, const float b)
{
	for (unsigned int i = 0; i < 9; ++i)
		a.mat_array[i] -= b;
	return a;
}
inline Matrix3& operator*=(Matrix3& a, const float b)
{
	for (unsigned int i = 0; i < 9; ++i)
		a.mat_array[i] *= b;
	return a;
}
inline Matrix3& operator/=(Matrix3& a, const float b)
{
	for (unsigned int i = 0; i < 9; ++i)
		a.mat_array[i] /= b;
	return a;
}

inline Matrix3 operator+(Matrix3& a, const float b)
{
	Matrix3 m;
	for (unsigned int i = 0; i < 9; ++i)
		m.mat_array[i] = a.mat_array[i] + b;
	return m;
}

inline Matrix3 operator-(const Matrix3& a, const float b)
{
	Matrix3 m;
	for (unsigned int i = 0; i < 9; ++i)
		m.mat_array[i] = a.mat_array[i] - b;
	return m;
}
inline Matrix3 operator*(const Matrix3& a, const float b)
{
	Matrix3 m;
	for (unsigned int i = 0; i < 9; ++i)
		m.mat_array[i] = a.mat_array[i] * b;
	return m;
}
inline Matrix3 operator/(const Matrix3& a, const float b)
{
	Matrix3 m;
	for (unsigned int i = 0; i < 9; ++i)
		m.mat_array[i] = a.mat_array[i] / b;
	return m;
}


inline Matrix3 operator*(const Matrix3& a, const Matrix3& b)
{
	Matrix3 out;

	out._11 = a._11 * b._11 + a._12 * b._21 + a._13 * b._31;
	out._12 = a._11 * b._12 + a._12 * b._22 + a._13 * b._32;
	out._13 = a._11 * b._13 + a._12 * b._23 + a._13 * b._33;

	out._21 = a._21 * b._11 + a._22 * b._21 + a._23 * b._31;
	out._22 = a._21 * b._12 + a._22 * b._22 + a._23 * b._32;
	out._23 = a._21 * b._13 + a._22 * b._23 + a._23 * b._33;

	out._31 = a._31 * b._11 + a._32 * b._21 + a._33 * b._31;
	out._32 = a._31 * b._12 + a._32 * b._22 + a._33 * b._32;
	out._33 = a._31 * b._13 + a._32 * b._23 + a._33 * b._33;

	return out;
}

inline void InplaceMatrix3MultMatrix3(Matrix3* out, const Matrix3& a, const Matrix3& b)
{
	out->_11 = a._11 * b._11 + a._12 * b._21 + a._13 * b._31;
	out->_12 = a._11 * b._12 + a._12 * b._22 + a._13 * b._32;
	out->_13 = a._11 * b._13 + a._12 * b._23 + a._13 * b._33;

	out->_21 = a._21 * b._11 + a._22 * b._21 + a._23 * b._31;
	out->_22 = a._21 * b._12 + a._22 * b._22 + a._23 * b._32;
	out->_23 = a._21 * b._13 + a._22 * b._23 + a._23 * b._33;

	out->_31 = a._31 * b._11 + a._32 * b._21 + a._33 * b._31;
	out->_32 = a._31 * b._12 + a._32 * b._22 + a._33 * b._32;
	out->_33 = a._31 * b._13 + a._32 * b._23 + a._33 * b._33;
}

inline void InplaceMatrix3MultMatrix3Additive(Matrix3* out, const Matrix3& a, const Matrix3& b)
{
	out->_11 += a._11 * b._11 + a._12 * b._21 + a._13 * b._31;
	out->_12 += a._11 * b._12 + a._12 * b._22 + a._13 * b._32;
	out->_13 += a._11 * b._13 + a._12 * b._23 + a._13 * b._33;

	out->_21 += a._21 * b._11 + a._22 * b._21 + a._23 * b._31;
	out->_22 += a._21 * b._12 + a._22 * b._22 + a._23 * b._32;
	out->_23 += a._21 * b._13 + a._22 * b._23 + a._23 * b._33;

	out->_31 += a._31 * b._11 + a._32 * b._21 + a._33 * b._31;
	out->_32 += a._31 * b._12 + a._32 * b._22 + a._33 * b._32;
	out->_33 += a._31 * b._13 + a._32 * b._23 + a._33 * b._33;
}

inline Vector3 operator*(const Matrix3& a, const Vector3& b)
{
	Vector3 out;

	out.x = a._11 * b.x;
	out.y = a._21 * b.x;
	out.z = a._31 * b.x;

	out.x += a._12 * b.y;
	out.y += a._22 * b.y;
	out.z += a._32 * b.y;

	out.x += a._13 * b.z;
	out.y += a._23 * b.z;
	out.z += a._33 * b.z;

	/*out.x = a._11 * b.x
		+ a._12 * b.y
		+ a._13 * b.z;

	out.y = a._21 * b.x
		+ a._22 * b.y
		+ a._23 * b.z;

	out.z = a._31 * b.x
		+ a._32 * b.y
		+ a._33 * b.z;*/

	return out;
}

inline void InplaceMatrix3MultVector3(Vector3* out, const Matrix3& a, const Vector3& b)
{
	out->x = a._11 * b.x;
	out->y = a._21 * b.x;
	out->z = a._31 * b.x;

	out->x += a._12 * b.y;
	out->y += a._22 * b.y;
	out->z += a._32 * b.y;

	out->x += a._13 * b.z;
	out->y += a._23 * b.z;
	out->z += a._33 * b.z;
}

inline void InplaceMatrix3MultVector3Additve(Vector3* out, const Matrix3& a, const Vector3& b)
{
	out->x += a._11 * b.x;
	out->y += a._21 * b.x;
	out->z += a._31 * b.x;

	out->x += a._12 * b.y;
	out->y += a._22 * b.y;
	out->z += a._32 * b.y;

	out->x += a._13 * b.z;
	out->y += a._23 * b.z;
	out->z += a._33 * b.z;
}

inline void InplaceMatrix3MultVector3Subtract(Vector3* out, const Matrix3& a, const Vector3& b)
{
	out->x -= a._11 * b.x;
	out->y -= a._21 * b.x;
	out->z -= a._31 * b.x;

	out->x -= a._12 * b.y;
	out->y -= a._22 * b.y;
	out->z -= a._32 * b.y;

	out->x -= a._13 * b.z;
	out->y -= a._23 * b.z;
	out->z -= a._33 * b.z;
}

#endif //MAT33_H