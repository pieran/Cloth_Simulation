
#pragma once
#include "Matrix3.h"
#include "Matrix4.h"






Matrix3& operator+=(Matrix3& a, const Matrix3& b)
{
	for (unsigned int i = 0; i < 9; ++i)
		a.mat_array[i] += b.mat_array[i];
	return a;
}

Matrix3& operator-=(Matrix3& a, const Matrix3& b)
{
	for (unsigned int i = 0; i < 9; ++i)
		a.mat_array[i] -= b.mat_array[i];
	return a;
}

Matrix3 operator+(const Matrix3& a, const Matrix3& b)
{
	Matrix3 m;
	for (unsigned int i = 0; i < 9; ++i)
		m.mat_array[i] = a.mat_array[i] + b.mat_array[i];
	return m;
}

Matrix3 operator-(const Matrix3& a, const Matrix3& b)
{
	Matrix3 m;
	for (unsigned int i = 0; i < 9; ++i)
		m.mat_array[i] = a.mat_array[i] - b.mat_array[i];
	return m;
}

Matrix3 operator*(const Matrix3& a, const Matrix3& b)
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

Matrix3& operator+=(Matrix3& a, const float b)
{
	for (unsigned int i = 0; i < 9; ++i)
		a.mat_array[i] += b;
	return a;
}

Matrix3& operator-=(Matrix3& a, const float b)
{
	for (unsigned int i = 0; i < 9; ++i)
		a.mat_array[i] -= b;
	return a;
}
Matrix3& operator*=(Matrix3& a, const float b)
{
	for (unsigned int i = 0; i < 9; ++i)
		a.mat_array[i] *= b;
	return a;
}
Matrix3& operator/=(Matrix3& a, const float b)
{
	for (unsigned int i = 0; i < 9; ++i)
		a.mat_array[i] /= b;
	return a;
}

Matrix3 operator+(Matrix3& a, const float b)
{
	Matrix3 m;
	for (unsigned int i = 0; i < 9; ++i)
		m.mat_array[i] = a.mat_array[i] + b;
	return m;
}

Matrix3 operator-(const Matrix3& a, const float b)
{
	Matrix3 m;
	for (unsigned int i = 0; i < 9; ++i)
		m.mat_array[i] = a.mat_array[i] - b;
	return m;
}
Matrix3 operator*(const Matrix3& a, const float b)
{
	Matrix3 m;
	for (unsigned int i = 0; i < 9; ++i)
		m.mat_array[i] = a.mat_array[i] * b;
	return m;
}
Matrix3 operator/(const Matrix3& a, const float b)
{
	Matrix3 m;
	for (unsigned int i = 0; i < 9; ++i)
		m.mat_array[i] = a.mat_array[i] / b;
	return m;
}

Vector3 operator*(const Matrix3& a, const Vector3& b)
{
	Vector3 out;

	out.x = a._11 * b.x
	+ a._12 * b.y
	+ a._13 * b.z;

	out.y = a._21 * b.x
	+ a._22 * b.y
	+ a._23 * b.z;

	out.z = a._31 * b.x
	+ a._32 * b.y
	+ a._33 * b.z;

	return out;
}
