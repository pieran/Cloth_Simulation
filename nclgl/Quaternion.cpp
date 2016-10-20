#include "Quaternion.h"

Quaternion::Quaternion(void)
{
	x = y = z = 0.0f;
	w = 1.0f;
}

Quaternion::Quaternion(float x, float y, float z, float w)
{
	this->x = x;
	this->y = y;
	this->z = z;
	this->w = w;
}

Quaternion::~Quaternion(void)
{
}

float Quaternion::Dot(const Quaternion &a, const Quaternion &b){
	return (a.x * b.x) + (a.y * b.y) + (a.z * b.z) + (a.w * b.w);
}

void Quaternion::Normalise(){
	float magnitude = sqrt(Dot(*this, *this));

	if (magnitude > 0.0f){
		float t = 1.0f / magnitude;

		x *= t;
		y *= t;
		z *= t;
		w *= t;
	}
}

Vector3 Quaternion::RotatePointByQuaternion(const Quaternion &q, const Vector3 &point)
{
	//Expansion of q . point . q' (where q' is the conjugate)
	Vector3 vector;
	float num12 = q.x + q.x;
	float num2 = q.y + q.y;
	float num = q.z + q.z;
	float num11 = q.w * num12;
	float num10 = q.w * num2;
	float num9 = q.w * num;
	float num8 = q.x * num12;
	float num7 = q.x * num2;
	float num6 = q.x * num;
	float num5 = q.y * num2;
	float num4 = q.y * num;
	float num3 = q.z * num;

	vector.x = ((point.x * ((1.0f - num5) - num3)) + (point.y * (num7 - num9))) + (point.z * (num6 + num10));
	vector.y = ((point.x * (num7 + num9)) + (point.y * ((1.0f - num8) - num3))) + (point.z * (num4 - num11));
	vector.z = ((point.x * (num6 - num10)) + (point.y * (num4 + num11))) + (point.z * ((1.0f - num8) - num5));

	/*vector.x =
		point.x * (1.0f - 2.0f * q.y * q.y - 2.0f * q.z * q.z) +
		point.y * (2.0f * q.x * q.y - 2.0f * q.w * q.z) +
		point.z * (2.0f * q.x * q.z + 2.0f * q.w * q.y);

	vector.y =
		point.x * (2.0f * q.x * q.y + 2.0f * q.w * q.z) +
		point.y * (1.0f - 2.0f * q.x * q.x - 2.0f * q.z * q.z) +
		point.z * (2.0f * q.y * q.z - 2.0f * q.w * q.x);

	vector.z =
		point.x * (2.0f * q.x * q.z - 2.0f * q.w * q.y) +
		point.y * (2.0f * q.y * q.z + 2.0f * q.w * q.x) +
		point.z * (1.0f - 2.0f * q.x * q.x - 2.0f * q.y * q.y);


	vector.x = point.x*(q.x*q.x + q.w*q.w - q.y*q.y - q.z*q.z) + point.y*(2 * q.x*q.y - 2 * q.w*q.z) + point.z*(2 * q.x*q.z + 2 * q.w*q.y);
	vector.y = point.x*(2 * q.w*q.z + 2 * q.x*q.y) + point.y*(q.w*q.w - q.x*q.x + q.y*q.y - q.z*q.z) + point.z*(-2 * q.w*q.x + 2 * q.y*q.z);
	vector.z = point.x*(-2 * q.w*q.y + 2 * q.x*q.z) + point.y*(2 * q.w*q.x + 2 * q.y*q.z) + point.z*(q.w*q.w - q.x*q.x - q.y*q.y + q.z*q.z);*/

	return vector;

}


Quaternion Quaternion::operator *(const Quaternion &b) const{
	Quaternion ans;

	ans.w = (w * b.w) - (x * b.x) - (y * b.y) - (z * b.z);
	ans.x = (x * b.w) + (w * b.x) + (y * b.z) - (z * b.y);
	ans.y = (y * b.w) + (w * b.y) + (z * b.x) - (x * b.z);
	ans.z = (z * b.w) + (w * b.z) + (x * b.y) - (y * b.x);

	return ans;
}

Quaternion Quaternion::operator *(const Vector3 &b) const{
	/*Quaternion ans;

	ans.w = -(x * b.x) - (y * b.y) - (z * b.z);

	ans.x =  (w * b.x) + (y * b.z) - (z * b.y);
	ans.y =  (w * b.y) + (z * b.x) - (x * b.z);
	ans.z =  (w * b.z) + (x * b.y) - (y * b.x);

	return ans;*/

	//This (^) is equiv to q * b, where the below is equiv to b * q (needed for physics)

	Quaternion ans;

	ans.w = -(x * b.x) - (y * b.y) - (z * b.z);

	ans.x = (w * b.x) + (b.y * z) - (b.z * y);
	ans.y = (w * b.y) + (b.z * x) - (b.x * z);
	ans.z = (w * b.z) + (b.x * y) - (b.y * x);

	return ans;
}

Matrix4 Quaternion::ToMatrix4() const{
	Matrix4 mat;

	float yy = y*y;
	float zz = z*z;
	float xy = x*y;
	float zw = z*w;
	float xz = x*z;
	float yw = y*w;
	float xx = x*x;
	float yz = y*z;
	float xw = x*w;

	mat.values[0] = 1 - 2 * yy - 2 * zz;
	mat.values[1] = 2 * xy + 2 * zw;
	mat.values[2] = 2 * xz - 2 * yw;

	mat.values[4] = 2 * xy - 2 * zw;
	mat.values[5] = 1 - 2 * xx - 2 * zz;
	mat.values[6] = 2 * yz + 2 * xw;

	mat.values[8] = 2 * xz + 2 * yw;
	mat.values[9] = 2 * yz - 2 * xw;
	mat.values[10] = 1 - 2 * xx - 2 * yy;

	return mat;
}

Matrix3 Quaternion::ToMatrix3() const{
	Matrix3 mat;

	float yy = y*y;
	float zz = z*z;
	float xy = x*y;
	float zw = z*w;
	float xz = x*z;
	float yw = y*w;
	float xx = x*x;
	float yz = y*z;
	float xw = x*w;

	mat.mat_array[0] = 1 - 2 * yy - 2 * zz;
	mat.mat_array[1] = 2 * xy + 2 * zw;
	mat.mat_array[2] = 2 * xz - 2 * yw;

	mat.mat_array[3] = 2 * xy - 2 * zw;
	mat.mat_array[4] = 1 - 2 * xx - 2 * zz;
	mat.mat_array[5] = 2 * yz + 2 * xw;

	mat.mat_array[6] = 2 * xz + 2 * yw;
	mat.mat_array[7] = 2 * yz - 2 * xw;
	mat.mat_array[8] = 1 - 2 * xx - 2 * yy;

	return mat;
}

Quaternion Quaternion::EulerAnglesToQuaternion(float pitch, float yaw, float roll)	{
	float y2 = (float)DegToRad(yaw / 2.0f);
	float p2 = (float)DegToRad(pitch / 2.0f);
	float r2 = (float)DegToRad(roll / 2.0f);


	float cosy = (float)cos(y2);
	float cosp = (float)cos(p2);
	float cosr = (float)cos(r2);

	float siny = (float)sin(y2);
	float sinp = (float)sin(p2);
	float sinr = (float)sin(r2);

	Quaternion q;


	q.x = cosr * sinp * cosy + sinr * cosp * siny;
	q.y = cosr * cosp * siny - sinr * sinp * cosy;
	q.z = sinr * cosp * cosy - cosr * sinp * siny;
	q.w = cosr * cosp * cosy + sinr * sinp * siny;

	return q;
};

Quaternion Quaternion::AxisAngleToQuaterion(const Vector3& vector, float degrees)	{
	float theta = (float)DegToRad(degrees);
	float result = (float)sin(theta / 2.0f);

	return Quaternion((float)(vector.x * result), (float)(vector.y * result), (float)(vector.z * result), (float)cos(theta / 2.0f));
}

void Quaternion::GenerateW()	{
	w = 1.0f - (x*x) - (y*y) - (z*z);
	if (w < 0.0f) {
		w = 0.0f;
	}
	else{
		w = -sqrt(w);
	}
}

Quaternion Quaternion::Conjugate() const
{
	return Quaternion(-x, -y, -z, w);
}

Quaternion Quaternion::FromMatrix(const Matrix4 &m)	{
	Quaternion q;

	q.w = sqrt(max(0.0f, (1.0f + m.values[0] + m.values[5] + m.values[10]))) / 2;
	q.x = sqrt(max(0.0f, (1.0f + m.values[0] - m.values[5] - m.values[10]))) / 2;
	q.y = sqrt(max(0.0f, (1.0f - m.values[0] + m.values[5] - m.values[10]))) / 2;
	q.z = sqrt(max(0.0f, (1.0f - m.values[0] - m.values[5] + m.values[10]))) / 2;

	q.x = (float)_copysign(q.x, m.values[9] - m.values[6]);
	q.y = (float)_copysign(q.y, m.values[2] - m.values[8]);
	q.z = (float)_copysign(q.z, m.values[4] - m.values[1]);

	return q;
}

Quaternion Quaternion::FromMatrix(const Matrix3 &m) {
	Quaternion q;

	q.w = sqrt(max(0.0f, (1.0f + m._11 + m._22 + m._33))) / 2.0;
	q.x = sqrt(max(0.0f, (1.0f + m._11 - m._22 - m._33))) / 2.0;
	q.y = sqrt(max(0.0f, (1.0f - m._11 + m._22 - m._33))) / 2.0;
	q.z = sqrt(max(0.0f, (1.0f - m._11 - m._22 + m._33))) / 2.0;

	q.x = (float)_copysign(q.x, m._23 - m._32);
	q.y = (float)_copysign(q.y, m._31 - m._13);
	q.z = (float)_copysign(q.z, m._12 - m._21);

	return q;
}

Quaternion Quaternion::Interpolate(const Quaternion& pStart, const Quaternion& pEnd, float pFactor)
{
	// calc cosine theta
	float cosom = pStart.x * pEnd.x + pStart.y * pEnd.y + pStart.z * pEnd.z + pStart.w * pEnd.w;

	// adjust signs (if necessary)
	Quaternion end = pEnd;
	if (cosom < 0.0f)
	{
		cosom = -cosom;
		end.x = -end.x;   // Reverse all signs
		end.y = -end.y;
		end.z = -end.z;
		end.w = -end.w;
	}

	// Calculate coefficients
	float sclp, sclq;
	if ((1.0f - cosom) > 0.0001f) // 0.0001 -> some epsillon
	{
		// Standard case (slerp)
		float omega, sinom;
		omega = acos(cosom); // extract theta from dot product's cos theta
		sinom = sin(omega);
		sclp = sin((1.0f - pFactor) * omega) / sinom;
		sclq = sin(pFactor * omega) / sinom;
	}
	else
	{
		// Very close, do linear interp (because it's faster)
		sclp = 1.0f - pFactor;
		sclq = pFactor;
	}

	Quaternion out;
	out.x = sclp * pStart.x + sclq * end.x;
	out.y = sclp * pStart.y + sclq * end.y;
	out.z = sclp * pStart.z + sclq * end.z;
	out.w = sclp * pStart.w + sclq * end.w;

	return out;
}