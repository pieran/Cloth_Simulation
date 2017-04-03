#include "utils.h"

Vector3 hsv2rgb(const Vector3& c) {
	const Vector3 K = Vector3(1.0, 2.0 / 3.0, 1.0 / 3.0);
	float intpart;
	Vector3 p = Vector3(
		abs(modf(c.x + K.x, &intpart) * 6.0 - 3.0),
		abs(modf(c.x + K.y, &intpart) * 6.0 - 3.0),
		abs(modf(c.x + K.z, &intpart) * 6.0 - 3.0)
		);
	return Vector3::InterpolateLinear(Vector3(1, 1, 1), Vector3::ClampVars(p - Vector3(1, 1, 1), 0.0, 1.0), c.y) * c.z;
}