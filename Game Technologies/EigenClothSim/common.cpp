#include "common.h"

void hsv2rgb(Vector3& rgb, const Vector3& hsv)
{
	float      hh, p, q, t, ff;
	long        i;

	hh = hsv.x;
	if (hh >= 360.0f) hh = 0.0f;
	hh /= 60.0f;
	i = (long)hh;
	ff = hh - i;
	p = hsv.z * (1.0f - hsv.y);
	q = hsv.z * (1.0f - (hsv.y * ff));
	t = hsv.z * (1.0f - (hsv.y * (1.0f - ff)));

	switch (i) {
	case 0:
		rgb.x = hsv.z;
		rgb.y = t;
		rgb.z = p;
		break;
	case 1:
		rgb.x = q;
		rgb.y = hsv.z;
		rgb.z = p;
		break;
	case 2:
		rgb.x = p;
		rgb.y = hsv.z;
		rgb.z = t;
		break;

	case 3:
		rgb.x = p;
		rgb.y = q;
		rgb.z = hsv.z;
		break;
	case 4:
		rgb.x = t;
		rgb.y = p;
		rgb.z = hsv.z;
		break;
	case 5:
	default:
		rgb.x = hsv.z;
		rgb.y = p;
		rgb.z = q;
		break;
	}
}