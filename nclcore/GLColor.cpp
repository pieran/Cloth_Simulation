#include "GLColor.h"



const GLColor GLColor::PredefinedColors::White	= GLColor(1.0f, 1.0f, 1.0f, 1.0f);
const GLColor GLColor::PredefinedColors::Black	= GLColor(0.0f, 0.0f, 0.0f, 1.0f);
const GLColor GLColor::PredefinedColors::Grey	= GLColor(0.5f, 0.5f, 0.5f, 1.0f);
const GLColor GLColor::PredefinedColors::Red	= GLColor(1.0f, 0.0f, 0.0f, 1.0f);
const GLColor GLColor::PredefinedColors::Green	= GLColor(0.0f, 1.0f, 0.0f, 1.0f);
const GLColor GLColor::PredefinedColors::Blue	= GLColor(0.0f, 0.0f, 1.0f, 1.0f);







GLColor GLColor::FromRGBA(uint rgba)
{
	char r   = (rgba >> 24) & 0xFF;
	char g = (rgba >> 16) & 0xFF;
	char b  = (rgba >> 8) & 0xFF;
	char a = rgba & 0xFF;

	return FromRGBA(r, g, b, a);
}

GLColor GLColor::FromRGBA(char red, char green, char blue, char alpha)
{
	return GLColor(
		((float)red) / 255.0f,
		((float)green) / 255.0f,
		((float)blue) / 255.0f,
		((float)alpha) / 255.0f);
}

//Taken from http://stackoverflow.com/questions/3018313/algorithm-to-convert-rgb-to-hsv-and-hsv-to-rgb-in-range-0-255-for-both
GLColor GLColor::FromHSV(float hue, float saturation, float value, float alpha)
{
	float      hh, p, q, t, ff;
	int        i;
	GLColor out(value, value, value, alpha);

	if (saturation <= 0.0) {
		return out;
	}

	hh = hue;
	if (hh >= 360.0f) hh = 0.0f;
	hh /= 60.0f;

	i = (int)hh;
	ff = hh - i;

	p = value * (1.0f - saturation);
	q = value * (1.0f - (saturation * ff));
	t = value * (1.0f - (saturation * (1.0f - ff)));

	switch (i) {
	case 0:
		out.r = value;
		out.g = t;
		out.b = p;
		break;
	case 1:
		out.r = q;
		out.g = value;
		out.b = p;
		break;
	case 2:
		out.r = p;
		out.g = value;
		out.b = t;
		break;

	case 3:
		out.r = p;
		out.g = q;
		out.b = value;
		break;
	case 4:
		out.r = t;
		out.g = p;
		out.b = value;
		break;
	case 5:
	default:
		out.r = value;
		out.g = p;
		out.b = q;
		break;
	}
	return out;
}