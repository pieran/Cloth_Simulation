#pragma once
#include "common.h"

class GLColor
{
public:
	static const struct PredefinedColors
	{
		static const GLColor White;
		static const GLColor Black;
		static const GLColor Grey;
		static const GLColor Red;
		static const GLColor Green;
		static const GLColor Blue;
	};



	static GLColor FromRGBA(uint rgba);
	static GLColor FromRGBA(char red, char green, char blue, char alpha);

	static GLColor FromHSV(float hue, float saturation, float value, float alpha = 1.0f);





public:
	GLColor() : r(0.0f), g(0.0f), b(0.0f), a(1.0f) {}
	GLColor(float grey) : r(grey), g(grey), b(grey), a(1.0f) {}
	GLColor(float red, float green, float blue) : r(red), g(green), b(blue), a(1.0f) {}
	GLColor(float red, float green, float blue, float alpha) : r(red), g(green), b(blue), a(alpha) {}
	~GLColor() {}
	

	union
	{
		float col_array[4];
		struct {
			float r;
			float g;
			float b;
			float a;
		};
	};
};