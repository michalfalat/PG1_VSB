#include "stdafx.h"
#include "structs.h"
#include "vector3.h"

Color4f::Color4f(const float * v)
{
	assert(v != NULL);

	r = v[0];
	g = v[1];
	b = v[2];
	a = v[4];
}

Color4f operator*(const Color4f &c, const  float a) {
	return Color4f{ a * c.r, a * c.g, a * c.b, a * c.a };
}

Color4f operator*(const float a, const  Color4f &c) {
	return Color4f{ a * c.r, a * c.g, a * c.b, a * c.a };
}

Color4f operator*(const Vector3 &v, const  Color4f &c) {
	return Color4f{ v.x * c.r, v.y * c.g, v.z * c.b, c.a };
}

Color4f operator*(const Color4f &u, const  Color4f &v) {
	return Color4f{ u.r * v.r, u.g * v.g, u.b * v.b, u.a * v.a };
}

Color4f operator+(const Color4f &c1, const  Color4f &c2) {
	return Color4f{ c1.r + c2.r, c1.g + c2.g, c1.b + c2.b, c1.a + c2.a };
}

Color4f operator+(const Color4f &c1, const  Vector3 &c2) {
	return Color4f{ c1.r + c2.x, c1.g + c2.y, c1.b + c2.z, c1.a  };
}
Color4f operator/(const Color4f &c, const  float a) {
	return Color4f{ c.r / a, c.g / a, c.b / a, c.a / a };
}

void operator+=(Color4f &c1, const Color4f &c2) {
	c1.r += c2.r;
	c1.g += c2.g;
	c1.b += c2.b;
	c1.a += c2.a;
}