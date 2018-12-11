#pragma once
#define IOR_AIR 1.000293f
#include "vector3.h"

struct Vertex3f { float x, y, z; }; // a single vertex position structure matching certain format

using Normal3f = Vertex3f; // a single vertex normal structure matching certain format

struct Coord2f { float u, v; }; // texture coord structure

struct Triangle3ui { unsigned int v0, v1, v2; }; // indicies of a single triangle, the struct must match certain format, e.g. RTC_FORMAT_UINT3

struct RTC_ALIGN( 16 ) Color4f
{
	struct { float r, g, b, a; }; // a = 1 means that the pixel is opaque
	Color4f(const float r, const float g, const float b, const float a) : r(r), g(g), b(b), a(a) { }
	Color4f(const float* v);

	friend Color4f operator*(const Color4f &c, const  float a);
	friend Color4f operator*(const float a, const  Color4f &c);
	friend Color4f operator*(const Vector3 &v, const  Color4f &c);
	friend Color4f operator*(const Color4f &u, const  Color4f &v);
	friend Color4f operator+(const Color4f &c1, const  Color4f &c2);
	friend Color4f operator+(const Color4f &c1, const  Vector3 &c2);
	friend Color4f operator/(const Color4f &c, const  float a);
	friend void operator+=(Color4f &c1, const Color4f &c2);
};

struct Color3f { float r, g, b; };

struct RTCRayHitWithIor {
	RTCRayHit ray_hit;
	float ior = IOR_AIR;
};

inline void reorient_against(Normal3f & n, const float v_x, const float v_y, const float v_z) {
	if ((n.x * v_x + n.y * v_y + n.z * v_z) > 0.0f) {
		n.x *= -1;
		n.y *= -1;
		n.z *= -1;
	}
}


