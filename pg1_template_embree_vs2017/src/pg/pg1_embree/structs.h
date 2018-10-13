#pragma once
#define IOR_AIR 1.000293f

struct Vertex3f { float x, y, z; }; // a single vertex position structure matching certain format

using Normal3f = Vertex3f; // a single vertex normal structure matching certain format

struct Coord2f { float u, v; }; // texture coord structure

struct Triangle3ui { unsigned int v0, v1, v2; }; // indicies of a single triangle, the struct must match certain format, e.g. RTC_FORMAT_UINT3

struct RTC_ALIGN( 16 ) Color4f
{
	struct { float r, g, b, a; }; // a = 1 means that the pixel is opaque
};

struct Color3f { float r, g, b; };

struct MyRTCRayHit {
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