#include "stdafx.h"
#include "SphericalBackground.h"
#define _USE_MATH_DEFINES
#include <math.h>



SphericalBackground::SphericalBackground(const char filename[])
{
	texture = new Texture(filename);
}



Color3f SphericalBackground::GetBackground(const float x, const float y, const float z) {
	const float theta = acosf(z);
	const float phi = atan2f(y, x) + float(M_PI);

	const float u = 1.0f - phi * 0.5f * float(M_1_PI);
	const float v = theta  * float(M_1_PI);

	return  texture->get_texel(u, v);
}


SphericalBackground::~SphericalBackground()
{
}
