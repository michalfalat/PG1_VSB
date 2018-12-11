#pragma once
#include "simpleguidx11.h"
#include "surface.h"
#include "camera.h"
#include "structs.h"
#include "Background.h"

/*! \class Raytracer
\brief General ray tracer class.

\author Tomáš Fabián
\version 0.1
\date 2018
*/
class Raytracer : public SimpleGuiDX11
{
public:
	Raytracer( const int width, const int height, 
		const float fov_y, const Vector3 view_from, const Vector3 view_at,
		const char * config = "threads=0,verbose=3" );
	~Raytracer();

	int InitDeviceAndScene( const char * config );

	int ReleaseDeviceAndScene();

	void LoadScene( const std::string file_name );

	Color4f get_pixel( const int x, const int y, const float t = 0.0f ) override;

	Color4f trace_ray(RTCRayHitWithIor ray, int depth);

	float trace_shadow_ray(const Vector3 & p, const Vector3 & l_d, const float dist, RTCIntersectContext context);
	float linearToSrgb(float color);
	float getGeometryTerm(Vector3 omegaI, RTCIntersectContext context, Vector3 vectorToLight, Vector3 intersectionPoint, Vector3 normal);
	float  castShadowRay(RTCIntersectContext context, Vector3 vectorToLight, float dstToLight, Vector3 intersectionPoint, Vector3 normal);

	Vector3 sampleHemisphere(Vector3 normal);

	Vector3 getInterpolatedPoint(RTCRay ray);
	int Ui();

private:
	std::vector<Surface *> surfaces_;
	std::vector<Material *> materials_;

	RTCDevice device_;
	RTCScene scene_;
	Camera camera_;
	Background background_;
};
