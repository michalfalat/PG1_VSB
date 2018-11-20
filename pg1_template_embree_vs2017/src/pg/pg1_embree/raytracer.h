#pragma once
#include "simpleguidx11.h"
#include "surface.h"
#include "camera.h"
#include "structs.h"
#include "Sphericalbackground.h"

/*! \class Raytracer
\brief General ray tracer class.

\author Tom� Fabi�n
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

	Color4f trace_ray(MyRTCRayHit ray, int depth,  const float t = 0.0f);

	float trace_shadow_ray(const Vector3 & p, const Vector3 & l_d, const float dist);
	float linearToSrgb(float color);
	Vector3 sampleHemisphere(Vector3 normal);

	Vector3 getInterpolatedPoint(RTCRay ray);
	int Ui();

private:
	std::vector<Surface *> surfaces_;
	std::vector<Material *> materials_;

	RTCDevice device_;
	RTCScene scene_;
	Camera camera_;
	SphericalBackground background_;
};
