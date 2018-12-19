#include "stdafx.h"
#include "raytracer.h"
#include "objloader.h"
#include "tutorials.h"
#include "material.h"
#include "background.h"
#include "utils.h"
#define _USE_MATH_DEFINES
#include <math.h>


float SQR(float r) {
	return r * r;
}

inline Vector3 reflect(const Vector3 & v, const Vector3 & n) {
	return 2.0f * (n.DotProduct(v))* n - v;
}

Raytracer::Raytracer(const int width, const int height,
	const float fov_y, const Vector3 view_from, const Vector3 view_at,
	const char * config) : SimpleGuiDX11(width, height)
{
	InitDeviceAndScene(config);

	camera_ = Camera(width, height, fov_y, view_from, view_at);
	background_ = Background("../../../data/background.jpg");
}




Raytracer::~Raytracer()
{
	ReleaseDeviceAndScene();
}

int Raytracer::InitDeviceAndScene(const char * config)
{
	device_ = rtcNewDevice(config);
	error_handler(nullptr, rtcGetDeviceError(device_), "Unable to create a new device.\n");
	rtcSetDeviceErrorFunction(device_, error_handler, nullptr);

	ssize_t triangle_supported = rtcGetDeviceProperty(device_, RTC_DEVICE_PROPERTY_TRIANGLE_GEOMETRY_SUPPORTED);

	// create a new scene bound to the specified device
	scene_ = rtcNewScene(device_);

	return S_OK;
}

int Raytracer::ReleaseDeviceAndScene()
{
	rtcReleaseScene(scene_);
	rtcReleaseDevice(device_);

	return S_OK;
}


float Raytracer::trace_shadow_ray(const Vector3 & p, const Vector3 & l_d, const float dist, RTCIntersectContext contedxt) {
	RTCHit hit;
	hit.geomID = RTC_INVALID_GEOMETRY_ID;
	hit.primID = RTC_INVALID_GEOMETRY_ID;

	RTCRay ray = RTCRay();
	ray.org_x = p.x;
	ray.org_y = p.y;
	ray.org_z = p.z;

	ray.dir_x = l_d.x;
	ray.dir_y = l_d.y;
	ray.dir_z = l_d.z;

	ray.tnear = 0.1f;
	ray.tfar = dist;

	ray.time = 0.0f;

	ray.mask = 0; // can be used to mask out some geometries for some rays
	ray.id = 0; // identify a ray inside a callback function
	ray.flags = 0; // reserved

	RTCIntersectContext context;
	rtcInitIntersectContext(&context);
	rtcOccluded1(scene_, &context, &ray);

	if (ray.tfar < dist) {
		return 0.00f;
	}
	else {
		return 1.0f;
	}
}

void Raytracer::LoadScene(const std::string file_name)
{
	const int no_surfaces = LoadOBJ(file_name.c_str(), surfaces_, materials_);

	// surfaces loop
	for (auto surface : surfaces_)
	{
		RTCGeometry mesh = rtcNewGeometry(device_, RTC_GEOMETRY_TYPE_TRIANGLE);

		Vertex3f * vertices = (Vertex3f *)rtcSetNewGeometryBuffer(
			mesh, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3,
			sizeof(Vertex3f), 3 * surface->no_triangles());

		Triangle3ui * triangles = (Triangle3ui *)rtcSetNewGeometryBuffer(
			mesh, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3,
			sizeof(Triangle3ui), surface->no_triangles());

		rtcSetGeometryUserData(mesh, (void*)(surface->get_material()));

		rtcSetGeometryVertexAttributeCount(mesh, 2);

		Normal3f * normals = (Normal3f *)rtcSetNewGeometryBuffer(
			mesh, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, RTC_FORMAT_FLOAT3,
			sizeof(Normal3f), 3 * surface->no_triangles());

		Coord2f * tex_coords = (Coord2f *)rtcSetNewGeometryBuffer(
			mesh, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 1, RTC_FORMAT_FLOAT2,
			sizeof(Coord2f), 3 * surface->no_triangles());

		// triangles loop
		for (int i = 0, k = 0; i < surface->no_triangles(); ++i)
		{
			Triangle & triangle = surface->get_triangle(i);

			// vertices loop
			for (int j = 0; j < 3; ++j, ++k)
			{
				const Vertex & vertex = triangle.vertex(j);

				vertices[k].x = vertex.position.x;
				vertices[k].y = vertex.position.y;
				vertices[k].z = vertex.position.z;

				normals[k].x = vertex.normal.x;
				normals[k].y = vertex.normal.y;
				normals[k].z = vertex.normal.z;

				tex_coords[k].u = vertex.texture_coords[0].u;
				tex_coords[k].v = vertex.texture_coords[0].v;
			}

			triangles[i].v0 = k - 3;
			triangles[i].v1 = k - 2;
			triangles[i].v2 = k - 1;
		}

		rtcCommitGeometry(mesh);
		unsigned int geom_id = rtcAttachGeometry(scene_, mesh);
		rtcReleaseGeometry(mesh);
	} // end of surfaces loop

	rtcCommitScene(scene_);
}

Color4f Raytracer::get_pixel(const int x, const int y, const float t)
{

	RTCRayHitWithIor my_ray_hit;

	// Uniform supersampling
	Color4f colorSum = Color4f(0.0f, 0.0f, 0.0f, 1.0f);
	int size = 1;

	float offsetX = -0.5f;
	float offsetY = -0.5f;
	float offsetAddition = 1.0f / size;

	/*int samples = 15;

	#pragma omp parallel for schedule(dynamic, 5) shared(colorSum)
		for (int i = 0; i < samples; i++) {
			float offsetX = (x + Random());
			float offsetY = (y + Random());
	
			my_ray_hit.ray_hit.ray = camera_.GenerateRay(offsetX, offsetY);
			my_ray_hit.ray_hit.hit = createEmptyHit();
			my_ray_hit.ior = IOR_AIR;
			Color4f traced = trace_ray(my_ray_hit, 5);
			colorSum += traced;
		}
	
		return colorSum / static_cast<float>(samples);*/
	offsetX = (x + Random());
	offsetY = (y + Random());

	my_ray_hit.ray_hit.ray = camera_.GenerateRay(offsetX, offsetY);
	my_ray_hit.ray_hit.hit = createEmptyHit();
	my_ray_hit.ior = IOR_AIR;
	Color4f traced = trace_ray(my_ray_hit, 4);

	return traced;
}

Color4f Raytracer::trace_ray(RTCRayHitWithIor my_ray_hit, int depth) {
	// TODO generate primary ray and perform ray cast on the scene
	// setup a hit

	// intersect ray with the scene
	RTCIntersectContext context;
	rtcInitIntersectContext(&context);
	rtcIntersect1(scene_, &context, &my_ray_hit.ray_hit);

	if (my_ray_hit.ray_hit.hit.geomID != RTC_INVALID_GEOMETRY_ID)
	{
		// we hit something
		RTCGeometry geometry = rtcGetGeometry(scene_, my_ray_hit.ray_hit.hit.geomID);
		Normal3f normal;
		// get interpolated normal
		rtcInterpolate0(geometry, my_ray_hit.ray_hit.hit.primID, my_ray_hit.ray_hit.hit.u, my_ray_hit.ray_hit.hit.v,
			RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, &normal.x, 3);

		//reorient_against(normal, my_ray_hit.ray_hit.ray.dir_x, my_ray_hit.ray_hit.ray.dir_y, my_ray_hit.ray_hit.ray.dir_z);

		// and texture coordinates
		Coord2f tex_coord;
		rtcInterpolate0(geometry, my_ray_hit.ray_hit.hit.primID, my_ray_hit.ray_hit.hit.u, my_ray_hit.ray_hit.hit.v,
			RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 1, &tex_coord.u, 2);

		tex_coord.v = 1.0f - tex_coord.v;

		Material * material = (Material *)(rtcGetGeometryUserData(geometry));

		//const Triangle & triangle = surfaces_[ray_hit]
		Vector3 l_position = Vector3(50, -50, 300);

		Vector3 p = getInterpolatedPoint(my_ray_hit.ray_hit.ray);
		Vector3 l_d = l_position - p;
		l_d.Normalize();

		Vector3 rd = Vector3(my_ray_hit.ray_hit.ray.dir_x, my_ray_hit.ray_hit.ray.dir_y, my_ray_hit.ray_hit.ray.dir_z);
		Vector3 normal_v = Vector3(normal.x, normal.y, normal.z);
		if (rd.DotProduct(normal_v) > 0) {
			normal_v *= -1;
		}

		if (depth <= 0) {
			Color4f background = changeGamma(background_.GetBackground(my_ray_hit.ray_hit.ray.dir_x, my_ray_hit.ray_hit.ray.dir_y, my_ray_hit.ray_hit.ray.dir_z));



			//return bck;
			//return Color4f(1.0f,bck.g, bck.b,1.0f);
			return Color4f(getSRGBColorValueForComponent(background.r), getSRGBColorValueForComponent(background.g), getSRGBColorValueForComponent(background.b), 1.0f);
		}
		// u,v tex_coord

		switch (material->shader_)
		{
		case Shader::NORMAL:
		{
			return Color4f(normal.x * 0.5f + 0.5f, normal.y * 0.5f + 0.5f, normal.z * 0.5f + 0.5f, 1.0f);
			break;
		}
		case Shader::LAMBERT:
		{
			Vector3 diffuse = material->doDiffuse(&tex_coord);
			float dot = l_d.DotProduct(normal_v);
			Vector3 lambert_color = max(0, dot) * diffuse;
			return Color4f(lambert_color.x, lambert_color.y, lambert_color.z, 1.0f);
			break;
		}
		case Shader::PHONG:
		{
			Vector3 v = Vector3(-my_ray_hit.ray_hit.ray.dir_x, -my_ray_hit.ray_hit.ray.dir_y, -my_ray_hit.ray_hit.ray.dir_z);
			Vector3 l_r = reflect(l_d, normal_v);

			float normal_dotProduct_l_d = normal_v.DotProduct(l_d);
			// get diffuse
			Vector3 diffuse = material->doDiffuse(&tex_coord);

			const float enlight = trace_shadow_ray(p, l_d, l_d.L2Norm(), context);
			Color4f final_color = Color4f{
				(material->ambient.x + enlight * ((diffuse.x * normal_dotProduct_l_d) + pow(material->specular.x * v.DotProduct(l_r), material->shininess))),
				(material->ambient.y + enlight * ((diffuse.y * normal_dotProduct_l_d) + pow(material->specular.y * v.DotProduct(l_r), material->shininess))),
				(material->ambient.z + enlight * ((diffuse.z * normal_dotProduct_l_d) + pow(material->specular.z * v.DotProduct(l_r), material->shininess))),
				1 } * material->reflectivity;
			return final_color;

			break;
		}
		case Shader::GLASS:
		{
			RTCRayHitWithIor refracted_ray_hit, reflected_ray_hit;

			Vector3 diffuse = material->diffuse;
			Vector3 rv = Vector3(-rd.x, -rd.y, -rd.z);
			Vector3 vector = getInterpolatedPoint(my_ray_hit.ray_hit.ray);

			float n1 = my_ray_hit.ior;
			float n2 = ((n1 == IOR_AIR) ? material->ior : IOR_AIR);

			float n_divided = n1 / n2;
			float cos_01 = (normal_v.DotProduct(rv));

			Vector3 rr = (2.0f * (normal_v.DotProduct(rv))) * normal_v - rv;

			float refractComponent = 1.0f - SQR(n_divided) * (1.0f - SQR(cos_01));

			
			reflected_ray_hit = createRayWithEmptyHitAndIor(vector, rr, FLT_MAX, 0.001f, n2);

			if (refractComponent > 0) {
				float cos_02 = sqrt(refractComponent);
				Vector3 rl = (n_divided * rd) + ((n_divided * cos_01 - cos_02) * normal_v);
				// Fresnel
				float Rs = SQR((n2 * cos_02 - n1 * cos_01) / (n2 * cos_02 + n1 * cos_01));
				float Rp = SQR((n2 * cos_01 - n1 * cos_02) / (n2 * cos_01 + n1 * cos_02));
				float part_reflect = 0.5f * (Rs + Rp);

				// Calculate coefficients
				float part_refract = 1.0f - part_reflect;
							   
				// refracted ray
				refracted_ray_hit = createRayWithEmptyHitAndIor(vector, rl, FLT_MAX, 0.001f, n2);

				return (diffuse * trace_ray(reflected_ray_hit, depth - 1) * part_reflect) + (diffuse * trace_ray(refracted_ray_hit, depth - 1) * part_refract);

				//FOR DEBUG
				//return diffuse * trace_ray(myReflectedRTCRayHit, depth - 1);
				//return diffuse * trace_ray(myRefractedRTCRayHit, depth - 1) * coefRefract;
			}
			else {
				return diffuse * trace_ray(reflected_ray_hit, depth - 1);
			}

		}
		case Shader::PATHTRACER:
		{
			Color4f emmision = Color4f{ material->emission.x, material->emission.y, material->emission.z, 1 };
			if (emmision.r != 0 && emmision.g != 0 && emmision.b != 0) {
				return emmision;
			}

			Vector3 r_d = Vector3(my_ray_hit.ray_hit.ray.dir_x, my_ray_hit.ray_hit.ray.dir_y, my_ray_hit.ray_hit.ray.dir_z);
			Vector3 r_v = Vector3(-r_d.x, -r_d.y, -r_d.z);


			Vector3 omegaI = sampleHemisphere(normal_v);
			float pdf = 1 / (2 * M_PI);

			Color4f l_i = trace_ray(createRayWithEmptyHitAndIor(getInterpolatedPoint(my_ray_hit.ray_hit.ray), omegaI, FLT_MAX, 0.001f, IOR_AIR), depth - 1);
			Vector3 fR = material->diffuse / M_PI;

			Color4f final_color = fR * l_i * (normal_v.DotProduct(omegaI) / pdf);

			return final_color;
			break;
		}
		case Shader::MIRROR:
		{
			RTCRayHitWithIor reflected_ray_hit;

			Vector3 diffuse = material->diffuse;
			Vector3 rd = Vector3(my_ray_hit.ray_hit.ray.dir_x, my_ray_hit.ray_hit.ray.dir_y, my_ray_hit.ray_hit.ray.dir_z);
			Vector3 rv = Vector3(-rd.x, -rd.y, -rd.z);

			float n1 = my_ray_hit.ior;
			float n2 = ((n1 == IOR_AIR) ? material->ior : IOR_AIR);

			Vector3 rr = (2.0f * (normal_v.DotProduct(rv))) * normal_v - rv;
			Vector3 vector = getInterpolatedPoint(my_ray_hit.ray_hit.ray);


			reflected_ray_hit = createRayWithEmptyHitAndIor(vector, rr, FLT_MAX, 0.001f, n2);

			return diffuse * trace_ray(reflected_ray_hit, depth - 1);
		}
		case Shader::CLEAR_GLASS:
		{
			RTCRayHitWithIor refracted_ray_hit;

			Vector3 diffuse = material->diffuse;
			Vector3 rv = Vector3(-rd.x, -rd.y, -rd.z);
			Vector3 vector = getInterpolatedPoint(my_ray_hit.ray_hit.ray);

			float n1 = my_ray_hit.ior;
			float n2 = ((n1 == IOR_AIR) ? material->ior : IOR_AIR);

			float n_divided = n1 / n2;
			float cos_01 = (normal_v.DotProduct(rv));

			Vector3 rr = (2.0f * (normal_v.DotProduct(rv))) * normal_v - rv;

			float refractComponent = 1.0f - SQR(n_divided) * (1.0f - SQR(cos_01));
			if (refractComponent > 0) {
				float cos_02 = sqrt(refractComponent);
				Vector3 rl = (n_divided * rd) + ((n_divided * cos_01 - cos_02) * normal_v);
				// Fresnel
				float Rs = SQR((n2 * cos_02 - n1 * cos_01) / (n2 * cos_02 + n1 * cos_01));
				float Rp = SQR((n2 * cos_01 - n1 * cos_02) / (n2 * cos_01 + n1 * cos_02));
				float part_reflect = 0.5f * (Rs + Rp);

				// Calculate coefficients
				float part_refract = 1.0f - part_reflect;

				// Generate refracted ray
				refracted_ray_hit = createRayWithEmptyHitAndIor(vector, rl, FLT_MAX, 0.001f, n2);

				return (diffuse * trace_ray(refracted_ray_hit, depth - 1) );

				//FOR DEBUG
				//return diffuse * trace_ray(myReflectedRTCRayHit, depth - 1);
				//return diffuse * trace_ray(myRefractedRTCRayHit, depth - 1) * coefRefract;
			}
			else {
				Color4f background = background_.GetBackground(my_ray_hit.ray_hit.ray.dir_x, my_ray_hit.ray_hit.ray.dir_y, my_ray_hit.ray_hit.ray.dir_z);

				return Color4f(getSRGBColorValueForComponent(background.r), getSRGBColorValueForComponent(background.g), getSRGBColorValueForComponent(background.b), 1.0f);
			}
		}
		default:
		{
			Vector3 diff = material->doDiffuse(&tex_coord);
			float dot = l_d.DotProduct(normal_v);
			Vector3 temp = max(0, dot) * diff;
			return Color4f(temp.x, temp.y, temp.z, 1.0f);
			break;
		}
		}

		//return Color4f(0.0f, 0.0f, 0.0f, 1.0f);
	}
	else {
		Color4f background = changeGamma(background_.GetBackground(my_ray_hit.ray_hit.ray.dir_x, my_ray_hit.ray_hit.ray.dir_y, my_ray_hit.ray_hit.ray.dir_z));

		return Color4f(getSRGBColorValueForComponent(background.r), getSRGBColorValueForComponent(background.g), getSRGBColorValueForComponent(background.b), 1.0f);
	}
}



//float Raytracer::getGeometryTerm(Vector3 omegaI, RTCIntersectContext context, Vector3 vectorToLight, Vector3 intersectionPoint, Vector3 normal) {
//	float randomX = Random() * 100 - 50; //generate in interval -50 and 50
//	float randomY = Random() * 100 - 50; //generate in interval -50 and 50
//	Vector3 sourcePoint(randomX, randomY, 489); //489 heaight of source
//
//	RTCRayHitWithIor rayToSource = createRayWithEmptyHitAndIor(intersectionPoint, (sourcePoint - intersectionPoint).Normalize(), FLT_MAX, 0.1f, IOR_AIR);
//	rtcIntersect1(scene_, &context, &(rayToSource.ray_hit));
//
//	double geometryTerm = 0.0;
//
//	vectorToLight = (sourcePoint - intersectionPoint);
//	float dstToLight = (vectorToLight).L2Norm();
//	(vectorToLight).Normalize();
//
//	double visibilityTerm = castShadowRay(context, vectorToLight, dstToLight, intersectionPoint, normal);
//
//	if (visibilityTerm == 1)
//	{
//		Vector3 sourceNormal{ 0, 0, -1 };
//
//		geometryTerm = ((normal.DotProduct(omegaI) *  sourceNormal.DotProduct(-omegaI))
//			/ SQR((intersectionPoint - sourcePoint).L2Norm())) * 10000;
//	}
//	return geometryTerm;
//}

float  Raytracer::castShadowRay(RTCIntersectContext context, Vector3 vectorToLight, float dstToLight, Vector3 intersectionPoint, Vector3 normal) {
	RTCRay rayFromIntersectPointToLight = createRay(intersectionPoint, vectorToLight, dstToLight, 0.1f);
	rtcOccluded1(scene_, &context, &rayFromIntersectPointToLight);
	return rayFromIntersectPointToLight.tfar < dstToLight ? 0.0f : 1.0f;
}

Vector3 Raytracer::getInterpolatedPoint(RTCRay ray) {
	return Vector3{
			ray.org_x + ray.tfar * ray.dir_x,
			ray.org_y + ray.tfar * ray.dir_y,
			ray.org_z + ray.tfar * ray.dir_z
	};
}


Vector3 Raytracer::sampleHemisphere(Vector3 normal) {
	float randomU = Random();
	float randomV = Random();

	float x = 2 * cosf(2 * M_PI * randomU) * sqrt(randomV * (1 - randomV));
	float y = 2 * sinf(2 * M_PI * randomU) * sqrt(randomV * (1 - randomV));
	float z = 1 - 2 * randomV;

	Vector3 omegaI = Vector3{ x, y, z };
	omegaI.Normalize();

	if (omegaI.DotProduct(normal) < 0) {
		omegaI *= -1;
	}
	return omegaI;
}

//float Raytracer::linearToSrgb(float color) {
//
//	/*if (color <= 0.0f)
//		return 0.0f;
//	else if (color >= 1.0f)
//		return 1.0f;
//	else */if (color < 0.0031308f)
//		return color * 12.92f;
//	else
//		return std::powf(color, 1.0f / 2.4f) * 1.055f - 0.055f;
//
//}


int Raytracer::Ui()
{
	static float f = 0.0f;
	static int counter = 0;

	// we use a Begin/End pair to created a named window
	ImGui::Begin("Ray Tracer Params");

	ImGui::Text("Surfaces = %d", surfaces_.size());
	ImGui::Text("Materials = %d", materials_.size());
	ImGui::Separator();
	ImGui::Checkbox("Vsync", &vsync_);

	//ImGui::Checkbox( "Demo Window", &show_demo_window );      // Edit bools storing our window open/close state
	//ImGui::Checkbox( "Another Window", &show_another_window );

	ImGui::SliderFloat("float", &f, 0.0f, 1.0f);            // Edit 1 float using a slider from 0.0f to 1.0f    
	//ImGui::ColorEdit3( "clear color", ( float* )&clear_color ); // Edit 3 floats representing a color

	if (ImGui::Button("Button"))                            // Buttons return true when clicked (most widgets return true when edited/activated)
		counter++;
	ImGui::SameLine();
	ImGui::Text("counter = %d", counter);

	ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
	ImGui::End();

	// 3. Show another simple window.
	/*if ( show_another_window )
	{
	ImGui::Begin( "Another Window", &show_another_window );   // Pass a pointer to our bool variable (the window will have a closing button that will clear the bool when clicked)
	ImGui::Text( "Hello from another window!" );
	if ( ImGui::Button( "Close Me" ) )
	show_another_window = false;
	ImGui::End();
	}*/

	return 0;
}