#include "stdafx.h"
#include "raytracer.h"
#include "objloader.h"
#include "tutorials.h"
#include "material.h"
#include "SphericalBackground.h"

Raytracer::Raytracer(const int width, const int height,
	const float fov_y, const Vector3 view_from, const Vector3 view_at,
	const char * config) : SimpleGuiDX11(width, height)
{
	InitDeviceAndScene(config);

	camera_ = Camera(width, height, fov_y, view_from, view_at);
	background_ = SphericalBackground("../../../data/background.jpg");
}

Vector3 get_hit_point(const RTCRay &ray)
{
	return Vector3(
		ray.org_x + ray.tfar * ray.dir_x,
		ray.org_y + ray.tfar * ray.dir_y,
		ray.org_z + ray.tfar * ray.dir_z);
}


float SQR(float r) {
	return r * r;
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

inline Vector3 reflect(const Vector3 & v, const Vector3 & n) {
	return (2.0f*(v.DotProduct(n))) *n - v;
}

float Raytracer::trace_shadow_ray(const Vector3 & p, const Vector3 & l_d, const float dist) {
	RTCHit hit;
	hit.geomID = RTC_INVALID_GEOMETRY_ID;
	hit.primID = RTC_INVALID_GEOMETRY_ID;

	RTCRay ray = RTCRay();
	ray.org_x = p.x; // ray origin
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
		return 0.0;
	}
	else {
		return 1.0;
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
			} // end of vertices loop

			triangles[i].v0 = k - 3;
			triangles[i].v1 = k - 2;
			triangles[i].v2 = k - 1;
		} // end of triangles loop

		rtcCommitGeometry(mesh);
		unsigned int geom_id = rtcAttachGeometry(scene_, mesh);
		rtcReleaseGeometry(mesh);
	} // end of surfaces loop

	rtcCommitScene(scene_);
}

Color4f Raytracer::get_pixel(const int x, const int y, const float t)
{
	RTCHit hit;
	hit.geomID = RTC_INVALID_GEOMETRY_ID;
	hit.primID = RTC_INVALID_GEOMETRY_ID;
	hit.Ng_x = 0.0f;
	hit.Ng_y = 0.0f;
	hit.Ng_z = 0.0f;

	// merge ray and hit structures
	MyRTCRayHit my_ray_hit;
	my_ray_hit.ray_hit.hit = hit;

	// Uniform supersampling
	Color4f colorSum = Color4f(0.0f, 0.0f, 0.0f, 0.0f);
	int size = 1;

	float offsetX = -0.5f;
	float offsetY = -0.5f;
	float offsetAddition = 1.0f / size;

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			// Generate ray
			my_ray_hit.ray_hit.ray = camera_.GenerateRay(x + offsetX, y + offsetY);
			Color4f traced = trace_ray(my_ray_hit, 4, t);
			colorSum = colorSum + traced;

			// Increase offset
			offsetX += offsetAddition;
		}

		// Reset offset
		offsetX = -0.5f;
		// Increase row offset
		offsetY += offsetAddition;
	}
	return colorSum / SQR(size);
}

Color4f Raytracer::trace_ray(MyRTCRayHit my_ray_hit, int depth, const float t) {

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

		reorient_against(normal, my_ray_hit.ray_hit.ray.dir_x, my_ray_hit.ray_hit.ray.dir_y, my_ray_hit.ray_hit.ray.dir_z);


		// and texture coordinates
		Coord2f tex_coord;
		rtcInterpolate0(geometry, my_ray_hit.ray_hit.hit.primID, my_ray_hit.ray_hit.hit.u, my_ray_hit.ray_hit.hit.v,
			RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 1, &tex_coord.u, 2);

		/*printf("normal = (%0.3f, %0.3f, %0.3f)\n", normal.x, normal.y, normal.z);
		printf("tex_coord = (%0.3f, %0.3f)\n", tex_coord.u, tex_coord.v);*/
		Material * material = (Material *)(rtcGetGeometryUserData(geometry));

		//const Triangle & triangle = surfaces_[ray_hit]
		Vector3 light_pos = Vector3(100, -100, 400);
		Vector3 light_color = Vector3(1,1,1);

		Vector3 p = get_hit_point(my_ray_hit.ray_hit.ray);
		Vector3 l_d = light_pos - p;
		l_d.Normalize();

		Vector3 normal_v = Vector3(normal.x, normal.y, normal.z);
		normal_v.Normalize();

		if (depth <= 0) {
			Color3f bck = background_.GetBackground(my_ray_hit.ray_hit.ray.dir_x, my_ray_hit.ray_hit.ray.dir_y, my_ray_hit.ray_hit.ray.dir_z);
			return Color4f(bck.r, bck.g, bck.b, 1.0f);
			/*Vector3 n = Vector3(normal.x, normal.y, normal.z);
			Vector3 v = Vector3(-my_ray_hit.ray_hit.ray.dir_x, -my_ray_hit.ray_hit.ray.dir_y, -my_ray_hit.ray_hit.ray.dir_y);
			Vector3 l_r = reflect(l_d, n);

			const float enlight = trace_shadow_ray(p, l_d, l_d.L2Norm());

			Vector3 C = material->ambient*light_color +
				enlight * material->diffuse * (max(0.0f, n.DotProduct(l_d))) * light_color +
				enlight * material->specular * powf(max(0.0f, l_r.DotProduct(v)), 1);


			Color4f myColor = Color4f(C.x, C.y, C.z, 1.0f);
			Color4f temp = myColor * material->reflectivity;
			Color4f temp2 = material->specular * temp;
			Color4f temp3 = temp2 + material->ambient;
			return temp3;
			return Color4f(1.0f, 0.0f, 0.0f, 1.0f);*/
		}

		switch (material->shader_)
		{
			case Shader::NORMAL:
			{
				return Color4f(normal.x * 0.5f + 0.5f, normal.y*0.5f + 0.5f, normal.z*0.5f + 0.5f, 1.0f);
				break;
			}
			case Shader::LAMBERT:
			{
				return Color4f(l_d.PosDotProduct(normal_v), l_d.PosDotProduct(normal_v), l_d.PosDotProduct(normal_v), 1.0f);
				break;
			}
			case Shader::PHONG:
			{
				Vector3 n = Vector3(normal.x, normal.y, normal.z);
				Vector3 v = Vector3(-my_ray_hit.ray_hit.ray.dir_x, -my_ray_hit.ray_hit.ray.dir_y, -my_ray_hit.ray_hit.ray.dir_z);
				Vector3 l_r = reflect(l_d, n);

				const float enlight = trace_shadow_ray(p, l_d, l_d.L2Norm());

				Vector3 C =
					0.2f* material->ambient*light_color +
					 enlight * material->diffuse * (max(0.0f, n.DotProduct(l_d))) * light_color +
					 enlight * material->specular * powf(max(0.0f, l_r.DotProduct(v)), 1);
				return Color4f(C.x, C.y, C.z, 1.0f);

				break;
			}
			case Shader::GLASS:
			{
				Color4f refractFragment = Color4f(0.0f, 0.0f, 0.0f, 1.0f);
				float coefReflect = 1.0f;
				normal_v = Vector3(0.0f, 1.0f, 0.0f);
				Vector3 diffuse = material->diffuse;
				Vector3 rd = Vector3(my_ray_hit.ray_hit.ray.dir_x, my_ray_hit.ray_hit.ray.dir_y, my_ray_hit.ray_hit.ray.dir_z);
				//Vector3 rd = Vector3(-0.429f, -0.903f, 0.0f);
				Vector3 rv = -rd;

				float n2 = my_ray_hit.ior;
				float n1 = (n2 == IOR_AIR) ? IOR_GLASS : IOR_AIR;

				float n_divided = n1 / n2;
				float cos_01 = (normal_v.DotProduct(rv));


				Vector3 rr = (2.0f * (normal_v.DotProduct(rv))) * normal_v - rv;

				float tmp = 1.0f - SQR(n_divided) * (1.0f - SQR(cos_01));
				if (tmp > 0) {
						float cos_02 = sqrt(tmp);
						if ( cos_02 < 0) {
							normal_v = -normal_v;
							cos_02 = (-normal_v).DotProduct(rd);
						}
						Vector3 rl = (n_divided * rd) + ((n_divided * cos_01 - cos_02) * normal_v);
						// Fresnel
						float Rs = SQR(((n2 * cos_02) - (n1 * cos_01)) / ((n2 * cos_02) + (n1 * cos_01)));
						float Rp = SQR(((n2 * cos_01) - (n1 * cos_02)) / ((n2 * cos_01) + (n1 * cos_02)));
						float R = 0.5f * (Rs + Rp);

						// Calculate coefficients
						coefReflect = R;
						float coefRefract = 1.0f - coefReflect;

						// Generate refracted ray
						Vector3 vector = get_hit_point(my_ray_hit.ray_hit.ray);
						RTCRay refractedRay = camera_.GenerateRay(0,0);
						refractedRay.dir_x = rl.x;
						refractedRay.dir_y = rl.y;
						refractedRay.dir_z = rl.z;
						refractedRay.tnear = 0.01f;
						refractedRay.org_x = vector.x;
						refractedRay.org_y = vector.y;
						refractedRay.org_z = vector.z;

						MyRTCRayHit myRefractedRTCRayHit = MyRTCRayHit{ refractedRay };

						//myRefractedRTCRayHit.ior = n2;
						Color4f temp2 = trace_ray(myRefractedRTCRayHit, depth - 1, t);
						Color4f temp4 = temp2 * coefRefract;
						//refractFragment = Color4f(1.0f, 0.0f, 0.0f, 1.0f);
						refractFragment =  diffuse * temp4;

				}	


				

			

				// Generate reflected ray
				Vector3 vector = get_hit_point(my_ray_hit.ray_hit.ray);
				RTCRay reflectedRay = camera_.GenerateRay(0, 0);
				reflectedRay.dir_x = rr.x;
				reflectedRay.dir_y = rr.y;
				reflectedRay.dir_z = rr.z;
				reflectedRay.tnear = 0.01f;
				reflectedRay.org_x = vector.x;
				reflectedRay.org_y = vector.y;
				reflectedRay.org_z = vector.z;

				MyRTCRayHit myReflectedRTCRayHit = MyRTCRayHit{ reflectedRay };

			


				//// Set IOR
				//myReflectedRTCRayHit.ior == n1 ? n2 : n1;

				Color4f temp1 = trace_ray(myReflectedRTCRayHit, depth - 1, t);
				Color4f temp3 = temp1 * coefReflect;
				Color4f temp5 = diffuse * temp3;


				return temp5;

				//TODO   fix  black refractFragment  condition if(tmp > 0) never triggered !!!
				//return temp5 +  refractFragment;

				//FOR TEST
				// return refractFragment;
			}
			default:
			{
				return Color4f(l_d.PosDotProduct(normal_v), l_d.PosDotProduct(normal_v), l_d.PosDotProduct(normal_v), 1.0f);
				break;
			}
		}

		/*Color3f bck = background_.GetBackground(my_ray_hit.ray_hit.ray.dir_x, my_ray_hit.ray_hit.ray.dir_y, my_ray_hit.ray_hit.ray.dir_z);
		return Color4f(bck.r, bck.g, bck.b, 1.0f);
		return Color4f(0.0f, 0.0f, 0.0f, 1.0f);*/
	}
	else {
		Color3f bck = background_.GetBackground(my_ray_hit.ray_hit.ray.dir_x, my_ray_hit.ray_hit.ray.dir_y, my_ray_hit.ray_hit.ray.dir_z);
		return Color4f(bck.r, bck.g, bck.b, 1.0f);
	}
}


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