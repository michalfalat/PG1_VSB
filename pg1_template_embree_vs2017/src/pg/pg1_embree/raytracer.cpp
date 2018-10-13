#include "stdafx.h"
#include "raytracer.h"
#include "objloader.h"
#include "tutorials.h"
#include "material.h"

Raytracer::Raytracer(const int width, const int height,
	const float fov_y, const Vector3 view_from, const Vector3 view_at,
	const char * config) : SimpleGuiDX11(width, height)
{
	InitDeviceAndScene(config);

	camera_ = Camera(width, height, fov_y, view_from, view_at);
}

Vector3 get_hit_point(const RTCRay &ray)
{
	return Vector3(ray.dir_x + ray.tfar * ray.dir_x,
		ray.dir_y + ray.tfar * ray.dir_y,
		ray.dir_z + ray.tfar * ray.dir_z);
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
	return trace_ray(x, y, t);
	// TODO generate primary ray and perform ray cast on the scene
	// setup a hit
	//RTCHit hit;
	//hit.geomID = RTC_INVALID_GEOMETRY_ID;
	//hit.primID = RTC_INVALID_GEOMETRY_ID;
	//hit.Ng_x = 0.0f;
	//hit.Ng_y = 0.0f;
	//hit.Ng_z = 0.0f;

	//// merge ray and hit structures
	//RTCRayHit ray_hit;
	//ray_hit.ray = camera_.GenerateRay(x + 0.5f, y + 0.5f);
	//ray_hit.hit = hit;

	//// intersect ray with the scene
	//RTCIntersectContext context;
	//rtcInitIntersectContext(&context);
	//rtcIntersect1(scene_, &context, &ray_hit);

	//if (ray_hit.hit.geomID != RTC_INVALID_GEOMETRY_ID)
	//{
	//	// we hit something
	//	RTCGeometry geometry = rtcGetGeometry(scene_, ray_hit.hit.geomID);
	//	Normal3f normal;
	//	// get interpolated normal
	//	rtcInterpolate0(geometry, ray_hit.hit.primID, ray_hit.hit.u, ray_hit.hit.v,
	//		RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, &normal.x, 3);

	//	reorient_against(normal, ray_hit.ray.dir_x, ray_hit.ray.dir_y, ray_hit.ray.dir_z);
	//	// and texture coordinates
	//	Coord2f tex_coord;
	//	rtcInterpolate0(geometry, ray_hit.hit.primID, ray_hit.hit.u, ray_hit.hit.v,
	//		RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 1, &tex_coord.u, 2);

	//	/*printf("normal = (%0.3f, %0.3f, %0.3f)\n", normal.x, normal.y, normal.z);
	//	printf("tex_coord = (%0.3f, %0.3f)\n", tex_coord.u, tex_coord.v);*/
	//	Material * material = (Material *)(rtcGetGeometryUserData(geometry));

	//	//const Triangle & triangle = surfaces_[ray_hit]
	//	Vector3 light_pos = Vector3(0,0,300);
	//	Vector3 light_color = Vector3(1,1,1);

	//	Vector3 p = Vector3(ray_hit.ray.dir_x + ray_hit.ray.tfar * ray_hit.ray.dir_x,
	//		ray_hit.ray.dir_y + ray_hit.ray.tfar * ray_hit.ray.dir_y,
	//		ray_hit.ray.dir_z + ray_hit.ray.tfar * ray_hit.ray.dir_z);

	//	Vector3 l_d  = light_pos - p;
	//	l_d.Normalize();

	//	Vector3 n = Vector3(normal.x, normal.y, normal.z);
	//	Vector3 v = Vector3(-ray_hit.ray.dir_x, -ray_hit.ray.dir_y, -ray_hit.ray.dir_y);
	//	Vector3 l_r = reflect(l_d, n);

	//	Vector3 C = material->ambient*light_color +
	//		1 * material->diffuse * (max(0.0f, n.DotProduct(l_d))) * light_color +
	//		1 * material->specular * powf(max(0.0f, l_r.DotProduct(v)), 1);
	//	//return Color4f{ material->diffuse.x, material->diffuse.y, material->diffuse.z, 1.0f };
	//	return Color4f{ C.x, C.y, C.z, 1.0f };
	//}
	//else
	//{
	//	return Color4f{ 0.0f, 0.0f, 0.0f, 1.0f };
	//}
}

Color4f Raytracer::trace_ray(const int x, const int y, const float t) {
	// TODO generate primary ray and perform ray cast on the scene
	// setup a hit
	RTCHit hit;
	hit.geomID = RTC_INVALID_GEOMETRY_ID;
	hit.primID = RTC_INVALID_GEOMETRY_ID;
	hit.Ng_x = 0.0f;
	hit.Ng_y = 0.0f;
	hit.Ng_z = 0.0f;

	// merge ray and hit structures
	MyRTCRayHit my_ray_hit;
	my_ray_hit.ray_hit.ray = camera_.GenerateRay(x + 0.5f, y + 0.5f);
	my_ray_hit.ray_hit.hit = hit;

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
		Vector3 light_pos = Vector3(0, 0, 300);
		Vector3 light_color = Vector3(1, 1, 1);


		switch (material->shader_)
		{
		case Shader::NORMAL:
		{
			return Color4f{ normal.x * 0.5f + 0.5f, normal.y*0.5f + 0.5f, normal.z*0.5f + 0.5f, 1.0f };
			break; 
		}
		case Shader::LAMBERT:
		{
			return Color4f{ 0.0f, 0.0f, 0.0f, 1.0f };
			break;
		}
		case Shader::PHONG:
		{
			Vector3 p = get_hit_point(my_ray_hit.ray_hit.ray);
			Vector3 l_d = light_pos - p;
			l_d.Normalize();
			p.Normalize();

			Vector3 n = Vector3(normal.x, normal.y, normal.z);
			Vector3 v = Vector3(-my_ray_hit.ray_hit.ray.dir_x, -my_ray_hit.ray_hit.ray.dir_y, -my_ray_hit.ray_hit.ray.dir_y);
			Vector3 l_r = reflect(l_d, n);
		
			Vector3 C = material->ambient*light_color +
				1 * material->diffuse * (max(0.0f, n.DotProduct(l_d))) * light_color +
				1 * material->specular * powf(max(0.0f, l_r.DotProduct(v)), 1);
			//return Color4f{ material->diffuse.x, material->diffuse.y, material->diffuse.z, 1.0f };
			return Color4f{ C.x, C.y, C.z, 1.0f };

			break;
		}
		case Shader::GLASS:
		{
			Vector3 rd = Vector3(my_ray_hit.ray_hit.ray.dir_x, my_ray_hit.ray_hit.ray.dir_y, my_ray_hit.ray_hit.ray.dir_z);
			rd.Normalize();

			float n1 = my_ray_hit.ior;
			float n2 = (n1 == IOR_AIR) ? material->ior : n1;

			//Vector3 cos_01 = n * v;

			Vector3 normal_v = Vector3(normal.x, normal.y, normal.z);

			// Calc cos_02
			float cos_02 = (-normal_v).DotProduct(rd);
			if (cos_02 < 0) {
				normal_v = -normal_v;
				cos_02 = (-normal_v).DotProduct(rd);
			}

			// Vector rs
			Vector3 rs = rd - (2 * normal_v.DotProduct(rd)) * normal_v;
			rs.Normalize();

			// Generate reflected ray
			Vector3 vector = get_hit_point(my_ray_hit.ray_hit.ray);
			RTCRay reflectedRay = my_ray_hit.ray_hit.ray;
			reflectedRay.dir_x = rs.x;
			reflectedRay.dir_y = rs.y;
			reflectedRay.dir_z = rs.z;
			reflectedRay.tnear = 0.01f;
			reflectedRay.org_x = vector.x;
			reflectedRay.org_y = vector.y;
			reflectedRay.org_z = vector.z;

			MyRTCRayHit myReflectedRTCRayHit = MyRTCRayHit{ reflectedRay };

			// Calc cos_01
			float n_d = n1 / n2;
			float sqrt_d = (1 - (n_d * n_d) * (1 - (cos_02 * cos_02)));

			// Absolute reflection
			//TODO finish
			/*if (sqrt_d < 0.0f) {
				return trace_ray(x, depth - 1) * 1.0f * diffuse;
			}*/

			float cos_01 = sqrt(sqrt_d);

			// Vector rr
			Vector3 rr = -n_d * rd - (n_d * cos_02 + cos_01) * normal_v;
			rr.Normalize();

			// Vector lr
			Vector3 lr = rr - (2 * normal_v.DotProduct(rr)) * normal_v;
			lr.Normalize();
			lr = -lr; // l => lr

			// Fresnel
			float Rs = SQR((n1 * cos_02 - n2 * cos_01) / (n1 * cos_02 + n2 * cos_01));
			float Rp = SQR((n1 * cos_01 - n2 * cos_02) / (n1 * cos_01 + n2 * cos_02));
			float R = 0.5f * (Rs + Rp);

			// Calculate coefficients
			float coefReflect = R;
			float coefRefract = 1.0f - coefReflect;

			// Generate refracted ray
			Vector3 vector = get_hit_point(my_ray_hit.ray_hit.ray);
			RTCRay refractedRay = my_ray_hit.ray_hit.ray;
			refractedRay.dir_x = lr.x;
			refractedRay.dir_y = lr.y;
			refractedRay.dir_z = lr.z;
			refractedRay.tnear = 0.01f;
			refractedRay.org_x = vector.x;
			refractedRay.org_y = vector.y;
			refractedRay.org_z = vector.z;

			MyRTCRayHit myRefractedRTCRayHit = MyRTCRayHit{ refractedRay };

			// Set IOR
			myRefractedRTCRayHit.ior = n2;
			myReflectedRTCRayHit.ior = n1;

			Vector3 C = material->ambient*light_color +
				1 * material->diffuse * (max(0.0f, n.DotProduct(l_d))) * light_color +
				1 * material->specular * powf(max(0.0f, l_r.DotProduct(v)), 1);
			//return Color4f{ material->diffuse.x, material->diffuse.y, material->diffuse.z, 1.0f };
			return Color4f{ C.x, C.y, C.z, 1.0f };

			return this->light_pos.ambient + calcColor(ray, mat, normal) + trace(reflectedRay, depth - 1) * mat->reflectivity * specular;
			break;
		}
		default:
		{
			return Color4f{ 0.0f, 0.0f, 0.0f, 1.0f };
			break;
		}
		}
		////ray. .....
		//RTCIntersectContext context;
		//rtcInitIntersectContext(context);
		//rtcOf


		return Color4f{ 0.0f, 0.0f, 0.0f, 1.0f };
	}
	else {
		return Color4f{ 0.0f, 0.0f, 0.0f, 1.0f };
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
