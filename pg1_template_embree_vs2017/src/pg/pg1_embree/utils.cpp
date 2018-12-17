#include "stdafx.h"
#include "vector3.h"
#include "structs.h"

using std::mt19937;
using std::uniform_real_distribution;

typedef mt19937                                     Engine;
typedef uniform_real_distribution<float>            Distribution;

auto uniform_generator = std::bind(Distribution(0.0f, 1.0f), Engine(1));


float Random(const float range_min, const float range_max)
{
	float ksi;

	//#pragma omp critical ( random ) 
	{
		//ksi = static_cast<float>( rand() ) / ( RAND_MAX + 1 );		
		ksi = static_cast<float>(uniform_generator());

		/*static float randoms[] = { 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f };
		static int next = 0;
		ksi = randoms[next];
#pragma omp atomic
		next++;
		if ( next > 10 ) next = 0;*/
	}

	return ksi * (range_max - range_min) + range_min;
}

long long GetFileSize64(const char * file_name)
{
	FILE * file = fopen(file_name, "rb");

	if (file != NULL)
	{
		_fseeki64(file, 0, SEEK_END); // pøesun na konec souboru
		long long file_size = _ftelli64(file); // zjištìní aktuální pozice
		_fseeki64(file, 0, SEEK_SET); // pøesun zpìt na zaèátek
		fclose(file);
		file = NULL;

		return file_size;
	}

	return 0;
}

void PrintTime(double t, char * buffer)
{
	// rozklad èasu
	int days = (int)(t / (24.0 * 60.0 * 60.0));
	int hours = (int)((t - days * 24.0 * 60.0 * 60.0) / (60.0 * 60.0));
	int minutes = (int)((t - days * 24.0 * 60.0 * 60.0 - hours * 60.0 * 60.0) / 60.0);
	double seconds = t - days * 24.0 * 60.0 * 60.0 - hours * 60.0 * 60.0 - minutes * 60.0;

	// ošetøení chybných stavù jako napø. 1m60s a pøevedení na korektní zápis 2m00s
	if (seconds >= 59.5)
	{
		seconds = 0.0;
		minutes++;
		if (minutes == 60)
		{
			minutes = 0;
			hours++;
			if (hours == 24)
			{
				hours = 0;
				days++;
			}
		}
	}

	// tisk èasu s odfiltrováním nulových poèáteèních hodnot 0d0h10m14s => 10m14s
	if (days == 0)
	{
		if (hours == 0)
		{
			if (minutes == 0)
			{
				if (seconds < 1.0)
				{
					sprintf(buffer, "%0.0f ms", seconds*1e+3);
				}
				else
				{
					sprintf(buffer, "%0.1f s", seconds);
				}
			}
			else
			{
				sprintf(buffer, "%d m%02.0f s", minutes, seconds);
			}
		}
		else
		{
			sprintf(buffer, "%d h%02d m%02.0f s", hours, minutes, seconds);
		}
	}
	else
	{
		sprintf(buffer, "%d d%02d h%02d m%02.0f s", days, hours, minutes, seconds);
	}
}

std::string TimeToString(const double t)
{
	// rozklad èasu
	int days = static_cast<int>(t / (24.0 * 60.0 * 60.0));
	int hours = static_cast<int>((t - days * 24.0 * 60.0 * 60.0) / (60.0 * 60.0));
	int minutes = static_cast<int>((t - days * 24.0 * 60.0 * 60.0 - hours * 60.0 * 60.0) / 60.0);
	double seconds = t - days * 24.0 * 60.0 * 60.0 - hours * 60.0 * 60.0 - minutes * 60.0;

	// ošetøení chybných stavù jako napø. 1m60s a pøevedení na korektní zápis 2m00s
	if (seconds >= 59.5)
	{
		seconds = 0.0;
		++minutes;

		if (minutes == 60)
		{
			minutes = 0;
			++hours;

			if (hours == 24)
			{
				hours = 0;
				++days;
			}
		}
	}

	char buffer[32] = { 0 };

	// tisk èasu s odfiltrováním nulových poèáteèních hodnot 0d0h10m14s => 10m14s
	if (days == 0)
	{
		if (hours == 0)
		{
			if (minutes == 0)
			{
				if (seconds < 10)
				{
					if (seconds < 1)
					{
						sprintf(buffer, "%0.1fms", seconds * 1e+3);
					}
					else
					{
						sprintf(buffer, "%0.1fs", seconds);
					}
				}
				else
				{
					sprintf(buffer, "%0.0fs", seconds);
				}
			}
			else
			{
				sprintf(buffer, "%dm%02.0fs", minutes, seconds);
			}
		}
		else
		{
			sprintf(buffer, "%dh%02dm%02.0fs", hours, minutes, seconds);
		}
	}
	else
	{
		sprintf(buffer, "%dd%02dh%02dm%02.0fs", days, hours, minutes, seconds);
	}

	return std::string(buffer);
}

char * LTrim(char * s)
{
	while (isspace(*s) || (*s == 0))
	{
		++s;
	};

	return s;
}

char * RTrim(char * s)
{
	char * back = s + strlen(s);

	while (isspace(*--back));

	*(back + 1) = '\0';

	return s;
}

char * Trim(char * s)
{
	return RTrim(LTrim(s));
}

float getLRGBColorValueForComponent(float colorComponent) {
	float colorFactor = 1 / 12.92f;
	return colorComponent <= 0.04045 ? colorComponent * colorFactor : pow((colorComponent + 0.055f) / (1 + 0.055f), 2.4f);
}

float getSRGBColorValueForComponent(float colorComponent) {
	float colorFactor = 1 / 2.4f;
	return colorComponent <= 0.0031308f ? 12.92f * colorComponent : (1 + 0.055f) * pow(colorComponent, colorFactor) - 0.055f;
}

float changeGamma(float value, float coef ) {
	return value * coef;
}

Color4f changeGamma(Color4f value, float coef ) {
	return value * coef;
}

RTCHit createEmptyHit() {
	RTCHit hit;
	hit.geomID = RTC_INVALID_GEOMETRY_ID;
	hit.primID = RTC_INVALID_GEOMETRY_ID;
	hit.Ng_x = 0.0f; // geometry normal
	hit.Ng_y = 0.0f;
	hit.Ng_z = 0.0f;
	return hit;
}

RTCRay createRay(Vector3 origin, Vector3 dir, float tfar, float tnear) {
	RTCRay ray;

	ray.tnear = tnear; // start of ray segment

	ray.org_x = origin.x;
	ray.org_y = origin.y;
	ray.org_z = origin.z;

	ray.dir_x = dir.x; // ray direction
	ray.dir_y = dir.y;
	ray.dir_z = dir.z;
	ray.time = 0.0f; // time of this ray for motion blur
	ray.tfar = tfar; // end of ray segment (set to hit distance)

	ray.mask = 0; // can be used to mask out some geometries for some rays
	ray.id = 0; // identify a ray inside a callback function
	ray.flags = 0; // reserved
	return ray;
}

RTCRayHitWithIor createRayWithEmptyHitAndIor(Vector3 origin, Vector3 dir, float tfar, float tnear, float ior) {
	RTCRayHitWithIor myRtcRay;
	myRtcRay.ray_hit.ray = createRay(origin, dir, tfar, tnear);
	myRtcRay.ray_hit.hit = createEmptyHit();
	myRtcRay.ior = ior;
	return myRtcRay;
}
