#include "stdafx.h"
#include "material.h"
#include "utils.h"

const char Material::kDiffuseMapSlot = 0;
const char Material::kSpecularMapSlot = 1;
const char Material::kNormalMapSlot = 2;
const char Material::kOpacityMapSlot = 3;

Material::Material()
{
	// defaultní materiál
	ambient = Vector3( 0.1f, 0.1f, 0.1f );
	diffuse = Vector3( 0.4f, 0.4f, 0.4f );
	specular = Vector3( 0.8f, 0.8f, 0.8f );	

	emission = Vector3( 0.0f, 0.0f, 0.0f );	

	reflectivity = static_cast<float>( 0.99 );
	shininess = 1;

	ior = -1;

	memset( textures_, 0, sizeof( *textures_ ) * NO_TEXTURES );

	name_ = "default";
}

Material::Material( std::string & name, const Vector3 & ambient, const Vector3 & diffuse,
	const Vector3 & specular, const Vector3 & emission, const float reflectivity, 
	const float shininess, const float ior, Texture ** textures, const int no_textures )
{
	name_ = name;

	this->ambient = ambient;
	this->diffuse = diffuse;
	this->specular = specular;

	this->emission = emission;

	this->reflectivity = reflectivity;
	this->shininess = shininess;	

	this->ior = ior;

	if ( textures )
	{
		memcpy( textures_, textures, sizeof( textures ) * no_textures );
	}
}

Material::~Material()
{
	for ( int i = 0; i < NO_TEXTURES; ++i )
	{
		if ( textures_[i] )
		{
			delete[] textures_[i];
			textures_[i] = nullptr;
		};
	}
}

void Material::set_name( const char * name )
{
	name_ = std::string( name );
}

std::string Material::get_name() const
{
	return name_;
}

void Material::set_texture( const int slot, Texture * texture )
{
	textures_[slot] = texture;
}

void Material::set_shader(Shader shader)
{
	this->shader_ = shader;
}

Texture * Material::get_texture( const int slot ) const
{
	return textures_[slot];
}

Vector3 Material::doDiffuse(const Coord2f * tex_coord) const
{
	if (tex_coord)
	{
		Texture * texture = textures_[kDiffuseMapSlot];

		if (texture) {
			Color4f texel = texture->get_texel(tex_coord->u, tex_coord->v);
			return Vector3{ texel.r, texel.g, texel.b };
			//return Vector3(getSRGBColorValueForComponent(texel.r), getSRGBColorValueForComponent(texel.g), getSRGBColorValueForComponent(texel.b));
		}

		return this->diffuse;
	}
}

//Vector3 Material::doDiffuse(Coord2f tex_coord) {
//
//	if (textures_[kDiffuseMapSlot]) {
//		Color4f color = textures_[kDiffuseMapSlot]->get_texel(tex_coord.u, 1 - tex_coord.v);
//		return Vector3{ color.r, color.g, color.b };
//	}
//
//	return this->diffuse;
//}

//Vector3 Material::getSpecular(Coord2f tex_coord) {
//
//	if (textures_[kSpecularMapSlot]) {
//		Color4f color = textures_[kDiffuseMapSlot]->get_texel(tex_coord.u, 1 - tex_coord.v);
//		return Vector3{ color.r, color.g, color.b };
//	}
//
//	return this->specular;
//}
