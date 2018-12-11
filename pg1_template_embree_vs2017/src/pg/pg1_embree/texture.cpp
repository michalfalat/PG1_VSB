#include "stdafx.h"
#include "texture.h"
#include "utils.h"

Texture::Texture( const char * file_name )
{
	// image format
	FREE_IMAGE_FORMAT fif = FIF_UNKNOWN;
	// pointer to the image, once loaded
	FIBITMAP * dib =  nullptr;
	// pointer to the image data
	BYTE * bits = nullptr;

	// check the file signature and deduce its format
	fif = FreeImage_GetFileType( file_name, 0 );
	// if still unknown, try to guess the file format from the file extension
	if ( fif == FIF_UNKNOWN )
	{
		fif = FreeImage_GetFIFFromFilename( file_name );
	}
	// if known
	if ( fif != FIF_UNKNOWN )
	{
		// check that the plugin has reading capabilities and load the file
		if ( FreeImage_FIFSupportsReading( fif ) )
		{
			dib = FreeImage_Load( fif, file_name );
		}
		// if the image loaded
		if ( dib )
		{
			// retrieve the image data
			bits = FreeImage_GetBits( dib );
			//FreeImage_ConvertToRawBits()
			// get the image width and height
			width_ = int( FreeImage_GetWidth( dib ) );
			height_ = int( FreeImage_GetHeight( dib ) );

			// if each of these is ok
			if ( ( bits != 0 ) && ( width_ != 0 ) && ( height_ != 0 ) )
			{				
				// texture loaded
				scan_width_ = FreeImage_GetPitch( dib ); // in bytes
				pixel_size_ = FreeImage_GetBPP( dib ) / 8; // in bytes

				data_ = new BYTE[scan_width_ * height_]; // BGR(A) format
				FreeImage_ConvertToRawBits( data_, dib, scan_width_, pixel_size_ * 8, FI_RGBA_RED_MASK, FI_RGBA_GREEN_MASK, FI_RGBA_BLUE_MASK, TRUE );
			}

			FreeImage_Unload( dib );
			bits = nullptr;
		}
	}	
}

Texture::~Texture()
{	
	if ( data_ )
	{
		// free FreeImage's copy of the data
		delete[] data_;
		data_ = nullptr;
		
		width_ = 0;
		height_ = 0;
	}
}

Color4f Texture::get_texel( const float u, const float v ) const
{
	//assert( ( u >= 0.0f && u <= 1.0f ) && ( v >= 0.0f && v <= 1.0f ) );	
	/*
	const int x = max( 0, min( width_ - 1, int( u * width_ ) ) );
	const int y = max( 0, min( height_ - 1, int( v * height_ ) ) );

	const int offset = y * scan_width_ + x * pixel_size_;
	const float b = data_[offset] / 255.0f;
	const float g = data_[offset + 1] / 255.0f;
	const float r = data_[offset + 2] / 255.0f;
	
	return Color4f{ r, g, b, 1 };*/

	/*const int x = max(0, min(width_ - 1, int(u * width_)));
	const int y = max(0, min(height_ - 1, int(v * height_)));
	int offset = y * scan_width_ + x * pixel_size_;
	float b = data_[offset] / 255.0f;
	float g = data_[offset + 1] / 255.0f;
	float r = data_[offset + 2] / 255.0f;

	return Color4f{ getLRGBColorValueForComponent(r), getLRGBColorValueForComponent(g), getLRGBColorValueForComponent(b), 1 };*/

	const float x = max(0, min(width_ - 1, u * width_));
	const float y = max(0, min(height_ - 1, v * height_));

	const int x0 = static_cast<int>(floor(x));
	const int y0 = static_cast<int>(floor(y));

	const int x1 = min(width_ - 1, x0 + 1);
	const int y1 = min(height_ - 1, y0 + 1);

	unsigned char * p1 = &data_[x0 * pixel_size_ + y0 * scan_width_];
	unsigned char * p2 = &data_[x1 * pixel_size_ + y0 * scan_width_];
	unsigned char * p3 = &data_[x1 * pixel_size_ + y1 * scan_width_];
	unsigned char * p4 = &data_[x0 * pixel_size_ + y1 * scan_width_];

	const float kx = x - x0;
	const float ky = y - y0;

	Color4f final_color =  (Color4f(p1[0], p1[1], p1[2], p1[3]) * (1 - kx) * (1 - ky) +
		Color4f(p2[0], p2[1], p2[2], p2[3]) * kx * (1 - ky) +
		Color4f(p3[0], p3[1], p3[2], p3[3]) * kx * ky +
		Color4f(p4[0], p4[1], p4[2], p4[3]) * (1 - kx) * ky) *
		static_cast<float>(1.0 / 255.0);

	return Color4f(getLRGBColorValueForComponent(final_color.b), getLRGBColorValueForComponent(final_color.g), getLRGBColorValueForComponent(final_color.r), 1);
}

int Texture::width() const
{
	return width_;
}

int Texture::height() const
{
	return height_;
}
