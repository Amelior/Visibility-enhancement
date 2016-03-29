#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "image.h"

#pragma warning(disable:4996)

extern "C"
{
#include <jpeglib.h>
#include <jerror.h>
}

using namespace std;

class JPGReader
{
public:
	JPGReader()
	{
		dinfo.err = jpeg_std_error(&jerr);
		cinfo.err = jpeg_std_error(&jerr);
	}
	
	image open(char filename[256])
	{
		FILE *inpfile;
		if ((inpfile = fopen(filename, "rb")) == NULL)
		{
			fprintf(stderr,"can't open %s\n",filename);
			exit(1);
		}

		jpeg_create_decompress(&dinfo);
		jpeg_stdio_src(&dinfo,inpfile);
		jpeg_read_header(&dinfo,true);
		jpeg_start_decompress(&dinfo);

		image I(dinfo.output_height,dinfo.output_width);
		I.width = dinfo.output_width;
		I.height = dinfo.output_height;

		int buffer_height = 1,
			row_stride = I.width*dinfo.output_components;
		//physical row width in output buffer
		JSAMPARRAY buff = (JSAMPARRAY)malloc(sizeof(JSAMPROW));
		buff[0] = (JSAMPROW)malloc(sizeof(JSAMPLE)*row_stride);

		
		unsigned char *image = new unsigned char[I.width*I.height*dinfo.output_components];
		long counter = 0;
		int i = 0;
		while (dinfo.output_scanline < dinfo.output_height)
		{
			i++;
			jpeg_read_scanlines(&dinfo,buff,1);
			memcpy(image+counter,buff[0],row_stride);
			counter+=row_stride;
		}

		for (int i = 0; i<I.height; i++)
			for (int j = 0; j<I.width; j++)
			{
				I[i][j].R = (int)*image++;
				I[i][j].G = (int)*image++;
				I[i][j].B = (int)*image++;
			}

		jpeg_finish_decompress(&dinfo);
		jpeg_destroy_decompress(&dinfo);
		fclose(inpfile);

		return I;
	}

	void save(char filename[256], image O, int quality, bool grayscale)
	{
		FILE *f;

		if ((f = fopen(filename, "wb")) == NULL)
		{
			fprintf(stderr,"can't open %s\n",filename);
			exit(1);
		}

		jpeg_create_compress(&cinfo);
		jpeg_stdio_dest(&cinfo,f);
		cinfo.image_height = O.height;
		cinfo.image_width = O.width;
		cinfo.input_components = 3;
		cinfo.in_color_space = JCS_RGB;
		jpeg_set_defaults(&cinfo);
		cinfo.num_components = 3;
		cinfo.dct_method = JDCT_FLOAT;
		jpeg_set_quality(&cinfo,quality,true);

		jpeg_start_compress(&cinfo,true);

		JSAMPROW row_pointer[1];
		
		unsigned char *img = new unsigned char[O.width*O.height*cinfo.input_components];
		for (int i = 0; i<O.height; i++)
			for (int j = 0; j<O.width; j++)
			{
				if (grayscale)
				{
					img[i*O.width*3+j*3]   = (char)O[i][j].R;
					img[i*O.width*3+j*3 + 1] = (char)O[i][j].R;
					img[i*O.width*3+j*3 + 2] = (char)O[i][j].R;
				}
				else
				{
					img[i*O.width*3+j*3]   = (char)O[i][j].R;
					img[i*O.width*3+j*3 + 1] = (char)O[i][j].G;
					img[i*O.width*3+j*3 + 2] = (char)O[i][j].B;
				}
			}

			int l = 0;
		while (cinfo.next_scanline < cinfo.image_height)
		{
			l++;
			row_pointer[0] = &img[cinfo.next_scanline*O.width*3];
			jpeg_write_scanlines(&cinfo,row_pointer,1);
		}
		jpeg_finish_compress(&cinfo);
		jpeg_destroy_compress(&cinfo);
		fclose(f);
	}
	
private:
	struct jpeg_decompress_struct dinfo;
	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;
};
