#include <iostream>

struct pixel
{
	int R;
	int G;
	int B;
};

class image
{
public:
	image()
	{

	}
	image(int h, int w)
	{
		width = w;
		height = h;
		colormap = new pixel*[height];
		for (int i = 0; i<height; i++)
			colormap[i] = new pixel[width];
	}
	int width;
	int height;
	friend class Just4index;
	class Just4index
	{
		image &img;
		int y;
		Just4index(image &img_input, int y_input):
		img(img_input),y(y_input){}
		friend class image;
	public:
		pixel &operator[](int x)
		{
			return img.colormap[y][x];
		}

	};
	Just4index operator[](int y)
	{
		return Just4index(*this,y);
	}

//private:
	pixel **colormap;	
};
