#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "image reading.h"
#include <conio.h>

#include <time.h>
#include <limits.h>
#include <math.h>
#include "omp.h"
//#ifdef _DEBUG
//#include <crtdbg.h>
//#define _CRTDBG_MAP_ALLOC
//#endif
//
//#ifdef _CRTDBG_MAP_ALLOC
//  inline void* __cdecl operator new(unsigned int s)
//     { return ::operator new(s, _NORMAL_BLOCK, __FILE__, __LINE__); }
//  #endif /* _CRTDBG_MAP_ALLOC */ 


static const double M_PI = 3.14;

class visibility_enhancer
{
public:
	image I, O, YUV;
	visibility_enhancer(JPGReader R, image I_input)
	{
		I = I_input;
		O = I_input;

		r = R;//temp
	}

	void enhance()
	{
		int it = 0;
		cout<<"Processing image.."<<endl;

		calculate_gystogramm();
		calculate_borders();

		//while (define_fragment(it++))
		//{
		//define_fragment(0);
			cout<<"Area "<<it<<"/49"<<":("<<x0<<","<<y0<<")"<<"->("<<x1<<","<<y1<<")"<<endl;

			calculate_Q();

			for (int i = y0; i<y1; i++)
				for (int j = x0; j<x1; j++)
				{
					//O.colormap[i][j].R = simple_hermite(O.colormap[i][j].R);
					//O.colormap[i][j].G = simple_hermite(O.colormap[i][j].G);
					//O.colormap[i][j].B = simple_hermite(O.colormap[i][j].B);
					O.colormap[i][j].R = hermite_function_via_step(O.colormap[i][j].R);
					O.colormap[i][j].G = hermite_function_via_step(O.colormap[i][j].G);
					O.colormap[i][j].B = hermite_function_via_step(O.colormap[i][j].B);
				}
		//}

	}

private:
	int** Grayscale;
	int *G; //Gystogramm
	int x0, x1, y0, y1; //Fragment coords
	int H0,	H1, C0, C1, Outlining, //Hermit coeffs
		D0, D1;
	double denom, a, b, Q;

	double P0, P1, Q0x, Q0y, Q1x, Q1y;

	void calculate_gystogramm()
	{
		Grayscale = new int*[I.height];
		for (int i = 0; i < I.height; i++)
		{
			Grayscale[i] = new int[I.width];
			for (int j = 0; j < I.width; j++)
			{
				int R = I[i][j].R,
					G = I[i][j].G,
					B = I[i][j].B;
				Grayscale[i][j] = (int)(0.299*R + 0.587*G + 0.114*B);
			}
		}

		G = new int[256];
		for (int i = 0; i<256; i++)
			G[i] = 0;
		for (int i = 0; i<I.height; i++)
			for (int j = 0; j<I.width; j++)
				G[Grayscale[i][j]]++;
	}


	void calculate_borders()
	{
		x0 = 0;
		x1 = I.width;
		y0 = 0;
		y1 = I.height;
		int Max = 0;
		for (int i = 0; i<256; i++)
			Max = max(Max,G[i]);

		int n = 0;
		int *peaks = find_peaks(n, Max);
			
		C0 = 0;
		C1 = 0;
		Outlining = 40;
		D0 = Outlining;
		D1 = 255 - 40;
		P0 = peaks[0]; 
		P1 = 0;

		for (int i = 0; i < 256; i++)
			C0 += G[i];
		C0 /= 100;
		H0 = C0/2;
		C1 = C0;
		H1 = H0;

		for (int i = 0; i<D0; i++)
		{
			if (G[i] >= H0)
			{
				P0 = i;
				break;
			}
			int sum = 0;
			for (int j = 0; j<i; j++)
				sum+=G[j];
			
			if (sum > C0)
			{
				P0 = i;
				break;
			}
		}
		P0 = min((int)P0,D0);

		for (int i = 255; i>= 0; i--)
		{
			if (G[i] >= H1)
			{
				P1 = i;
				break;
			}
			int sum = 0;
			for (int j = 255; j>= D1; j--)
				sum+=G[j];
			
			if (sum > C1)
			{
				P1 = i;
				break;
			}
		}
		P1 = max((int)P1, D1);
		P0/=255;
		P1/=255;
		delete[] peaks;
	}

	bool define_fragment(int it)
	{
		if (it == 0)
		{
			x0 = 0;
			x1 = I.width/4; 
			y0 = 0;
			y1 = I.height/4;
		}
		else
		{
			x0 += I.width / 8;
			x1 += I.width / 8;
			
			if (x1 > I.width)
			{
				x0 = 0;
				x1 = I.width/4;
				y0 += I.height/8;
				y1 += I.height/8;
				if (y1 > I.height)
					return false;				
			}
		}
		//if (it == 0)
		//{
		//	x0 = 0;
		//	x1 = I.width/4 - 1;
		//	y0 = 0;
		//	y1 = I.height/4;
		//}
		//else
		//{
		//	x0 += I.width / 8;
		//	x1 += I.width / 8;
		//	
		//	if (x1 > I.width)
		//	{
		//		x0 = 0;
		//		x1 = I.width/4;
		//		y0 += I.height/8;
		//		y1 += I.height/8;
		//		if (y1 > I.height)
		//			return false;				
		//	}
		//}
		return true;
	}

	void convert_YUV_into_RGB()
	{
		//pixel p = O.colormap[0][0];
		for (int i = y0; i < y1; i++)
		{
			for (int j = x0; j < x1; j++)
			{
				int Y = O[i][j].R,
					U = O[i][j].G,
					V = O[i][j].B;
				O.colormap[i][j].R = (int)(Y + 1.13983*(V - 128));
				O.colormap[i][j].G = (int)(Y - 0.39465*(U - 128) - 0.58060*(V - 128));
				O.colormap[i][j].B = (int)(Y + 2.03211*(U - 128));

			}
		}
	}

	void calculate_Q()
	{
		double k = otsu();

		if (k>0.5)
		{
			Q0x = 1 + 3*(k - 0.5);
			Q0y = 0;
			Q1x = 1;
			Q1y = 0;			
		}
		else
		{
			Q0x = 1;
			Q0y = 0;
			Q1x = 1 + 3*(0.5 - k);
			Q1y = 0;
		}
	}
	
	double otsu()
	{
		double actual_disp = 0,
			max_disp = 0,
			w1 = 0, w2 = 0, m1 = 0, m2 = 0,
			barrier = 0,
			total = O.height*O.width,
			E = 0,
			total_E = 0;

		for (int i = 1; i<256; i++)


		for (int i = 0; i<256; i++)
		{
			w1 += G[i];
			w2 = total - w1;
			E += i*G[i];
			
			m1 = E/w1;
			m2 = (total_E - E)/w2;


			actual_disp = w1*w2*(m1-m2)*(m1-m2);
			if (actual_disp > max_disp)
			{
				max_disp = actual_disp;
				barrier = i;
			}
		}
		
		return ((double)barrier/255);
	}

	double hermite_function_via_step(int X)
	{
		double x = (double)X/255,
			epsilon = 0.005,
			T = -1;
		for (double t = 0; t<=1; t+=0.01)
		{
			//double t = (double)i;
			double value = (2*pow(t,3)-t*t+1)*P0+(3*t*t-2*pow(t,3))*P1+
			(pow(t,3)-t*t+t)*Q0x+(-t*t+pow(t,3))*Q1x;
			if ((value < 0)||(value > 255))
			{
				value = 0;
			}
			if ((value >= x - epsilon)&&(value <= x + epsilon))
			{
				T = t;
				break;
			}			
		}

		if (T == -1.)
		{
			T = 0;
		}
		return(255*(3*T*T-2*pow(T,3)));
	}

	double hermite_function(int X)
	{
		double x = (double)X/255,
			c = (P0 - x) / denom,
			R = (2*pow(a,3) - 9*a*b + 27*c)/54,
			value = 0,
			x1 = -1, x2 = -1, x3 = -1;

		//if ((P0 > x)||(x > P1))
		//	return X;
		if (P0 > x)
			return 0;
		if (x > P1)
			return 255;

		if (R*R < pow(Q,3))
		{
			a/=3;
			double temp = acos(R/sqrt(Q*Q*Q)),
				x1 = -2*sqrt(Q)*cos(temp/3)-a,
				x2 = -2*sqrt(Q)*cos((temp+2*M_PI)/3)-a,
				x3 = -2*sqrt(Q)*cos((temp-2*M_PI)/3)-a;
		}
		else
		{
			if (R*R>=Q*Q*Q)
			{
				double A = -pow(abs(R)+sqrt(R*R-Q*Q*Q),1./3.),
					B;
				if (R < 0)
					A*=-1;
				if (A!=0)
					B = Q/A;
				else
					B = 0;
				x1 = (A+B)-a/3;
			}
		}
		double t;
		if ((x1 >= 0)&&(x1 <= 1))
			t = x1;
		else
		{
			if ((x2 >= 0)&&(x2 <= 1))
				t = x2;
			else
			{
				if ((x3 >= 0)&&(x3 <= 1))
					t = x3;
				else
				{
					t = 0;
				}
			}
		}

		return((3*t*t - 2*pow(t,3)+(t-2*t*t+pow(t,3))*Q0y+(t*t-pow(t,3))*Q1y)*255);
	}

	double simple_hermite(int t)
	{
		double x = (double)t/255;
		return(((2*pow(x,3)-3*x*x+1)*P0+(3*x*x-2*pow(x,3))*P1+
			(pow(x,3)-2*x*x+x)*Q0x+(x*x-pow(x,3))*Q1x)*255);
		//return(((2*pow(x,3)-x*x+1)*P0+(3*x*x-2*pow(x,3))*P1+
		//	(pow(x,3)-x*x+x)*Q0x+(-x*x+pow(x,3))*Q1x)*255);
	}

	int* find_peaks(int &k, int Max)
	{
		int MinDist = 5;
		int *MaybeP = new int[150];
		bool Up = true;

		//0 
		if ((G[0] > G[1])||((G[0] == G[1])&&(check_right(0) == 1)))
				MaybeP[k++] = 0;

		//1 - 254
		for (int i = 1; i<255; i++)
		{
			if ((G[i] >= G[i-1])&&(G[i] >= G[i+1]))
			{
				int left, right;
				if (G[i] == G[i-1])
					left = check_left(i - 1);
				else
					left = 1;

				if (G[i] == G[i+1])
					right = check_right(i + 1);
				else
					right = 1;
				if (left + right > 0) //possible fallout with static gystogramm ( 0 + 0 !> 0)
				{
					MaybeP[k++] = i;
					bool irrelevant = false;
					for (int j = 0; j < MinDist; j++)
					{
						if (i-j >= 0)
							if (G[MaybeP[k-1]] < G[i-j])
							{
								irrelevant = true;
								break;
							}
						if (i+j <= 255)
							if (G[MaybeP[k-1]] < G[i+j])
							{
								irrelevant = true;
								break;
							}
					}

					if (!irrelevant)
						if (k > 1)
							if (MaybeP[k - 1] - MaybeP[k - 2] < MinDist)
								irrelevant = true;

					if (!irrelevant)
						if (G[MaybeP[k-1]] < Max/5)
							irrelevant = true;

					if (irrelevant)
						k--;

				}
			}
		}

		//255
		if ((G[255] > G[254])||((G[255] == G[254])&&(check_left(254) == 1)))
			MaybeP[k++] = 255;

		int *peaks = new int[k];
		for (int i = 0; i<k; i++)
			peaks[i] = MaybeP[i];
		delete[] MaybeP;

		return peaks;
	}

	int check_left(int j)
	{
		bool border = false,
			fits = false;
		while (!border)
		{	
			if (G[j] == G[j-1])
			{
				if (j > 1)
					j--;
				else
					return 0;
			}
			else
			{
				border = true;
				if (G[j] > G[j-1])
					return 1;

				if (G[j] < G[j-1])
					return -1;
			}
		}
		return -1;
	}

	int check_right(int j)
	{
		bool border = false,
			fits = false;
		while (!border)
		{	
			if (G[j] == G[j+1])
			{
				if (j < 254)
					j++;
				else
					return 0;
			}
			else
			{
				border = true;
				if (G[j] > G[j+1])
					return 1;

				if (G[j] < G[j+11])
					return -1;
			}
		}
		return -1;
	}

	//temporary
	JPGReader r;
	void show_gystogramm()
	{
		int Max = 0,
			Height = 180;
		double Scale = 0;
		for (int i = 0; i<256; i++)
		{
			Max = max(Max,G[i]);
		}
		Scale = (double)Height/(double)Max;
		image GImage(Height + 1,256);
		for (int i = Height; i>0; i--)
			for (int j = 0; j<256; j++)
			{
				if (i > Height - Scale*G[j])
					GImage[i][j].R = 0;
				else
					GImage[i][j].R = 255;
			}
		for (int j = 0; j<256;j++)
			GImage[0][j].R = 255;
		r.save("C:\\IP\\Gystogramm.jpg", GImage, 100, true);
	}
};

class omp_visibility_enhancer
{
public:
	image I, O;
	omp_visibility_enhancer(JPGReader R, image I_input)
	{
		I = I_input;
		O = I_input;

		r = R;//temp
	}

	void enhance()
	{
		int it = 0;
		cout<<"Processing image.."<<endl;
		while (define_fragment(it++))
		{
			cout<<"Area "<<it<<"/49"<<":("<<x0<<","<<y0<<")"<<"->("<<x1<<","<<y1<<")"<<endl;
			convert_RGB_into_YUV();
			calculate_gystogramm();
			analyze_gystogramm();
			convert_YUV_into_RGB();
			calculate_function_coeffients();
			for (int i = y0; i<y1; i++)
				for (int j = x0; j<x1; j++)
				{
					O.colormap[i][j].R = simple_hermite(O.colormap[i][j].R);
					O.colormap[i][j].G = simple_hermite(O.colormap[i][j].G);
					O.colormap[i][j].B = simple_hermite(O.colormap[i][j].B);
				}
			//show_gystogramm();
			//delete [] G;
			//for (int i = 0; i<y0 - y1; i++)
			//	delete [] Grayscale[i];
			//delete [] Grayscale;
		}

		//x0 = 0;
		//x1 = I.width;
		//y0 = 0;
		//y1 = I.height;
		//cout<<"Area "<<it<<"/50"<<":("<<x0<<","<<y0<<")"<<"->("<<x1<<","<<y1<<")"<<endl;
		//convert_RGB_into_YUV();
		//calculate_gystogramm();
		//analyze_gystogramm();
		//convert_YUV_into_RGB();
		//for (int i = y0; i<y1; i++)
		//	for (int j = x0; j<x1; j++)
		//	{
		//		O.colormap[i][j].R = hermit_function(O.colormap[i][j].R);
		//		O.colormap[i][j].G = hermit_function(O.colormap[i][j].G);
		//		O.colormap[i][j].B = hermit_function(O.colormap[i][j].B);
		//	}
	}

private:
	int** Grayscale;
	int *G; //Gystogramm
	int x0, x1, y0, y1; //Fragment coords
	int H0,	H1, C0, C1, Outlining, //Hermit coeffs
		D0, D1;
	double denom, a, b, Q;

	double P0, P1, Q0x, Q0y, Q1x, Q1y;

	bool define_fragment(int it)
	{
		if (it == 0)
		{
			x0 = 0;
			x1 = I.width/4; 
			y0 = 0;
			y1 = I.height/4;
		}
		else
		{
			x0 += I.width / 8;
			x1 += I.width / 8;
			
			if (x1 > I.width)
			{
				x0 = 0;
				x1 = I.width/4;
				y0 += I.height/8;
				y1 += I.height/8;
				if (y1 > I.height)
					return false;				
			}
		}
		//if (it == 0)
		//{
		//	x0 = 0;
		//	x1 = I.width/4 - 1;
		//	y0 = 0;
		//	y1 = I.height/4;
		//}
		//else
		//{
		//	x0 += I.width / 8;
		//	x1 += I.width / 8;
		//	
		//	if (x1 > I.width)
		//	{
		//		x0 = 0;
		//		x1 = I.width/4;
		//		y0 += I.height/8;
		//		y1 += I.height/8;
		//		if (y1 > I.height)
		//			return false;				
		//	}
		//}
		return true;
	}

	void convert_RGB_into_YUV()
	{
		Grayscale = new int*[y1 - y0];	
			double start = omp_get_wtime();
			#pragma omp parallel for
			for (int i = y0; i < y1; i++)
			{
				Grayscale[i-y0] = new int[x1 - x0];
				#pragma omp parallel for
				for (int j = x0; j < x1; j++)
				{
					int R = O[i][j].R,
						G = O[i][j].G,
						B = O[i][j].B;
					O[i][j].R = (int)(0.299*R + 0.587*G + 0.114*B);
					O[i][j].G = (int)(-0.14713*R - 0.28886*G + 0.436*B + 128);
					O[i][j].B = (int)(0.615*R - 0.51499*G - 0.10001*B + 128);
					Grayscale[i-y0][j-x0] = O[i][j].R;
				}
			}
			double end = omp_get_wtime(),
				exec_time = end - start;		
	}

	void convert_YUV_into_RGB()
	{
		//pixel p = O.colormap[0][0];
		#pragma omp parallel for
		for (int i = y0; i < y1; i++)
		{
			#pragma omp parallel for
			for (int j = x0; j < x1; j++)
			{
				int Y = O[i][j].R,
					U = O[i][j].G,
					V = O[i][j].B;
				O.colormap[i][j].R = (int)(Y + 1.13983*(V - 128));
				O.colormap[i][j].G = (int)(Y - 0.39465*(U - 128) - 0.58060*(V - 128));
				O.colormap[i][j].B = (int)(Y + 2.03211*(U - 128));

			}
		}
	}

	void calculate_gystogramm()
	{
		G = new int[256];
		#pragma omp parallel for
		for (int i = 0; i<256; i++)
			G[i] = 0;
		#pragma omp parallel for
		for (int i = y0; i<y1; i++)
		{
			#pragma omp parallel for
			for (int j = x0; j<x1; j++)
				G[O.colormap[i][j].R]++;
		}
	}

	void analyze_gystogramm()
	{
		int Max = 0;
		for (int i = 0; i<256; i++)
			Max = max(Max,G[i]);

		int n = 0;
		int *peaks = find_peaks(n, Max);
			
		C0 = 0;
		C1 = 0;
		Outlining = 40;
		D0 = Outlining;
		D1 = 255 - 40;
		P0 = peaks[0]; 
		P1 = 0;

		for (int i = 0; i < 256; i++)
			C0 += G[i];
		C0 /= 100;
		H0 = C0/2;
		C1 = C0;
		H1 = H0;

		for (int i = 0; i<D1; i++)
		{
			if (G[i] >= H0)
			{
				P0 = i;
				break;
			}
			int sum = 0;
			for (int j = 0; j<i; j++)
				sum+=G[j];
			
			if (sum > C0)
			{
				P0 = i;
				break;
			}
		}
		P0 = min((int)P0,D0);

		for (int i = 255; i>= 0; i--)
		{
			if (G[i] >= H1)
			{
				P1 = i;
				break;
			}
			int sum = 0;
			for (int j = 255; j>= D1; j--)
				sum+=G[j];
			
			if (sum > C1)
			{
				P1 = i;
				break;
			}
		}
		P1 = max((int)P1, D1);
		
		double k = otsu();

		if (k>0.5)
		{
			Q0x = 1 + 3*(k - 0.5);
			Q0y = 0;
			Q1x = 1;
			Q1y = 0;			
		}
		else
		{
			Q0x = 1;
			Q0y = 0;
			Q1x = 1 + 3*(0.5 - k);
			Q1y = 0;
		}
		delete[] peaks;
		P0/=255;
		P1/=255;
		Q0x/=255;
		Q1x/=255;
		Q0y/=255;
		Q1y/=255;
	}
	
	double otsu()
	{
		double actual_disp = 0,
			max_disp = 0,
			w1 = 0, w2 = 0, m1 = 0, m2 = 0,
			barrier = 0,
			total = O.height*O.width,
			E = 0,
			total_E = 0;

		for (int i = 1; i<256; i++)


		for (int i = 0; i<256; i++)
		{
			w1 += G[i];
			w2 = total - w1;
			E += i*G[i];
			
			m1 = E/w1;
			m2 = (total_E - E)/w2;


			actual_disp = w1*w2*(m1-m2)*(m1-m2);
			if (actual_disp > max_disp)
			{
				max_disp = actual_disp;
				barrier = i;
			}
		}
		
		return ((double)barrier/255);
	}

	void calculate_function_coeffients()
	{
		denom = 2*P0 - 2*P1 + Q0x - Q1x;
		a = (3*P1 - 3*P0 - Q0x + Q1x) / denom;
		b = Q0x / denom;
		Q = (a*a - 3*b)/9;
	}


	double hermite_function(int X)
	{
		double x = (double)X/255,
			c = (P0 - x) / denom,
			R = (2*pow(a,3) - 9*a*b + 27*c)/54,
			value = 0,
			x1 = -1, x2 = -1, x3 = -1;

		//if ((P0 > x)||(x > P1))
		//	return X;
		if (P0 > x)
			return 0;
		if (x > P1)
			return 255;

		if (R*R < pow(Q,3))
		{
			a/=3;
			double temp = acos(R/sqrt(Q*Q*Q)),
				x1 = -2*sqrt(Q)*cos(temp/3)-a,
				x2 = -2*sqrt(Q)*cos((temp+2*M_PI)/3)-a,
				x3 = -2*sqrt(Q)*cos((temp-2*M_PI)/3)-a;
		}
		else
		{
			if (R*R>=Q*Q*Q)
			{
				double A = -pow(abs(R)+sqrt(R*R-Q*Q*Q),1./3.),
					B;
				if (R < 0)
					A*=-1;
				if (A!=0)
					B = Q/A;
				else
					B = 0;
				x1 = (A+B)-a/3;
			}
		}
		double t;
		if ((x1 >= 0)&&(x1 <= 1))
			t = x1;
		else
		{
			if ((x2 >= 0)&&(x2 <= 1))
				t = x2;
			else
			{
				if ((x3 >= 0)&&(x3 <= 1))
					t = x3;
				else
				{
					t = 0;
				}
			}
		}

		return((3*t*t - 2*pow(t,3)+(t-2*t*t+pow(t,3))*Q0y+(t*t-pow(t,3))*Q1y)*255);
	}

	double simple_hermite(int t)
	{
		double x = (double)t/255;
		return(((2*pow(x,3)-3*x*x+1)*P0+(3*x*x-2*pow(x,3))*P1+
			(pow(x,3)-2*x*x+x)*Q0x+(x*x-pow(x,3))*Q1x)*255);
		//return(((2*pow(x,3)-x*x+1)*P0+(3*x*x-2*pow(x,3))*P1+
		//	(pow(x,3)-x*x+x)*Q0x+(-x*x+pow(x,3))*Q1x)*255);
	}

	int* find_peaks(int &k, int Max)
	{
		int MinDist = 5;
		int *MaybeP = new int[150];
		bool Up = true;

		//0 
		if ((G[0] > G[1])||((G[0] == G[1])&&(check_right(0) == 1)))
				MaybeP[k++] = 0;

		//1 - 254
		for (int i = 1; i<255; i++)
		{
			if ((G[i] >= G[i-1])&&(G[i] >= G[i+1]))
			{
				int left, right;
				if (G[i] == G[i-1])
					left = check_left(i - 1);
				else
					left = 1;

				if (G[i] == G[i+1])
					right = check_right(i + 1);
				else
					right = 1;
				if (left + right > 0) //possible fallout with static gystogramm ( 0 + 0 !> 0)
				{
					MaybeP[k++] = i;
					bool irrelevant = false;
					for (int j = 0; j < MinDist; j++)
					{
						if (i-j >= 0)
							if (G[MaybeP[k-1]] < G[i-j])
							{
								irrelevant = true;
								break;
							}
						if (i+j <= 255)
							if (G[MaybeP[k-1]] < G[i+j])
							{
								irrelevant = true;
								break;
							}
					}

					if (!irrelevant)
						if (k > 1)
							if (MaybeP[k - 1] - MaybeP[k - 2] < MinDist)
								irrelevant = true;

					if (!irrelevant)
						if (G[MaybeP[k-1]] < Max/5)
							irrelevant = true;

					if (irrelevant)
						k--;

				}
			}
		}

		//255
		if ((G[255] > G[254])||((G[255] == G[254])&&(check_left(254) == 1)))
			MaybeP[k++] = 255;

		int *peaks = new int[k];
		for (int i = 0; i<k; i++)
			peaks[i] = MaybeP[i];
		delete[] MaybeP;

		return peaks;
	}

	int check_left(int j)
	{
		bool border = false,
			fits = false;
		while (!border)
		{	
			if (G[j] == G[j-1])
			{
				if (j > 1)
					j--;
				else
					return 0;
			}
			else
			{
				border = true;
				if (G[j] > G[j-1])
					return 1;

				if (G[j] < G[j-1])
					return -1;
			}
		}
		return -1;
	}

int check_right(int j)
	{
		bool border = false,
			fits = false;
		while (!border)
		{	
			if (G[j] == G[j+1])
			{
				if (j < 254)
					j++;
				else
					return 0;
			}
			else
			{
				border = true;
				if (G[j] > G[j+1])
					return 1;

				if (G[j] < G[j+11])
					return -1;
			}
		}
		return -1;
	}

	//temporary
	JPGReader r;
	void show_gystogramm()
	{
		int Max = 0,
			Height = 180;
		double Scale = 0;
		for (int i = 0; i<256; i++)
		{
			Max = max(Max,G[i]);
		}
		Scale = (double)Height/(double)Max;
		image GImage(Height + 1,256);
		for (int i = Height; i>0; i--)
			for (int j = 0; j<256; j++)
			{
				if (i > Height - Scale*G[j])
					GImage[i][j].R = 0;
				else
					GImage[i][j].R = 255;
			}
		for (int j = 0; j<256;j++)
			GImage[0][j].R = 255;
		r.save("C:\\IP\\Gystogramm.jpg", GImage, 100, true);
	}
};

int main()
{	
	char inputfile[256] = "C:\\IP\\squeezed.jpg",
		outputfile[256] = "C:\\IP\\parts.jpg";

	JPGReader reader;
	visibility_enhancer E(reader, reader.open(inputfile));
	E.enhance();
	//quality - [0,..,100]
	reader.save(outputfile,E.O,100, false);
	
	return 0;
}