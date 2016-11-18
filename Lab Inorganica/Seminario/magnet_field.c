#include <stdio.h>
#include <stdlib.h>

int L = 10, m, N, A0 = 2;
double l = 1.5, d = 1.5, h = 0.05;

int pos(int i, int j);

double *init(int x0, int x1, int y0, int y1, double *array);
double *init2(int x00, int x01, int x10, int x11, int y0, int y1, double *array);


int main(int argc, char **argv)
{
	m = L/h; // discretization
	int up, down, left, right, i = 0, j = 0, n = 0;
	double average, N = 2*m*m;
	
	double* space = malloc(m*m*sizeof(double)); // spatial array
	
	int min = (L - l)/(2*h), max = (L + l)/(2*h); // singe magnet coordinates
	
	init(min, max, min, max, space);
	
    /*
        starts iteration for the first magnet
    */
	while (n < N)
	{
		for(i = 1; i < m-1; i++)
		{
			for(j = 1; j < m-1; j++)
			{
				up = pos(i-1, j);
				down = pos(i+1, j);
				left = pos(i, j-1);
				right = pos(i, j+1);
				
				if (!(j >= min && j <= max && i >= min && i <= max))
				{
					average = (space[up] + space[down] + space[left] + space[right])/4;
					space[pos(i,j)] = average;
				}
			}
		}
		n += 1;
	}
	for(i = 0; i < m*m; i++)
	{
		printf("%f ", space[i]); // prints results with spaces inbetween
	}
	
    /*
        second magnet variables relate
    */
	n = 0;
	
	int x00, x01, x10, x11, y0, y1;
	
	x00 = (L/2.0 - l)/(2*h);
	x01 = (L/2.0 + l)/(2*h);
	
	x10 = (3*L/2.0 - l)/(2*h);
	x11 = (3*L/2.0 + l)/(2*h);
	
	y0 = (L - l)/(2*h);
	y1 = (L + l)/(2*h);
	
	init2(x00, x01, x10, x11, y0, y1, space);
	
	while (n < N)
	{
		for(i = 1; i < m-1; i++)
		{
			for(j = 1; j < m-1; j++)
			{
				up = pos(i-1, j);
				down = pos(i+1, j);
				left = pos(i, j-1);
				right = pos(i, j+1);
				
				if (!(j >= x00 && j <= x01 && i >= y0 && i <= y1) && !(j >= x10 && j <= x11 && i >= y0 && i <= y1))
				{
					average = (space[up] + space[down] + space[left] + space[right])/4;
					space[pos(i,j)] = average;
				}
			}
		}
		n += 1;
	}
	
	for(i = 0; i < m*m; i++)
	{
		printf("%f ", space[i]);
	}
	
	return 0;
}

int pos(int i, int j)
{
	return i*m + j;
}

double *init(int x0, int x1, int y0, int y1, double *array)
{	
	int a, b;
	for(a = x0; a <= x1; a++)
	{
		for(b = y0; b <= y0+(y1-y0)/2; b++)
		{
			array[pos(b, a)] = A0/2;
		}
		for(b = y0+(y1-y0)/2; b <= y1; b++)
		{
			array[pos(b, a)] = -A0/2;
		}
	}
	return array;
}

double *init2(int x00, int x01, int x10, int x11, int y0, int y1, double *array)
{	
	int a, b;
	
	for(b = y0; b <= y1; b++)
	{
		for(a = x00; a <= x00 + (x01-x00)/2; a++)
		{
			array[pos(b, a)] = -A0/2;
		}
		for(a = x00 + (x01-x00)/2; a <= x01; a++)
		{
			array[pos(b, a)] = A0/2;
		}
		
		for(a = x10; a <= x10 + (x11-x10)/2; a++)
		{
			array[pos(b, a)] = -A0/2;
		}
		for(a = x10 + (x11-x10)/2; a <= x11; a++)
		{
			array[pos(b, a)] = A0/2;
		}
	}
	return array;
}
