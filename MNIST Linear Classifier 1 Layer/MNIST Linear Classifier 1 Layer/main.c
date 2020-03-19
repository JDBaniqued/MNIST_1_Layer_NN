#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrixlib.h"

int main()
{
	int i, j,n;
	FILE* outw = fopen("outw.csv", "w");
	matrix x, y, z, w;
	FILE* test = fopen("Test.txt", "r");
	matrix a = makemat(10, 1);
	int real, prediction, errornum;
	errornum = 0;
	//printf("julius");
	double max;
	x = makemat(785, 1);
	y = makemat(10, 1);
	z = makemat(10, 1);
	//initialize w
	w = makemat(10, 785);
	w = wupdate(x, y, z, w);
	
	for (n = 0; n < 10000; n++)
	{
		fscanf(test, "%d", &real);
		x.mat[0][0] = 1;
		for (i = 1; i < 785; i++)
		{
			fscanf(test, "%lf", &x.mat[i][0]);
			x.mat[i][0] = x.mat[i][0] / 254;
		}
		z = multmat(w, x);
		max = 0;
		for (i = 0; i < a.row; i++)
		{
			if (z.mat[i][0] < -13)
			{
				a.mat[i][0] = 0;
			}
			else if (z.mat[i][0] > 13)
			{
				a.mat[i][0] = 1;
			}
			else
			{
				a.mat[i][0] = pow(1 + exp(-1 * z.mat[i][0]), -1);
			}
			if (a.mat[i][0] > max)
			{
				max = a.mat[i][0];
				prediction = i;
			}
		}
		
		if (prediction != real)
		{
			//printf("oops\n");
			errornum += 1;
		}

		/*if (n % 100 == 0)
		{
			printf("%d\n", n);
		}*/
	}
	printf("%d", errornum);
}