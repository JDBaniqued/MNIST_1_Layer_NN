#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrixlib.h"

matrix wupdate(matrix x, matrix y, matrix z, matrix w)
{
	int i, n, j;
	double error, alpha;
	char junk[20];
	matrix a = makemat(z.row, z.col);
	FILE* train;
	train = fopen("Train.txt", "r");

	for (n = 0; n < 60000; n++)
	{
		//input training data
		fscanf(train, "%s", &junk);
		for (i = 0; i < y.row; i++)
		{
			fscanf(train, "%lf", &y.mat[i][0]);
		}
		//printf("y\n");
		//printmat(y);
		x.mat[0][0] = 1;
		fscanf(train, "%s", &junk);
		for (i = 1; i < x.row; i++)
		{
			fscanf(train, "%lf", &x.mat[i][0]);
			x.mat[i][0] = x.mat[i][0] / 254;
		}
		//printf("x\n");
		//printmat(x);

		//compute predicted output
		z = multmat(w, x);
		//printf("z\n");
		//printmat(z);
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
		}
		//printmat(a);
		//printf("\n");

		//compute error
		error = 0;
		for (i = 0; i < a.row; i++)
		{
			error += pow(a.mat[i][0] - y.mat[i][0], 2);
		}
		//printf("%d error = \t%lf\n", n, error);

		//adjust w
		alpha = 0.3;
		for (i = 0; i < w.row; i++)
		{
			for (j = 0; j < w.col; j++)
			{
				w.mat[i][j] -= alpha * 2 * (a.mat[i][0] - y.mat[i][0]) * pow(1 + exp(-1 * z.mat[i][0]), -2) * exp(-1 * z.mat[i][0]) * x.mat[j][0];
			}
		}
		//printf("w\n");
		//printmat(w);
	}
	fclose(train);

	return w;
}
