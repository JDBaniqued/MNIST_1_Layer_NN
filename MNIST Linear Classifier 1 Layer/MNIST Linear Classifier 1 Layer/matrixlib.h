#include <stdio.h>
#include <stdlib.h>
#include <math.h>


struct Matrix
{
    int row;
    int col;
    double** mat;
};
typedef struct Matrix matrix;

matrix makemat(int row, int col)
{
    matrix A;
    int i, j;
    A.row = row;
    A.col = col;

    A.mat = (double**)malloc(sizeof(double*) * A.row);
    for (i = 0;i < A.row;i++)
    {
        A.mat[i] = (double*)malloc(sizeof(double) * A.col);
    }

    for (i = 0;i < A.row;i++)
        for (j = 0;j < A.col;j++)
        {
            A.mat[i][j] = 0;
        }
    return A;
}

void printmat(matrix A)
{
    int i, j;
    for (i = 0; i < A.row; i++)
    {
        for (j = 0; j < A.col; j++)
        {
            printf("%lf\t", A.mat[i][j]);
        }
        printf("\n");
    }
}

matrix augmentmat(matrix A, matrix B)
{
    if (A.row != B.row)
        printf("Cannot augment matrices. A.row != B.row");

    int i, j;
    matrix C;
    C = makemat(A.row, A.col + B.col);

    for (i = 0; i < C.row; i++)
    {
        for (j = 0; j < C.col; j++)
        {
            if (j < A.col)
                C.mat[i][j] = A.mat[i][j];
            else
                C.mat[i][j] = B.mat[i][j - A.col];
        }
    }

    return C;
}

//rows i and j will interchange
matrix rowpivot(matrix A, int i, int j)
{
    int p;
    double temp;
    for (p = 0; p < A.col; p++)
    {
        temp = A.mat[i][p];
        A.mat[i][p] = A.mat[j][p];
        A.mat[j][p] = temp;
    }

    return A;
}

matrix gaussjordan(matrix A, matrix B)
{
    int i = 0, j, k, max;
    matrix C, D;

    //Gauss Elimination
    C = augmentmat(A, B);
    for (i = 0; i < C.row - 1; i++)
    {
        max = i;
        for (j = i + 1; j < C.row; j++)
        {
            if (C.mat[max][i] < C.mat[j][i])
                max = j;
        }
        C = rowpivot(C, i, max);
        for (j = i + 1; j < C.col; j++)
        {
            C.mat[i][j] = C.mat[i][j] / C.mat[i][i];
        }
        C.mat[i][i] = 1;
        for (j = i + 1; j < C.row; j++)
        {
            for (k = i + 1; k < C.col; k++)
            {
                C.mat[j][k] = C.mat[j][k] - C.mat[i][k] * C.mat[j][i];
            }
            C.mat[j][i] = 0;
        }
    }
    //backsubstitution
    C.mat[i][i + 1] = C.mat[i][i + 1] / C.mat[i][i];
    C.mat[i][i] = 1;
    for (j = C.col - 2; j > 0; j--)
    {
        for (i = 0; i < j; i++)
        {
            C.mat[i][C.col - 1] = C.mat[i][C.col - 1] - C.mat[i][j] * C.mat[j][C.col - 1];
            C.mat[i][j] = 0;
        }
    }

    D = makemat(C.row, 1);
    for (i = 0; i < C.row; i++)
    {
        D.mat[i][0] = C.mat[i][C.col - 1];
    }
    return D;

}

matrix addmat(matrix A, matrix B)
{
    int i, j;
    matrix C;

    C = makemat(A.row, A.col);

    for (i = 0; i < A.row; i++)
    {
        for (j = 0; j < A.col; j++)
        {
            C.mat[i][j] = A.mat[i][j] + B.mat[i][j];
        }
    }

    return C;
}

matrix submat(matrix A, matrix B)
{
    int i, j;
    matrix C;

    C = makemat(A.row, A.col);

    for (i = 0; i < A.row; i++)
    {
        for (j = 0; j < A.col; j++)
        {
            C.mat[i][j] = A.mat[i][j] - B.mat[i][j];
        }
    }

    return C;
}

matrix multmat(matrix A, matrix B)
{
    matrix C;
    int i, j, k;

    C = makemat(A.row, B.col);

    for (i = 0; i < C.row; i++)
    {
        for (j = 0; j < C.col; j++)
        {
            for (k = 0; k < A.col; k++)
            {
                C.mat[i][j] += A.mat[i][k] * B.mat[k][j];
            }
        }
    }

    return C;
}

matrix transmat(matrix A)
{
    matrix B;
    int i, j;

    B = makemat(A.row, A.col);

    for (i = 0; i < B.row; i++)
    {
        for (j = 0; j < B.col; j++)
        {
            B.mat[i][j] = A.mat[j][i];
        }
    }

    return B;
}

matrix scamult(matrix A, double a)
{
    int i, j;

    for (i = 0; i < A.row; i++)
    {
        for (j = 0; j < A.col; j++)
        {
            A.mat[i][j] = A.mat[i][j] * a;
        }
    }

    return A;
}

matrix CG(matrix A, matrix b)
{
    int k;
    double p;
    matrix x, r, w, pf, pw;

    x = makemat(b.row, 1);
    r = makemat(b.row, 1);
    w = makemat(b.row, 1);
    pf = makemat(1, 1);
    pw = makemat(1, 1);

    for (k = 0; k < x.row; k++)
    {
        x.mat[k][0] = 0;
    }

    for (k = 0; k < 1000; k++)
    {
        r = submat(b, multmat(A, x));
        if (k == 0)
        {
            w = r;
            pf = multmat(transmat(w), w);
            pf.mat[0][0] = pf.mat[0][0] * pf.mat[0][0];
            pw = multmat(multmat(transmat(w), A), w);
            x = addmat(x, scamult(w, pf.mat[0][0] * pow(pw.mat[0][0], -1)));
        }
        else
        {
            p = pf.mat[0][0];
            pf = multmat(transmat(r), r);
            pf.mat[0][0] = pf.mat[0][0] * pf.mat[0][0];
            w = addmat(r, scamult(w, pow(p, -1) * pf.mat[0][0]));
            pw = multmat(multmat(transmat(w), A), w);
            x = addmat(x, scamult(w, pf.mat[0][0] * pow(pw.mat[0][0], -1)));
        }

        if (sqrt(p) < 0.00001)
            break;
    }

    return x;
}

void printanumber(int n)
{
    printf("%d", n);
}

matrix gaussseidel(matrix A, matrix b)
{
    int i, j, k;
    matrix x, s;
    double e = 0;

    x = makemat(b.row, b.col);
    s = makemat(x.row, x.col);
    for (i = 0; i < x.row; i++)
    {
        x.mat[i][0] = 0;
    }

    for (i = 0; i < 1000; i++)
    {
        e = 0;
        for (j = 0; j < x.row; j++)
        {
            s.mat[j][0] = x.mat[j][0];
        }
        for (j = 0; j < x.row; j++)
        {
            x.mat[j][0] = 0;
            for (k = 0; k < A.col; k++)
            {
                if (j != k)
                    x.mat[j][0] += -1 * A.mat[j][k] * x.mat[k][0] * pow(A.mat[j][j], -1);
            }
            x.mat[j][0] += b.mat[j][0] * pow(A.mat[j][j], -1);
        }
        for (j = 0; j < x.row; j++)
        {
            if (fabs(s.mat[j][0] - x.mat[j][0]) > e)
            {
                e = x.mat[j][0];
            }
        }
        if (e < 0.00001)
            break;
    }

    return x;
}

matrix thomas(matrix A, matrix b)
{
    int k;
    double m;
    matrix x;

    x = makemat(b.row, b.col);

    for (k = 1; k < A.row; k++)
    {
        m = A.mat[k][k - 1] * pow(A.mat[k - 1][k - 1], -1);
        A.mat[k][k] = A.mat[k][k] - m * A.mat[k - 1][k];
        b.mat[k][0] = b.mat[k][0] - m * b.mat[k - 1][0];
    }

    x.mat[x.row - 1][0] = b.mat[b.row - 1][0] * pow(A.mat[A.row - 1][A.col - 1], -1);

    for (k = x.row - 2; k >= 0; k--)
    {
        x.mat[k][0] = (b.mat[k][0] - A.mat[k][k + 1] * x.mat[k + 1][0]) * pow(A.mat[k][k], -1);
    }

    return x;
}

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
        //if (n % 100 == 0)
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
