#include "../inc/gauss.h"

void Gauss(int N, double* A, double* B)
{
	//分配内存空间
	double** a;
	int i, j, k, n;
	a = new double* [N];
	for(i = 0; i != N; i++)
	{
		a[i] = new double [N];
	}
	double* b;
	b = new double [N];
	double* x;
	x = new double [N];

	//赋值
	for(i = 0; i != N; i++)
	{
		for(j = 0; j != N; j++)
		{
			a[i][j] = A[i*N + j];
		}
	}
	
	for(j = 0; j != N; j++)
	{
		b[j] = B[j];
	}


	//选列主元Gauss消去法求解线性方程组
	int num;
	double major,temp,sum;
	double* m;
	m = new double [N];

	for(k = 0; k != N-1; ++k)
	{
		//step-1:选主元
		num = k;
		major = a[k][k];
		for(i = k+1; i != N; i++)
		{
			if( fabs(a[i][k]) > fabs(major) )
			{
				num = i;
				major = a[i][k];
			}
		}

		//step-2:行交换
		for(j = k; j != N; j++)
		{
			temp = a[k][j];
			a[k][j] = a[num][j];
			a[num][j] = temp;
		}
		temp = b[k];
		b[k] = b[num];
		b[num] = temp;

		//step-3:Gauss消元法
		for(i = k+1; i != N; i++)
		{
			m[i] = a[i][k]/a[k][k];

			for(j = k+1; j != N; j++)
			{
				a[i][j] = a[i][j] - m[i] * a[k][j];
			}

			b[i] = b[i] - m[i] * b[k];
		}

	}
	
	//step-4: 回代过程
	x[N-1] = b[N-1]/a[N-1][N-1];
	for(n = N-2; n != -1; n--)
	{
		sum = 0.0;
		for(j = n+1; j != N; j++)
		{
			sum = sum + a[n][j] * x[j];
		}

		x[n] = (b[n] - sum)/a[n][n];
	}

	//返回解向量
	for(i = 0; i != N; i++)
	{
		B[i] = x[i];
	}

	delete[] m;
	delete[] b;
	delete[] x;
	for(i = 0; i != N; i++)
	{
		delete[] a[i];
	}
	delete[] a;
	
}