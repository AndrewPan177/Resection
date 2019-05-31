// Resection.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "math.h"
#include "malloc.h"
#include "iostream"
#include "Eigen\Dense"

using namespace Eigen;
using namespace std;

#define m  50000


int n;			//已知点的对数
int counter = 0;	//迭代次数

double Xs0;			
double Ys0;
double Zs0;
double φ0, ω0, κ0;


//获取旋转矩阵R的值
void GetR(double φ,double ω,double κ,double *R) {
	double a1, a2, a3, b1, b2, b3, c1, c2, c3;

	a1 = cos(φ)*cos(κ) - sin(φ)*sin(ω)*sin(κ);
	a2 = -cos(φ)*sin(κ) - sin(φ)*sin(ω)*cos(κ);
	a3 = -sin(φ)*cos(ω);
	b1 = cos(ω)*sin(κ);
	b2 = cos(ω)*cos(κ);
	b3 = -sin(ω);
	c1 = sin(φ)*cos(κ)+cos(φ)*sin(ω)*sin(κ);
	c2 = -sin(φ)*sin(κ) + cos(φ)*sin(ω)*cos(κ);
	c3 = cos(φ)*cos(ω);

	R[0] = a1;
	R[1] = a2;
	R[2] = a3;
	R[3] = b1;
	R[4] = b2;
	R[5] = b3;
	R[6] = c1;
	R[7] = c2;
	R[8] = c3;

}

//计算lx和ly，以及(x)和(y)
void GetCoordinate(double *X, double *Y, double *Z, double Xs, double Ys, double Zs, 
					double *L, double f, double *R, int k, double *x, double *y) {

	double X_ = -f*(R[0] * (X[k] - Xs) + R[3] * (Y[k] - Ys) + R[6] * (Z[k] - Zs))
		/ (R[2] * (X[k] - Xs) + R[5] * (Y[k] - Ys) + R[8] * (Z[k] - Zs));

	double Y_ = -f*(R[1] * (X[k] - Xs) + R[4] * (Y[k] - Ys) + R[7] * (Z[k] - Zs))
		/ (R[2] * (X[k] - Xs) + R[5] * (Y[k] - Ys) + R[8] * (Z[k] - Zs));

	L[2 * k + 0] = x[k] - X_;
	L[2 * k + 1] = y[k] - Y_;
}

//对每一个已知控制点计算其A
void GetNorEquation(double *X, double *Y, double *Z, double Xs, double Ys, double Zs, double *R, double f,
					double φ, double ω, double κ, int k, double *Atemp, double *x, double*y) {

	double X_bar = 0, Y_bar = 0, Z_bar = 0;
	double a11, a12, a13, a14, a15, a16, a21, a22, a23, a24, a25, a26;
	double t_x = 0, t_y = 0;

	Z_bar = R[2] * (X[k] - Xs) + R[5] * (Y[k] - Ys) + R[8] * (Z[k] - Zs);

	t_x = x[k];
	t_y = y[k];

	a11 = (R[0] * f + R[2] * t_x) / Z_bar;
	a12 = (R[3] * f + R[5] * t_x) / Z_bar;
	a13 = (R[6] * f + R[8] * t_x) / Z_bar;
	a21 = (R[1] * f + R[2] * t_y) / Z_bar;
	a22 = (R[4] * f + R[5] * t_y) / Z_bar;
	a23 = (R[7] * f + R[8] * t_y) / Z_bar;

	a14 = t_y*sin(ω) - ((t_x / f)*(t_x*cos(κ) - t_y*sin(κ)) + f*cos(κ))*cos(ω);
	a15 = -f*sin(κ) - (t_x / f)*(t_x*sin(κ) + t_y*cos(κ));
	a16 = t_y;
	a24 = -t_x*sin(ω) - ((t_y / f)*(t_x*cos(κ) - t_y*sin(κ)) - f*sin(κ))*cos(ω);
	a25 = -f*cos(κ) - (t_y / f)*(t_x*sin(κ) + t_y*cos(κ));
	a26 = -t_x;

	Atemp[0] = a11;	Atemp[1] = a12;	Atemp[2] = a13;	Atemp[3] = a14;	Atemp[4] = a15;		Atemp[5] = a16;
	Atemp[6] = a21;	Atemp[7] = a22;	Atemp[8] = a23;	Atemp[9] = a24;	Atemp[10] = a25;	Atemp[11] = a26;
}


//主要功能
void Function(double *X, double *Y, double *Z, double *x, double *y, double f) {

	counter++;	//迭代次数+1

	//变量声明

	double R[9];
	double *A = (double*)malloc(sizeof(double) * 2 * n * 6);
	double *Atemp = (double*)malloc(sizeof(double) * 2 * 6);
	double *L = (double*)malloc(sizeof(double) * 2 * n);

	MatrixXd MatA(4 * 2, 6);
	MatrixXd MatAtemp;
	MatrixXd MatL(2 * n, 1);
	MatrixXd MatX;
	MatrixXd MatV;

	double mx0 = 0, my0 = 0;

	//组成旋转矩阵

	GetR(φ0, ω0, κ0, R);

	//逐点计算系数矩阵L和A的值
	for (int k = 0; k < n; k++) {

		//计算(x)、(y)和lx、ly
		GetCoordinate(X, Y, Z, Xs0, Ys0, Zs0, L, f, R, k, x, y);

		//逐点组成误差方程式并法化
		GetNorEquation(X, Y, Z, Xs0, Ys0, Zs0, R, f, φ0, ω0, κ0, k, Atemp, x, y);

		//将Atemp的值赋值到A中
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 6; j++) {
				A[(i + 2 * k) * 6 + j] = Atemp[i * 6 + j];
			}
		}
	}

	//解法方程

	//将法方程的各个矩阵用Matrix来表示以便易于相乘
	for (int i = 0; i < 2 * n; i++) 
	{
		for (int j = 0; j < 6; j++){
			MatA(i, j) = A[i * 6 + j];	//矩阵A
		}
		MatL(i, 0) = L[i];				//矩阵L
	}
	
	MatrixXd B = MatA.transpose()*MatA;
	MatX = B.inverse()*(MatA.transpose())*MatL;

	if (MatX(3,0)>0.00001|| MatX(4, 0)>0.00001|| MatX(5, 0)>0000.1) {
		Xs0 += MatX(0, 0);
		Ys0 += MatX(1, 0);
		Zs0 += MatX(2, 0);
		φ0 += MatX(3, 0);
		ω0 += MatX(4, 0);
		κ0 += MatX(5, 0);

		Function(X, Y, Z, x, y, f);
	}
	else {
		printf("迭代了%d次。\n\n", counter);
		printf("\n计算结果为：\n\n%lf\t%lf\t%lf\n\n", Xs0, Ys0, Zs0);
		printf("\n旋转矩阵R为：\n\n");

		//打印R
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				printf(" %.5lf\t", R[i * 3 + j]);
			}
			printf("\n");
		}

		MatV = MatA*MatX - MatL;
		//cout << MatV;

		MatrixXd m0;
		m0 = MatV.transpose()*MatV;

		double m_ = 0; 	
		m_ = sqrt(m0(0, 0) / (2 * n - 6));
		cout <<"\n\n中误差为：\n"<< endl << m_;
	}

}

int main()
{
	n = 4;

	double X[4] = { 36589.41,37631.08,39100.97,40426.54 };
	double Y[4] = { 25273.32,31324.51,24934.98,30319.81 };
	double Z[4] = { 2195.17,728.69,2386.5,757.31 };

	double x[4] = { -86.15,-53.40,-14.78,10.46 };
	double y[4] = { -68.99,82.21,-76.63,64.43 };

	for (int i = 0; i < n; i++) {
		x[i] /= 1000;
		y[i] /= 1000; 
	}
	
	double f = 153.24;
	f /= 1000;
	double x0 = 0, y0 = 0;

	Xs0 = 0;	Ys0 = 0;	Zs0 = 0;
	φ0 = 0;	ω0 = 0;	κ0 = 0;
	
	//计算和确定初值

	for (int i = 0; i < n; i++) {
		Xs0 += X[i];
		Ys0 += Y[i];
	}

	Xs0 /= 4;
	Ys0 /= 4;
	Zs0 = m*f;

	printf("初始值为：\n\n");
	printf("Xs0 = %.2lf\t  Ys0 = %.2lf\t  Zs0 = %.2lf\n", Xs0, Ys0, Zs0);
	printf("φ0 = %.0lf\t          ω0 = %.0lf\t          κ0 = %.0lf\n\n\n", φ0, ω0, κ0);


	Function(X, Y, Z, x, y, f);

	printf("\n\n\n\n");
	return 0;
}

