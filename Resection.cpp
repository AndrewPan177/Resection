// Resection.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"
#include "math.h"
#include "malloc.h"
#include "iostream"
#include "Eigen\Dense"

using namespace Eigen;
using namespace std;

#define m  50000


int n;			//��֪��Ķ���
int counter = 0;	//��������

double Xs0;			
double Ys0;
double Zs0;
double ��0, ��0, ��0;


//��ȡ��ת����R��ֵ
void GetR(double ��,double ��,double ��,double *R) {
	double a1, a2, a3, b1, b2, b3, c1, c2, c3;

	a1 = cos(��)*cos(��) - sin(��)*sin(��)*sin(��);
	a2 = -cos(��)*sin(��) - sin(��)*sin(��)*cos(��);
	a3 = -sin(��)*cos(��);
	b1 = cos(��)*sin(��);
	b2 = cos(��)*cos(��);
	b3 = -sin(��);
	c1 = sin(��)*cos(��)+cos(��)*sin(��)*sin(��);
	c2 = -sin(��)*sin(��) + cos(��)*sin(��)*cos(��);
	c3 = cos(��)*cos(��);

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

//����lx��ly���Լ�(x)��(y)
void GetCoordinate(double *X, double *Y, double *Z, double Xs, double Ys, double Zs, 
					double *L, double f, double *R, int k, double *x, double *y) {

	double X_ = -f*(R[0] * (X[k] - Xs) + R[3] * (Y[k] - Ys) + R[6] * (Z[k] - Zs))
		/ (R[2] * (X[k] - Xs) + R[5] * (Y[k] - Ys) + R[8] * (Z[k] - Zs));

	double Y_ = -f*(R[1] * (X[k] - Xs) + R[4] * (Y[k] - Ys) + R[7] * (Z[k] - Zs))
		/ (R[2] * (X[k] - Xs) + R[5] * (Y[k] - Ys) + R[8] * (Z[k] - Zs));

	L[2 * k + 0] = x[k] - X_;
	L[2 * k + 1] = y[k] - Y_;
}

//��ÿһ����֪���Ƶ������A
void GetNorEquation(double *X, double *Y, double *Z, double Xs, double Ys, double Zs, double *R, double f,
					double ��, double ��, double ��, int k, double *Atemp, double *x, double*y) {

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

	a14 = t_y*sin(��) - ((t_x / f)*(t_x*cos(��) - t_y*sin(��)) + f*cos(��))*cos(��);
	a15 = -f*sin(��) - (t_x / f)*(t_x*sin(��) + t_y*cos(��));
	a16 = t_y;
	a24 = -t_x*sin(��) - ((t_y / f)*(t_x*cos(��) - t_y*sin(��)) - f*sin(��))*cos(��);
	a25 = -f*cos(��) - (t_y / f)*(t_x*sin(��) + t_y*cos(��));
	a26 = -t_x;

	Atemp[0] = a11;	Atemp[1] = a12;	Atemp[2] = a13;	Atemp[3] = a14;	Atemp[4] = a15;		Atemp[5] = a16;
	Atemp[6] = a21;	Atemp[7] = a22;	Atemp[8] = a23;	Atemp[9] = a24;	Atemp[10] = a25;	Atemp[11] = a26;
}


//��Ҫ����
void Function(double *X, double *Y, double *Z, double *x, double *y, double f) {

	counter++;	//��������+1

	//��������

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

	//�����ת����

	GetR(��0, ��0, ��0, R);

	//������ϵ������L��A��ֵ
	for (int k = 0; k < n; k++) {

		//����(x)��(y)��lx��ly
		GetCoordinate(X, Y, Z, Xs0, Ys0, Zs0, L, f, R, k, x, y);

		//����������ʽ������
		GetNorEquation(X, Y, Z, Xs0, Ys0, Zs0, R, f, ��0, ��0, ��0, k, Atemp, x, y);

		//��Atemp��ֵ��ֵ��A��
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 6; j++) {
				A[(i + 2 * k) * 6 + j] = Atemp[i * 6 + j];
			}
		}
	}

	//�ⷨ����

	//�������̵ĸ���������Matrix����ʾ�Ա��������
	for (int i = 0; i < 2 * n; i++) 
	{
		for (int j = 0; j < 6; j++){
			MatA(i, j) = A[i * 6 + j];	//����A
		}
		MatL(i, 0) = L[i];				//����L
	}
	
	MatrixXd B = MatA.transpose()*MatA;
	MatX = B.inverse()*(MatA.transpose())*MatL;

	if (MatX(3,0)>0.00001|| MatX(4, 0)>0.00001|| MatX(5, 0)>0000.1) {
		Xs0 += MatX(0, 0);
		Ys0 += MatX(1, 0);
		Zs0 += MatX(2, 0);
		��0 += MatX(3, 0);
		��0 += MatX(4, 0);
		��0 += MatX(5, 0);

		Function(X, Y, Z, x, y, f);
	}
	else {
		printf("������%d�Ρ�\n\n", counter);
		printf("\n������Ϊ��\n\n%lf\t%lf\t%lf\n\n", Xs0, Ys0, Zs0);
		printf("\n��ת����RΪ��\n\n");

		//��ӡR
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
		cout <<"\n\n�����Ϊ��\n"<< endl << m_;
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
	��0 = 0;	��0 = 0;	��0 = 0;
	
	//�����ȷ����ֵ

	for (int i = 0; i < n; i++) {
		Xs0 += X[i];
		Ys0 += Y[i];
	}

	Xs0 /= 4;
	Ys0 /= 4;
	Zs0 = m*f;

	printf("��ʼֵΪ��\n\n");
	printf("Xs0 = %.2lf\t  Ys0 = %.2lf\t  Zs0 = %.2lf\n", Xs0, Ys0, Zs0);
	printf("��0 = %.0lf\t          ��0 = %.0lf\t          ��0 = %.0lf\n\n\n", ��0, ��0, ��0);


	Function(X, Y, Z, x, y, f);

	printf("\n\n\n\n");
	return 0;
}

