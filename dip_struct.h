#pragma once
#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
using namespace std;

struct Arz {//���������� ������ �������
	int col, row;
	double value;
	Arz(const int& c, const int& r, const double& v) : col(c), row(r), value(v) {};
};

int 							Nx, Ny, NN; 			//����������� �������
float							a1 = 1, a2 = 1;			//������ � ��� ������� B
double							hx, hy;					//��� �������
double							x0, y0_;				//�������
double							x_end, y_end;		
double							norma = 0, norma2 = 0;
vector <Arz>					Globrz;
vector <double> 				F, X, X2; 				//������ �����, ������ X, ������. �������
vector <double>					lyam;					//������ (����)
vector <double>					OX, OY;					//�� ������� ����� 


