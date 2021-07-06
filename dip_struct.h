#pragma once
#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
using namespace std;

struct Arz {//разряженый формат матрицы
	int col, row;
	double value;
	Arz(const int& c, const int& r, const double& v) : col(c), row(r), value(v) {};
};

int 							Nx, Ny, NN; 			//размерность матрицы
float							a1 = 1, a2 = 1;			//вектор а для матрицы B
double							hx, hy;					//шаг матрицы
double							x0, y0_;				//нулевые
double							x_end, y_end;		
double							norma = 0, norma2 = 0;
vector <Arz>					Globrz;
vector <double> 				F, X, X2; 				//правая часть, вектор X, аналит. решение
vector <double>					lyam;					//лямбда (коэф)
vector <double>					OX, OY;					//на сколько узлов 


