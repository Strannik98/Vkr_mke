#include "dip_struct.h"

//double FuncU(double x, double y) { return (exp(x * y)); };
//double FuncF(double x, double y) { return -(2 * exp(x * y)); };

double FuncU(double x, double y) { return (x * (1 - x) * y * (1 - y)); };
double FuncF(double x, double y) { return (-(2 * ((x * (x - 1)) + (y * (y - 1))))); };

void initialization(int r, int c) {			//зануление значений (подготовка) глоб. матрицы
	int i = 0, s = 1;
	while (s < Ny) {
		for (; i < (Nx * s) - 1; i++) 
		{								
			F.push_back(0);											//вектор правой части	
			if (s > 1) {											
				if (i % Nx != 0) Globrz.push_back(Arz(i + 1, i - Nx, 0));	
				Globrz.push_back(Arz(i + 1, i + 1 - Nx, 0));				
				Globrz.push_back(Arz(i + 1, i + 2 - Nx, 0));				
			}		
			if (i % Nx != 0) Globrz.push_back(Arz(i + 1, i, 0));				
			Globrz.push_back(Arz(i + 1, i + 1, 0)); 					
			Globrz.push_back(Arz(i + 1, i + 2, 0));					
			if (i % Nx != 0) Globrz.push_back(Arz(i + 1, i + Nx, 0));	
			Globrz.push_back(Arz(i + 1, i + 1 + Nx, 0));				
			Globrz.push_back(Arz(i + 1, i + 1 + Nx + 1, 0));			
		}
		F.push_back(0);
		if (s > 1) {
			if (i % Nx != 0) Globrz.push_back(Arz(i + 1, i - Nx, 0));	
			Globrz.push_back(Arz(i + 1, i + 1 - Nx, 0));				
		}
		Globrz.push_back(Arz(i + 1, i, 0));							
		Globrz.push_back(Arz(i + 1, i + 1, 0));						
		Globrz.push_back(Arz(i + 1, i + Nx, 0));						
		Globrz.push_back(Arz(i + 1, i + 1 + Nx, 0));
		i++;	s++;
	}
	for (; i < Nx * Ny - 1; i++) {									//последний ряд
		F.push_back(0);
		if (s > 1) {
			if (i % Nx != 0) Globrz.push_back(Arz(i + 1, i - Nx, 0));	
			Globrz.push_back(Arz(i + 1, i + 1 - Nx, 0));				
			Globrz.push_back(Arz(i + 1, i + 2 - Nx, 0));				
		}
		if (i % Nx != 0) Globrz.push_back(Arz(i + 1, i, 0));		
		Globrz.push_back(Arz(i + 1, i + 1, 0));					
		Globrz.push_back(Arz(i + 1, i + 2, 0));					
	}
	F.push_back(0);												//последний кэ
	if (i % Nx != 0) Globrz.push_back(Arz(i + 1, i - Nx, 0));	
	Globrz.push_back(Arz(i + 1, i + 1 - Nx, 0));
	if (i % Nx != 0) Globrz.push_back(Arz(i + 1, i, 0));
	Globrz.push_back(Arz(i + 1, i + 1, 0));					
}

void lyambda() {
	for (int i = 0; i < NN; i++) 
		lyam[i] = 0.1;
}

int find(int i, int j) {										//поиск элемента глобальной матрицы
	int k = 0;
	if (i == NN && j == NN) return (Globrz.size() - 1);
	for (; (Globrz[k].col != i) && (Globrz.size() > k); k++);
	for (; (Globrz[k].row != j) && (Globrz[k].col == i) && (Globrz.size()-1 > k); k++);	
	if (Globrz[k].col == i && Globrz[k].row == j)
		return k;
	else	{
		return -1;
	}
}

int find(int i) {	//по координате диагонального элемента ищу место в глобальной матрице
	int k = 0;
	for (; (Globrz[k].col != i) && (Globrz.size()-1 > k); k++);
	if (Globrz[k].col == i)
		return k;
	else return -1;
}

//уравнения
	double Ann(double x1, double x2, double y1, double y2) {	/*11=22,33,44*/
		return ((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)) / (3 * (x1 - x2) * (y1 - y2));
	}
	double A12(double x1, double x2, double y1, double y2) {/*12=34*/
		return ((x1 - x2) / (6 * (y1 - y2))) + ((y2 - y1) / (3 * (x1 - x2)));
	}
	double A13(double x1, double x2, double y1, double y2) {/*13=24*/
		return (((x2 - x1) / (3 * (y1 - y2))) + ((y1 - y2) / (6 * (x1 - x2))));
	}
	double A14(double x1, double x2, double y1, double y2) {/*14=23*/
		return -(((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)) / (6 * (x1 - x2) * (y1 - y2)));
	}

	double B11(double x1, double x2, double y1, double y2, float a1, float a2) {
		return -(a1 * (y2 - y1) + a2 * (x2 - x1)) / 6;
	}
	double B12(double x1, double x2, double y1, double y2, float a1, float a2) {
		return (a1 * (y2 - y1) + a2 * (x2 - x1) / 2) / 6;
	}
	double B13(double x1, double x2, double y1, double y2, float a1, float a2) {
		return (a1 * -(y2 - y1) / 2 + a2 * (x2 - x1)) / 6 ;
	}
	double B14(double x1, double x2, double y1, double y2, float a1, float a2) {
		return (a1 * (y2 - y1) / 2 + a2 * (x2 - x1) / 2) / 6;
	}

	double B21(double x1, double x2, double y1, double y2, float a1, float a2) {
		return -(a1 * (y2 - y1) + a2 * -(x2 - x1) / 2) / 6;
	}
	double B22(double x1, double x2, double y1, double y2, float a1, float a2) {
		return (a1 * (y2 - y1) + a2 * -(x2 - x1)) / 6;
	}
	double B23(double x1, double x2, double y1, double y2, float a1, float a2) {
		return (a1 * -(y2 - y1) / 2 + a2 * (x2 - x1) / 2) / 6 ;
	}
	double B24(double x1, double x2, double y1, double y2, float a1, float a2) {
		return (a1 * (y2 - y1) / 2 + a2 * (x2 - x1)) / 6;
	}

	double B31(double x1, double x2, double y1, double y2, float a1, float a2) {
		return -(a1 * (y2 - y1) / 2 + a2 * (x2 - x1)) / 6;
	}
	double B32(double x1, double x2, double y1, double y2, float a1, float a2) {
		return (a1 * (y2 - y1) / 2 + a2 * -(x2 - x1) / 2) / 6;
	}
	double B33(double x1, double x2, double y1, double y2, float a1, float a2) {
		return (a1 * -(y2 - y1) + a2 * (x2 - x1)) / 6;
	}
	double B34(double x1, double x2, double y1, double y2, float a1, float a2) {
		return (a1 * (y2 - y1) + a2 * (x2 - x1) / 2) / 6;
	}

	double B41(double x1, double x2, double y1, double y2, float a1, float a2) {
		return -(a1 * (y2 - y1) / 2 + a2 * (x2 - x1) / 2) / 6;
	}
	double B42(double x1, double x2, double y1, double y2, float a1, float a2) {
		return (a1 * (y2 - y1) / 2 + a2 * -(x2 - x1)) / 6;
	}
	double B43(double x1, double x2, double y1, double y2, float a1, float a2) {
		return (a1 * -(y2 - y1) + a2 * (x2 - x1)) / 6;
	}
	double B44(double x1, double x2, double y1, double y2, float a1, float a2) {
		return (a1 * (y2 - y1) + a2 * (x2 - x1)) / 6;
	}

void read() {
	ifstream fin("input.txt");
	fin >> x0 >> x_end >> hx >> y0_ >> y_end >> hy;
	fin.close();
}

void OX_OY(double a, double b, double hx, double c, double d, double hy,
	vector <double> &OX, vector <double> &OY)	// заполняю значения на осях
{
	double x = a;
	while (x < b)
	{
		OX.push_back(x);
		x = x + hx;
	}
	if (x >= b)	OX.push_back(b);

	double y = c;
	while (y < d)
	{
		OY.push_back(y);
		y = y + hy;
	}
	if(y >= d)	OY.push_back(d);
}

void add_elem(int j, int l, int I, double val, vector <Arz>& Globrz) {
	auto begin = Globrz.begin();
	int  k = j;

	for (; find(l, k) < 0 && k <= NN; k++); //ищем следующее значение 

	if (k <= NN)
		Globrz.insert(begin + find(l, k), (Arz(l, j, -val)));
	else
		if (l + 1 <= NN) 
			Globrz.insert(begin + find(l + 1), (Arz(l, j, -val)));
}

void right_side(vector <double>& F) {				//вектор правой части и учёт краевых условий глобальной матриц
	
	//for (int j = 0; j < Ny - 1; j++)
	//{
	//	for (int i = 0; i < Nx - 1; i++)
	//	{
	//		int uz1 = i + j * Nx;
	//		int uz2 = uz1 + 1;
	//		int uz3 = uz1 + Nx;
	//		int uz4 = uz3 + 1;
	//		F[uz1] += ((hx * hy) / 36) * (4 * FuncF(OX[i], OY[j]) + 2 * FuncF(OX[i + 1], OY[j]) + 2 * FuncF(OX[i], OY[j + 1]) + FuncF(OX[i + 1], OY[j + 1]));
	//		F[uz2] += ((hx * hy) / 36) * (2 * FuncF(OX[i], OY[j]) + 4 * FuncF(OX[i + 1], OY[j]) + FuncF(OX[i], OY[j + 1]) + 2 * FuncF(OX[i + 1], OY[j + 1]));
	//		F[uz3] += ((hx * hy) / 36) * (2 * FuncF(OX[i], OY[j]) + FuncF(OX[i + 1], OY[j]) + 4 * FuncF(OX[i], OY[j + 1]) + 2 * FuncF(OX[i + 1], OY[j + 1]));
	//		F[uz4] += ((hx * hy) / 36) * (FuncF(OX[i], OY[j]) + 2 * FuncF(OX[i + 1], OY[j]) + 2 * FuncF(OX[i], OY[j + 1]) + 4 * FuncF(OX[i + 1], OY[j + 1]));
	//	}
	//}

	//int n = 1; //первые краевые
	//for (int j = 0; j < Ny; j++)
	//{
	//	for (int i = 0; i < Nx; i++)
	//	{
	//		if (i == 0 || j == 0 || i == Nx - 1 || j == Ny - 1)
	//		{
	//			F[n - 1] = FuncU(OX[i], OY[j]);
	//			for (int n2 = 1; n2 < NN; n2++)
	//			{
	//				if (n == n2)
	//				{
	//					Globrz[find(n, n2)].value = 1;
	//				}
	//				else
	//				{
	//					if (find(n, n2) >= 0)
	//						/*Globrz.erase()*/
	//						Globrz[find(n, n2)].value = 0;
	//				}
	//			}
	//		}
	//		n++;
	//	}
	//}
	
	//термоусловия слева справа:
	int i = 1; //строка
	for (int n = 0; n < Ny; n++)
	{
		F[n * Ny] = 50;			//левые
		F[n * Ny + Nx -1] = 10;	//правые	
		for (int j = 0; j < NN; j++)		//для левых 
		{
			if (i - 1 == j)
			{
				Globrz[find(i , j + 1)].value = 1;
			}
			else
			{
				if (find(i, j + 1) > 0)
					Globrz[find(i, j + 1)].value = 0;
			}
		}
		i += Nx - 1;
		for (int j = 0; j < NN; j++)		//для правых
		{
			if (i - 1 == j)
			{
				Globrz[find(i, j + 1)].value = 1;
			}
			else
			{
				if (find(i, j + 1) > 0)
					Globrz[find(i, j + 1)].value = 0;
			}
		}
		i++;
	}	
}

void gauss(vector <Arz> &Globrz, vector <double> &F, vector <double> &X)	//i номер верхней строки | l перебор строк ниже | j перебор столбцов 
{
	double coefficient;
	for (int i = 1; i <= NN-1; i++)	//берем строку
	{
		for (int l = i + 1; l <= NN; l++)	//берем строки ниже нее
		{
			if (find(l, i) >= 0)	//если поддиагональное значение существует, то изменяем строку
			{
				coefficient = Globrz[find(l, i)].value / Globrz[find(i, i)].value;	//
				if (coefficient != 0)
				{
					for (int j = i; j <= NN; j++)	//начиная с поддиагонального элемента
					{	
						int I = find(i, j);
						int L = find(l, j);		
						double val = coefficient * Globrz[I].value;
						if (L >= 0)
						{
							if (I >= 0) //если значения существуют 
								if (val != 0)
									Globrz[L].value -= val;
						}
						else //если нет элемента, под тем который есть, то добавляем
						{
							if (I >= 0 && val != 0) {				//найти значение в которое положить	
								add_elem(j, l, I, val, Globrz);
							}
						}
					}
					F[l - 1] -= coefficient * F[i - 1];
				}
			}
		}
	}
	for (int i = NN - 1; i >= 0; i--) {
		for (int j = NN - 1; j > i; j--)
		{
			int IJ = find(i + 1, j + 1);
			if (IJ >= 0)
				F[i] -= Globrz[IJ].value * X[j];
		}
		X[i] = F[i] / Globrz[find(i + 1, i + 1)].value;
	}
}

void result_F() {
	ofstream fout;
	fout.open("output_F.txt");
	fout << "Правая часть:" << endl;
	for (int i = 0; i < F.size(); i++) {
		fout << i + 1 << "\t";
		fout << F[i] << endl;
	}
	fout.close();
}
void result_G(string nof) {
	ofstream fout;
	fout.open(nof);
	fout << "col\trow\tvalue" << endl;
	for (int i = 0; i < Globrz.size(); i++) {
		fout << Globrz[i].col << "\t";
		fout << Globrz[i].row << "\t";
		fout << Globrz[i].value << endl;
	}
	fout << Nx << "*" << Ny << endl;
	fout << x0 << " step: " << hx << "\n" << y0_ << " step: " << hy << endl;
	fout.close();
}

void result_X_1st() {
	//проверка X первая строка
	ofstream fout;
	fout.open("X_1st.txt");
	//fout << "X:" << endl;
	for (int i = 0; i < Nx; i++) {
		fout << i + 1 << "\t";
		fout << X[i] << endl;
	}
	fout.close();
}
void result_X() {
	//проверка X
	ofstream fout;
	fout.open("X.txt");
	//fout << "X:" << endl;
	for (int i = 0; i < NN; i++) {
		fout << i + 1 << "\t";
		fout << X[i] << endl;
	}
	fout.close();
}
void result_Xkp() {
	//проверка X по строкам
	ofstream fout;
	fout.open("Xkp.txt");
	//fout << "X:" << endl;
	int i = 0, l = 0;
	for (; i < X.size(); i++, l++) {
		fout << l + 1 << "\t";
		fout << X[i] << endl;
		if (l == Nx - 1) l = -1;
	}
	fout.close();
}