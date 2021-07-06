#include "dip_func.cpp"

int main() {

	read();

	Nx = ceil(((x_end - x0) / hx) + 1);
	Ny = ceil(((y_end - y0_) / hy) + 1);

	NN = Nx * Ny;
	
	initialization(Nx, Ny);		

	lyam.resize(NN);
	lyambda();
	
	OX_OY(x0, x_end, hx, y0_, y_end, hy, OX, OY);

	//Заполнение глобальной матрицы
	for (int j = 1; j < Ny; j++)
	{
		for (int i = 1; i < Nx; i++)
		{
			Globrz[find(i + (j - 1) * Nx, i + (j - 1) * Nx)].value		+= Ann(OX[i - 1], OX[i], OY[j - 1], OY[j]) * lyam[i] + B11(OX[i - 1], OX[i], OY[j - 1], OY[j], a1, a2);
			Globrz[find(i + (j - 1) * Nx, i + 1 + (j - 1) * Nx)].value	+= A12(OX[i - 1], OX[i], OY[j - 1], OY[j]) * lyam[i] + B12(OX[i - 1], OX[i], OY[j - 1], OY[j], a1, a2);
			Globrz[find(i + (j - 1) * Nx, i + j * Nx)].value			+= A13(OX[i - 1], OX[i], OY[j - 1], OY[j]) * lyam[i] + B13(OX[i - 1], OX[i], OY[j - 1], OY[j], a1, a2);
			Globrz[find(i + (j - 1) * Nx, i + 1 + j * Nx)].value		+= A14(OX[i - 1], OX[i], OY[j - 1], OY[j]) * lyam[i] + B14(OX[i - 1], OX[i], OY[j - 1], OY[j], a1, a2);

			Globrz[find(i + 1 + (j - 1) * Nx, i + (j - 1) * Nx)].value		+= A12(OX[i - 1], OX[i], OY[j - 1], OY[j]) * lyam[i] + B21(OX[i - 1], OX[i], OY[j - 1], OY[j], a1, a2);
			Globrz[find(i + 1 + (j - 1) * Nx, i + 1 + (j - 1) * Nx)].value	+= Ann(OX[i - 1], OX[i], OY[j - 1], OY[j]) * lyam[i] + B22(OX[i - 1], OX[i], OY[j - 1], OY[j], a1, a2);
			Globrz[find(i + 1 + (j - 1) * Nx, i + j * Nx)].value			+= A14(OX[i - 1], OX[i], OY[j - 1], OY[j]) * lyam[i] + B23(OX[i - 1], OX[i], OY[j - 1], OY[j], a1, a2);
			Globrz[find(i + 1 + (j - 1) * Nx, i + 1 + j * Nx)].value		+= A13(OX[i - 1], OX[i], OY[j - 1], OY[j]) * lyam[i] + B24(OX[i - 1], OX[i], OY[j - 1], OY[j], a1, a2);

			Globrz[find(i + j * Nx, i + (j - 1) * Nx)].value			+= A13(OX[i - 1], OX[i], OY[j - 1], OY[j]) * lyam[i] + B31(OX[i - 1], OX[i], OY[j - 1], OY[j], a1, a2);
			Globrz[find(i + j * Nx, i + 1 + (j - 1) * Nx)].value		+= A14(OX[i - 1], OX[i], OY[j - 1], OY[j]) * lyam[i] + B32(OX[i - 1], OX[i], OY[j - 1], OY[j], a1, a2);
			Globrz[find(i + j * Nx, i + j * Nx)].value					+= Ann(OX[i - 1], OX[i], OY[j - 1], OY[j]) * lyam[i] + B33(OX[i - 1], OX[i], OY[j - 1], OY[j], a1, a2);
			Globrz[find(i + j * Nx, i + 1 + j * Nx)].value				+= A12(OX[i - 1], OX[i], OY[j - 1], OY[j]) * lyam[i] + B34(OX[i - 1], OX[i], OY[j - 1], OY[j], a1, a2);

			Globrz[find(i + 1 + j * Nx, i + (j - 1) * Nx)].value		+= A14(OX[i - 1], OX[i], OY[j - 1], OY[j]) * lyam[i] + B41(OX[i - 1], OX[i], OY[j - 1], OY[j], a1, a2);
			Globrz[find(i + 1 + j * Nx, i + 1 + (j - 1) * Nx)].value	+= A13(OX[i - 1], OX[i], OY[j - 1], OY[j]) * lyam[i] + B42(OX[i - 1], OX[i], OY[j - 1], OY[j], a1, a2);
			Globrz[find(i + 1 + j * Nx, i + j * Nx)].value				+= A12(OX[i - 1], OX[i], OY[j - 1], OY[j]) * lyam[i] + B43(OX[i - 1], OX[i], OY[j - 1], OY[j], a1, a2);
			Globrz[find(i + 1 + j * Nx, i + 1 + j * Nx)].value			+= Ann(OX[i - 1], OX[i], OY[j - 1], OY[j]) * lyam[i] + B44(OX[i - 1], OX[i], OY[j - 1], OY[j], a1, a2);
		}
	}

	right_side(F);

	//вывод F
	result_F();

	//проверка матрицы после учета краевых условий
	result_G("output_G_F");

	X.resize(F.size());

	gauss(Globrz, F, X);

	//проверка матрицы на правильность приведения к треугольному виду
	result_G("output_G");

	//вывод X первая строка 
	result_X_1st();

	//вывод X  
	result_X();

	//вывод X по строкам
	result_Xkp();

	//аналит. решение
	X2.resize(NN);
	ofstream fout;
	fout.open("X2.txt");
	for (int i = 0, y = 0; y < Ny; y++)
	{
		for (int x = 0; x < Nx; x++)
		{
			X2[i] = FuncU(OX[x], OY[y]);
			fout << X2[i] << endl;
			i++;
		}
	}

	fout << endl;
	for (int i = 0; i < NN; i++)
	{
		fout << abs(X[i] - X2[i]) << endl;
	}
	fout.close();

	double norma = 0, norma2 = 0;

	// вычисление нормы
	for (int i = 0; i < NN; i++)		
	{
		norma += pow((X[i] - X2[i]), 2);
		norma2 += abs(pow(X2[i], 2));
	}
	cout << (sqrt(norma)) / (pow(NN, 2)) << endl;

	cout << NN << endl;
	cout << Nx << " * " << Ny << endl;
	cout << x0 << "__" << x_end << " step: " << hx << endl;
	cout << y0_ << "__" << y_end << " step: " << hy << endl;
	for (int i = 0; OX.size() > i; i++) { cout << OX[i] << "  "; } cout << endl;
	for (int i = 0; OY.size() > i; i++) { cout << OY[i] << "  "; } cout << endl;
	cout << "The end." << endl;
	return 0;
}
