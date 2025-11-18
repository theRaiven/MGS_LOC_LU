// MGS_LOC_LU.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include "SparseMatrix.h"


int main()
{
	setlocale(LC_ALL, "rus");
	SparseMatrix A;
	A.ReadFromFiles();
	A.PrintDense();
	double* x = new double[A.n];
	A.SolveMSG(x, "diag");
	for (int i = 0; i < A.n; i++)
	{
		cout << x[i] << ' ';
	}
}
