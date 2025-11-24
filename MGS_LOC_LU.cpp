// MGS_LOC_LU.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include "SparseMatrix.h"


int main()
{
	setlocale(LC_ALL, "rus");
	SparseMatrix A;
	A.ReadFromFiles();
	A.PrintDense();
	double* x1 = new double[A.n];
	double* x2 = new double[A.n];
	double* x3 = new double[A.n];
	

	/*A.SolveMSG(x1, "diag");
	for (int i = 0; i < A.n; i++)
	{
		cout << x1[i] << ' ';
	}
	cout << "\n\n\n";
	A.SolveMSG(x2, "ilu");
	for (int i = 0; i < A.n; i++)
	{
		cout << x2[i] << ' ';
	}*/
	cout << "\n\n\n";
	A.SolveMsgLUFactorization(x3);
	for (int i = 0; i < A.n; i++)
	{
		cout << x3[i] << ' ';
	}
	A.PrintDense();
}
