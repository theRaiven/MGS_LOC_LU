// MGS_LOC_LU.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include "SparseMatrix.h"


int main()
{
	setlocale(LC_ALL, "rus");
	SparseMatrix A;
	A.ReadFromFiles();
	A.PrintDense();

	double* x1 = new double[A.n] {0,0,0,0,0,0,0,0,0,0,0,0};
	double* x2 = new double[A.n] {0,0,0,0,0,0,0,0,0,0,0,0};
	double* x3 = new double[A.n] {0,0,0,0,0,0,0,0,0,0,0,0};
	double* xTrue = new double[A.n] {1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12};

	/*A.SolveMSG(x1, xTrue, "diag");
	cout << "\n\n"; x1 = new double[A.n] {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	A.SolveMSG(x2, xTrue, "ilu");
	cout << "\n\n"; x2 = new double[A.n] {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	A.SolveMSG(x3, xTrue); 
	cout << "\n\n"; x3 = new double[A.n] {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};*/


	A.SolveLOS(x1, xTrue, "diag");
	cout << "\n\n";
	A.SolveLOS(x2, xTrue, "ilu");
	cout << "\n\n";
	A.SolveLOS(x3, xTrue);
}
