// main.cpp

#include <iostream>
#include <iomanip>
#include <string>
#include <chrono>
#include <algorithm> 
#include "SparseMatrix.h"
void RunSolverTests(const SparseMatrix& A, real* x, real* xTrue, const string& methodName)
{
    string preconds[] = { "diag", "ilu", "none" };

    for (const string& prec : preconds)
    {
        for (int i = 0; i < A.n; i++)
        {
            x[i] = 0.0;
        }

        cout << "=== " << methodName << " (" << prec << ") ===" << endl;

        auto start_time = std::chrono::high_resolution_clock::now();

        bool success = false;
        if (methodName == "MSG") 
        {
            success = A.SolveMSG(x, xTrue, prec);
        }
        else if (methodName == "LOS")
        {
            success = A.SolveLOS(x, xTrue, prec);
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<real> elapsed_time = end_time - start_time;

        if (success)
        {
            cout << "Time: " << fixed << setprecision(6) << elapsed_time.count() << " sec" << endl;
        }
        else 
        {
            cout << "Time: " << fixed << setprecision(6) << elapsed_time.count() << " sec (Not Converged)" << endl;
        }
        cout << "\n";
    }

}

void WorkingWithHilbert()
{
    cout << "\n======================================================\n";
    cout << "      ЗАПУСК ТЕСТОВ НА МАТРИЦАХ ГИЛЬБЕРТА \n";
    cout << "======================================================\n";

    SparseMatrix A;

    A.maxIter = 100000;
    A.eps = 1e-12;

    for (int n = 5; n <= 10; n+=5)
    {
        cout << ">>> Matrix Hilbert " << n << "x" << n << " <<<" << endl;

        A.GenerateHilbertMatrix(n);

        A.PrintDense();

        real* x = new real[n];
        real* xTrue = new real[n];

        for (int i = 0; i < n; i++) xTrue[i] = i + 1.0;

        if (n <= 3) A.PrintDense();

        RunSolverTests(A, x, xTrue, "MSG");
        RunSolverTests(A, x, xTrue, "LOS");

        delete[] x;
        delete[] xTrue;

        cout << "------------------------------------------------------" << endl;
    }
}

int main()
{
    setlocale(LC_ALL, "rus");

    // Выбираем режим работы:
    int mode;
    cout << "1 - Чтение из файлов (kuslau)\n2 - Тест Гильберта (5 и 10)\nВыбор: ";
    cin >> mode;

    if (mode == 1)
    {
        SparseMatrix A;
        if (A.ReadFromFiles())
        {
            real* x = new real[A.n];
            real* xTrue = new real[A.n];
            // Тут логика для файлового режима (или заглушка)
            for (int i = 0; i < A.n; i++) { x[i] = 0; xTrue[i] = i + 1; } // заглушка xTrue
            RunSolverTests(A, x, xTrue, "MSG");
            delete[] x; delete[] xTrue;
        }
    }
    else if (mode == 2)
    {
        WorkingWithHilbert();
    }

    return 0;
}