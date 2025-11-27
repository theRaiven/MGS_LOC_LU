// main.cpp

#include <iostream>
#include <iomanip>
#include <string>
#include <chrono>
#include <algorithm> 
#include "SparseMatrix.h"
void RunSolverTests(const SparseMatrix& A, double* x, double* xTrue, const string& methodName)
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
        std::chrono::duration<double> elapsed_time = end_time - start_time;

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

int main()
{
    setlocale(LC_ALL, "rus");

    SparseMatrix A;
    
    if (!A.ReadFromFiles())
    {
        return 1;
    }
    double* x = new double[A.n];
    double* xTrue = new double[A.n];

    for (int i = 0; i < A.n; i++)
    {
        xTrue[i] = i + 1.0;
        x[i] = 0.0;
    }
    A.PrintDense();
    RunSolverTests(A, x, xTrue, "MSG");

    cout << "----------------------------------------------------" << endl;
    cout << endl;

    RunSolverTests(A, x, xTrue, "LOS");

    delete[] xTrue;
    return 0;
}