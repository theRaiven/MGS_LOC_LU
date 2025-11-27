// SparseMatrix.h

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>
using namespace std;

class SparseMatrix
{
private:

    double* Atb;    // A^T b
    double* r;
    double* z;
    double* p;
    double* Ap;
    double* tmp;
    double* invDiag;   // если предобусловливание diag
    double* diagAtA;

    double* inxL;
    double* inxU;

    // ---------- Основные операции ----------
    double Dot(const double* a, const double* c) const;
    double Norm2(const double* a) const;

    void Multiply(const double* x, double* y) const;
    void MultiplyTranspose(const double* v, double* y) const;

    void ComputeDiagAtA(double* diagAtA) const;
    bool InitialSolution(double* x) const;

    void PreparePreconditioner();
    void ApplyDiagPreconditioner(const double* r, double* z) const;
    void ApplyILUPreconditioner(const double* r, double* z) const;
    void SolveL(const double* y, double* x) const;
    void SolveU(const double* y, double* x) const;
    void SolveLTranspose(const double* y, double* x) const;
    void SolveUTranspose(const double* y, double* x) const;

    void WriteSolution(const string& filename, const double* x) const;
    double GetRelativeResidual(const double* x) const;
public:
    int n;
    int sizeGGL{ 0 }, sizeGGU{ 0 };
    int* ig;      //  n+1
    int* jg;      //  jgSize = ig[n]
    double* ggl;     //  jgSize (lower)
    double* ggu;     //  jgSize (upper)
    double* di;      //  n
    double* b;

    int maxIter;
    double eps;

    SparseMatrix();
    ~SparseMatrix();
    bool AllocateWorkVectors();
    // ---------- Чтение ----------
    bool ReadFromFiles(const string& kuslauFile = "kuslau.txt");
    bool ReadIntArray(const string& file, int*& masData, int size);
    bool ReadDoubleArray(const string& file, double*& masData, int size);
    bool ReadDoubleArray(const string& file, double*& masData, int size, int& GetRealSize);

    // МГС
    bool SolveMSG(double* x, double* xTrue, const string& prec = "") const;
    bool SolveLOS(double* x, double* xTrue, const string& prec = "") const;

    // Парсер для протново представления в консольке(для проверки корректности входных данных наглядно
    void PrintDense() const;
    void PrintResults(int iter, double normR, const double* x, const double* xTrue, int n) const;
};
