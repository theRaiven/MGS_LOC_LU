// SparseMatrix.h

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>
#include <limits>
#include <cmath>
using namespace std;
using real = double;

template<typename T>
constexpr T relativeEPS();

template<>
constexpr float relativeEPS<float>() { return 1e-34f; }

template<>
constexpr double relativeEPS<double>() { return 1e-300; }

template<>
constexpr long double relativeEPS<long double>() { return 1e-3000; }

const real EPS = relativeEPS<real>();
class SparseMatrix
{
private:

    real* Atb;    // A^T b
    real* r;
    real* z;
    real* p;
    real* Ap;
    real* tmp;
    real* invDiag;   // если предобусловливание diag
    real* diagAtA;

    real* inxL;
    real* inxU;

    // ---------- Основные операции ----------
    real Dot(const real* a, const real* c) const;
    real Norm2(const real* a) const;

    void Multiply(const real* x, real* y) const;
    void MultiplyTranspose(const real* v, real* y) const;

    void ComputeDiagAtA(real* diagAtA) const;
    bool InitialSolution(real* x) const;

    void PreparePreconditioner();
    void ApplyDiagPreconditioner(const real* r, real* z) const;
    void ApplyILUPreconditioner(const real* r, real* z) const;
    void SolveL(const real* y, real* x) const;
    void SolveU(const real* y, real* x) const;
    void SolveLTranspose(const real* y, real* x) const;
    void SolveUTranspose(const real* y, real* x) const;

    void WriteSolution(const string& filename, const real* x) const;
    real GetRelativeResidual(const real* x) const;
public:
    int n;
    int sizeGGL{ 0 }, sizeGGU{ 0 };
    int* ig;      //  n+1
    int* jg;      //  jgSize = ig[n]
    real* ggl;     //  jgSize (lower)
    real* ggu;     //  jgSize (upper)
    real* di;      //  n
    real* b;

    int maxIter;
    real eps;

    SparseMatrix();
    ~SparseMatrix();
    bool AllocateWorkVectors();
    // ---------- Чтение ----------
    bool ReadFromFiles(const string& kuslauFile = "kuslau.txt");
    bool ReadIntArray(const string& file, int*& masData, int size);
    bool ReadrealArray(const string& file, real*& masData, int size);
    bool ReadrealArray(const string& file, real*& masData, int size, int& GetRealSize);

    // Генератор матриц гильберта 
    void GenerateHilbertMatrix(int Size);

    // МГС и ЛОС
    bool SolveMSG(real* x, real* xTrue, const string& prec = "") const;
    bool SolveLOS(real* x, real* xTrue, const string& prec = "") const;

    // Парсер для протново представления в консольке(для проверки корректности входных данных наглядно
    void PrintDense() const;
    void PrintResults(int iter, real normR, const real* x, const real* xTrue, int n) const;
};
