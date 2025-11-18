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

    double Dot(const double* a, const double* c) const;
    double Norm2(const double* a) const;
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
    // ---------- Основные операции ----------
    void Multiply(const double* x, double* y) const;
    void MultiplyTranspose(const double* v, double* y) const;
    // МГС
    bool SolveMSG(double* x, const string& prec) const;


    // Парсер для протново представления в консольке(для проверки корректности входных данных наглядно
    void PrintDense() const;
    
};
