#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>
using namespace std;

class SparseMatrix
{
private:

    double* Atb;       // A^T b
    double* r;
    double* z;
    double* p;
    double* Ap;
    double* tmp;
    double* tmp2;
    double* invDiag;   // если предобусловливание diag

    // Массивы для ILU разложения
    double* ilu_ggl{ nullptr };
    double* ilu_ggu{ nullptr };
    double* ilu_di{ nullptr };
    
    int* ggl_row_start{ nullptr };
    int* ggu_row_start{ nullptr };

    double Dot(const double* a, const double* b) const;

    double Norm2(const double* a) const;

    void BuildILU();
    // Решение L * y = b
    void Solve_Ly_b(const double* b, double* y) const;

    // Решение U * x = y
    void Solve_Ux_y(const double* y, double* x) const;

    // Решение L^T * x = y (Для МСГ несимм)
    void SolveLTranspose(const double* y, double* x) const;

    // Решение U^T * x = y (Для МСГ несимм)
    void SolveUTranspose(const double* y, double* x) const;
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
    bool SolveMsgLUFactorization(double* x);


    // Парсер для протново представления в консольке(для проверки корректности входных данных наглядно
    void PrintDense() const;
    
};
