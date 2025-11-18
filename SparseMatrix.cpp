#include "SparseMatrix.h"
#include <vector>
SparseMatrix::SparseMatrix() : n(0), ig(nullptr), jg(nullptr), ggl(nullptr), ggu(nullptr),
                                di(nullptr), b(nullptr),
                                Atb(nullptr), r(nullptr), z(nullptr), p(nullptr),
                                Ap(nullptr), tmp(nullptr), invDiag(nullptr)
{ }

bool SparseMatrix::AllocateWorkVectors()
{
    if (n <= 0) return false;

    delete[] Atb; delete[] r; delete[] z; delete[] p; delete[] Ap; delete[] tmp; delete[] invDiag;

    Atb = new double[n];
    r = new double[n];
    z = new double[n];
    p = new double[n];
    Ap = new double[n];
    tmp = new double[n];
    invDiag = new double[n];

    for (int i = 0; i < n; ++i) {
        Atb[i] = r[i] = z[i] = p[i] = Ap[i] = tmp[i] = invDiag[i] = 0.0;
    }

    return true;
}

SparseMatrix::~SparseMatrix()
{
    delete[] ig;  ig = nullptr;
    delete[] jg;  jg = nullptr;
    delete[] ggl; ggl = nullptr;
    delete[] ggu; ggu = nullptr;
    delete[] di;  di = nullptr;
    delete[] b;   b = nullptr;
}

bool SparseMatrix::ReadIntArray(const string& file, int*& masData, int size)
{
    if (size < 0)
    {
        cerr << "Размерность массива не действиетльна: " << endl;
        return false;
    }
    ifstream fmas(file);
    if (!fmas)
    {
        cerr << "Невозможно открыть: " << file << endl;
        return false;
    }
    masData = new int[size];
    for (int i = 0; i < size; i++)
    {
        if (!(fmas >> masData[i]))
        {
            cerr << "Ошибка чтения: " << file << endl;
            delete[] masData;
            masData = nullptr;
            return false;
        }
    }
    return true;
}
bool SparseMatrix::ReadDoubleArray(const string& file, double*& masData, int size)
{
    if (size < 0)
    {
        cerr << "Размерность массива не действиетльна: " << endl;
        return false;
    }
    ifstream fmas(file);
    if (!fmas)
    {
        cerr << "Невозможно открыть: " << file << endl;
        return false;
    }
    masData = new double[size];
    for (int i = 0; i < size; i++)
    {
        if (!(fmas >> masData[i]))
        {
            cerr << "Ошибка чтения: " << file << endl;
            delete[] masData;
            masData = nullptr;
            return false;
        }
    }
    return true;
}
bool SparseMatrix::ReadDoubleArray(const string& file, double*& masData, int size, int& GetRealSize)
{
    if (size < 0)
    {
        cerr << "Размерность массива не действиетльна: " << endl;
        return false;
    }
    ifstream fmas(file);
    if (!fmas)
    {
        cerr << "Невозможно открыть: " << file << endl;
        return false;
    }
    masData = new double[size];
    for (int i = 0; i < size; i++)
    {
        if (!(fmas >> masData[i]))
        {
            if (fmas.eof())
            {
                break;
            }
            else
            {
                cerr << "Ошибка чтения: " << file << " на элементе " << i << endl;
                delete[] masData;
                masData = nullptr;
                return false;
            }
        }
        GetRealSize++;
    }
    return true;
}
bool SparseMatrix::ReadFromFiles(const string& kuslauFile)
{
    ifstream fkus(kuslauFile.c_str());
    if (!fkus)
    {
        cerr << "Невозможно открыть " << kuslauFile << "\n";
        return false;
    }
    fkus >> n >> maxIter >> eps;
    fkus.close();

    if (!ReadIntArray("ig.txt", ig, n + 1)) return false;
    int jgSize = ig[n];
    if (jgSize <= 0)
    {
        cerr << "Неверный jgSize (ig[n]) = " << jgSize << endl;
        return false;
    }
    if (!ReadIntArray("jg.txt", jg, jgSize)) return false;
    if (!ReadDoubleArray("ggl.txt", ggl, jgSize, sizeGGL)) return false;
    if (!ReadDoubleArray("ggu.txt", ggu, jgSize, sizeGGU)) return false;
    if (!ReadDoubleArray("di.txt", di, n)) return false;
    if (!ReadDoubleArray("pr.txt", b, n)) return false;
    AllocateWorkVectors();
    return true;
}
void SparseMatrix::Multiply(const double* x, double* y) const
{
    for (int i = 0; i < n; i++)
    {
        y[i] = di[i] * x[i];
    }

    int ggl_k = 0;
    int ggu_k = 0;

    for (int i = 0; i < n; i++)
    {
        int i0 = ig[i], i1 = ig[i + 1];
        for (int p = i0; p < i1; p++)
        {
            int j = jg[p];
            if (j < i)
            {
                if (ggl_k < sizeGGL)
                {
                    y[i] += ggl[ggl_k++] * x[j];
                }
            }
            else if (j > i)
            {
                if (ggu_k < sizeGGU)
                {
                    y[i] += ggu[ggu_k++] * x[j];
                }
            }
        }
    }
}

void SparseMatrix::MultiplyTranspose(const double* v, double* y) const
{
    for (int i = 0; i < n; i++) 
    {
        y[i] = di[i] * v[i];
    }

    int ggl_k = 0;
    int ggu_k = 0;

    for (int i = 0; i < n; i++)
    {
        int i0 = ig[i], i1 = ig[i + 1];
        for (int p = i0; p < i1; p++)
        {
            int j = jg[p];
            if (j < i)
            { // A[i][j]
                if (ggl_k < sizeGGL)
                {
                    y[j] += ggl[ggl_k++] * v[i];
                }
            }
            else if (j > i)
            { // A[i][j]
                if (ggu_k < sizeGGU)
                {
                    y[j] += ggu[ggu_k++] * v[i];
                }
            }
        }
    }
}


void SparseMatrix::PrintDense() const
{
    cout << "Плотная матрица " << n << "x" << n << ":\n";
    int lowerUpperSize = (ig[n] - n) / 2;

    double** dense = new double* [n];

    for (int i = 0; i < n; i++)
    {
        dense[i] = new double[n];
        for (int j = 0; j < n; j++)
        {
            dense[i][j] = 0.0;
        }
    }

    for (int i = 0; i < n; i++)
    {
        dense[i][i] = di[i];
    }

    int ggl_index = 0;
    int ggu_index = 0;

    for (int i = 0; i < n; i++)
    {
        int start = ig[i];
        int end = ig[i + 1];

        for (int p = start; p < end; p++)
        {
            int j = jg[p];

            if (i > j)
            { 
                if (ggl_index < sizeGGL)
                {
                    dense[i][j] = ggl[ggl_index++];
                }
            }
            else if (i < j)
            {
                if (ggu_index < sizeGGU)
                {
                    dense[i][j] = ggu[ggu_index++];
                }
            }
        }
    }

    // Печать матрицы
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << setw(8) << setprecision(3) << dense[i][j] << " ";
        }
        cout << endl;
    }
}

double SparseMatrix::Dot(const double* a, const double* b) const
{
    double result = 0.0;
    for (int i = 0; i < n; ++i)
    {
        result += a[i] * b[i];
    }
    return result;
}

double SparseMatrix::Norm2(const double* a) const
{
    return sqrt(Dot(a, a));
}

bool SparseMatrix::SolveMSG(double* x, const string& prec) const
{
    MultiplyTranspose(b, Atb);

    for (int i = 0; i < n; ++i)
    {
        x[i] = 0.0;
        r[i] = Atb[i];
    }

    for (int i = 0; i < n; ++i) z[i] = r[i];

    for (int i = 0; i < n; ++i) p[i] = z[i];

    double rz_old = Dot(r, z);
    double normAtb = Norm2(Atb);
    if (normAtb < 1e-16) normAtb = 1.0;

    for (int iter = 0; iter < maxIter; ++iter)
    {
        Multiply(p, tmp);
        MultiplyTranspose(tmp, Ap);

        double pAp = Dot(p, Ap);
        if (fabs(pAp) < 1e-30)
        {
            cout << "Решение сошлось (pAp ~ 0) на итерации " << iter << endl;
            return true;
        }

        double alpha = rz_old / pAp;

        for (int i = 0; i < n; ++i)
        {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        double relRes = Norm2(r) / normAtb;
        if (relRes < eps)
        {
            cout << "Решение сошлось за " << iter + 1 << " итераций.\n";
            return true;
        }

        for (int i = 0; i < n; ++i) z[i] = r[i];

        double rz_new = Dot(r, z);

        double beta = rz_new / rz_old;
        rz_old = rz_new;

        for (int i = 0; i < n; ++i)
        {
            p[i] = z[i] + beta * p[i];
        }
    }

    cerr << "Решение не сошлось за " << maxIter << " итераций.\n";
    return false;
}